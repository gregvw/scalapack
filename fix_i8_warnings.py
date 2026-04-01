#!/usr/bin/env python3
"""
Conservative fixer for printf-family format warnings.

Reads GCC/Clang warning output from stdin or a file, finds format mismatches,
and rewrites only the warned conversion specifiers to use a project macro
(default: SLINT_FMT).

Typical usage:
  make 2> build.log
  python3 tools/fix_format_warnings.py build.log > format_fixes.diff
  git apply --check format_fixes.diff
  git apply format_fixes.diff

Or apply directly:
  python3 tools/fix_format_warnings.py build.log --apply

Intended use case:
  - ILP64 / _I8 migration where %d is wrong for long long / wide integer fields
  - targeted replacements driven by compiler diagnostics

Conservative limitations:
  - only fixes printf/fprintf/sprintf/snprintf
  - only when the format argument is a plain or concatenated string literal
  - skips positional parameters, %* width/precision, and non-literal format args
  - skips calls it cannot parse safely
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


PRINTF_SPECS = {"d", "i", "u", "o", "x", "X"}
TARGET_FUNCS = {
    "printf": 0,
    "fprintf": 1,
    "sprintf": 1,
    "snprintf": 2,
}
STRING_LITERAL_RE = re.compile(r'"(?:\\.|[^"\\])*"', re.DOTALL)

# GCC/Clang-style warning lines, e.g.:
# foo.c:123:45: warning: format '%d' expects argument of type 'int', but argument 4 has type 'long long int'
WARNING_RE = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s*warning:
    \s*format
    \s+
    (?P<fmtquote>'[^']+'|`[^`]+`)
    \s+
    expects\ argument\ of\ type\ .+?
    ,\ but\ argument\ (?P<argnum>\d+) \s+ has\ type\ \s+ (?P<argtype>.+?)
    (?:
      \s*\[-Wformat[^\]]*\]
    )?
    \s*$
    """,
    re.VERBOSE,
)

# Alternate wording sometimes seen from Clang.
WARNING_RE2 = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s*warning:
    .+?
    format\ specifies\ type\ .+?
    \ but\ the\ argument\ has\ type\ \s+ (?P<argtype>.+?)
    \s*
    (?:
      \[-Wformat[^\]]*\]
    )?
    \s*$
    """,
    re.VERBOSE,
)

NOTE_RE = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s*note:
    \s*format\ string\ is\ defined\ here
    """,
    re.VERBOSE,
)


@dataclass
class WarningSite:
    file: str
    line: int
    col: int
    argnum: Optional[int]
    raw_spec: Optional[str]
    argtype: str
    raw_line: str


@dataclass
class CallSite:
    func: str
    start_idx: int
    end_idx: int
    text: str
    args: List[str]
    format_arg_index: int
    first_value_arg_index: int


@dataclass
class SpecInfo:
    spec_index: int          # counts actual conversion specs excluding %%
    start: int               # start index in format-content
    end: int                 # end index in format-content (exclusive)
    conv: str                # d, i, u, ...
    raw: str                 # e.g. "%d"
    flags_width_prec_len: str


def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def write_text(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(text)


def parse_warning_lines(lines: List[str]) -> List[WarningSite]:
    sites: List[WarningSite] = []
    for line in lines:
        m = WARNING_RE.match(line.rstrip("\n"))
        if m:
            raw_spec = m.group("fmtquote")[1:-1]
            argnum = int(m.group("argnum"))
            sites.append(
                WarningSite(
                    file=m.group("file"),
                    line=int(m.group("line")),
                    col=int(m.group("col")),
                    argnum=argnum,
                    raw_spec=raw_spec,
                    argtype=m.group("argtype").strip(),
                    raw_line=line.rstrip("\n"),
                )
            )
            continue

        m2 = WARNING_RE2.match(line.rstrip("\n"))
        if m2:
            # Less useful because it lacks argument number and specific spec.
            sites.append(
                WarningSite(
                    file=m2.group("file"),
                    line=int(m2.group("line")),
                    col=int(m2.group("col")),
                    argnum=None,
                    raw_spec=None,
                    argtype=m2.group("argtype").strip(),
                    raw_line=line.rstrip("\n"),
                )
            )
    return sites


def line_to_index(text: str, line_no: int) -> int:
    if line_no <= 1:
        return 0
    pos = 0
    cur = 1
    while cur < line_no:
        nxt = text.find("\n", pos)
        if nxt == -1:
            return len(text)
        pos = nxt + 1
        cur += 1
    return pos


def find_call_around(text: str, line_no: int) -> Optional[CallSite]:
    """
    Find a printf-family call that encloses or is near the given line.
    Conservative heuristic:
      - search backward a limited window for target func + '('
      - from there, parse until matching ')' and terminating ';'
    """
    anchor = line_to_index(text, line_no)
    window_start = max(0, anchor - 4000)
    snippet = text[window_start: min(len(text), anchor + 4000)]

    candidates = []
    for func in TARGET_FUNCS:
        for m in re.finditer(r"\b" + re.escape(func) + r"\s*\(", snippet):
            candidates.append((m.start(), func))
    if not candidates:
        return None

    # Prefer the last candidate before or near anchor.
    abs_candidates = [(window_start + s, func) for s, func in candidates if window_start + s <= anchor + 500]
    if not abs_candidates:
        return None
    abs_candidates.sort(key=lambda x: x[0], reverse=True)

    for start, func in abs_candidates:
        call = parse_call_at(text, start, func)
        if call and call.start_idx <= anchor <= call.end_idx + 1:
            return call
    # fallback: return nearest parsable one above the line
    for start, func in abs_candidates:
        call = parse_call_at(text, start, func)
        if call:
            return call
    return None


def parse_call_at(text: str, start_idx: int, func: str) -> Optional[CallSite]:
    m = re.match(r"\b" + re.escape(func) + r"\s*\(", text[start_idx:])
    if not m:
        return None
    open_paren = start_idx + m.end() - 1

    end_paren = find_matching_paren(text, open_paren)
    if end_paren is None:
        return None

    semi = skip_ws(text, end_paren + 1)
    if semi >= len(text) or text[semi] != ";":
        return None

    inner = text[open_paren + 1:end_paren]
    args = split_top_level_args(inner)
    if args is None:
        return None

    fmt_idx = TARGET_FUNCS[func]
    if len(args) <= fmt_idx:
        return None

    value_start = fmt_idx + 1
    return CallSite(
        func=func,
        start_idx=start_idx,
        end_idx=semi,
        text=text[start_idx:semi + 1],
        args=args,
        format_arg_index=fmt_idx,
        first_value_arg_index=value_start,
    )


def find_matching_paren(text: str, open_idx: int) -> Optional[int]:
    depth = 0
    i = open_idx
    n = len(text)
    in_str = False
    in_chr = False
    in_line_comment = False
    in_block_comment = False
    escape = False

    while i < n:
        ch = text[i]
        nxt = text[i + 1] if i + 1 < n else ""

        if in_line_comment:
            if ch == "\n":
                in_line_comment = False
            i += 1
            continue

        if in_block_comment:
            if ch == "*" and nxt == "/":
                in_block_comment = False
                i += 2
                continue
            i += 1
            continue

        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_str = False
            i += 1
            continue

        if in_chr:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == "'":
                in_chr = False
            i += 1
            continue

        if ch == "/" and nxt == "/":
            in_line_comment = True
            i += 2
            continue
        if ch == "/" and nxt == "*":
            in_block_comment = True
            i += 2
            continue
        if ch == '"':
            in_str = True
            i += 1
            continue
        if ch == "'":
            in_chr = True
            i += 1
            continue

        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return None


def skip_ws(text: str, idx: int) -> int:
    while idx < len(text) and text[idx].isspace():
        idx += 1
    return idx


def split_top_level_args(s: str) -> Optional[List[str]]:
    args: List[str] = []
    start = 0
    depth_paren = depth_brack = depth_brace = 0
    in_str = in_chr = False
    in_line_comment = in_block_comment = False
    escape = False
    i = 0
    n = len(s)

    while i < n:
        ch = s[i]
        nxt = s[i + 1] if i + 1 < n else ""

        if in_line_comment:
            if ch == "\n":
                in_line_comment = False
            i += 1
            continue

        if in_block_comment:
            if ch == "*" and nxt == "/":
                in_block_comment = False
                i += 2
                continue
            i += 1
            continue

        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_str = False
            i += 1
            continue

        if in_chr:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == "'":
                in_chr = False
            i += 1
            continue

        if ch == "/" and nxt == "/":
            in_line_comment = True
            i += 2
            continue
        if ch == "/" and nxt == "*":
            in_block_comment = True
            i += 2
            continue
        if ch == '"':
            in_str = True
            i += 1
            continue
        if ch == "'":
            in_chr = True
            i += 1
            continue

        if ch == "(":
            depth_paren += 1
        elif ch == ")":
            depth_paren -= 1
        elif ch == "[":
            depth_brack += 1
        elif ch == "]":
            depth_brack -= 1
        elif ch == "{":
            depth_brace += 1
        elif ch == "}":
            depth_brace -= 1
        elif ch == "," and depth_paren == depth_brack == depth_brace == 0:
            args.append(s[start:i].strip())
            start = i + 1

        i += 1

    tail = s[start:].strip()
    if tail:
        args.append(tail)
    return args


def extract_concatenated_string_literal(expr: str) -> Optional[Tuple[str, List[Tuple[int, int, str]]]]:
    """
    Returns:
      (decoded_content, pieces)
    where pieces is a list of (start_idx, end_idx, raw_literal_text) relative to expr.
    Only succeeds if expr consists solely of string literals plus whitespace/comments.
    """
    pieces: List[Tuple[int, int, str]] = []
    i = 0
    n = len(expr)

    def skip_junk(j: int) -> int:
        while j < n:
            if expr[j].isspace():
                j += 1
                continue
            if expr.startswith("//", j):
                nl = expr.find("\n", j)
                return n if nl == -1 else nl + 1
            if expr.startswith("/*", j):
                k = expr.find("*/", j + 2)
                return n if k == -1 else k + 2
            break
        return j

    i = skip_junk(i)
    while i < n:
        m = STRING_LITERAL_RE.match(expr, i)
        if not m:
            return None
        pieces.append((m.start(), m.end(), m.group(0)))
        i = skip_junk(m.end())

    if not pieces:
        return None

    content = "".join(decode_c_string_literal(raw) for _, _, raw in pieces)
    return content, pieces


def decode_c_string_literal(raw: str) -> str:
    assert raw.startswith('"') and raw.endswith('"')
    body = raw[1:-1]
    # Conservative enough for diagnostics strings; unicode escapes not required here.
    return bytes(body, "utf-8").decode("unicode_escape")


def encode_c_string_literal(content: str) -> str:
    content = content.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
    return '"' + content + '"'


def parse_format_specs(fmt: str) -> Optional[List[SpecInfo]]:
    """
    Parse C printf conversion specs in a conservative way.
    Skips if it sees positional args or '*' width/precision.
    """
    specs: List[SpecInfo] = []
    i = 0
    spec_index = 0
    n = len(fmt)

    while i < n:
        if fmt[i] != "%":
            i += 1
            continue
        if i + 1 < n and fmt[i + 1] == "%":
            i += 2
            continue

        start = i
        i += 1

        # positional args unsupported
        m = re.match(r"\d+\$", fmt[i:])
        if m:
            return None

        # flags
        while i < n and fmt[i] in "-+ #0'":
            i += 1

        # width
        if i < n and fmt[i] == "*":
            return None
        while i < n and fmt[i].isdigit():
            i += 1

        # precision
        if i < n and fmt[i] == ".":
            i += 1
            if i < n and fmt[i] == "*":
                return None
            while i < n and fmt[i].isdigit():
                i += 1

        # length modifiers
        while i < n and fmt[i] in "hljztLq":
            # handle doubled letters like ll or hh naturally by loop
            i += 1

        if i >= n:
            return None

        conv = fmt[i]
        raw = fmt[start:i + 1]
        if conv.isalpha():
            spec_index += 1
            specs.append(
                SpecInfo(
                    spec_index=spec_index,
                    start=start,
                    end=i + 1,
                    conv=conv,
                    raw=raw,
                    flags_width_prec_len=fmt[start + 1:i],
                )
            )
            i += 1
            continue

        return None

    return specs


def replacement_for_spec(spec: SpecInfo, macro_name: str) -> Optional[str]:
    if spec.conv not in PRINTF_SPECS:
        return None
    # Only replace integer-like conversions. Caller decides based on warning type.
    return macro_name


def build_rewritten_format(fmt: str, target_spec_indices: Dict[int, str]) -> Optional[str]:
    specs = parse_format_specs(fmt)
    if specs is None:
        return None

    pieces: List[str] = []
    cursor = 0

    for spec in specs:
        replacement = target_spec_indices.get(spec.spec_index)
        if replacement is None:
            continue
        pieces.append(encode_c_string_literal(fmt[cursor:spec.start]))
        pieces.append(replacement)
        cursor = spec.end

    pieces.append(encode_c_string_literal(fmt[cursor:]))

    # Remove empty string literals where possible.
    cleaned = [p for p in pieces if p != '""']

    if not cleaned:
        return '""'

    # Adjacent string literals/macros are valid C concatenation.
    return " ".join(cleaned)


def argnum_to_spec_index(call: CallSite, warning_argnum: int) -> Optional[int]:
    """
    Map compiler's 1-based argument index in the whole call to
    1-based format-spec index among actual conversion specs.
    """
    local_index = warning_argnum - 1
    value_arg_index = local_index - call.first_value_arg_index
    if value_arg_index < 0:
        return None
    # No support for '*' width/precision consuming arguments, since parser skips those.
    return value_arg_index + 1


def choose_sites(sites: List[WarningSite]) -> Dict[Tuple[str, int], List[WarningSite]]:
    grouped: Dict[Tuple[str, int], List[WarningSite]] = {}
    for s in sites:
        grouped.setdefault((s.file, s.line), []).append(s)
    return grouped


def apply_edits(
    source_text: str,
    file_path: str,
    line_groups: Dict[int, List[WarningSite]],
    macro_name: str,
    verbose: bool = False,
) -> Tuple[str, List[str]]:
    """
    Apply edits from bottom to top to preserve offsets.
    """
    messages: List[str] = []
    candidates: List[Tuple[int, int, str]] = []

    for line_no, warnings in sorted(line_groups.items(), reverse=True):
        call = find_call_around(source_text, line_no)
        if call is None:
            messages.append(f"skip {file_path}:{line_no}: could not locate printf-family call")
            continue

        fmt_expr = call.args[call.format_arg_index]
        extracted = extract_concatenated_string_literal(fmt_expr)
        if extracted is None:
            messages.append(f"skip {file_path}:{line_no}: format arg is not a plain string literal")
            continue

        fmt_content, _pieces = extracted
        specs = parse_format_specs(fmt_content)
        if specs is None:
            messages.append(f"skip {file_path}:{line_no}: unsupported format string features")
            continue

        spec_replacements: Dict[int, str] = {}

        for w in warnings:
            if w.argnum is None:
                messages.append(f"skip {file_path}:{line_no}: warning lacks argument number")
                continue
            spec_index = argnum_to_spec_index(call, w.argnum)
            if spec_index is None:
                messages.append(f"skip {file_path}:{line_no}: could not map argument {w.argnum} to format spec")
                continue
            spec = next((s for s in specs if s.spec_index == spec_index), None)
            if spec is None:
                messages.append(f"skip {file_path}:{line_no}: no matching format spec for argument {w.argnum}")
                continue
            repl = replacement_for_spec(spec, macro_name)
            if repl is None:
                messages.append(
                    f"skip {file_path}:{line_no}: spec {spec.raw!r} is not an integer conversion we rewrite"
                )
                continue
            spec_replacements[spec_index] = repl

        if not spec_replacements:
            continue

        new_fmt_expr = build_rewritten_format(fmt_content, spec_replacements)
        if new_fmt_expr is None:
            messages.append(f"skip {file_path}:{line_no}: failed to rebuild format string")
            continue

        new_args = list(call.args)
        new_args[call.format_arg_index] = new_fmt_expr

        open_paren = call.text.find("(")
        close_paren = call.text.rfind(")")
        if open_paren == -1 or close_paren == -1 or close_paren < open_paren:
            messages.append(f"skip {file_path}:{line_no}: malformed call text")
            continue

        new_inner = ", ".join(new_args)
        new_call = call.text[:open_paren + 1] + new_inner + call.text[close_paren:]

        candidates.append((call.start_idx, call.end_idx + 1, new_call))
        messages.append(f"fix  {file_path}:{line_no}: rewrote {call.func} format string")

    # Apply replacements back-to-front.
    new_text = source_text
    for start, end, repl in sorted(candidates, key=lambda t: t[0], reverse=True):
        new_text = new_text[:start] + repl + new_text[end:]

    return new_text, messages


def build_unified_diff(old: str, new: str, path: str) -> str:
    return "".join(
        difflib.unified_diff(
            old.splitlines(keepends=True),
            new.splitlines(keepends=True),
            fromfile=path,
            tofile=path,
        )
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("warning_log", nargs="?", help="Compiler warning log file. If omitted, read stdin.")
    ap.add_argument("--apply", action="store_true", help="Write changes in place instead of printing diff.")
    ap.add_argument("--macro", default="SLINT_FMT", help="Format macro to splice in. Default: SLINT_FMT")
    ap.add_argument("--verbose", action="store_true", help="Print skip/fix messages to stderr")
    args = ap.parse_args()

    if args.warning_log:
        warn_text = read_text(args.warning_log)
    else:
        warn_text = sys.stdin.read()

    warning_sites = parse_warning_lines(warn_text.splitlines())
    if not warning_sites:
        print("No matching format warnings found.", file=sys.stderr)
        return 1

    # Group by file then line.
    file_map: Dict[str, Dict[int, List[WarningSite]]] = {}
    for site in warning_sites:
        file_map.setdefault(site.file, {}).setdefault(site.line, []).append(site)

    any_changes = False
    all_diffs: List[str] = []

    for path, line_groups in sorted(file_map.items()):
        if not os.path.exists(path):
            if args.verbose:
                print(f"skip {path}: file does not exist", file=sys.stderr)
            continue

        old = read_text(path)
        new, messages = apply_edits(old, path, line_groups, args.macro, verbose=args.verbose)

        if args.verbose:
            for msg in messages:
                print(msg, file=sys.stderr)

        if new != old:
            any_changes = True
            if args.apply:
                write_text(path, new)
            else:
                all_diffs.append(build_unified_diff(old, new, path))

    if not any_changes:
        print("No safe rewrites were found.", file=sys.stderr)
        return 2

    if not args.apply:
        sys.stdout.write("".join(all_diffs))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
