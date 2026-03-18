#!/usr/bin/env python3
"""Compare sibling Fortran variants after precision-aware normalization.

This tool is intended for families like LATRD/LARFG where:
- S vs D should be almost entirely mechanical
- C vs Z should be almost entirely mechanical
- real vs complex differences should remain visible

It groups fixed-form Fortran into logical statements, normalizes precision-only
differences, and reports where variants still diverge semantically.
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DIR = REPO_ROOT / "SRC"

VARIANT_ORDER = ["ps", "pd", "pc", "pz"]
FIXED_FORM_SUFFIXES = {".f", ".F"}


@dataclass
class Statement:
    start_line: int
    end_line: int
    text: str
    normalized: str


NORMALIZE_PATTERNS = [
    (re.compile(r"\bDOUBLE\s+PRECISION\b"), "REAL_T"),
    (re.compile(r"\bREAL\b"), "REAL_T"),
    (re.compile(r"\bDOUBLE\s+COMPLEX\b"), "COMPLEX_T"),
    (re.compile(r"\bCOMPLEX\*16\b"), "COMPLEX_T"),
    (re.compile(r"\bCOMPLEX\b"), "COMPLEX_T"),
    (re.compile(r"\bINTEGER\*8\b"), "INTEGER_T"),
    (re.compile(r"\bINTEGER\b"), "INTEGER_T"),
    (re.compile(r"\bDBLE\s*\("), "TO_REAL_T("),
    (re.compile(r"\bREAL\s*\("), "TO_REAL_T("),
    (re.compile(r"\bDCMPLX\s*\("), "TO_COMPLEX_T("),
    (re.compile(r"\bCMPLX\s*\("), "TO_COMPLEX_T("),
    (re.compile(r"\bP[SDCZ]([A-Z0-9_]{3,})\b"), r"P@\1"),
    (re.compile(r"\b[SDCZ]LARFG\b"), "TLARFG"),
    (re.compile(r"\bP[SD]SYMV(_I8)?\b"), r"P@SYHEMV\1"),
    (re.compile(r"\bP[CZ]HEMV(_I8)?\b"), r"P@SYHEMV\1"),
    (re.compile(r"\b[SD]SYMV\b"), "TSYHEMV"),
    (re.compile(r"\b[CZ]HEMV\b"), "TSYHEMV"),
    (re.compile(r"\bP[SD]SYR2K(_I8)?\b"), r"P@SYHER2K\1"),
    (re.compile(r"\bP[CZ]HER2K(_I8)?\b"), r"P@SYHER2K\1"),
    (re.compile(r"\bP[SD]SYR2(_I8)?\b"), r"P@SYHER2\1"),
    (re.compile(r"\bP[CZ]HER2(_I8)?\b"), r"P@SYHER2\1"),
    (re.compile(r"\b[SD]SYR2K\b"), "TSYHER2K"),
    (re.compile(r"\b[CZ]HER2K\b"), "TSYHER2K"),
    (re.compile(r"\b[SD]SYR2\b"), "TSYHER2"),
    (re.compile(r"\b[CZ]HER2\b"), "TSYHER2"),
    (re.compile(r"\bP[SD]DOT(_I8)?\b"), r"P@DOT\1"),
    (re.compile(r"\bP[CZ]DOTC(_I8)?\b"), r"P@DOT\1"),
    (re.compile(r"\b[SD]DOT\b"), "TDOT"),
    (re.compile(r"\b[CZ]DOTC\b"), "TDOT"),
    (re.compile(r"\bP[SDCZ]AXPY(_I8)?\b"), r"P@AXPY\1"),
    (re.compile(r"\b[SDCZ]AXPY\b"), "TAXPY"),
    (re.compile(r"\bP[SDCZ]SCAL(_I8)?\b"), r"P@SCAL\1"),
    (re.compile(r"\b[SDCZ]SCAL\b"), "TSCAL"),
    (re.compile(r"\b[SDCZ]GEMV\b"), "TGEMV"),
    (re.compile(r"\bP[SDCZ]GEMV(_I8)?\b"), r"P@GEMV\1"),
    (re.compile(r"\bP[SD]SYNTRD(_I8)?\b"), r"P@REDTRD\1"),
    (re.compile(r"\bP[CZ]HENTRD(_I8)?\b"), r"P@REDTRD\1"),
    (re.compile(r"\bP[SD]SYTD2(_I8)?\b"), r"P@TD2\1"),
    (re.compile(r"\bP[CZ]HETD2(_I8)?\b"), r"P@TD2\1"),
    (re.compile(r"\bP[SD]SYTTRD(_I8)?\b"), r"P@TTRD\1"),
    (re.compile(r"\bP[CZ]HETTRD(_I8)?\b"), r"P@TTRD\1"),
    (re.compile(r"\b[SD]SYTRD\b"), "TREDTRD"),
    (re.compile(r"\b[CZ]HETRD\b"), "TREDTRD"),
    (re.compile(r"\b[0-9]+(?:\.[0-9]*)?[DE][+-]?[0-9]+\b"), "FLOAT_LIT"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare sibling Fortran variants after precision-aware normalization."
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="Explicit files to compare. If omitted, use --stem discovery.",
    )
    parser.add_argument(
        "--stem",
        help="Family stem to discover, e.g. latrd or latrd_i8.",
    )
    parser.add_argument(
        "--dir",
        default=str(DEFAULT_DIR),
        help="Directory to search when using --stem. Defaults to SRC.",
    )
    parser.add_argument(
        "--baseline",
        choices=VARIANT_ORDER,
        default="pd",
        help="Baseline variant prefix. Defaults to pd.",
    )
    parser.add_argument(
        "--max-hunks",
        type=int,
        default=8,
        help="Maximum diff hunks to show per comparison.",
    )
    parser.add_argument(
        "--context",
        type=int,
        default=2,
        help="Context statements to show around each hunk.",
    )
    return parser.parse_args()


def resolve_path(raw: str) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path


def discover_files(stem: str, directory: Path) -> list[Path]:
    matches: list[Path] = []
    for prefix in VARIANT_ORDER:
        candidate = directory / f"{prefix}{stem}.f"
        if candidate.exists():
            matches.append(candidate)
            continue
        candidate_upper = directory / f"{prefix}{stem}.F"
        if candidate_upper.exists():
            matches.append(candidate_upper)
    return matches


def is_comment(line: str) -> bool:
    if not line:
        return False
    first = line[0]
    return first in {"c", "C", "*", "!"}


def is_continuation(line: str) -> bool:
    if len(line) >= 6 and line[5] not in {" ", "0"}:
        return True
    return line.lstrip().startswith("$") or line.rstrip().endswith("&")


def fixed_form_payload(line: str, continuation: bool) -> str:
    if len(line) >= 6:
        payload = line[6:] if continuation else line[6:]
    else:
        payload = line
    if continuation:
        payload = payload.lstrip("$& ").rstrip()
    return payload.rstrip()


def normalize_statement(text: str) -> str:
    normalized = text.upper()
    for pattern, replacement in NORMALIZE_PATTERNS:
        normalized = pattern.sub(replacement, normalized)
    normalized = re.sub(r"\s+", " ", normalized).strip()
    return normalized


def parse_statements(path: Path) -> list[Statement]:
    raw_lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    statements: list[Statement] = []
    buffer: list[str] = []
    start_line = 0
    end_line = 0

    def flush() -> None:
        nonlocal buffer, start_line, end_line
        if not buffer:
            return
        text = " ".join(part for part in buffer if part).strip()
        if text:
            statements.append(
                Statement(
                    start_line=start_line,
                    end_line=end_line,
                    text=text,
                    normalized=normalize_statement(text),
                )
            )
        buffer = []
        start_line = 0
        end_line = 0

    for idx, raw in enumerate(raw_lines, start=1):
        if is_comment(raw):
            flush()
            continue
        if not raw.strip():
            flush()
            continue

        continuation = is_continuation(raw)
        payload = fixed_form_payload(raw, continuation).strip()
        if not payload:
            continue

        if continuation and buffer:
            buffer.append(payload)
            end_line = idx
            continue

        flush()
        buffer = [payload]
        start_line = idx
        end_line = idx

    flush()
    return statements


def variant_from_name(path: Path) -> str:
    name = path.stem.lower()
    for prefix in VARIANT_ORDER:
        if name.startswith(prefix):
            return prefix
    return name[:2]


def similarity(a: list[str], b: list[str]) -> float:
    return difflib.SequenceMatcher(None, a, b).ratio()


def render_hunk(
    base_name: str,
    base_statements: list[Statement],
    other_name: str,
    other_statements: list[Statement],
    opcode: tuple[str, int, int, int, int],
    context: int,
) -> list[str]:
    tag, i1, i2, j1, j2 = opcode
    lines: list[str] = []
    base_start = max(0, i1 - context)
    base_end = min(len(base_statements), i2 + context)
    other_start = max(0, j1 - context)
    other_end = min(len(other_statements), j2 + context)

    lines.append(f"  hunk [{tag}]")
    lines.append(f"    {base_name}:")
    for idx, stmt in enumerate(base_statements[base_start:base_end], start=base_start):
        marker = ">" if i1 <= idx < i2 else " "
        lines.append(
            f"      {marker} L{stmt.start_line}-{stmt.end_line}: {stmt.text}"
        )
    lines.append(f"    {other_name}:")
    for idx, stmt in enumerate(other_statements[other_start:other_end], start=other_start):
        marker = ">" if j1 <= idx < j2 else " "
        lines.append(
            f"      {marker} L{stmt.start_line}-{stmt.end_line}: {stmt.text}"
        )
    return lines


def render_comparison(
    baseline_path: Path,
    baseline: list[Statement],
    other_path: Path,
    other: list[Statement],
    max_hunks: int,
    context: int,
) -> str:
    base_norm = [stmt.normalized for stmt in baseline]
    other_norm = [stmt.normalized for stmt in other]
    matcher = difflib.SequenceMatcher(None, base_norm, other_norm)
    opcodes = [op for op in matcher.get_opcodes() if op[0] != "equal"]
    lines: list[str] = []
    lines.append(
        f"{baseline_path.name} vs {other_path.name}: "
        f"normalized_similarity={similarity(base_norm, other_norm):.3f} "
        f"statement_count={len(base_norm)}/{len(other_norm)} "
        f"diff_hunks={len(opcodes)}"
    )
    for opcode in opcodes[:max_hunks]:
        lines.extend(
            render_hunk(
                baseline_path.name,
                baseline,
                other_path.name,
                other,
                opcode,
                context,
            )
        )
    if len(opcodes) > max_hunks:
        lines.append(f"  ... {len(opcodes) - max_hunks} more hunk(s) omitted")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()

    if not args.files and not args.stem:
        print("Provide explicit files or use --stem.", file=sys.stderr)
        return 2

    if args.files:
        files = [resolve_path(raw) for raw in args.files]
    else:
        directory = resolve_path(args.dir)
        files = discover_files(args.stem, directory)

    files = [path for path in files if path.exists() and path.suffix in FIXED_FORM_SUFFIXES]
    if len(files) < 2:
        print("Need at least two Fortran variant files to compare.", file=sys.stderr)
        return 2

    deduped: dict[str, Path] = {}
    for path in files:
        key = variant_from_name(path)
        existing = deduped.get(key)
        if existing is None or existing.suffix == ".F":
            deduped[key] = path

    files = sorted(
        deduped.values(),
        key=lambda path: VARIANT_ORDER.index(variant_from_name(path))
        if variant_from_name(path) in VARIANT_ORDER
        else 99,
    )
    by_variant = {variant_from_name(path): path for path in files}
    baseline_path = by_variant.get(args.baseline, files[0])

    parsed = {path: parse_statements(path) for path in files}

    lines: list[str] = []
    lines.append("Fortran variant diff")
    lines.append(f"baseline: {baseline_path.relative_to(REPO_ROOT)}")
    lines.append("variants:")
    for path in files:
        lines.append(f"  {variant_from_name(path)} -> {path.relative_to(REPO_ROOT)}")
    lines.append("")

    for path in files:
        if path == baseline_path:
            continue
        lines.append(
            render_comparison(
                baseline_path,
                parsed[baseline_path],
                path,
                parsed[path],
                args.max_hunks,
                args.context,
            )
        )
        lines.append("")

    print("\n".join(lines).rstrip() + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
