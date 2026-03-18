#!/usr/bin/env python3
"""Inventory PBLAS routine signatures and direct dependencies.

This is a lightweight static inventory, not a full parser. Its goal is to make
Phase 5 planning more data-driven by extracting:

- routine entry points found under PBLAS/SRC
- public-family grouping for the top-level PBLAS wrappers
- descriptor and integer-like arguments
- direct call sites from C and Fortran routines
- in-tree call resolution against other routines in PBLAS/SRC
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ROOT = REPO_ROOT / "PBLAS" / "SRC"
SOURCE_SUFFIXES = {".c", ".f", ".F", ".f90", ".F90"}

C_CALL_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\(")
C_FUNC_RE = re.compile(r"\bvoid\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(")
FORTRAN_ENTRY_RE = re.compile(
    r"^\s*(SUBROUTINE|[A-Z0-9* ]+\bFUNCTION)\s+([A-Z0-9_]+)\s*\(",
    re.IGNORECASE,
)
FORTRAN_CALL_RE = re.compile(r"\bCALL\s+([A-Z0-9_]+)\s*\(", re.IGNORECASE)
FORTRAN_EXTERNAL_RE = re.compile(r"^\s*EXTERNAL\s+(.+)$", re.IGNORECASE)

C_KEYWORDS = {
    "if",
    "for",
    "while",
    "switch",
    "return",
    "sizeof",
    "free",
    "malloc",
    "calloc",
    "realloc",
}

BLACS_PREFIXES = ("Cblacs_", "BI_", "Cd", "Cs", "Cz", "Cc", "igam", "sgam", "dgam")
PB_PREFIXES = ("PB_",)


@dataclass
class Routine:
    name: str
    path: str
    language: str
    area: str
    family: str
    variant: str
    arg_count: int
    args: list[str]
    descriptor_args: list[str]
    int_like_args: list[str]
    externals: list[str]
    direct_calls: list[str]
    resolved_calls: list[str]
    call_categories: dict[str, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inventory PBLAS/SRC routine signatures and direct dependencies."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=[str(DEFAULT_ROOT)],
        help="Files or directories to scan. Defaults to PBLAS/SRC.",
    )
    parser.add_argument(
        "--root",
        default=str(REPO_ROOT),
        help="Repository root. Defaults to the current script's repo.",
    )
    parser.add_argument(
        "--focus",
        help="Regex used to filter detailed output by routine name, family, or path.",
    )
    parser.add_argument(
        "--area",
        choices=("all", "public", "pbblas", "ptzblas"),
        default="all",
        help="Restrict detailed output to one source area.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format.",
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="Show detailed per-routine output for the filtered set.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top families to show in summary output.",
    )
    return parser.parse_args()


def resolve_targets(root: Path, raw_paths: Iterable[str]) -> list[Path]:
    targets: list[Path] = []
    for raw in raw_paths:
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = root / candidate
        if candidate.exists():
            targets.append(candidate)
    return targets


def iter_source_files(targets: Iterable[Path]) -> Iterable[Path]:
    for target in targets:
        if target.is_file():
            if target.suffix in SOURCE_SUFFIXES:
                yield target
            continue
        for path in sorted(target.rglob("*")):
            if path.is_file() and path.suffix in SOURCE_SUFFIXES:
                yield path


def strip_c_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*", "", text)
    return text


def classify_area(path: Path) -> str:
    parts = path.parts
    if "PTZBLAS" in parts:
        return "ptzblas"
    if "PBBLAS" in parts:
        return "pbblas"
    return "public"


def classify_public_family(name: str) -> tuple[str, str]:
    stem = name.rstrip("_").lower()
    if len(stem) >= 3 and stem[0] == "p" and stem[1] in "sdczi":
        return stem[2:], stem[1]
    return stem, ""


def classify_call(name: str) -> str:
    if name.startswith(PB_PREFIXES):
        return "pb_internal"
    if name.startswith(BLACS_PREFIXES) or "blacs" in name.lower():
        return "blacs_mpi"
    if name.endswith("_"):
        if name.lower().startswith("p"):
            return "pblas_or_scalapack"
        return "dense_blas_lapack"
    return "other"


def unique_sorted(items: Iterable[str]) -> list[str]:
    return sorted({item for item in items if item})


def parse_c_args(arg_text: str) -> tuple[list[str], list[str], list[str]]:
    raw_args = [part.strip() for part in arg_text.replace("\n", " ").split(",")]
    args: list[str] = []
    descriptor_args: list[str] = []
    int_like_args: list[str] = []

    for raw in raw_args:
        if not raw or raw == "void":
            continue
        tokens = raw.replace("*", " ").split()
        if not tokens:
            continue
        name = tokens[-1].strip("[]")
        args.append(name)
        if name.upper().startswith("DESC"):
            descriptor_args.append(name)
        if "Int" in raw or "F_CHAR_T" in raw:
            int_like_args.append(name)
    return args, descriptor_args, int_like_args


def extract_c_routine(path: Path, text: str) -> Routine | None:
    uncommented = strip_c_comments(text)
    match = C_FUNC_RE.search(uncommented)
    if not match:
        return None

    name = match.group(1)
    start = match.end()
    depth = 1
    idx = start
    while idx < len(uncommented) and depth > 0:
        char = uncommented[idx]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        idx += 1
    arg_text = uncommented[start : idx - 1]
    args, descriptor_args, int_like_args = parse_c_args(arg_text)

    body_start = uncommented.find("{", idx)
    body = uncommented[body_start:] if body_start != -1 else ""
    calls = [
        call
        for call in C_CALL_RE.findall(body)
        if call not in C_KEYWORDS and call != name
    ]

    family, variant = classify_public_family(name)
    categories = Counter(classify_call(call) for call in calls)
    return Routine(
        name=name,
        path=str(path.relative_to(REPO_ROOT)),
        language="c",
        area=classify_area(path),
        family=family,
        variant=variant,
        arg_count=len(args),
        args=args,
        descriptor_args=unique_sorted(descriptor_args),
        int_like_args=unique_sorted(int_like_args),
        externals=[],
        direct_calls=unique_sorted(calls),
        resolved_calls=[],
        call_categories=dict(sorted(categories.items())),
    )


def preprocess_fortran_lines(text: str) -> list[str]:
    lines: list[str] = []
    for raw in text.splitlines():
        if not raw:
            lines.append("")
            continue
        first = raw[0]
        if first in {"c", "C", "*", "!"}:
            lines.append("")
            continue
        lines.append(raw.rstrip())
    return lines


def parse_fortran_args(lines: list[str], start_index: int) -> list[str]:
    buffer = lines[start_index].strip()
    idx = start_index
    while ")" not in buffer and idx + 1 < len(lines):
        idx += 1
        cont = lines[idx].replace("$", " ").strip()
        buffer += " " + cont
    open_paren = buffer.find("(")
    close_paren = buffer.rfind(")")
    if open_paren == -1 or close_paren == -1 or close_paren <= open_paren:
        return []
    arg_text = buffer[open_paren + 1 : close_paren]
    return [part.strip() for part in arg_text.split(",") if part.strip()]


def extract_fortran_externals(lines: list[str]) -> list[str]:
    externals: list[str] = []
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        match = FORTRAN_EXTERNAL_RE.match(line)
        if not match:
            idx += 1
            continue
        chunk = match.group(1).strip()
        idx += 1
        while idx < len(lines):
            continuation = lines[idx]
            if not continuation.lstrip().startswith("$"):
                break
            chunk += " " + continuation.replace("$", " ").strip()
            idx += 1
        for name in chunk.replace(",", " ").split():
            upper = name.strip().upper()
            if upper:
                externals.append(upper)
    return unique_sorted(externals)


def extract_fortran_routine(path: Path, text: str) -> Routine | None:
    lines = preprocess_fortran_lines(text)
    entry_match = None
    entry_index = -1
    for idx, line in enumerate(lines):
        match = FORTRAN_ENTRY_RE.match(line)
        if match:
            entry_match = match
            entry_index = idx
            break
    if entry_match is None:
        return None

    name = entry_match.group(2).upper()
    args = parse_fortran_args(lines, entry_index)
    descriptor_args = [arg for arg in args if arg.upper().startswith("DESC")]
    int_like_args = [
        arg for arg in args if re.match(r"^(I|J|K|M|N|L|INC|DESC|RSRC|CSRC)", arg.upper())
    ]
    calls = [match.group(1).upper() for line in lines for match in FORTRAN_CALL_RE.finditer(line)]
    externals = extract_fortran_externals(lines)
    family, variant = classify_public_family(name)
    categories = Counter(classify_call(call) for call in calls)

    return Routine(
        name=name,
        path=str(path.relative_to(REPO_ROOT)),
        language="fortran",
        area=classify_area(path),
        family=family,
        variant=variant,
        arg_count=len(args),
        args=args,
        descriptor_args=unique_sorted(descriptor_args),
        int_like_args=unique_sorted(int_like_args),
        externals=externals,
        direct_calls=unique_sorted(calls),
        resolved_calls=[],
        call_categories=dict(sorted(categories.items())),
    )


def extract_routine(path: Path) -> Routine | None:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None

    if path.suffix == ".c":
        return extract_c_routine(path, text)
    return extract_fortran_routine(path, text)


def resolve_calls(routines: list[Routine]) -> None:
    names = {routine.name.lower() for routine in routines}
    for routine in routines:
        routine.resolved_calls = unique_sorted(
            call for call in routine.direct_calls if call.lower() in names
        )


def public_family_summary(routines: list[Routine]) -> list[dict[str, object]]:
    grouped: dict[str, list[Routine]] = defaultdict(list)
    for routine in routines:
        if routine.area == "public":
            grouped[routine.family].append(routine)

    summary: list[dict[str, object]] = []
    for family, members in grouped.items():
        variants = sorted(r.variant for r in members if r.variant)
        summary.append(
            {
                "family": family,
                "variant_count": len(set(variants)),
                "variants": variants,
                "member_count": len(members),
                "files": sorted(r.path for r in members),
            }
        )
    summary.sort(key=lambda item: (-int(item["variant_count"]), item["family"]))
    return summary


def filter_routines(routines: list[Routine], pattern: str | None, area: str) -> list[Routine]:
    filtered = routines
    if area != "all":
        filtered = [routine for routine in filtered if routine.area == area]
    if pattern is None:
        return filtered
    regex = re.compile(pattern, re.IGNORECASE)
    return [
        routine
        for routine in filtered
        if regex.search(routine.name)
        or regex.search(routine.family)
        or regex.search(routine.path)
        or any(regex.search(call) for call in routine.direct_calls)
    ]


def render_text(routines: list[Routine], filtered: list[Routine], top_n: int, details: bool) -> str:
    lines: list[str] = []
    lines.append(
        f"Scanned {len(routines)} routine(s) under {DEFAULT_ROOT.relative_to(REPO_ROOT)}."
    )
    by_area = Counter(r.area for r in routines)
    lines.append(
        "Area counts: "
        + ", ".join(f"{area}={by_area.get(area, 0)}" for area in ("public", "pbblas", "ptzblas"))
    )
    by_lang = Counter(r.language for r in routines)
    lines.append(
        "Language counts: "
        + ", ".join(f"{lang}={by_lang.get(lang, 0)}" for lang in ("c", "fortran"))
    )
    lines.append("")

    lines.append("Top public families:")
    for item in public_family_summary(routines)[:top_n]:
        variants = ",".join(item["variants"]) if item["variants"] else "-"
        lines.append(
            f"  {item['family']}: variants={variants} members={item['member_count']}"
        )
    lines.append("")

    if not details:
        lines.append(
            f"Filtered routines available for detail: {len(filtered)}. Use --details or --format json."
        )
        return "\n".join(lines) + "\n"

    lines.append(f"Detailed routines: {len(filtered)}")
    for routine in sorted(filtered, key=lambda r: (r.area, r.family, r.name)):
        lines.append(f"- {routine.name} [{routine.language}, {routine.area}]")
        lines.append(f"  path: {routine.path}")
        lines.append(f"  family: {routine.family}  variant: {routine.variant or '-'}")
        lines.append(f"  arg_count: {routine.arg_count}")
        lines.append(
            f"  descriptor_args: {', '.join(routine.descriptor_args) or '-'}"
        )
        lines.append(f"  int_like_args: {', '.join(routine.int_like_args) or '-'}")
        if routine.externals:
            lines.append(f"  externals: {', '.join(routine.externals)}")
        lines.append(f"  direct_calls({len(routine.direct_calls)}): {', '.join(routine.direct_calls) or '-'}")
        lines.append(
            f"  resolved_calls({len(routine.resolved_calls)}): {', '.join(routine.resolved_calls) or '-'}"
        )
        if routine.call_categories:
            cat_text = ", ".join(
                f"{name}={count}" for name, count in sorted(routine.call_categories.items())
            )
            lines.append(f"  call_categories: {cat_text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    targets = resolve_targets(root, args.paths)
    routines = [routine for path in iter_source_files(targets) if (routine := extract_routine(path))]
    resolve_calls(routines)
    filtered = filter_routines(routines, args.focus, args.area)

    if args.format == "json":
        payload = {
            "summary": {
                "routine_count": len(routines),
                "filtered_count": len(filtered),
                "by_area": Counter(r.area for r in routines),
                "by_language": Counter(r.language for r in routines),
                "public_families": public_family_summary(routines),
            },
            "routines": [asdict(routine) for routine in filtered],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(render_text(routines, filtered, args.top, args.details), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
