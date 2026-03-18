#!/usr/bin/env python3
"""Report legacy/I8 variant completeness for a ScaLAPACK family stem.

This script is intended as a lightweight phase tracker for families like:
  latrd, larfg, syntrd, hentrd

For each of the standard p{s,d,c,z} variants it reports:
- legacy file presence
- _i8 file presence
- tests that call the legacy or _i8 routine
- whether those tests appear in a CMakeLists.txt
- _i8 callers that still reference the legacy routine name
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CODE_DIRS = ("SRC", "TOOLS", "PBLAS/SRC", "REDIST/SRC")
DEFAULT_TEST_DIRS = ("SRC", "TOOLS", "REDIST/TESTING", "PBLAS/TESTING")
VARIANTS = ("ps", "pd", "pc", "pz")
SOURCE_SUFFIXES = {".f", ".F", ".f90", ".F90", ".c", ".h"}
TEST_PREFIX = "test_"


@dataclass
class VariantRow:
    variant: str
    legacy_name: str
    i8_name: str
    legacy_file: str | None
    i8_file: str | None
    legacy_tests: list[str] = field(default_factory=list)
    i8_tests: list[str] = field(default_factory=list)
    legacy_tests_wired: list[str] = field(default_factory=list)
    i8_tests_wired: list[str] = field(default_factory=list)
    i8_callers_using_legacy: list[str] = field(default_factory=list)
    i8_callers_using_i8: list[str] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report variant completeness and migration status for a family stem."
    )
    parser.add_argument("stem", help="Family stem, e.g. latrd, larfg, syntrd, hentrd.")
    parser.add_argument(
        "--root",
        default=str(REPO_ROOT),
        help="Repository root. Defaults to the current script's repo.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format.",
    )
    parser.add_argument(
        "--code-dirs",
        nargs="*",
        default=list(DEFAULT_CODE_DIRS),
        help="Directories to scan for source files and callers.",
    )
    parser.add_argument(
        "--test-dirs",
        nargs="*",
        default=list(DEFAULT_TEST_DIRS),
        help="Directories to scan for tests.",
    )
    return parser.parse_args()


def resolve_targets(root: Path, raw_paths: Iterable[str]) -> list[Path]:
    targets: list[Path] = []
    for raw in raw_paths:
        path = Path(raw)
        if not path.is_absolute():
            path = root / path
        if path.exists():
            targets.append(path)
    return targets


def iter_files(targets: Iterable[Path]) -> Iterable[Path]:
    for target in targets:
        if target.is_file():
            if target.suffix in SOURCE_SUFFIXES:
                yield target
            continue
        for path in sorted(target.rglob("*")):
            if path.is_file() and path.suffix in SOURCE_SUFFIXES:
                yield path


def slurp(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def discover_family_files(code_files: list[Path], stem: str) -> dict[str, dict[str, Path]]:
    mapping: dict[str, dict[str, Path]] = {variant: {} for variant in VARIANTS}
    pattern = re.compile(rf"^(p[sdcz]){re.escape(stem)}(_i8)?$", re.IGNORECASE)

    for path in code_files:
        match = pattern.match(path.stem)
        if not match:
            continue
        variant = match.group(1).lower()
        flavor = "i8" if match.group(2) else "legacy"
        existing = mapping[variant].get(flavor)
        if existing is None or existing.suffix == ".F":
            mapping[variant][flavor] = path
    return mapping


def build_call_regex(symbol: str) -> re.Pattern[str]:
    return re.compile(rf"\b(?:CALL\s+)?{re.escape(symbol)}\s*\(", re.IGNORECASE)


def build_definition_regex(symbol: str) -> re.Pattern[str]:
    return re.compile(
        rf"^\s*(SUBROUTINE|[A-Z0-9* ]+\bFUNCTION|void)\s+{re.escape(symbol)}\s*\(",
        re.IGNORECASE | re.MULTILINE,
    )


def find_callers(
    files: list[Path],
    symbol: str,
    *,
    only_i8_files: bool = False,
    exclude_tests: bool = False,
) -> list[str]:
    regex = build_call_regex(symbol)
    def_regex = build_definition_regex(symbol)
    callers: list[str] = []
    for path in files:
        if exclude_tests and path.name.lower().startswith(TEST_PREFIX):
            continue
        if only_i8_files and "_i8" not in path.stem.lower():
            continue
        if path.stem.lower() == symbol.lower():
            continue
        try:
            text = slurp(path)
        except OSError:
            continue
        if def_regex.search(text):
            # Skip definition-only self matches in files whose stem does not
            # exactly match the symbol name.
            text = def_regex.sub("", text)
        if regex.search(text):
            callers.append(rel(path))
    return sorted(callers)


def find_tests_calling(test_files: list[Path], symbol: str) -> list[Path]:
    regex = build_call_regex(symbol)
    hits: list[Path] = []
    for path in test_files:
        if not path.name.lower().startswith(TEST_PREFIX):
            continue
        try:
            text = slurp(path)
        except OSError:
            continue
        if regex.search(text):
            hits.append(path)
    return sorted(hits)


def wired_tests(test_paths: list[Path], root: Path) -> list[str]:
    cmake_files = sorted(root.rglob("CMakeLists.txt"))
    results: list[str] = []
    for test_path in test_paths:
        basename = test_path.name
        for cmake in cmake_files:
            try:
                text = slurp(cmake)
            except OSError:
                continue
            if basename in text:
                results.append(rel(test_path))
                break
    return sorted(set(results))


def derive_status(rows: list[VariantRow]) -> str:
    all_i8 = all(row.i8_file for row in rows)
    all_tested = all(row.i8_tests for row in rows if row.i8_file)
    any_legacy_in_i8 = any(row.i8_callers_using_legacy for row in rows)

    if all_i8 and all_tested and not any_legacy_in_i8:
        return "complete"
    if any(row.i8_file for row in rows):
        return "partial"
    return "blocked"


def render_text(stem: str, rows: list[VariantRow], status: str) -> str:
    lines: list[str] = []
    lines.append(f"family: {stem}")
    lines.append(f"status: {status}")
    lines.append("")
    lines.append(
        "variant  legacy  i8  legacy_tests  i8_tests  legacy_wired  i8_wired  i8->legacy  i8->i8"
    )
    for row in rows:
        lines.append(
            f"{row.variant:<7} "
            f"{'yes' if row.legacy_file else 'no ':<6} "
            f"{'yes' if row.i8_file else 'no ':<3} "
            f"{len(row.legacy_tests):<12} "
            f"{len(row.i8_tests):<8} "
            f"{len(row.legacy_tests_wired):<12} "
            f"{len(row.i8_tests_wired):<8} "
            f"{len(row.i8_callers_using_legacy):<10} "
            f"{len(row.i8_callers_using_i8):<6}"
        )
    lines.append("")

    for row in rows:
        lines.append(f"[{row.variant}] {row.legacy_name} / {row.i8_name}")
        lines.append(f"  legacy_file: {row.legacy_file or '-'}")
        lines.append(f"  i8_file: {row.i8_file or '-'}")
        lines.append(f"  legacy_tests: {', '.join(row.legacy_tests) or '-'}")
        lines.append(f"  i8_tests: {', '.join(row.i8_tests) or '-'}")
        lines.append(f"  legacy_tests_wired: {', '.join(row.legacy_tests_wired) or '-'}")
        lines.append(f"  i8_tests_wired: {', '.join(row.i8_tests_wired) or '-'}")
        lines.append(
            f"  i8_callers_using_legacy: {', '.join(row.i8_callers_using_legacy) or '-'}"
        )
        lines.append(
            f"  i8_callers_using_i8: {', '.join(row.i8_callers_using_i8) or '-'}"
        )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    code_dirs = resolve_targets(root, args.code_dirs)
    test_dirs = resolve_targets(root, args.test_dirs)
    code_files = list(iter_files(code_dirs))
    test_files = list(iter_files(test_dirs))

    family_files = discover_family_files(code_files, args.stem.lower())
    rows: list[VariantRow] = []

    for variant in VARIANTS:
        legacy_name = f"{variant}{args.stem}".lower()
        i8_name = f"{legacy_name}_i8"
        legacy_file = family_files[variant].get("legacy")
        i8_file = family_files[variant].get("i8")

        legacy_tests = [rel(path) for path in find_tests_calling(test_files, legacy_name)]
        i8_tests = [rel(path) for path in find_tests_calling(test_files, i8_name)]

        rows.append(
            VariantRow(
                variant=variant,
                legacy_name=legacy_name,
                i8_name=i8_name,
                legacy_file=rel(legacy_file) if legacy_file else None,
                i8_file=rel(i8_file) if i8_file else None,
                legacy_tests=legacy_tests,
                i8_tests=i8_tests,
                legacy_tests_wired=wired_tests([root / path for path in legacy_tests], root),
                i8_tests_wired=wired_tests([root / path for path in i8_tests], root),
                i8_callers_using_legacy=find_callers(
                    code_files, legacy_name, only_i8_files=True, exclude_tests=True
                ),
                i8_callers_using_i8=find_callers(
                    code_files, i8_name, only_i8_files=True, exclude_tests=True
                ),
            )
        )

    status = derive_status(rows)

    if args.format == "json":
        payload = {
            "family": args.stem.lower(),
            "status": status,
            "rows": [asdict(row) for row in rows],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(render_text(args.stem.lower(), rows, status), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
