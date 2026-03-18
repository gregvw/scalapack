#!/usr/bin/env python3
"""Lightweight audit for I8 width-boundary hazards.

This intentionally favors high-signal checks over full parsing. It is designed
to catch the kinds of issues that have already appeared during the refactor:

- REAL(...) used to transport workspace/count-like INTEGER*8 values
- hard-coded MPI_INT usage in C/MPI code
- obvious INT(...) narrowing at Fortran call boundaries
- DBLE(...) transport of integer-like values (less severe, still worth review)
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_DIRS = ("BLACS", "PBLAS", "REDIST", "SRC", "TOOLS")
SOURCE_SUFFIXES = {".c", ".h", ".f", ".F", ".f90", ".F90", ".inc"}
SEVERITY_ORDER = {"info": 0, "medium": 1, "high": 2}

REPO_ROOT = Path(__file__).resolve().parents[3]

REAL_OR_DBLE_RE = re.compile(r"\b(REAL|DBLE)\s*\(\s*([A-Z][A-Z0-9_]*)", re.IGNORECASE)
MPI_INT_RE = re.compile(r"\bMPI_INT\b")
INT_CALL_RE = re.compile(r"\bINT\s*\(", re.IGNORECASE)

TRANSPORT_EXACT = {
    "COUNT",
    "INTMAX",
    "LEN",
    "LLRWORK",
    "LLWORK",
    "LRWORK",
    "LRWMIN",
    "LWORK",
    "LWMIN",
    "ONEPMIN",
    "ONEPRMIN",
    "SIZE",
    "TTLRWMIN",
    "TTLWMIN",
}


@dataclass(frozen=True)
class Finding:
    severity: str
    category: str
    path: Path
    line_no: int
    message: str
    line: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit the ScaLAPACK tree for high-signal I8 width-boundary hazards."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=list(DEFAULT_DIRS),
        help="Directories or files to audit, relative to the repo root by default.",
    )
    parser.add_argument(
        "--root",
        default=str(REPO_ROOT),
        help="Repository root. Defaults to the current script's repo.",
    )
    parser.add_argument(
        "--fail-on",
        choices=("high", "medium", "info"),
        help="Exit nonzero if any finding at or above this severity is found.",
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
                if "build-" in path.parts or ".git" in path.parts:
                    continue
                yield path


def is_comment_or_blank(path: Path, line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return True

    if path.suffix.lower() in {".f", ".f90", ".inc"}:
        if line and line[0] in {"c", "C", "*", "!"}:
            return True
        return False

    return stripped.startswith(("//", "/*", "*", "*/"))


def is_transport_integerish(name: str) -> bool:
    upper = name.upper()
    if upper in TRANSPORT_EXACT:
        return True
    return upper.endswith(("LWORK", "RWORK", "WMIN", "COUNT", "SIZE", "LEN"))


def int_calls(line: str) -> Iterable[str]:
    for match in INT_CALL_RE.finditer(line):
        start = match.end()
        end = line.find(")", start)
        if end == -1:
            continue
        yield line[start:end]


def is_width_widening(expr: str) -> bool:
    compact = expr.replace(" ", "")
    return compact.endswith(",8") or ",8," in compact


def expr_looks_integerish(expr: str) -> bool:
    upper = expr.upper()
    if re.search(r"\b[A-Z0-9_]+8\b", upper):
        return True
    return any(
        token in upper
        for token in (
            "DESCA",
            "DESCB",
            "DESCW",
            "LWORK",
            "LRWORK",
            "LLWORK",
            "LLRWORK",
            "TTLWMIN",
            "TTLRWMIN",
            "INTMAX",
            "COUNT",
            "SIZE",
            "NPROW",
            "NPCOL",
            "MYROW",
            "MYCOL",
            "RSRC",
            "CSRC",
            "LLD",
            " IA",
            " JA",
            " IB",
            " JB",
            " NP",
            " NQ",
            " NB",
            " MB",
            " N-",
            " M-",
        )
    )


def audit_file(path: Path) -> list[Finding]:
    findings: list[Finding] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return findings

    for line_no, line in enumerate(lines, start=1):
        if is_comment_or_blank(path, line):
            continue

        for match in REAL_OR_DBLE_RE.finditer(line):
            func = match.group(1).upper()
            name = match.group(2).upper()
            if not is_transport_integerish(name):
                continue
            severity = "high" if func == "REAL" else "medium"
            findings.append(
                Finding(
                    severity=severity,
                    category=f"{func.lower()}-integer-transport",
                    path=path,
                    line_no=line_no,
                    message=f"{func}() applied to integer-like value '{name}'",
                    line=line.rstrip(),
                )
            )

        if MPI_INT_RE.search(line):
            findings.append(
                Finding(
                    severity="high",
                    category="mpi-int",
                    path=path,
                    line_no=line_no,
                    message="Hard-coded MPI_INT; verify API width and ILP64 safety",
                    line=line.rstrip(),
                )
            )

        for expr in int_calls(line):
            if is_width_widening(expr):
                continue
            if not expr_looks_integerish(expr):
                continue
            findings.append(
                Finding(
                    severity="medium",
                    category="int-narrowing",
                    path=path,
                    line_no=line_no,
                    message=f"Possible narrowing INT(...) boundary on '{expr.strip()}'",
                    line=line.rstrip(),
                )
            )

    return findings


def format_report(root: Path, findings: list[Finding]) -> str:
    lines: list[str] = []
    lines.append(
        f"Found {len(findings)} finding(s) across {len({f.path for f in findings})} file(s)."
    )

    counts = {severity: 0 for severity in SEVERITY_ORDER}
    for finding in findings:
        counts[finding.severity] += 1
    lines.append(
        "Severity counts: "
        + ", ".join(f"{name}={counts[name]}" for name in ("high", "medium", "info"))
    )
    lines.append("")

    for finding in sorted(
        findings,
        key=lambda f: (
            -SEVERITY_ORDER[f.severity],
            str(f.path.relative_to(root)),
            f.line_no,
            f.category,
        ),
    ):
        rel = finding.path.relative_to(root)
        lines.append(f"[{finding.severity}] {finding.category}: {rel}:{finding.line_no}")
        lines.append(f"  {finding.message}")
        lines.append(f"  {finding.line.strip()}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def should_fail(findings: list[Finding], threshold: str | None) -> bool:
    if threshold is None:
        return False
    min_level = SEVERITY_ORDER[threshold]
    return any(SEVERITY_ORDER[finding.severity] >= min_level for finding in findings)


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    targets = resolve_targets(root, args.paths)
    if not targets:
        print("No valid files or directories to audit.", file=sys.stderr)
        return 2

    findings: list[Finding] = []
    for path in iter_source_files(targets):
        findings.extend(audit_file(path))

    print(format_report(root, findings), end="")
    return 1 if should_fail(findings, args.fail_on) else 0


if __name__ == "__main__":
    raise SystemExit(main())
