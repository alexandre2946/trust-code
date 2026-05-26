#!/usr/bin/env python3
"""
expand_env.py — replace every @@VAR@@ token in a text stream with the
corresponding environment variable value.

Examples
--------
    expand_env.py < in.txt > out.txt          # stdin → stdout
    expand_env.py in.txt                      # in.txt → stdout
    expand_env.py in.txt -o out.txt           # in.txt → out.txt
    expand_env.py -i in.txt                   # rewrite in.txt in place
    cat in.txt | expand_env.py -o out.txt     # stdin → out.txt

Exit codes
----------
    0  every @@VAR@@ was resolved
    1  one or more variables were missing (left unexpanded so the gap is
       visible), or an I/O error occurred
    2  bad CLI usage
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

# @@VAR_NAME@@ — letters, digits, underscores; must start with letter or _
_ENV_VAR_RE = re.compile(r"@@([A-Za-z_][A-Za-z0-9_]*)@@")


def expand_env_vars(text: str) -> tuple[str, list[str]]:
    """Return (expanded_text, unique_missing_variable_names).

    Missing variables are left unexpanded in the output so the gap is visible.
    """
    missing: list[str] = []
    seen: set[str] = set()

    def replacer(match: re.Match) -> str:
        name = match.group(1)
        value = os.environ.get(name)
        if value is None:
            if name not in seen:
                seen.add(name)
                missing.append(name)
            return match.group(0)
        return value

    return _ENV_VAR_RE.sub(replacer, text), missing


def _read(source: Path | None) -> str:
    if source is None:
        return sys.stdin.read()
    return source.read_text(encoding="utf-8")


def _write(destination: Path | None, data: str) -> None:
    if destination is None:
        sys.stdout.write(data)
    else:
        destination.write_text(data, encoding="utf-8")


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="expand_env.py",
        description="Expand @@VAR@@ tokens using the current environment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Use '-' as INPUT or with -o to mean stdin/stdout explicitly.",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="input file (default: stdin, or '-')",
    )
    sink = parser.add_mutually_exclusive_group()
    sink.add_argument(
        "-o",
        "--output",
        default="-",
        help="output file (default: stdout, or '-')",
    )
    sink.add_argument(
        "-i",
        "--in-place",
        action="store_true",
        help="rewrite INPUT in place (requires a real input file)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    input_path: Path | None = None if args.input == "-" else Path(args.input)

    if args.in_place:
        if input_path is None:
            print("error: --in-place requires a real input file (not stdin)", file=sys.stderr)
            return 2
        output_path: Path | None = input_path
    else:
        output_path = None if args.output == "-" else Path(args.output)

    try:
        source = _read(input_path)
    except FileNotFoundError:
        print(f"error: input file not found: {input_path}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"error: cannot read {input_path}: {exc}", file=sys.stderr)
        return 1

    expanded, missing = expand_env_vars(source)

    if missing:
        label = str(input_path) if input_path is not None else "<stdin>"
        print(
            f"warning: {label}: the following environment variables are not set "
            f"and were left unexpanded:\n  " + "\n  ".join(missing),
            file=sys.stderr,
        )

    try:
        _write(output_path, expanded)
    except OSError as exc:
        print(f"error: cannot write {output_path}: {exc}", file=sys.stderr)
        return 1

    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
