#!/usr/bin/env python3
"""Generate (or validate) the root index.html from index.html.in.

The version number is read from ``$TRUST_ROOT/VERSION`` and injected into the
template located beside this script (``index.html.in``). Substitution uses the
usual ``${key}`` Python template syntax.

Two keys are filled:

* ``${VERSION}``           the version string, e.g. ``1.9.8_beta``
* ``${READTHEDOCS_LABEL}`` the readthedocs slug used in documentation URLs:
                           ``latest`` for a *_beta version, otherwise the
                           release tag ``v<x.x.x>`` (e.g. ``v1.9.8``).

By default the script runs in dry-run mode: it renders the expected index.html
and compares it against the existing ``$TRUST_ROOT/index.html``, reporting any
problem (missing file or diff) and exiting non-zero. Pass ``--apply`` to write
the file.
"""
from __future__ import annotations

import argparse
import difflib
import os
import sys
from pathlib import Path
from string import Template

TEMPLATE_PATH = Path(__file__).resolve().parent / "index.html.in"


def read_version(version_file: Path) -> str:
    """Return the stripped contents of the VERSION file."""
    return version_file.read_text(encoding="utf-8").strip()


def compute_substitutions(raw_version: str) -> dict[str, str]:
    """Map a raw VERSION string to the template substitution values."""
    core = raw_version[1:] if raw_version.startswith("v") else raw_version
    is_beta = core.endswith("_beta")
    return {
        "VERSION": core,
        "READTHEDOCS_LABEL": "next" if is_beta else f"v{core}",
    }


def render(template_text: str, substitutions: dict[str, str]) -> str:
    """Render the template, raising on any unknown/malformed ${placeholder}."""
    return Template(template_text).substitute(substitutions)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write the generated index.html (default: dry-run validation only)",
    )
    args = parser.parse_args(argv)

    trust_root = os.environ.get("TRUST_ROOT")
    if not trust_root:
        print("error: TRUST_ROOT is not set", file=sys.stderr)
        return 1
    trust_root = Path(trust_root)

    version_file = trust_root / "VERSION"
    if not version_file.is_file():
        print(f"error: version file not found: {version_file}", file=sys.stderr)
        return 1

    if not TEMPLATE_PATH.is_file():
        print(f"error: template not found: {TEMPLATE_PATH}", file=sys.stderr)
        return 1

    substitutions = compute_substitutions(read_version(version_file))
    try:
        rendered = render(TEMPLATE_PATH.read_text(encoding="utf-8"), substitutions)
    except (KeyError, ValueError) as exc:
        print(f"error: failed to render {TEMPLATE_PATH}: {exc}", file=sys.stderr)
        return 1

    index_path = trust_root / "index.html"

    if args.apply:
        index_path.write_text(rendered, encoding="utf-8")
        print(
            f"wrote {index_path} "
            f"(VERSION={substitutions['VERSION']}, "
            f"READTHEDOCS_LABEL={substitutions['READTHEDOCS_LABEL']})"
        )
        return 0

    # Dry-run: validate the existing index.html against the rendered template.
    if not index_path.is_file():
        print(
            f"error: {index_path} does not exist (run with --apply to create it)",
            file=sys.stderr,
        )
        return 1

    current = index_path.read_text(encoding="utf-8")
    if current == rendered:
        print(f"OK: {index_path} is up to date (VERSION={substitutions['VERSION']})")
        return 0

    diff = difflib.unified_diff(
        current.splitlines(keepends=True),
        rendered.splitlines(keepends=True),
        fromfile=f"{index_path} (current)",
        tofile="expected (from index.html.in)",
    )
    sys.stdout.writelines(diff)
    print(
        f"error: {index_path} is out of date (run with --apply to update it)",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
