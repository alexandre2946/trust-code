#!/usr/bin/env python3
"""Build a PDF of the trustify-generated keyword reference.

Consumes the output of `trustify generate_markdown --out <dir>` and
emits a single PDF via doxygen's LaTeX backend.  Standalone tool --
no coupling to docs/Makefile.

Smoke test from a TRUST checkout (the docs/ HTML build already
produces a trustify markdown tree we can reuse):

    python3 Outils/trustify/pdf/keyword_reference_to_pdf.py \\
        docs/build/generated_content/kw_ref \\
        --output /tmp/keyword_reference.pdf
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import NoReturn


SCRIPT_DIR = Path(__file__).resolve().parent
DOXYFILE = SCRIPT_DIR / "Doxyfile.pdf"


def die(msg: str, code: int = 1) -> NoReturn:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def validate_input_dir(raw: str) -> Path:
    """Resolve and sanity-check a `trustify generate_markdown` output dir.

    Returns the absolute path.  Calls die() on any hard failure.
    """
    p = Path(raw).resolve()
    # `is_dir()` calls os.stat under the hood and can raise on sandboxed
    # / restricted filesystems where a non-existent path under a
    # locked-down parent gives EACCES instead of "no such file".  Treat
    # any OSError the same way as "not a directory" -- the diagnostic
    # ("run trustify generate_markdown first") still applies.
    try:
        is_dir = p.is_dir()
    except OSError:
        is_dir = False
    if not is_dir:
        die(f"{p} is not a directory.\n"
            f"       Run `trustify generate_markdown --out {p}` first.")
    if not (p / "keyword_reference.md").is_file():
        die(f"{p}/keyword_reference.md is missing.\n"
            f"       This directory does not look like a "
            f"`trustify generate_markdown` output.")
    if not any(p.glob("kw_*.md")):
        die(f"No kw_*.md files found in {p}.\n"
            f"       This directory does not look like a "
            f"`trustify generate_markdown` output.")
    if not (p / "figures").is_dir():
        warn(f"{p}/figures/ is absent -- images in the markdown "
             f"will not resolve.")
    return p


def resolve_doxygen(explicit: str | None) -> Path:
    """Find a doxygen binary.

    Resolution order, first match wins:
      1. ``explicit`` (the --doxygen CLI flag)
      2. ``$TRUST_DOXYGEN_BINARY`` env var
      3. ``doxygen`` on PATH

    The chosen binary is probed with ``--version``; we reject anything
    whose first stdout line does not start with a version-shaped token
    (e.g. "1.16.1") to catch unrelated binaries that happen to be
    named ``doxygen``.
    """
    candidate: str | None = None
    if explicit:
        candidate = explicit
    elif os.environ.get("TRUST_DOXYGEN_BINARY"):
        candidate = os.environ["TRUST_DOXYGEN_BINARY"]
    else:
        candidate = shutil.which("doxygen")

    if not candidate:
        die("doxygen not found.  Resolution order:\n"
            "       1. --doxygen <path>\n"
            "       2. $TRUST_DOXYGEN_BINARY\n"
            "       3. `doxygen` on PATH")

    path = Path(candidate)
    try:
        out = subprocess.run([str(path), "--version"],
                             capture_output=True, text=True, check=False)
    except OSError as exc:
        die(f"{path} could not be executed: {exc}")

    first = (out.stdout.splitlines() or [""])[0].strip()
    if not re.match(r"^\d+\.\d+", first):
        die(f"{path} does not look like a doxygen binary "
            f"(`--version` produced: {first!r}).")
    return path


def make_build_dir(input_dir: Path, keep_build: bool) -> Path:
    """Create (or reuse) the scratch build directory.

    With --keep-build, the dir lives at ``<input_dir>/_pdf_build`` and
    is reused across runs (doxygen overwrites its outputs; pdflatex
    regenerates refman.pdf), which makes iteration faster than wiping.

    Without --keep-build, we use a one-shot ``tempfile.mkdtemp`` that
    the caller is expected to remove in a finally clause.
    """
    if keep_build:
        bd = input_dir / "_pdf_build"
        bd.mkdir(exist_ok=True)
        return bd
    return Path(tempfile.mkdtemp(prefix="trustify-pdf-"))


def run_doxygen(doxygen: Path, input_dir: Path, build_dir: Path) -> None:
    """Run doxygen on the trustify markdown.

    The Doxyfile is path-free and reads ``$(KW_INPUT_DIR)`` and
    ``$(KW_OUTPUT_DIR)`` from the environment -- same trick the shared
    ``doxyfile_common`` uses with ``$(TRUST_ROOT)``.
    """
    env = os.environ.copy()
    env["KW_INPUT_DIR"] = str(input_dir)
    env["KW_OUTPUT_DIR"] = str(build_dir)
    try:
        subprocess.run([str(doxygen), str(DOXYFILE)],
                       cwd=build_dir, env=env, check=True)
    except subprocess.CalledProcessError as exc:
        warnlog = build_dir / "doxygen_warnings.log"
        hint = f" (see {warnlog})" if warnlog.exists() else ""
        die(f"doxygen exited with status {exc.returncode}{hint}.")


def run_latex(build_dir: Path) -> None:
    """Run ``make -C <build_dir>/latex`` to drive pdflatex.

    Doxygen emits a Makefile that handles the multi-pass
    pdflatex/makeindex dance, so we don't reimplement it.

    On failure we surface make's stdout/stderr AND the tail of
    refman.log -- LaTeX errors are buried there and ``make``'s own
    output only says "Error 1".
    """
    latex = build_dir / "latex"
    # Drop a stale refman.pdf so finalize_pdf cannot falsely succeed
    # if pdflatex errors out before producing a fresh one.
    stale = latex / "refman.pdf"
    if stale.exists():
        stale.unlink()

    proc = subprocess.run(
        ["make", "-C", str(latex)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        log = latex / "refman.log"
        if log.exists():
            tail = log.read_text(errors="replace").splitlines()[-80:]
            sys.stderr.write("\n----- last 80 lines of refman.log -----\n")
            sys.stderr.write("\n".join(tail) + "\n")
        die(f"`make -C {latex}` failed "
            f"(rerun with --keep-build for the full log).")


def finalize_pdf(build_dir: Path, output: Path) -> None:
    """Copy refman.pdf to the requested output path."""
    src = build_dir / "latex" / "refman.pdf"
    if not src.is_file():
        die(f"{src} was not produced -- check "
            f"{build_dir}/latex/refman.log.")
    parent = output.parent
    if not parent.is_dir():
        # Don't auto-mkdir: a typo in --output would silently land the
        # PDF in an unintended location.
        die(f"parent of --output does not exist: {parent}")
    shutil.copy2(src, output)


def probe_latex_toolchain() -> None:
    """Fail fast if pdflatex or make is missing.

    Called BEFORE doxygen so the user doesn't sit through a long
    doxygen pass only to fail at the LaTeX step.
    """
    missing: list[tuple[str, str]] = []
    if not shutil.which("pdflatex"):
        missing.append((
            "pdflatex",
            "install texlive-latex-base + texlive-latex-recommended "
            "(Debian/Ubuntu) or the MacTeX/MiKTeX equivalent",
        ))
    if not shutil.which("make"):
        missing.append(("make", "install GNU make"))
    if missing:
        for name, hint in missing:
            print(f"ERROR: {name} not found in PATH "
                  f"({hint})", file=sys.stderr)
        sys.exit(1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="keyword_reference_to_pdf.py",
        description="Build a PDF from a `trustify generate_markdown` output.",
    )
    p.add_argument(
        "input_dir",
        help="Directory written by `trustify generate_markdown --out`.",
    )
    p.add_argument(
        "--output", "-o", default=None,
        help="PDF output path "
             "(default: <input-dir>/keyword_reference.pdf).",
    )
    p.add_argument(
        "--keep-build", action="store_true",
        help="Keep the scratch build dir for inspection.",
    )
    p.add_argument(
        "--doxygen", default=None,
        help="Path to doxygen binary "
             "(overrides $TRUST_DOXYGEN_BINARY and PATH).",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_dir = validate_input_dir(args.input_dir)
    output = (Path(args.output).resolve() if args.output
              else input_dir / "keyword_reference.pdf")
    doxygen = resolve_doxygen(args.doxygen)
    probe_latex_toolchain()
    build_dir = make_build_dir(input_dir, args.keep_build)
    try:
        run_doxygen(doxygen, input_dir, build_dir)
        run_latex(build_dir)
        finalize_pdf(build_dir, output)
    finally:
        if not args.keep_build:
            shutil.rmtree(build_dir, ignore_errors=True)
    print(f"[trustify-pdf] Wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
