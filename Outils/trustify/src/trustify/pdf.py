"""Build a PDF of the trustify-generated keyword reference.

Drives doxygen's LaTeX backend over a `trustify generate_markdown`
output to produce a single hyperlinked PDF. Used by
`api.generate_pdf` (and thus `trustify generate_pdf`).

stdlib-only by design: a near-identical standalone copy lives at
`pdf/keyword_reference_to_pdf.py` for users who cannot install
trustify but already have the markdown. The two are deliberately
duplicated; keep them in sync (see pdf/README.md). Functions here
raise exceptions instead of exiting so the CLI maps them to exit
codes; the standalone script keeps its own die()/sys.exit flow.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Embedded twin of pdf/Doxyfile.pdf. Kept as a Python string (not a
# package-data file) so it ships with this module under any build
# backend / install mode with zero packaging config. Path-free:
# doxygen expands $(KW_INPUT_DIR) / $(KW_OUTPUT_DIR) from the env that
# run_doxygen sets. If you edit this, mirror pdf/Doxyfile.pdf.
_DOXYFILE = r"""# Minimal Doxyfile for building a PDF of the trustify-generated
# keyword reference. KW_INPUT_DIR points at the
# `trustify generate_markdown --out` directory; KW_OUTPUT_DIR is a
# scratch build dir. Both are exported into doxygen's environment by
# trustify.pdf.run_doxygen and expanded via $(VAR) below.

PROJECT_NAME           = "TRUST Keywords Reference"
OUTPUT_DIRECTORY       = $(KW_OUTPUT_DIR)
INPUT                  = $(KW_INPUT_DIR)
IMAGE_PATH             = $(KW_INPUT_DIR)/figures
FILE_PATTERNS          = *.md
RECURSIVE              = NO
# keyword_reference.md already carries an @page directive; leaving
# USE_MDFILE_AS_MAINPAGE empty avoids a conflict.
USE_MDFILE_AS_MAINPAGE =

QUIET                  = YES
WARN_LOGFILE           = $(KW_OUTPUT_DIR)/doxygen_warnings.log

# Disable everything we don't need for a PDF of pure markdown.
GENERATE_HTML          = NO
GENERATE_XML           = NO
GENERATE_RTF           = NO
GENERATE_MAN           = NO
HAVE_DOT               = NO
SOURCE_BROWSER         = NO
ALPHABETICAL_INDEX     = NO

# LaTeX -> PDF (clickable cross-references via PDF_HYPERLINKS).
GENERATE_LATEX         = YES
LATEX_OUTPUT           = latex
USE_PDFLATEX           = YES
PDF_HYPERLINKS         = YES
LATEX_BATCHMODE        = YES
COMPACT_LATEX          = YES
PAPER_TYPE             = a4
# Trustify-generated markdown uses \text and \overline (amsmath) plus
# \langle / \rangle (amssymb) inside math mode; without these, pdflatex
# errors out on every formula even though it still emits a (broken) PDF.
EXTRA_PACKAGES         = amsmath amssymb
"""

DOXYFILE_NAME = "Doxyfile.pdf"


def validate_markdown_dir(raw: str | Path) -> Path:
    """Resolve and sanity-check a `generate_markdown` output dir.

    Returns the absolute path. Raises ValueError on any hard failure;
    a missing `figures/` is a stderr warning, not an error.
    """
    p = Path(raw).resolve()
    # is_dir() can raise on restricted filesystems (EACCES instead of
    # "no such file"); treat any OSError as "not a directory".
    try:
        is_dir = p.is_dir()
    except OSError:
        is_dir = False
    if not is_dir:
        raise ValueError(f"{p} is not a directory. Run `trustify generate_markdown --out {p}` first.")
    if not (p / "keyword_reference.md").is_file():
        raise ValueError(
            f"{p}/keyword_reference.md is missing; this does not look like a `trustify generate_markdown` output."
        )
    if not any(p.glob("kw_*.md")):
        raise ValueError(f"no kw_*.md files in {p}; this does not look like a `trustify generate_markdown` output.")
    if not (p / "figures").is_dir():
        print(
            f"trustify generate_pdf: warning: {p}/figures/ is absent; images in the markdown will not resolve.",
            file=sys.stderr,
        )
    return p


def resolve_doxygen(explicit: str | None = None) -> Path:
    """Find a doxygen binary.

    Resolution order, first match wins: ``explicit`` ->
    ``$TRUST_DOXYGEN_BINARY`` -> ``doxygen`` on PATH. The chosen binary
    is probed with ``--version``; anything whose first stdout line is
    not version-shaped (e.g. "1.16.1") is rejected. Raises RuntimeError
    on failure.
    """
    if explicit:
        candidate = explicit
    elif os.environ.get("TRUST_DOXYGEN_BINARY"):
        candidate = os.environ["TRUST_DOXYGEN_BINARY"]
    else:
        candidate = shutil.which("doxygen")

    if not candidate:
        raise RuntimeError(
            "doxygen not found. Resolution order: --doxygen <path>, $TRUST_DOXYGEN_BINARY, then `doxygen` on PATH."
        )

    path = Path(candidate)
    try:
        out = subprocess.run([str(path), "--version"], capture_output=True, text=True, check=False)
    except OSError as exc:
        raise RuntimeError(f"{path} could not be executed: {exc}") from exc

    first = (out.stdout.splitlines() or [""])[0].strip()
    if not re.match(r"^\d+\.\d+", first):
        raise RuntimeError(f"{path} does not look like a doxygen binary (`--version` produced: {first!r}).")
    return path


def probe_latex_toolchain() -> None:
    """Raise RuntimeError if pdflatex or make is missing.

    Called BEFORE doxygen so a missing TeX Live surfaces immediately
    rather than after a multi-minute doxygen pass.
    """
    missing: list[str] = []
    if not shutil.which("pdflatex"):
        missing.append(
            "pdflatex (install texlive-latex-base + texlive-latex-recommended "
            "on Debian/Ubuntu, or the MacTeX/MiKTeX equivalent)"
        )
    if not shutil.which("make"):
        missing.append("make (install GNU make)")
    if missing:
        raise RuntimeError("missing LaTeX toolchain: " + "; ".join(missing))


def run_doxygen(doxygen: Path, input_dir: Path, build_dir: Path) -> None:
    """Write the embedded Doxyfile into build_dir and run doxygen on it.

    KW_INPUT_DIR / KW_OUTPUT_DIR are exported into doxygen's env (the
    Doxyfile expands them via $(VAR)). Raises RuntimeError on non-zero
    exit, pointing at the warnings log when present.
    """
    doxyfile = build_dir / DOXYFILE_NAME
    doxyfile.write_text(_DOXYFILE)
    env = os.environ.copy()
    env["KW_INPUT_DIR"] = str(input_dir)
    env["KW_OUTPUT_DIR"] = str(build_dir)
    try:
        subprocess.run([str(doxygen), str(doxyfile)], cwd=build_dir, env=env, check=True)
    except subprocess.CalledProcessError as exc:
        warnlog = build_dir / "doxygen_warnings.log"
        hint = f" (see {warnlog})" if warnlog.exists() else ""
        raise RuntimeError(f"doxygen exited with status {exc.returncode}{hint}.") from exc


def run_latex(build_dir: Path) -> None:
    """Run ``make -C <build_dir>/latex`` to drive pdflatex.

    Doxygen emits a Makefile that handles the multi-pass
    pdflatex/makeindex dance. The stale refman.pdf is dropped first so
    a later finalize cannot falsely succeed on an old artifact. On
    failure, make's output and the tail of refman.log are written to
    stderr (LaTeX errors live there; make only says "Error 1") and a
    RuntimeError is raised.
    """
    latex = build_dir / "latex"
    stale = latex / "refman.pdf"
    if stale.exists():
        stale.unlink()

    proc = subprocess.run(["make", "-C", str(latex)], capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        log = latex / "refman.log"
        if log.exists():
            tail = log.read_text(errors="replace").splitlines()[-80:]
            sys.stderr.write("\n----- last 80 lines of refman.log -----\n")
            sys.stderr.write("\n".join(tail) + "\n")
        raise RuntimeError(f"`make -C {latex}` failed (rerun with --keep-build for the full log).")


def finalize_pdf(build_dir: Path, output: Path) -> None:
    """Copy refman.pdf to the requested output path.

    Raises RuntimeError if the PDF was not produced or if the output's
    parent dir is missing (we do NOT auto-mkdir: a typo in the output
    path would silently land the PDF somewhere unintended).
    """
    src = build_dir / "latex" / "refman.pdf"
    if not src.is_file():
        raise RuntimeError(f"{src} was not produced -- check {build_dir}/latex/refman.log.")
    parent = Path(output).parent
    if not parent.is_dir():
        raise RuntimeError(f"parent of output does not exist: {parent}")
    shutil.copy2(src, output)
