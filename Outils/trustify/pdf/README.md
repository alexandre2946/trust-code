@page Trustify_KeywordsReferencePDF Keyword Reference Manual as PDF

`keyword_reference_to_pdf.py` turns the output of
@ref Trustify "`trustify generate_markdown`" into a single
hyperlinked PDF of the TRUST keyword reference manual.

## Relationship to `trustify generate_pdf`

The primary way to build this PDF is now the CLI:

```bash
trustify generate_pdf --out /tmp/keyword_reference.pdf
```

`trustify generate_pdf` runs the whole pipeline (schema → markdown →
doxygen → LaTeX → PDF) end-to-end, and `--from-markdown DIR` reuses an
existing `generate_markdown` output. Its orchestration lives in
`src/trustify/pdf.py`, with the Doxyfile embedded as a Python string.

This standalone `keyword_reference_to_pdf.py` (and its sibling
`Doxyfile.pdf`) is a **deliberate duplicate** kept for the exotic case of
a user who *cannot install trustify* but already has a markdown tree (e.g.
obtained from elsewhere) and only Python stdlib + doxygen + TeX Live. It
is stdlib-only and self-contained. It carries no automated tests and is
maintained on demand; the tested copy is `src/trustify/pdf.py`. If you
change one, mirror the other.

It runs doxygen (with a deliberately minimal Doxyfile — see *Files in
this directory* below) on the trustify-emitted markdown tree, then
drives `make -C latex` to produce `refman.pdf` via `pdflatex`. The
script is standalone — it does not depend on `$TRUST_ROOT/docs/Makefile`
or the shared `doxyfile_common`.

## Quick start

```bash
# 1. Produce the markdown tree (see `trustify generate_markdown --help`).
trustify generate_markdown --out /tmp/kw_ref

# 2. Build the PDF.  The script is on $PATH inside the TRUST env via the
#    bin/ symlink.
keyword_reference_to_pdf.py /tmp/kw_ref --output /tmp/keyword_reference.pdf
```

A successful run prints one line:

```
[trustify-pdf] Wrote /tmp/keyword_reference.pdf
```

The resulting PDF (~290 pages for the current TRUST schema) has:

- A table of contents and PDF bookmarks (outline) for every keyword
  family.
- Clickable cross-references — every `@ref` / `@subpage` in the
  markdown becomes a `/Link` annotation that jumps to the target.
- Embedded figures from `<input-dir>/figures/`.

## CLI

```text
keyword_reference_to_pdf.py <input-dir> [--output PATH]
                                        [--keep-build]
                                        [--doxygen PATH]
```

| Argument         | Meaning                                                                                                                                       |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `<input-dir>`    | Required. Directory written by `trustify generate_markdown --out`. Must contain `keyword_reference.md` and at least one `kw_*.md`.            |
| `--output PATH`  | PDF output path. Default: `<input-dir>/keyword_reference.pdf`. The parent directory must already exist — typos are surfaced rather than masked. |
| `--keep-build`   | Keep the scratch build directory at `<input-dir>/_pdf_build/` instead of using a one-shot `tempfile.mkdtemp` that is wiped on exit. Useful to inspect `latex/refman.log` after a failure, or to iterate faster (doxygen overwrites; re-running reuses the existing tree). |
| `--doxygen PATH` | Explicit doxygen binary. Resolution order otherwise: `$TRUST_DOXYGEN_BINARY` → `doxygen` on `$PATH`. The chosen binary is probed with `--version`; a non-doxygen file by that name is rejected. |

## Requirements

The tool fails fast (before invoking doxygen) when any of these is
missing:

- **doxygen** — any `1.x` release. The TRUST docs build pins
  `1.16.1` from `externalpackages/`; this script will use whatever
  `doxygen` resolves to (see `--doxygen` resolution above).
- **`pdflatex`** — TeX Live, e.g. `texlive-latex-base` +
  `texlive-latex-recommended` on Debian / Ubuntu, or the MacTeX /
  MiKTeX equivalent.
- **GNU `make`** — drives the multi-pass
  `pdflatex` / `makeindex` dance via the `Makefile` doxygen emits
  into the `latex/` output directory.

## Files in this directory

| File                          | Role                                                                                              |
| ----------------------------- | ------------------------------------------------------------------------------------------------- |
| `keyword_reference_to_pdf.py` | CLI entry point. Stdlib-only Python orchestrator (~270 LOC).                                      |
| `Doxyfile.pdf`                | Minimal Doxyfile template. Path-free; the script exports `KW_INPUT_DIR` and `KW_OUTPUT_DIR` into doxygen's environment, and the Doxyfile expands `$(KW_INPUT_DIR)` / `$(KW_OUTPUT_DIR)` (same trick `$TRUST_ROOT/docs/doxyfile_common` uses with `$(TRUST_ROOT)`). |

The Doxyfile turns OFF every output format and feature that doesn't
serve a markdown-only PDF (HTML, XML, graphs, source browser,
alphabetical index, …), turns ON `GENERATE_LATEX` + `USE_PDFLATEX` +
`PDF_HYPERLINKS`, and pulls in `amsmath` + `amssymb` via
`EXTRA_PACKAGES` (the trustify-generated markdown uses `\text`,
`\overline`, `\langle` / `\rangle` in math mode — without those
packages `pdflatex` exits non-zero on every formula).

## Pipeline

1. Validate `<input-dir>` — must exist, contain `keyword_reference.md`
   and at least one `kw_*.md`. Missing `figures/` is a warning, not a
   hard error.
2. Resolve the doxygen binary (`--doxygen` → `$TRUST_DOXYGEN_BINARY`
   → `$PATH`).
3. Probe `pdflatex` and `make` on `$PATH`. Fail fast — surfacing a
   missing LaTeX install after a multi-minute doxygen run would be
   unfriendly.
4. Create a scratch build directory (`mkdtemp` by default;
   `<input-dir>/_pdf_build/` under `--keep-build`).
5. Run `doxygen Doxyfile.pdf` with `KW_INPUT_DIR` / `KW_OUTPUT_DIR`
   in the environment.
6. Run `make -C <build>/latex`. On non-zero exit the script prints
   `make`'s output and the last 80 lines of `refman.log` to stderr —
   LaTeX errors live there, `make`'s own output only says `Error 1`.
7. Copy `<build>/latex/refman.pdf` to `--output`.
8. Remove the scratch directory unless `--keep-build` was passed.

## Troubleshooting

- **"`pdflatex` not found in PATH"** — install TeX Live (see
  *Requirements* above).
- **"`make -C ... failed`"** — re-run with `--keep-build` and inspect
  `<input-dir>/_pdf_build/latex/refman.log`. The script's tail-dump to
  stderr usually contains the actual error, but the full log carries
  the surrounding context. If the failure mentions an undefined
  control sequence, the trustify markdown probably grew a new
  math-mode macro that needs another LaTeX package — extend
  `EXTRA_PACKAGES` in `Doxyfile.pdf`.
- **PDF looks empty or stops mid-document** — check
  `<build>/doxygen_warnings.log` for `@page` / `@subpage` /
  `@ref` resolution failures in the markdown (doxygen will skip
  unresolved cross-references silently in the LaTeX output).
- **"`doxygen` not found"** — set `TRUST_DOXYGEN_BINARY` to the
  binary's absolute path, or pass `--doxygen /path/to/doxygen`. The
  bundled `1.16.1` lives at `$TRUST_ROOT/exec/doxygen/bin/doxygen`
  after the docs build has run.
