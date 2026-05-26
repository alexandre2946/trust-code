# trustutils API documentation generator

This directory contains `gen_doc_md.py`, the script that extracts docstrings
from the `trustutils` Python package and produces Doxygen-compatible Markdown
pages for the TRUST HTML documentation.

> **Two parallel workflows coexist in this directory:**
> - `make md` → `gen_doc_md.py` generates **Doxygen-compatible Markdown** pages,
>   consumed by the main TRUST HTML documentation build (`docs/`).
> - `make html` → **Sphinx** builds a standalone HTML site from the RST sources
>   in `source/`.  This is independent of Doxygen and useful for publishing the
>   `trustutils` Python API on its own (e.g. readthedocs).
>
> Both are described below.

---

## Pipeline position

```
Validation/Outils/trustutils/   ← Python source (docstrings are the source of truth)
        │
        ▼  gen_doc_md.py  (this directory)
        │
  build/md/trustutils_<module>.md    ← one Doxygen @page per module
  build/md/validation_form_api.md    ← index page (@page Dev_ValidationFormAPI)
        │
        ▼  make copy_trustutils  (docs/Makefile)
        │
  docs/content/dev_corner/validation_api/
  docs/content/dev_corner/validation_form_api.md
        │
        ▼  doxygen
        │
  html/  (browsable under "Developer Corner → Validation form API guide")
```

---

## How `gen_doc_md.py` works

The script uses Python's `ast` module to parse each trustutils source file
**statically** — it never imports the package at runtime.  This means it has
zero runtime dependencies: no TRUST environment, no optional Python packages
(numpy, matplotlib, ...) need to be installed.

For each module listed in the `MODULES` table:

1. **Parse** the `.py` source file with `ast`.
2. **Extract** module, class, method, and top-level function docstrings together
   with argument lists.  Private symbols (names starting with `_`) are omitted;
   `__init__` is included only when it carries its own docstring.
3. **Convert** NumPy-style docstring sections to Markdown:
   - `Parameters` / `Attributes` / `Raises` / `Warns` → GFM table (Name | Type | Description)
   - `Returns` / `Yields` → inline `*type* — description` pairs
   - `Notes` → block-quote
   - `Examples` → fenced Python code block
4. **Emit** a Doxygen `@page` Markdown file named `<page_id_lowercase>.md`, with
   anchored headings (`{#anchor}`) for every public symbol so that cross-references
   work from anywhere in the documentation.
5. **Emit** an index page `validation_form_api.md` that lists every module page
   via `@subpage`.

Output is written to `build/md/`.

---

## Documented modules

The `MODULES` list at the top of `gen_doc_md.py` controls which modules are
documented.  Current entries:

| Source path (relative to `trustutils/`) | Page ID | Display name |
|---|---|---|
| `files.py` | `TRUSTUTILS_files` | `trustutils.files` |
| `jupyter/run.py` | `TRUSTUTILS_jupyter_run` | `trustutils.jupyter.run` |
| `jupyter/plot.py` | `TRUSTUTILS_jupyter_plot` | `trustutils.jupyter.plot` |
| `jupyter/filelist.py` | `TRUSTUTILS_jupyter_filelist` | `trustutils.jupyter.filelist` |
| `jupyter/widget.py` | `TRUSTUTILS_jupyter_widget` | `trustutils.jupyter.widget` |
| `visitutils/tools_for_visit.py` | `TRUSTUTILS_visitutils` | `trustutils.visitutils.tools_for_visit` |
| `visitutils/message.py` | `TRUSTUTILS_visitutils_msg` | `trustutils.visitutils.message` |

---

## How to invoke

From this directory:

```bash
make md      # generate Markdown into build/md/ (default target)
make clean   # remove build/
```

Or from `docs/`:

```bash
make api                            # regenerate trustutils API docs and rebuild HTML
make trustutils_doc copy_trustutils # regenerate only, without rebuilding HTML
```

---

## Adding a new module

1. Add an entry to the `MODULES` list in `gen_doc_md.py`:

```python
MODULES = [
    # (path_relative_to_trustutils_root, doxygen_page_id, display_name)
    ("mymodule.py", "TRUSTUTILS_mymodule", "trustutils.mymodule"),
    ...
]
```

2. The page ID **must** start with `TRUSTUTILS_`.  This prefix is what
   `patch_filter.py` (in `docs/`) uses to route the page to the
   *Trustutils* search tab in the generated HTML (see the `RULES` table in
   `docs/scripts/patch_filter.py`).

3. Run `make md` and check that `build/md/trustutils_mymodule.md` was written.

4. The new page will appear in the Doxygen output after the next `make api` (or
   `make all`) in `docs/`.

---

## Editing module documentation

Edit the docstrings directly in the relevant `.py` source file under
`Validation/Outils/trustutils/`.  Use NumPy-style sections (`Parameters`,
`Returns`, `Examples`, ...) — `gen_doc_md.py` converts them automatically.
No special markup is needed: plain prose also renders correctly.

A missing source file produces a warning on stderr and is skipped; it does not
abort the build.

---

## Sphinx workflow (standalone HTML / readthedocs)

The `source/` directory contains RST sources for a **Sphinx**-based
documentation site.  This is independent of the Doxygen build and can be used
to publish the `trustutils` API on its own (e.g. on readthedocs).

Prerequisites:

```bash
pip install -U Sphinx
```

Targets:

```bash
make html       # build HTML site into build/html/
make latexpdf   # build PDF via LaTeX into build/latex/
make clean      # remove build/ (shared with the md target above)
```

The RST sources live in `source/` and are maintained separately from the
docstrings — changes to the Python source are **not** automatically reflected
in the Sphinx output (unlike `gen_doc_md.py`, which always reads the live
source).  Keep both in sync when updating public API documentation.
