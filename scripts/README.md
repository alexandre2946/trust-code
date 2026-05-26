# `scripts/`

Small, generic utilities used elsewhere in the tree.  Each script is
self-contained and should carry its own usage docs (module docstring +
`--help`).  This file is just the index: one short paragraph per script
plus a map of every place it is called from, so renaming or removing a
script is a grep-free decision.

When adding a new script:

1. Drop a `--help`-ready script in this directory.
2. Add an entry below with **Purpose** (one or two sentences) and
   **Used by** (every call site, file path + line, with a one-line
   note on what it does there).
3. Keep entries alphabetical.

---

## `expand_env.py`

**Purpose.** Substitute `@@VAR@@` tokens in a text stream with the
matching environment variable.  Supports stdin/stdout, an `-o` output
file, and `-i` in-place rewrite.  See `expand_env.py --help` for the
full CLI.

**Used by.**

- `bin/trust` (baltik project creation, `-init` branch) — expands
  three templates from `bin/baltik/templates/basic/` into the freshly
  created project:
  - `project.cfg`
  - `docs/Makefile`
  - `docs/content/main_page.md`
- `bin/baltik/templates/basic/docs/Makefile` (`update_doxyfile`
  target) — expands the baltik `Doxyfile` template at first `make`,
  injecting `PROJECT_NAME`. This Makefile gets copied into baltiks, so
  the script will be used implicitly by baltiks
