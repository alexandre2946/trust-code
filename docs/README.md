# TRUST Documentation — Build System

This directory contains the Doxygen configuration, Markdown source pages, and
satellite scripts that together produce the TRUST HTML documentation.

---

## Independence from TRUST's build system

The documentation pipeline is **fully orthogonal** to TRUST's classic
`./configure` + `make optim` build. To build the docs you do **not**
need to:

- run `./configure`
- compile any TRUST target
- have a working `$exec_opt` / `$exec_debug` binary
- source `env_TRUST.sh`

Conversely, building TRUST itself never invokes the doc pipeline — no
`./configure` flag turns it on, and no `make` target in the root
Makefile recurses into `docs/`. The two trees are decoupled.


You can build the documentation either:

- from inside a TRUST environment configured **with** conda:
    - everything is already installed; `make` uses the conda Python at
      `exec/python/` (it ships `pydantic` and `trustify`) and works offline.
    - that conda Python is picked up even if you have **not** sourced
      `env_TRUST.sh`.

- from a TRUST environment configured **`-without-conda`**, or from a clean
  clone with no `./configure` at all. Here `make` falls back to a bare system
  Python (the `exec/python/bin/python` symlink, or `python3` on `PATH`), and
  installing into a shared/system site is **opt-in** — pick one:
    - **(recommended)** create + activate a venv, then
      `pip install "pydantic>=2,<3"`, then run `make`; or
    - set `TRUSTIFY_USER_INSTALL=1` to let `make` do a `pip --user` install
      (you must still provide `pydantic`: `pip install --user "pydantic>=2,<3"`).
    - otherwise `make` stops with an actionable error rather than silently
      writing into a system Python.

- in all cases `externalpackages/` is optional — if it is missing, the doxygen
  binary and MathJax are downloaded from the internet (the only step that truly
  needs network access). See [Prerequisites](#prerequisites) for the full list.


---

## Quick start

From a fresh clone:

```bash
cd docs
make
```

The resulting HTML is written to `docs/build/html/`. Doxygen itself is
unpacked from `externalpackages/` into `exec/doxygen/bin/` on first use
so all users get the same version — the extraction is driven by the
`docs/Makefile` itself (no separate `./configure` or installer script).

If the bundled Linux binary is not portable to your machine, point the
build at any other doxygen installation by exporting
`TRUST_DOXYGEN_BINARY`:

```bash
TRUST_DOXYGEN_BINARY=$(command -v doxygen) make
```

If you don't have a working doxygen at all, see
[Building doxygen from source](#building-doxygen-from-source) below.

---

## Building doxygen from source

When the bundled Linux x86_64 binary does not work on your machine
(other OS, libc mismatch, aarch64, ...), build doxygen yourself and
export `TRUST_DOXYGEN_BINARY`. TRUST tracks **doxygen 1.16.1**; newer
versions usually work too but may produce slightly different HTML
chrome.

```bash
# Requires cmake >= 3.16, a C++17 compiler, flex, bison, python3
git clone --branch Release_1_16_1 --depth 1 https://github.com/doxygen/doxygen.git
cmake -S doxygen -B doxygen/build -G "Unix Makefiles"
cmake --build doxygen/build -j

# Point TRUST at the freshly built binary and run the doc build.
export TRUST_DOXYGEN_BINARY=$PWD/doxygen/build/bin/doxygen
cd $TRUST_ROOT/docs && make
```

Add `TRUST_DOXYGEN_BINARY` to your shell rc file if you want the
override to stick across sessions.

---

## Prerequisites

| Tool | Version | Notes |
|---|---|---|
| `python3` | `>=3.10` | Required by trustify and the satellite scripts |
| `pydantic` | `>=2,<3` | trustify's only runtime dep, installed with `--no-deps`, so it must already be importable from the Python `make` resolves. A conda TRUST env ships it; off-env, install it yourself (`pip install "pydantic>=2,<3"` in a venv, or `pip install --user "pydantic>=2,<3"` alongside `TRUSTIFY_USER_INSTALL=1`). Debian's `python3-pydantic` is v1 and too old. |
| `bash` | any | Hard-set as `SHELL` in `common.mk`; the `run_doxygen` recipe uses process substitution (`> >(...)`) which POSIX `/bin/sh` does not support. Pre-installed on Linux and macOS; **Alpine/busybox**: `apk add bash`. |
| `curl` | any | Used by `install_doxygen` / `install_mathjax` **only when the matching tarball is missing from `externalpackages/`**. If both tarballs are present (the default in a freshly cloned tree), `curl` is not needed at all. |
| `graphviz` | any | Optional — required for class/call graphs |
| `tar` | any | Used once to unpack the bundled doxygen tarball |

The `bash` requirement is checked at parse time in `common.mk`; a missing
`/bin/bash` fails the build immediately with a pointer to this table
rather than crashing inside a recipe.

Doxygen itself is **not** a prerequisite — it is shipped with TRUST in
`externalpackages/` and extracted automatically by the Makefile (or
swapped out via `TRUST_DOXYGEN_BINARY`).

---

## Makefile targets

Run `make help` for an up-to-date list. Here is a description of the various Makefile targets that can invoked from `docs` directory :

### Combo targets

| Target | What it does |
|---|---|
| `all` *(default)* | Regenerate everything (keyword ref + API docs) and build HTML |
| `kw` | Regenerate keyword reference only, then build HTML |
| `api` | Regenerate trustutils API docs only, then build HTML |
| `doxygen` | Build HTML from **existing** sources — no regeneration |

### Low-level targets

| Target | What it does | Defined in |
|---|---|---|
| `clean` | Removes `html/` and generated fragment dirs | `Makefile` |
| `clean_html` | Removes `html/` only | `Makefile` |
| `clean_doxygen` | Removes the extracted doxygen binary | `common.mk` |
| `install_doxygen` | Unpacks the bundled doxygen (skipped if already extracted, or validates `$TRUST_DOXYGEN_BINARY`) | `common.mk` |
| `install_trustify` | `pip install`s the trustify package if its CLI is not on PATH | `common.mk` |
| `trustify_doc` | Runs `trustify generate_markdown` straight into the doxygen source tree | `Makefile` |
| `trustutils_doc` | Generates trustutils Python API Markdown docs | `Makefile` |
| `copy_trustutils` | Copies trustutils `.md` files into the doxygen source tree | `Makefile` |

The `common.mk` rows are shared with every baltik's `docs/Makefile`
via `include` — see [Sharing the pipeline with baltiks](#sharing-the-pipeline-with-baltiks) below.

### Common workflows

```bash
# Full rebuild from scratch
make all

# Only rebuild HTML after editing a narrative .md page
make run_doxygen

# Keyword reference changed (new // XD tags added to C++ sources)
make kw

# trustutils API changed (docstrings updated in .py files)
make api

# Regenerate source fragments without rebuilding HTML yet
make trustify_doc
make trustutils_doc copy_trustutils
```

---

## Directory layout

```
docs/
├── Makefile                        ← TRUST-side build system (described above)
├── common.mk                       ← shared rules included by Makefile AND
│                                     by every baltik's docs/Makefile
├── baltik.mk                       ← rules used ONLY by baltiks
│                                     (Doxyfile generation, doxyfile_chain,
│                                     trustify_doc with -p)
├── README.md                       ← you are here
├── How_to_write_doc.md             ← author-facing guide
├── Doxyfile                        ← doxygen configuration
│
├── content/                        ← hand-written narrative pages (sources)
│   ├── main_page.md                ← @mainpage (root of the page tree)
│   ├── references.bib              ← BibTeX bibliography
│   ├── groups.md                   ← @defgroup declarations (KW, TRUSTUTILS)
│   ├── quick_start.md
│   ├── references.md
│   ├── figures/                    ← images (IMAGE_PATH points here; copied to html/)
│   ├── tutorials/
│   │   └── *.md                    ← user tutorials
│   ├── user_guide/
│   │   ├── index.md                ← @page UserGuide
│   │   ├── general_guide/
│   │   └── numerical_methods/
│   └── dev_corner/
│       ├── index.md                ← @page DevCorner
│       ├── cpp_api.md
│       ├── debug.md
│       ├── stats_profilers.md
│       ├── guidelines/
│       └── tutorials/
│
├── theme/                          ← Doxygen theming + page chrome
│   ├── DoxygenLayout.xml           ← sidebar tab layout
│   ├── doxygen-awesome.css         ← doxygen-awesome base theme
│   ├── doxygen-awesome-sidebar-only.css
│   ├── doxygen-awesome-sidebar-only-darkmode-toggle.css
│   ├── doxygen-awesome-darkmode-toggle.js
│   ├── custom.css                  ← project overrides (load last)
│   ├── header.html                 ← HTML page chrome (GitHub corner, dark-mode init)
│   ├── ISAS_logo.png               ← sidebar logo (copied to html/ via HTML_EXTRA_FILES)
│   └── favicon.ico
│
├── scripts/                        ← post-processors run after Doxygen
│   ├── flatten_navtree.py
│   └── patch_filter.py
│
└── build/                          ← ALL generated artefacts (gitignored)
    ├── generated_content/          ← markdown generated by the pipeline
    │   │                             (mirror of content/, but auto-built;
    │   │                              also scanned by Doxygen via INPUT)
    │   ├── keyword_reference.md    ← @page KeywordsReference (trustify)
    │   ├── kw_ref/                 ← per-family keyword pages (trustify)
    │   │   ├── kw_<parent>.md
    │   │   └── figures/            ← images referenced by kw pages
    │   ├── validation_form_api.md  ← @page Dev_ValidationFormAPI (trustutils)
    │   └── validation_api/         ← per-module trustutils pages
    │       └── trustutils_<module>.md
    └── html/                       ← Doxygen HTML output (open html/index.html)
```

**Never edit anything under `build/`** — `make clean` (= `make clean_build`)
wipes the whole directory, and every target writes its outputs there on
the next run.

---

## How the documentation is structured

### Page hierarchy and navtree

The documentation is organised as a tree of Doxygen **pages**. Each `.md` file
declares itself with `@page <ID> <Title>` and is attached to the tree by a
parent page using `@subpage <ID>`. The root is `main_page.md`, which uses the
special `@mainpage` directive.

```
@mainpage  (main_page.md)
├── @subpage QuickStart
├── @subpage UserGuide
│   ├── @subpage KeywordsReference
│   │   └── @subpage KW_<parent>   ← one per keyword family
│   └── ... (general guide, numerical methods, ...)
├── @subpage UserTutorials
├── @subpage DevCorner
│   ├── @subpage DevTutorials
│   ├── @subpage DevGuidelines
│   ├── @subpage Dev_Debug
│   ├── @subpage Dev_StatsProfilers
│   └── @subpage Dev_ValidationFormAPI
│       ├── @subpage TRUSTUTILS_files
│       ├── @subpage TRUSTUTILS_jupyter_run
│       ├── @subpage TRUSTUTILS_jupyter_plot
│       ├── @subpage TRUSTUTILS_jupyter_filelist
│       ├── @subpage TRUSTUTILS_jupyter_widget
│       ├── @subpage TRUSTUTILS_visitutils
│       └── @subpage TRUSTUTILS_visitutils_msg
└── @subpage References  (bibliography — links to citelist.html)
```

Doxygen renders this tree as the left-hand navigation panel (navtree). For a
step-by-step guide on declaring pages and wiring them into this tree, see
`How_to_write_doc.md` in this same directory.

> **Important:** do **not** add `@ingroup` to pages that are part of the
> `@subpage` tree. In recent Doxygen versions `@ingroup` can conflict with the
> `@subpage` relationship and silently remove the page from the navtree.
> Group membership for search-filter purposes is handled separately — see
> [Groups](#groups) below.

> **Important:** do **not** use Markdown `##` headings on `main_page.md` when
> `@subpage` links are present. Doxygen registers every `##` heading as a
> navtree entry and — due to a Doxygen quirk — inserts it as a child of the
> nearest `@subpage` node rather than as a sibling. Use a raw HTML tag to get
> the visual heading without a navtree entry:
> ```html
> <h2>My Section</h2>
> ```

### Doxygen groups

Two **module groups** are declared in `content/groups.md`:

| Group ID | Label | Purpose |
|---|---|---|
| `KW` | TRUST Keyword Reference | Groups all keyword reference pages |
| `TRUSTUTILS` | Trustutils Python API | Groups all trustutils API pages |

Groups drive the **search filter** tabs (see below) and provide an alternative
index entry point. Group membership is declared with `@defgroup` (once, in
`groups.md`) and referenced from C++ source via `@ingroup` inside `/** */`
comment blocks. Pages are **not** tagged with `@ingroup`; their association to
the Trustutils search tab is handled purely by URL prefix matching in
`patch_filter.py`.

### C++ API

The C++ source tree (`src/`) is scanned directly by Doxygen. Classes,
namespaces, and free functions are extracted from `/** */` doc blocks and
appear in the **C++ API** tab of the sidebar. This content is independent of
the Markdown page tree — it lives in the `usergroup0` section of
`DoxygenLayout.xml`.

### Sidebar layout

`DoxygenLayout.xml` controls which tabs appear in the sidebar and in what
order. See the dedicated [DoxygenLayout.xml](#doxygenLayoutxml) section below
for the full breakdown.

---

## DoxygenLayout.xml

`DoxygenLayout.xml` is Doxygen's layout configuration file. It has two
responsibilities:

1. **`<navindex>`** — defines the sidebar tab structure: which top-level tabs
   exist, in what order, and which sub-tabs appear beneath each one.
2. **Page section layouts** (`<class>`, `<namespace>`, `<file>`, `<group>`,
   `<directory>`) — controls which sections appear on each type of auto-generated
   page and in what order (member lists, description blocks, graphs, etc.).

### Visibility attribute

Every element accepts a `visible` attribute with three possible values:

| Value | Meaning |
|---|---|
| `"yes"` | Always shown |
| `"no"` | Always hidden |
| `"$VARIABLE"` | Read from the matching `VARIABLE` setting in `Doxyfile` |

The `$VARIABLE` form keeps the layout in sync with the Doxyfile without
duplicating settings. For example `visible="$SHOW_INCLUDE_FILES"` mirrors
the `SHOW_INCLUDE_FILES` Doxyfile option directly.

### `<navindex>` — sidebar tabs

The current tab tree is:

```
Home                        ← mainpage tab (always present)
C++ API                     ← usergroup tab (custom label)
  Namespaces
    (namespace list)
    (namespace members)
  Classes
    (class list)
    (alphabetical index)
    (class hierarchy)
    (class members)
  Structs
    (struct list)
    (alphabetical index)
  Files                     ← visible="no" (SHOW_FILES=NO in Doxyfile)
    (file list)
    (global symbols)
```

Narrative pages (Quick Start, User Guide, Keyword Reference, Developer
Corner, ...) do **not** have explicit tabs. They are subpages of `@mainpage`
and are promoted to top-level navtree entries by `flatten_navtree.py` at
build time. Adding explicit `<pages>` or `<modules>` tabs here would
duplicate those entries in the sidebar.

#### Why the Structs tab is hidden

Doxygen 1.9+ merges structs into the combined annotated class list
(`annotated.html`) and no longer generates a separate `structlist.html`.
The `<tab type="structs">` entry would therefore point to a non-existent page,
making the section appear empty.  Structs remain fully documented and accessible
via **Classes → Class List**.

#### Why the Files tab is hidden

`SHOW_FILES = NO` is set in the Doxyfile to suppress individual source-file
pages from the output. However, Doxygen's `SHOW_FILES` option only hides the
page content — the sidebar tab is controlled independently by this layout
file. The `<tab type="files">` entry and its two children (`filelist`,
`globals`) are therefore explicitly set to `visible="no"` here to remove the
tab from the C++ API section.

> **Rule:** whenever you toggle a `SHOW_*` option in the Doxyfile, check
> whether the corresponding tab in `DoxygenLayout.xml` also needs updating.

#### Customising the tab structure

Common operations:

```xml
<!-- Hide a tab entirely -->
<tab type="files" visible="no" .../>

<!-- Rename a tab -->
<tab type="classes" visible="yes" title="C++ Classes">

<!-- Add a custom external link tab at the top level -->
<tab type="user" url="https://example.com" title="External"/>

<!-- Reorder tabs: move them by cutting/pasting the XML elements -->
```

### Page section layouts

The lower part of the file controls what appears on each auto-generated page
type. The most relevant sections for this project are:

| Section | What it controls |
|---|---|
| `<class>` | Layout of every C++ class page (inheritance graph, member list, detailed description, ...) |
| `<namespace>` | Layout of every namespace page |
| `<file>` | Layout of every source-file page (only relevant when `SHOW_FILES = YES`) |
| `<group>` | Layout of Doxygen group pages (`@defgroup`) — used for KW and TRUSTUTILS groups |
| `<directory>` | Layout of directory pages |

Within each section, child elements correspond to rendered blocks.  Setting
`visible="no"` on any child hides that block from the page without affecting
others.  For example, to hide inheritance graphs globally without touching the
Doxyfile:

```xml
<class>
  ...
  <inheritancegraph visible="no"/>
  ...
</class>
```

---

## How the keyword reference is generated

Keywords are extracted from `// XD` tags in the C++ source code by the
**trustify** CLI.  The `trustify_doc` Makefile target invokes
`trustify generate_markdown --out content/user_guide/keyword_ref
--kw-ref-path content/user_guide/keyword_reference.md`, which writes one
`kw_<parent>.md` Doxygen page per keyword family, plus a `figures/` subdir
holding any images referenced via `\includeimage{{...}}`, directly into the
Doxygen source tree.

For the full pipeline description, the XD tag syntax, and how to edit keyword
descriptions, see **`Outils/trustify/README.md`** and the trustify CLI help
(`trustify generate_markdown --help`).

---

## How the trustutils API is generated

The `trustutils` Python package (`Validation/Outils/trustutils/`) is documented
by statically parsing its docstrings with `gen_doc_md.py` (no TRUST environment
required).  Output goes to `Validation/Outils/documentation/build/md/`, then
`make copy_trustutils` imports those files into the Doxygen source tree.

For the full pipeline description, the list of documented modules, and how to
add a new module, see **`Validation/Outils/documentation/README.md`**.

---

## Satellite Python scripts

Two Python scripts post-process the Doxygen output. They are both invoked
automatically by `make run_doxygen`; you should not need to run them manually.

### `flatten_navtree.py` — navtree root flattening and bibliography redirect

Location: `docs/scripts/flatten_navtree.py`
Invoked by: `make run_doxygen` (runs after Doxygen completes)

Doxygen wraps the entire navigation tree under a single root node labelled
with the project name. This script post-processes `navtreedata.js`,
`index.js`, and all `navtreeindex*.js` files to:

1. Keep the project-name node as a **leaf** pointing to the main page (so it
   is still clickable), with no children.
2. **Promote** its former children to top-level entries in the navtree.
3. **Fix breadcrumb indices** in `navtreeindex*.js` — Doxygen encodes
   navigation paths as integer arrays; after hoisting by one level, every
   first index must be incremented by 1.
4. **Remove** the `unshift(0)` call in `navtree.js` that Doxygen hard-codes to
   force traversal through the root node — without this removal, expanded
   navtree items would break after the restructure.
5. **Redirect the Bibliography navtree entry** from `References.html` to
   `citelist.html` (see [Bibliography navtree slot](#bibliography-navtree-slot)
   below).

The end result is a cleaner sidebar where the top-level entries are the major
documentation sections (Quick Start, User Guide, ...) rather than a single
collapsible "TRUST" root.

#### Bibliography navtree slot

Doxygen auto-generates `citelist.html` (the bibliography page) as a special
internal page that it **never includes in the navtree**, even when referenced
with `@subpage citelist`.  To work around this, `doc_md/references.md`
declares a normal `@page References Bibliography` that Doxygen *does* include
in the navtree.  `flatten_navtree.py` then rewrites every occurrence of
`"References.html"` in `index.js` and `navtreeindex*.js` with
`"citelist.html"`, so the sidebar entry points directly at the actual
bibliography content.

**Doxygen version compatibility:** in Doxygen ≤ 1.9.x the full navtree is
stored inline in `navtreedata.js` (no separate `index.js`).  In Doxygen
≥ 1.16, children are lazy-loaded from `index.js`.  The script rewrites
`navtreedata.js` unconditionally and then rewrites `index.js` only if it
exists, so the redirect works on both versions.

**Consequence for maintenance:** if you rename or remove `references.md`, the
Bibliography entry will disappear from the navtree.  If you add a new top-level
section between `DevCorner` and `References` in `main_page.md`, no manual
index update is needed — the breadcrumb rewrite in step 3 is automatic.

### `patch_filter.py` — search bar filter customisation

Location: `docs/scripts/patch_filter.py`
Invoked by: `make run_doxygen` (runs after `flatten_navtree.py`)

Doxygen generates a monolithic search index. This script replaces it with
**five categorised tabs** in the search bar:

| Tab label | Section ID | What it contains |
|---|---|---|
| All | `custom_all` | Every entry (union of all other tabs) |
| Documentation | `trust_docs` | Narrative pages — user guide, tutorials, dev corner |
| Keyword Reference | `trust_kwref` | Keyword reference pages |
| C++ API | `trust_cppapi` | Classes, namespaces, files, functions, ... |
| Trustutils | `trust_trustutils` | `trustutils` Python API pages |

#### How classification works

Each search entry in Doxygen's `searchdata.js` contains a URL. The script
extracts the bare filename (without path or fragment) and applies prefix-based
rules in order:

```python
RULES = [
    (("KeywordsReference", "KW_"),                                    KWREF),
    (("group__TRUSTUTILS", "TRUSTUTILS_", "Dev_ValidationFormAPI"),  TRUSTUTILS),
    (("usergroup0", "class", "namespace", "struct"),                 CPPAPI),
]
```

Entries that do not match any rule fall through to `DOCS` (Documentation tab).

Entries from Doxygen's symbol-level sections (`classes`, `namespaces`,
`files`, `functions`, etc.) are always routed to `CPPAPI` regardless of their
URL, because those sections are exclusively populated by C++ source scanning.

#### Output

The script writes per-section, per-first-character `.js` files
(`<section_id>_<hex_index>.js`) alongside the standard Doxygen search files,
then rewrites `searchdata.js` with updated `indexSectionsWithContent`,
`indexSectionNames`, and `indexSectionLabels` objects that the Doxygen
search widget reads to populate the filter dropdown.

---

## Adding or modifying documentation

| What you want to change | Where to edit |
|---|---|
| A keyword description or attribute | The `// XD` comment in the corresponding `.cpp` file |
| A manually written keyword entry | `src/Kernel/Utilitaires/generic.xd` |
| An embedded text block in a keyword description | `docs/trustify/extras/<name>.md` (referenced from a `// XD` comment via `\input{{<name>}}`) |
| A trustutils class or function doc | The docstring in the relevant `.py` file under `Validation/Outils/trustutils/` |
| Add a new trustutils module to the API docs | Add an entry to `MODULES` in `Validation/Outils/documentation/gen_doc_md.py` |
| A user guide narrative page | The appropriate `.md` file under `doc_md/` |
| The main landing page | `docs/content/main_page.md` |
| The sidebar tab structure | `docs/theme/DoxygenLayout.xml` |
| Search filter categories or rules | `docs/scripts/patch_filter.py` |
| Custom CSS (theming, dark mode overrides) | `docs/theme/custom.css` |
| HTML page structure (header, dark mode toggle) | `docs/theme/header.html` |
| Doxygen version or build options | `docs/common.mk` (`DOXY_VERSION`) and `docs/Doxyfile` |

---

## Sharing the pipeline with baltiks

Baltiks (built via `trust -baltik`) ship a documentation tree with the
same shape as `docs/`. To avoid maintaining two parallel pipelines,
the common bits are factored into two `*.mk` files that live here and
get `include`d from each baltik's `docs/Makefile`:

| File | What it provides | Loaded by |
|---|---|---|
| `docs/common.mk` | `DOXY_*` paths, `TRUST_DOXYGEN_BINARY`-aware `DOXYGEN`, `install_doxygen` / `clean_doxygen`, `install_trustify`, the `doxygen` recipe (run + post-process) | TRUST `docs/Makefile` **and** every baltik `docs/Makefile` |
| `docs/baltik.mk` | `Doxyfile` (generated once from the baltik template), `doxyfile_chain` (regenerated from `BALTIK_DEPENDENCIES`), `trustify_doc -p $project_directory` | every baltik `docs/Makefile` (**not** TRUST) |

Each file has a **contract** at the top — the includer must set named
variables before the `include` line, and the file `$(error ...)`s out
with a clear message if anything is missing. This means changes to the
shared pipeline (e.g. a new post-processor, a new doxygen flag) are
made once in `docs/` and picked up by every baltik on next `make`.

Baltik-only rules are deliberately kept out of `common.mk` so a bug in
`baltik.mk` cannot affect the TRUST doc build (e.g. accidentally
regenerating TRUST's own `Doxyfile` from the baltik template).

---

## Theming and dark mode

The documentation uses the [doxygen-awesome](https://github.com/jothepro/doxygen-awesome-css)
theme in sidebar-only mode with a dark mode toggle.

### CSS load order

Stylesheets are loaded in this order (later files override earlier ones at equal
specificity):

1. `tabs.css` — Doxygen default tab styling
2. `navtree.css` — Doxygen default sidebar/navtree styling (**hardcoded hex colours**)
3. `search/search.css` — Doxygen search widget
4. `doxygen.css` — Doxygen base page styles (**hardcoded hex colours**)
5. `doxygen-awesome.css` — doxygen-awesome main theme (CSS variables)
6. `doxygen-awesome-sidebar-only.css` — sidebar layout adjustments
7. `doxygen-awesome-sidebar-only-darkmode-toggle.css` — toggle button placement
8. `custom.css` — **project overrides** (loads last — always wins)

### Overriding hardcoded colours

Doxygen's own `navtree.css` and `doxygen.css` use hardcoded hex colours on
several elements that `doxygen-awesome.css` does not override.  `custom.css`
corrects these using the doxygen-awesome CSS variables:

| Element | Variable used |
|---|---|
| `#page-nav` (in-page TOC panel, Doxygen 1.16+) | `--side-nav-background`, `--separator-color` |
| `#page-nav-resize-handle` | `--separator-color`, `--page-secondary-foreground-color` |
| `ul.page-outline li.vis` (highlighted TOC row) | `--menu-selected-background` |
| `#nav-path li.navelem:after` (breadcrumb chevrons) | `--page-background-color` |
| `#nav-path li.navelem:hover` | `--menu-selected-background` |
| `div.nav-sync-icon` (sidebar sync button) | `--side-nav-background`, `--separator-color` |

**Rule:** whenever you upgrade Doxygen or doxygen-awesome and notice new white
elements in dark mode, add an override in `custom.css` using the appropriate
CSS variable from `doxygen-awesome.css`.  Never add `!important` unless the
target selector has higher specificity than your rule — the load order alone is
sufficient for same-specificity overrides.

### Dark mode toggle

The toggle is powered by `doxygen-awesome-darkmode-toggle.js` and initialised
in `header.html`.  User preference is stored in `localStorage`.

> **Known behaviour with `file://`:** In Firefox, `localStorage` is scoped
> per-file when pages are opened via the `file://` protocol, so the dark/light
> preference may not persist when navigating between pages locally.  This is a
> browser limitation that does not affect the deployed documentation served from
> a web server.

