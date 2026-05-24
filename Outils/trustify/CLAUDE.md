# trustify — TRUST Dataset Tooling

Python package that turns the **TRUST** CFD code's input-keyword
declarations (`// XD` C++ comments) into a pydantic model + parser, then
uses that schema to validate, format, and document `.data` input files.
The companion LSP server lives in the sibling `Outils/trustify-lsp/`
directory.

This file documents conventions and gotchas for agentic LLMs working on
trustify. The TRUST repo, trustify, and trustify-lsp all live in one git
repository rooted at `trust-code/`.

## XD tag syntax

XD tags are special C++ `//`-comments embedded in TRUST sources (and in
standalone `.xd` files). They declare blocks (TRUST classes) and their
attributes. The scanner lives in
`src/trustify/core/trad2_utilities.py` (`scanOneCppLine`, `_ParseXD`).

### Block declaration (header line)

```
// XD <name> <name_base> <syno> <brace_flag> <description...>
```

- `<name>` — the keyword/class name (lowercase).
- `<name_base>` — the parent class (must already be declared, or be
  `objet_u`). `<name> == <name_base>` is an error.
- `<syno>` — pipe-separated synonyms (e.g. `lire|read`); use the class
  name to indicate "no synonyms".
- `<brace_flag>` — one of:
  - `BRACE` (or legacy `1` / `-3`) — keyword is followed by `{ ... }`.
  - `NO_BRACE` (or `0` / `-2`) — no braces.
  - `INHERITS_BRACE` (or `-1`) — inherits the parent's brace policy.
- `<description...>` — free text used in error messages and docs.

`name_base == listobj` produces a `TRAD2BlockList` instead of a
`TRAD2Block` — the line takes two extra fields (`<itemtype>
<comma_flag>`) before the description.

### Attribute line

```
//   XD attr <name> <type> <syno> <opt_flag> <description...>
```

Lives directly under its owning block (TRAD2 ordering matters — the
scanner ties an attribute to whichever block was most recently opened).

- `<type>` — `entier`, `floattant`, `chaine`, `rien` (`flag` in
  shorthand), `dico`, `listentier`, `chaine(into=["X","Y"])`, or any
  declared block name.
- `<opt_flag>` — `REQ` (mandatory, legacy `0`) or `OPT` (legacy `1`).

### C++-only shortcut tags

Inside `.cpp` files the scanner additionally recognises tags that pull
data directly out of the `Param::ajouter*` boilerplate:

- `// XD_ADD_P <type> <description...>` — placed on a line that already
  contains a `.ajouter(...)` / `.ajouter_flag(...)` / etc. call. The
  attribute name is read from the first argument of `ajouter`.
- `// XD_ADD_DICO <name>` — placed on a line that contains a
  `.dictionnaire(...)` call. Appends `name` to the preceding
  `XD_ADD_P dico ...` attribute's allowed-values list (turns
  `chaine(into=[])` into `chaine(into=["name",])`).

### Multi-line declarations

Long XD lines may be continued with `// XD_CONT`:

```cpp
// XD foo bar baz BRACE this description happens to be very long
// XD_CONT and continues here for another sentence or two
```

`XD_CONT` must immediately follow the line it continues; otherwise the
scanner raises.

### Nested levels (2XD, 3XD)

Tags prefixed with `2` or `3` (`// 2XD ...`, `// 2XD_ADD_P ...`,
`// 2XD_ADD_DICO ...`, `// 2XD_CONT`) declare attributes on a
*secondary* block opened by an enclosing `Read_Sub_Param` /
`Param_Mediums` style C++ helper. The scanner keeps three independent
"current block" slots per file so the secondary attributes attach to
the right parent. `3XD` is the same but one level deeper.

### Other macros scanned (not strictly XD)

The scanner also picks up two `Process::` boilerplate macros to
maintain the synonyms table:

- `Implemente_instanciable(<class>, "<syno1>|<syno2>|...", ...)` —
  registers `<class>`'s synonyms.
- `<class>::add_synonym("<other_name>");` — adds one more synonym.

Both end up in `TRAD2Content.synos`, which `assemble()` folds into the
block's `synos` list.

## Main pipeline

### Generate the schema

`trustify.generate_schema(trust_root=..., projects=[...])` (or the
`trustify generate_schema` CLI) walks each project's `src/` directory,
parses every `// XD` tag into a `TRAD2Content` (the in-memory schema),
then code-generates:

- `TRAD2_trustify` — canonical flat-text schema.
- `source_locations.json` — `block_name → {project, path, line}` map
  (per-attr nested under each block).
- `trustify_gen_pyd.py` — pydantic data model.
- `trustify_gen.py` — parser module (consumes `trustify_gen_pyd`).
- `provenance.json` — at-generation paths + timestamp + version + per-step
  timings. Display-only metadata; never participates in the cache hash.
  Doubles as the **completion sentinel**: `is_cached` requires it, so
  a process killed between writing `trustify_gen.py` and the provenance
  file leaves a directory that future runs correctly treat as
  not-cached. `generate_schema` writes the whole entry into a sibling
  `.<hash>.tmp.<pid>/` directory and `os.replace`s it onto
  `<hash>/` only after every file (provenance included) is on disk —
  concurrent readers either see the previous entry or the complete
  new entry, never a torn write.

**Cross-process cache lock.** `cache.cache_lock(mode, root=...)` is a
`fcntl.flock` advisory lock on `<cache_root>/.lock`. SHARED for
readers (`api._resolve_schema_module` around `import_parser_module`),
EXCLUSIVE for writers (`api.generate_schema` around the write+trim
region in cache mode — override `out=` skips the lock since the user
owns coordination). The kernel auto-releases on process death, so
no stale-lock cleanup. `_resolve_schema_module` retries the import
once on `FileNotFoundError` to cover the narrow race window between
`generate_schema`'s EXCLUSIVE release and its own SHARED acquire (a
concurrent writer's trim may have rmtree'd the dir in that gap).

These five files land in `~/.cache/trustify/<hash>/` (or
`$TRUSTIFY_CACHE_DIR/<hash>/` when that env var is set — for HPC
sites with a quota'd / read-only `$HOME`, ephemeral CI runners, or
multi-user shared TRUST installs that need one shared cache).
Whitespace-only env values fall back to the default. The hash is
computed from:

1. The trustify version string.
2. The canonical TRAD2 text.
3. The canonical JSON text of `source_locations.json`.

Two clones of the same project set at different absolute paths produce
byte-identical schema files (no absolute paths in either input), so
they share a single cache entry. A baltik moving a source file within
itself changes `source_locations.json` content and therefore the
hash — the cache regenerates cleanly.

### Use the schema

`trustify.load_dataset(path, trust_root=...)` and
`trustify.check(...)` / `batch_check(...)` import the cached parser via
`trustify.core.misc_utilities.import_parser_module`. The returned
`SchemaModule` exposes `get_parser_class(name)` /
`get_pyd_class(name)` and a synonym-aware lookup table.

CLI surface (defined in `src/trustify/cli/__init__.py`):

```
trustify generate_schema [--out DIR]
trustify check DATA_FILE
trustify batch-check DATA_PATH [DATA_PATH ...] [--exclude PATTERN] [--follow-symlinks] [--jobs N]
trustify batch-check-projects [--only ROLE[,ROLE...]] [--subdir REL_PATH] [--exclude PATTERN] [--jobs N] ...
trustify format DATA_FILE [--out PATH | --in-place]
trustify batch-format DATA_PATH [DATA_PATH ...] [--exclude PATTERN] [--apply] [--follow-symlinks]
trustify modernize [--apply]
trustify init-config [DIR] [--force]
trustify generate_markdown --out DIR
trustify generate_keywords [--out DIR]
trustify generate_pdf --out PATH.pdf [--from-markdown DIR] [--keep-build] [--doxygen PATH]
trustify cache list|clean|copy|path ...
trustify install-completion [--shell SHELL] [--dest PATH] [--print] [--force]
```

Global project-resolution flags (`--projects`, `--trust-root`,
`--schema`) come BEFORE the subcommand. `--schema` is rejected for
commands that don't consume a schema (`generate_schema`, `format`,
`batch-format`, `init-config`, `cache`).

### Modernize

```
trustify modernize [--apply]
```

Rewrites legacy XD-tag forms in the current top-level project's
`<project>/src/` tree.

**Scope** — narrower than schema-consuming commands because modernize
rewrites sources in place. It deliberately opts OUT of the four-signal
resolution (`--projects` > `.trustify.json` > baltik auto-detect > env
var) that `check` / `generate_schema` / etc. use. Rules live in
`cli/__init__.py::_scope_modernize`; first match wins:

1. `--projects A [B ...]` — those exact projects, NO transitive
   `[dependencies]` expansion, NO `--trust-root` overlay. An explicit
   multi-project warning prints to stderr at the end (dry-run and
   `--apply` both).
2. Auto-detected baltik — `project.cfg` walk-up from CWD, or
   `$project_directory` env. Modernize that baltik only; `--trust-root`
   / `$TRUST_ROOT` are honoured for schema purposes elsewhere but
   ignored here.
3. Only `--trust-root` / `$TRUST_ROOT` resolves (no baltik, no
   `$project_directory`) — modernize the TRUST tree itself.

Failure modes (both via `parser.error`, exit code 2):

- `.trustify.json` contributed `workspace.projects` /
  `workspace.trust_root`. `.trustify.json` is for user-environment
  management (decoupling a dataset folder from whatever env is
  active); driving in-place source rewrites from it is the wrong
  tool. Get inside a baltik or pass `--projects`.
- Nothing resolves at all — error lists the four options.

Anchor-inside-trust_root behaviour: when CWD lives under `$TRUST_ROOT`
and an implicit baltik signal (auto-detected `project.cfg` or
`$project_directory` env) is active, modernize targets that baltik and
emits a WARNING — it no longer suppresses the signal in favour of the
whole TRUST tree. This mirrors `_resolve_project_args` (the
schema-resolving path) and is what lets you modernize a baltik whose
sources live inside TRUST (e.g. ICoCo). The warning is sterner here
because modernize rewrites sources in place. Crucially, modernize's
scope is owned entirely by `_scope_modernize`: `modernize._resolve_src_dirs`
does NOT re-apply the `$project_directory` / `$TRUST_ROOT` env
fallbacks, so `trust_root=None` means "exclude TRUST" even when
`$TRUST_ROOT` is set — otherwise a single-baltik rewrite would silently
fold in all of TRUST.

Rewrite rules applied per file:

- Numeric brace flags (`-3`/`-2`/`-1`/`0`/`1`) on `XD` block headers
  → `BRACE`/`NO_BRACE`/`INHERITS_BRACE`.
- Numeric opt flags (`0`/`1`) on `XD attr` lines → `REQ`/`OPT`.
- XD lines (block headers, `XD attr`, `2`/`3`-prefixed variants)
  exceeding 120 chars split into base + `XD_CONT` continuation(s)
  at word boundaries via a greedy O(words) packer.
- `XD_ADD_P` (any level) canonicalised unconditionally to a
  metadata-only opener line (`<cpp_prefix>// <lvl>XD_ADD_P <type>`)
  plus a following `XD_CONT` line carrying the description. Empty
  description → bare `// XD_CONT` placeholder. The opener may exceed
  120 chars (C++ prefix is acceptably long); not counted as
  `unsplittable`.
- Multi-space whitespace between structural tokens inside any XD
  comment (`// XD   attr ...` → `// XD attr ...`) collapses to a
  single space. The C++ prefix (if any) is preserved verbatim — only
  the `// ...` payload is normalised. Applied as a side effect of
  rebuilding lines via `" ".join(...)` in rules 3 and 4.

Dry-run by default — emits unified diff on stdout, summary on
stderr, exits 1 if any changes are pending (CI-friendly, matches
`ruff format --check`). `--apply` writes the rewrite in place and
exits 0.

Idempotent: existing `XD_CONT` blocks fold first, then re-split
canonically. Running modernize on its own output is a no-op. The
splits counter reports "net new" lines vs the input block — an
already-correctly-split block contributes zero even when Rule 3
still touched it internally.

`generate_schema` separately emits a WARNING for every modernizable
pattern it sees — long lines (> 120 chars), legacy brace flags,
legacy opt flags — each suggesting `trustify modernize` as the fix.
The long-line warning is gated on `modernize.is_unsplittable(line)`:
lines whose structured prefix alone exceeds 120 chars (e.g. `XD attr`
with a huge `chaine(into=[...])` type) are skipped — modernize cannot
fix them, so the warning would be unactionable noise.

## Project structure

```
src/trustify/
  api.py                # public API: generate_schema, check, batch_check, ...
  cache.py              # ~/.cache/trustify/ management (hash + inspection)
  formatter.py          # .data formatter (shared with the LSP)
  projects.py           # .trustify.json discovery + project list resolution
  cli/                  # argparse dispatcher over trustify.api
  core/
    base.py             # runtime parser base classes (Abstract_Parser, ...)
    comment_scanner.py  # single source of truth for #...# and /*...*/ states
    hacks.py            # tiny post-load patches on the generated module
    misc_utilities.py   # import_parser_module + UnitUtils test mixin
    naming.py           # ToPydName / ToParserName / valid_variable_name
    schema_module.py    # SchemaModule wrapper around the imported parser mod
    source_location.py  # SourceLocation + build_projects_dict + relativize
    trad2_pydantic.py   # code generator (writes trustify_gen.py / _pyd.py)
    trad2_utilities.py  # XD scanner + TRAD2Content model
    trust_parser.py     # token-stream parser for .data files
  doc/
    generator.py        # markdown manual generator
    origin.py           # OriginClassifier (project name → label)
    extras.py           # \input / \includeimage handling for baltik docs
    render.py           # markdown rendering primitives
tests/                  # see "Test categories" below
  corpus/               # formatter golden files (inputs/ + expected/)
```

## Install + test

`TRUST_ROOT` must point at a checked-out TRUST source tree for any
test that needs the live schema (most of `test_api.py`,
`test_rw_full_datasets.py`, `test_doc_generator.py`, etc.). The
trustify env sets it; `make test` errors out without it.

### `make install-dev` — the development install

```bash
make install-dev
```

Runs `pip install -e ".[dev]"`. Two things happen:

- **Editable (`-e`)** — `src/trustify/` is symlinked into the active
  Python's site-packages. Source edits take effect on the next
  `import trustify` with NO reinstall step. Crucial for the
  edit-test-edit loop in agentic workflows; without `-e`, every code
  change requires another `pip install` before tests pick it up.
- **`[dev]` extra** — pulls in the full developer toolchain via the
  meta-extra declared in `pyproject.toml:27` (`dev =
  ["trustify[lint,util]"]`):
    - `lint` → `ruff` (mandatory for `make lint`).
    - `optional` → `argcomplete` (powers `trustify install-completion`;
      configurer_env best-effort installs this in the bundled python).
  Adding a new dev-only dependency means extending one of these
  sub-extras AND `dev` resolves it automatically.

`make install` (no `-dev`) is the PRODUCTION install path —
`pip install --force-reinstall --no-deps .` — used by the TRUST
build pipeline. It is NOT a substitute for `install-dev`: it copies
the files into site-packages (no editable link, source edits don't
take effect), and it omits the `[dev]` extra (no ruff → `make lint`
errors out).

### Verifying the install before development

Before touching trustify code, confirm the active Python actually
sees an editable+dev install. Two one-liners:

```bash
pip show trustify | grep -E "^(Editable project location|Location)"
ruff --version
```

The first line should print an `Editable project location:` field
pointing at this directory (`.../Outils/trustify`). If you only see
`Location: .../site-packages/...` without the `Editable` line, the
install is not editable — your edits will be invisible to the test
suite. The second line confirms the `[dev]` extra resolved (ruff is
present). If either fails, re-run `make install-dev` from
`Outils/trustify/` before doing anything else.

### Running tests

```bash
make test                 # python -m unittest discover -s tests -t . -v
                          # (+ a Jupyter smoke test unless
                          # $TRUST_DISABLE_JUPYTER=1)
make test_debug           # same with TRUSTIFY_DEBUG=1
make clean                # wipe egg-info, build/, tests/generated/

python -m unittest tests.test_api -v           # focus on one file
python -m unittest discover -s tests -t . -v   # full suite directly
```

## Test categories

| File | What it covers |
|---|---|
| `test_api.py` | Public API: `generate_schema`, `check`, `batch_check`, `load_dataset`, `init_config`. Needs `$TRUST_ROOT`. |
| `test_cache.py` | `compute_content_hash`, `cache_dir_for`, `is_cached`, the cache-inspection API (`iter_entries`, `resolve_id`, `copy_entry`). Pure unit tests, no `$TRUST_ROOT` required. |
| `test_cli.py` | End-to-end CLI tests via subprocess; covers every subcommand. `cache` tests use `HOME=<tmpdir>` to redirect the cache root. |
| `test_comment_scanner.py` | The shared `#...#` / `/*...*/` scanner (also used by the LSP). |
| `test_corpus.py` | Formatter golden-file regression — every pair under `tests/corpus/inputs/` + `tests/corpus/expected/` is parametrised. |
| `test_doc_extras.py` / `test_doc_generator.py` / `test_doc_origin.py` | The markdown-manual generator + the origin classifier. |
| `test_formatter.py` | Direct unit tests for the `formatter.format_dataset` etc. |
| `test_hacks.py` | `core.hacks.apply` — runtime patches applied to the generated module. |
| `test_naming.py` | Name-mangling helpers. |
| `test_parse_error_locations.py` | Source ranges + line/col reporting on parse errors. Uses the **"simple"** schema slot (`tests/trad2/TRAD_2_adr_simple`). |
| `test_projects.py` | `.trustify.json` discovery + project-list resolution. |
| `test_rw_elementary.py` | Round-trip parsing of small datasets against the **"simple"** schema slot. |
| `test_rw_full_datasets.py` | Round-trip parsing of full TRUST datasets against the **"full"** schema slot (generated live from `$TRUST_ROOT`). |
| `test_schema_module.py` | `SchemaModule` accessors. |
| `test_source_location.py` | `SourceLocation` dataclass + project helpers. |
| `test_source_ranges.py` | Token-stream source-range tracking. |
| `test_trad2_utilities.py` | XD scanner + `TRAD2Content` reader/writer. |

`tests/datasets/` holds reference `.data` files. Three of them
(`Canal_perio_VEF_2D.data`, `diffusion_implicite_jdd6.data`,
`distance_paroi_jdd1.data`) are **relative symlinks** into the live
`$TRUST_ROOT/tests/...` tree — do not edit them in place.

### Mock baltik fixtures

`tests/_baltik_fixtures.py` is the single source for building mock
baltiks on disk in tests. Use it instead of hand-rolling `project.cfg`
+ `src/markers.xd` writes — the helper guarantees the shape stays
identical to a real baltik (same `project.cfg` layout, same `// XD`
comment syntax, same `src/` convention).

```python
from tests._baltik_fixtures import mock_workspace

with mock_workspace() as ws:
    dep = ws.baltik("dep")
    main = ws.baltik(
        "main",
        deps={"dep": dep},
        xd_blocks=["// XD foo objet_lecture foo NO_BRACE Foo keyword."],
    )
    # main now has project.cfg + [dependencies] + src/markers_0.xd
```

The helper supports nested baltiks (pass `parent=`), `[dependencies]`
(mapping of name → path / relative string / `` `pwd` ``-style
expression), and XD content blocks (each becomes a separate `.xd` file
under `src/`). It does NOT belong in production code — strictly under
`tests/`.

## Gotchas

### Multi-schema coexistence — `_schema_module` activation

The single most common source of mysterious test failures.

Each generated `trustify_gen.py` runs a footer that pushes itself onto
every `Abstract_Parser` subclass via the class-level
`_schema_module` attribute. The base parser classes (`Dataset_Parser`,
`Bloc_lecture_Parser`, ...) live in `core/base.py` and are
**singletons in the process** — whichever schema loaded last "wins"
the back-ref until something re-activates the earlier one.

Symptom: a test that loads the "full" slot fails with errors that
clearly come from the "simple" slot (e.g. *"Unexpected attribute
'domain' in keyword 'read_med'"* with a `Model:` line pointing at
`TRAD_2_adr_simple`). This means the simple schema was activated last
and the full slot is now mis-resolving classes.

Mitigation already in place:

- The footer in the generated parser exposes
  `_activate_schema()` — re-callable, sets `_schema_module` on every
  Abstract_Parser subclass present in the module's globals.
- `core/misc_utilities.import_parser_module` always calls
  `_activate_schema()` before returning a module — including on
  cache-hit returns from `sys.modules`.
- `UnitUtils.generate_python_and_import` also calls
  `_activate_schema()` when the slot is already in `_TRUG[slot]`.

What can re-break this:

- Adding a NEW back-ref attribute to base parser classes that lives at
  the class level. Same hazard. Either: thread it through
  `_activate_schema`, store it on the SchemaModule wrapper instead, or
  make it per-instance.
- A code path that reaches into `_TRUG[slot]._module` directly without
  going through `generate_python_and_import` or `import_parser_module`.
  It must call `_activate_schema()` itself.
- Tests that assume a single global schema (e.g. directly setting
  `Dataset_Parser._schema_module` once and never updating). Don't do
  that.

### `_infoMain` / `_infoAttr` shape

After the source-locations refactor these are
`(project, path, line)` **triples** — not the old `[path, line]`
2-element lists. Consumers:

- `core/base.py::_render_genErr_message` reads `_infoMain[2]`
  (line) and `_infoMain[1]` (path).
- `doc/generator.py::_origin_of` reads `_infoMain[0]`
  (project name) and passes it to `OriginClassifier.classify_project`.
- The LSP's `schema_service.get_source_location` returns
  `(triple[1], triple[2])`.

The empty-info default is `()` (empty tuple), not `[]`.

### `trustify_gen_pyd` module name carries the content digest

The pyd module name is content-digest-suffixed:
`trustify_gen_pyd_<digest>.py` under the cache directory, and the
generated parser emits `from trustify_gen_pyd_<digest> import *`. The
digest is the same 16-char content hash that names the cache directory
(computed via `cache.compute_content_hash`), so two clones at different
absolute paths still produce byte-identical generated files.

This isolates concurrent schemas in `sys.modules`: each schema lands in
a distinct `trustify_gen_pyd_<digest>` slot, so an overlay baltik's
pydantic classes can never be shadowed by an earlier schema's cached
module.

When changing the pyd-emission filename convention, three call sites
must stay in lockstep:

- `api.generate_schema` computes the digest and writes the file.
- `cache.is_cached` globs `trustify_gen_pyd_*.py` (used by every cache
  inspection / lookup path).
- The hand-written test fixture path (`UnitUtils.generate_python_and_import`
  in `core/misc_utilities.py`) uses the slot name as the suffix; the
  parser's `from ... import *` is derived from `out_pyd_filename.stem`
  inside `trad2_pydantic.generate_pyd_and_pars`, so any new digest
  scheme must flow through that filename.

Pre-existing caches written by older trustify carry the un-digested
`trustify_gen_pyd.py`; `is_cached` (now glob-based on the digested
form) treats them as not-cached, so the next `generate_schema` writes
a fresh entry. Run `trustify cache clean` to evict the orphans.

### Pyd objects carry their schema affinity

Every pyd class generated by trustify inherits a `_schema_module`
ClassVar from its module's `TRUSTBaseModel`. The parser-module footer
sets it once at import time (`TRUSTBaseModel._schema_module =
_this_module`) — and crucially points it at the PARSER module, not
the pyd module, because only the parser module's namespace contains
both pyd and parser classes (via `from trustify.core.base import *`
+ `from trustify_gen_pyd_<digest> import *`). Two pyd modules loaded
in the same process keep distinct affinity (each has its own
TRUSTBaseModel class object).

Consumed by `Abstract_Parser._parser_cls_for(pyd_obj)`, which resolves
a parser class through `type(pyd_obj)._schema_module` instead of
through the shared parser-base `_schema_module` slot. Use the helper
(NOT `self.__class__._parser_from_pyd(pyd.__class__)`) in any
toDatasetTokens-time code that constructs a fresh parser for a pyd
value whose `_parser` may be `None`. Sites currently using the helper:

- `Dataset_Parser.toDatasetTokens` — per-entry `_read_type` read.
- `ListOfBase_Parser.getItemTokens` — fresh list element.
- `Read_Parser.toDatasetTokens` — fresh `Read.obj`.
- `BaseCommon_Parser._toDatasetTokens_braces` attr loop — fresh attr.

Parse-time class resolution (`Dimension_Parser`, `Declaration_Parser`,
`_init_itemParserType`, `GetAllClassesFromSyno`) deliberately stays
on `cls._schema_module` — `_activate_schema()` keeps that slot
correct during parsing, and parse paths don't survive a schema switch
the way serialization paths can.

When adding a NEW toDatasetTokens-time site that constructs a parser
for a pyd value, use `self._parser_cls_for(pyd_value)`. The
`self.__class__._parser_from_pyd(pyd_value.__class__)` form aliases
through the most-recently activated schema and silently mis-resolves
under concurrent schemas.

### `TRUSTIFY NOT` opt-out tag — inline justification

The opt-out tag accepts an optional inline justification after a `:`:

```
# TRUSTIFY NOT: this is a sed-substitution template, not a standalone dataset #
```

Surfaced everywhere trustify mentions the opt-out: live `batch-check`
progress, summary block, `--summary-out` file, the `--no-skip`
obsolete-marker WARNING. The bare `# TRUSTIFY NOT #` form keeps
working — `CheckResult.marker_justification` is empty for it.

Whitespace inside the inline text is collapsed
(`' '.join(text.split())`) so multi-space and accidental-newline
inputs normalise to single-space output. `# TRUSTIFY NOT: #` and `#
TRUSTIFY NOT:    #` are equivalent to the bare form.

When adding a new opt-out, prefer the inline form — the justification
shows up in CI summary output without needing to grep the file
itself. The existing two-line convention (bare sentinel followed by
a free-text TRUST comment) keeps parsing, but trustify does NOT
scrape the follow-up comment.

### `$project_directory` env leak (test contamination)

`trustify.projects.effective_projects(None)` falls back to
`$project_directory` if the env var is set (see `projects.py` —
"lowest-priority `projects` source"). The CLI, the API, and the test
suite all flow through this fallback. **In a shell where
`$project_directory` was set earlier** (e.g. a `cd` into a baltik or
a `source env_TrioCFD.sh`), running `make test` from
`Outils/trustify/` silently overlays that baltik onto every
generate_schema call.

Symptom: the test suite reports failures that look like missing
classes or unexpected attributes coming from a project the test code
doesn't mention. Recent occurrences:

- `test_ds_ftd` failing on `Unexpected attribute 'critere_remaillage'
  in keyword 'bloc_lecture_remaillage|nul'`
- `test_init_config_no_args_writes_full_scaffold` failing
- `test_doc_generator.TestDocSmokeFullTrust.test_generation_*` failing

All four pass in a clean shell (`unset project_directory`). Before
filing or debugging an obscure trustify test failure, **check the
env**:

```bash
env | grep -E 'project_directory|TRUST_ROOT'
unset project_directory && make test     # clean repro
```

The CLI is robust to this — passing `--projects` or `--trust-root`
explicitly overrides the env fallback — but the test infrastructure
relies on the default code path.

### `obj.info` mid-conversion shape

`TRAD2Attr.BuildFromTab` / `TRAD2Block.BuildFromTab` initialise
`obj.info = [absolute_path, lineno]` — a 2-element list. Only after
`TRAD2Content._convert_infos_to_source_locations` runs is each `info`
rewritten into a `SourceLocation`. That conversion is triggered by
`BuildFromOrgAndSources` when `self.projects` is non-empty (the
synthesis branch makes that automatic when only `trust_root=` is
passed). Code that builds `TRAD2Content` manually and goes straight to
`serialize_source_locations` must call the conversion itself.

### Schemas without `source_locations.json` — officially supported

`tests/trad2/TRAD_2_adr_simple` is hand-written and ships without a
sibling `source_locations.json`. More generally, any TRAD2 file loaded
via `BuildContentFromTRAD2(trad2, source_locations=None)` (the
default) is allowed to lack the sibling JSON — the loader auto-detects
when the file is present and silently falls through when it is absent.
`info` stays empty for blocks/attrs that never got a `SourceLocation`,
and `_info_to_triple` produces `("", abs_path, line)` triples (empty
project name) so consumers that render C++-source backrefs degrade
gracefully (no backref shown rather than a crash). Pass an explicit
path to `source_locations=` only when you have a JSON to load — that
form raises `FileNotFoundError` if the path doesn't exist.

### Cache auto-trim on regeneration

`api.generate_schema` deletes other cache entries whose
`provenance.projects` dict matches the new entry's exactly, every
time it writes a fresh entry to the default cache root. Rationale:
during active TRUST / baltik development the same project set
produces a new cache entry on every source change (line numbers in
`source_locations.json` shift), so without this the cache would
balloon. Same-project-set dedup carries most of the load — only the
newest entry per project set survives.

Opt-out knobs:

- `.trustify.json` → `{"cache": {"auto_trim": false}}` (validated by
  `projects.load_config`).
- `TRUSTIFY_NO_AUTO_TRIM` env var set to any non-empty value
  (one-off / per-shell escape; follows the NO_COLOR convention — the
  literal value is not parsed, so `=0` disables too).
- Passing `out=<dir>` to `api.generate_schema` (or `--out` on the
  CLI) skips auto-trim — the user-chosen output directory is
  considered user-controlled territory.

A manual lever is also exposed: `trustify cache trim` walks every
entry, groups by `provenance.projects`, and deletes all but the
newest per group. Use it after an `auto_trim: false` window.

### Hash determinism

Two consecutive `generate_schema` runs on the same source MUST
produce the same hash — otherwise the cache is worthless. The
synonym list inside `assemble()` used to leak `set()` iteration
order into the TRAD2 text, breaking this; it's now `sorted(syn)`.
If a future change adds another step that produces a list from a
set (or a dict iteration that depends on hash randomization), make
sure the canonical-text output is order-stable.

### Cache invalidation when the schema format changes

The hash takes (version, TRAD2 text, `source_locations.json` text). A
format change to any of those naturally invalidates pre-existing
cache entries — but if you change the **emission format** in
`trad2_pydantic.py` (e.g. add a new class attribute) without bumping
anything that goes into the hash, the cache will reuse old generated
files that are missing the new field. When in doubt, `rm -rf
~/.cache/trustify/*` (or run `trustify cache clean --force`) before
re-testing.

### LSP and trustify ship together

They live in the same monorepo and are versioned in lockstep. An
older trustify-lsp cannot read a cache produced by a newer trustify
(and vice versa). Cross-version compatibility isn't invested in —
when the runtime format changes, bump both.

### `generate_pdf` duplicates the standalone `pdf/` script

`trustify generate_pdf` orchestration lives in `src/trustify/pdf.py`
(tested). It is a deliberate duplicate of the standalone
`pdf/keyword_reference_to_pdf.py`, which is kept stdlib-only and
untested for users who cannot install trustify. The Doxyfile is
embedded in `pdf.py` as the `_DOXYFILE` string (not shipped as
package-data — that would tie bundling to setuptools), a twin of
`pdf/Doxyfile.pdf`. When editing either copy, mirror the other.
`generate_pdf` is in `_SCHEMA_CONSUMERS` (accepts `--schema`) but NOT
in `_NEEDS_PROJECT_RESOLUTION` — it resolves the workspace lazily in
its dispatch branch, only when `--from-markdown` is absent.

## Development rules

- **Before any development work on trustify, confirm the active
  Python has an editable + `[dev]` install.** From
  `Outils/trustify/`:

  ```bash
  pip show trustify | grep -E "^(Editable project location|Location)"
  ruff --version
  ```

  The first command MUST print an `Editable project location:` line
  pointing at this directory. If it only prints `Location:
  .../site-packages/trustify`, edits to `src/` will NOT be picked up
  by the test runner — every change would silently test against the
  stale copy in site-packages, and you'd waste a debugging cycle
  before noticing. The second command MUST succeed; `make lint`
  cannot run without ruff. If either check fails, run `make
  install-dev` (see "Install + test" above) before doing anything
  else.

- **Tests before every commit that touches code, implementation, or
  test files.** From `Outils/trustify/`:

  ```bash
  python -m unittest discover -s tests -t .
  ```

  All passing, except for tests legitimately skipped (Front Tracking
  BALTIK classes absent from core TRUST will skip `test_ds_ftd` and
  one other in `test_rw_full_datasets.py`).

  If LSP-side code is touched too, also run:

  ```bash
  cd ../trustify-lsp/server && python -m pytest -q
  ```

  Tests are NOT required for changes that cannot affect runtime
  behavior — `*.md`, comments, `.gitignore`, CI config.

- **Lint must pass before every commit that touches
  `Outils/trustify/src/` or `Outils/trustify/tests/`.** From
  `Outils/trustify/`:

  ```bash
  make lint
  ```

  This runs `ruff check src/ tests/` and `ruff format --check src/
  tests/`. Fix any reported issue (either by hand for the `check` rules
  or by running `ruff format src/ tests/` for formatting drift) — don't
  commit a red lint. Same scope rule as for tests: lint is NOT required
  when the diff is confined to `*.md`, comments, `.gitignore`, or CI
  config files.

- **`make check` whenever the change can affect dataset parsing,
  schema extraction, or schema code generation.** From
  `Outils/trustify/`:

  ```bash
  make check
  ```

  This runs `trustify batch-check $TRUST_ROOT/tests --jobs $(nproc)`
  (capped at 16 — `nproc` reports the full physical node on HPC login
  machines which are often 128+ cores; audit 6.7) against every
  `.data` file shipped with TRUST, surfacing any regression the unit
  tests would miss. Override the worker count via `make check
  TRUSTIFY_CHECK_JOBS=N`. Required whenever the diff touches:

    - the parser / formatter runtime
      (`src/trustify/core/{base,trust_parser,hacks}.py`,
      `src/trustify/formatter.py`),
    - the schema extraction or code-generation pipeline
      (`src/trustify/core/{trad2_utilities,trad2_pydantic}.py`,
      `src/trustify/core/source_location.py`, the LSP-shared
      comment scanner),
    - any `// XD ...` declaration in the TRUST C++ sources — including
      adding, removing, retyping, renaming, or changing the parent of
      an existing tag, or
    - a TRUST `.data` file under `$TRUST_ROOT/tests/`.

  Expected outcome: zero `FAILED`; the small `SKIPPED` set
  (`TRUSTIFY NOT` opt-outs) is allowed. If a previously-passing dataset
  starts failing, fix it before committing — don't grow the skip list
  to paper over a regression.

- **Don't claim a fix works until tests verify it.** Confirm the
  failures you see are exactly the ones you expected (e.g. only the
  new red-then-green test in a TDD cycle), not collateral damage. If
  a failure existed before your change, leave a note explaining why
  it's pre-existing.

- **Commit messages prefix with `[trustify]`** (or `[trustify+lsp]`
  when the same commit touches the LSP). Subject under ~72 chars;
  body explains the *why*, not just the *what*. Recent style — keep
  consistent:

  ```
  [trustify] cache list: record + show schema-generation timings

  api.generate_schema now times the three meaningful phases ...
  ```

- **Commit without asking for confirmation** when the user has
  approved the work and tests pass. Don't insert "are you sure?"
  prompts.

- **No emoji** in code, comments, or commit messages unless the user
  explicitly asks.

- **Comments only when the WHY is non-obvious.** No restating the
  code in prose. Avoid `TODO Task N:` style breadcrumbs that
  reference an in-flight refactor — file them as issues / spec doc
  notes instead, since they rot once the refactor lands.

- **Don't add backwards-compatibility shims** for code paths the
  refactor obsoletes. Cross-version interop between old/new trustify
  isn't a requirement; just bump both halves.

- **The CLI is fully documented in `README.md` and kept in sync with
  `src/trustify/cli/__init__.py`.** `README.md` is the canonical CLI
  reference: every subcommand has an entry in the "Subcommands at a
  glance" table and its own `### trustify <command>` section listing
  options and exit codes, and the global flags are described under
  "CLI API". Any change that touches the CLI surface — adding or
  removing a subcommand, adding/renaming/removing a flag, changing
  exit-code semantics, or altering which commands accept `--schema` —
  MUST update `README.md` in the same commit. Keep the abridged CLI
  surface block in this file (under "Use the schema") consistent too.
  A CLI change with no `README.md` diff is incomplete.

- **Per-project configuration** lives in `.trustify.json` (shared
  with the LSP). Filename is fixed. Three top-level sections, all
  optional:
    - `workspace.projects`, `workspace.trust_root`,
      `workspace.auto_resolve_dependencies` — project resolution.
    - `cache.auto_trim` — cache management.
    - `lsp.schema_dir`, `lsp.enum_dedup_threshold` — forwarded to the
      LSP server.
  Adding a new top-level section or key is fine (silently kept under
  `TrustifyConfig.raw` for forward compat); adding a new `lsp.*` key
  requires bumping the LSP reader too.

  The full four-signal precedence (CLI flag > `.trustify.json` >
  baltik auto-detection via `project.cfg` > env var), the CLI-vs-API
  split (auto-discovery in CLI, opt-in via `trustify.discover_config`
  in the API), and the multi-path discovery rule (per-arg, error on
  disagreement) live in `docs/project-resolution.md`. Update that doc
  in lockstep with `trustify/projects.py` whenever the resolution
  rules change.
