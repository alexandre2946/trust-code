@page Trustify trustify

`trustify` is a Python package that turns the TRUST CFD code's input
keyword declarations (`// XD` C++ comments in `$TRUST_ROOT/src/`, plus
any baltik overlay you point it at) into a `pydantic` data model and a
matching `.data` file parser. It's the reference Python implementation
of the TRUST dataset grammar.

The same generated schema powers three workflows:

- **Validation** — round-trip a `.data` file through the schema to
  catch typos, unknown keywords, type errors, and ordering issues.
- **Formatting** — re-emit a `.data` file with canonical indentation
  and whitespace, independent of the editor that produced it.
- **Scripting** — load any TRUST dataset as a `pydantic` model in
  Python, modify it, write it back. The same model also drives
  third-party tools (notably the LSP server in
  `Outils/trustify-lsp/`).

Each TRUST keyword becomes a Python class whose attributes are the
keyword's attributes. For example `read_med` is a class with
attributes `fichier`, `mesh`, etc.

This README covers the **CLI** — the principal entry point for
baltik users. The programmatic API and the `.trustify.json` config
each live on their own page; see *Further reading* below.

<h2>Further reading</h2>

- @subpage Trustify_ForTrustDevs — how to annotate C++ sources with `// XD` tags so trustify can pick them up.
- @subpage Trustify_ForTrustifyDevs — internals of trustify itself (TRAD_2 format, parsing/writing logic, error checking).
- @subpage Trustify_ProgrammaticAPI — everything importable from the top-level `trustify` package.
- @subpage Trustify_ConfigReference — every key of the per-workspace `.trustify.json` file.
- @subpage Trustify_ProjectResolution — full precedence ladder for `--projects` / `--trust-root` resolution.
- @subpage Trustify_KeywordsReferencePDF — the standalone, install-free script behind `trustify generate_pdf` (render the keyword reference as a single hyperlinked PDF).

## Installation

trustify ships with the TRUST conda environment. After sourcing it
(`source $TRUST_ROOT/env_TRUST.sh`), the `trustify` command is on
your `PATH` and `from trustify import ...` works.

For a manual install:


- If building TRUST with conda environment (default `./configure`):
```bash
pip install --force-reinstall --no-deps $TRUST_ROOT/Outils/trustify

```

- If building TRUST without conda (`./configure -without-conda`):
```bash
# with a venv:
pip install $TRUST_ROOT/Outils/trustify/
# with system python
pip install --user $TRUST_ROOT/Outils/trustify/
# with a conda managed environment that is not the TRUST one:
pip install --force-reinstall --no-deps $TRUST_ROOT/Outils/trustify/
```

- For developing trustify itself (recommend using the conda TRUST build):
```bash

cd $TRUST_ROOT/Outils/trustify
make install-dev    
# which does
pip install -e ".[dev]"  (editable + all dev dependencies, ruff and coverage)
```

### Dependencies

- `pydantic>=2,<3` is a hard runtime dependency, the only one. The TRUST conda env comes with it 
    (which is why `--no-deps` is used)

- `argcomplete` is an optional dependency for shell completion.
    `pip install $TRUST_ROOT/Outils/trustify[optional]` to get all optional deps (may increase in future).

- `ruff` and `coverage` are dev dependencies. Installed with
    `pip install $TRUST_ROOT/Outils/trustify[dev]`. 
    Prefer editable install for development (`-e`, see above).

---

## CLI API

### Invocation shape

```text
trustify [--projects DIR]... [--trust-root DIR] [--schema DIR] [--no-auto-deps] <command> [<args>]
```

**Project-resolution flags are GLOBAL** and must come **before** the
subcommand name. They tell trustify which TRUST source tree and which
baltik overlays to scan when building the schema:

| Flag | Effect |
|---|---|
| `--projects DIR`, `-p DIR` | Baltik directory to overlay. Repeat for multiple (`-p A -p B`). Each must contain a `project.cfg`. |
| `--trust-root DIR` | The TRUST source tree. Defaults to `$TRUST_ROOT`. |
| `--schema DIR` | Reuse an already-generated schema directory. Mutually exclusive with the two above. Only valid for `check`, `batch-check`, `batch-check-projects`, `generate_markdown`, `generate_keywords`, `generate_pdf`. |
| `--no-auto-deps` | Do not auto-resolve transitive `[dependencies]` from each project's `project.cfg`. Equivalent to setting `workspace.auto_resolve_dependencies: false` in `.trustify.json`. |

When **none** of these are passed, the CLI walks upward from the
current path looking for a `.trustify.json` and/or a `project.cfg`,
falling back to `$TRUST_ROOT` / `$project_directory`. Full precedence
ladder in @ref Trustify_ConfigReference "the .trustify.json reference"
and in @ref Trustify_ProjectResolution "the project-resolution reference".

### Subcommands at a glance

| Command | Purpose |
|---|---|
| `generate_schema` | Build the pydantic model + parser for the current project set; cache the result. |
| `check` | Validate one `.data` file against the schema (parse + round-trip). |
| `batch-check` | Run `check` on many files; emit an aggregate summary. |
| `batch-check-projects` | Run `batch-check` on a subdirectory (default `tests/Reference`) of each resolved project, filtered by `--only` role. |
| `format` | Reformat one `.data` file (pure text, no schema needed). |
| `batch-format` | Reformat many `.data` files in place. |
| `modernize` | Rewrite legacy XD-tag forms (brace flags, opt flags, long lines) in the resolved project source trees. |
| `init-config` | Freeze the currently-resolved workspace into a `.trustify.json`. |
| `projects` | Dump the resolved project list for shell scripting. |
| `generate_markdown` | Generate the keyword reference manual (markdown). |
| `generate_pdf` | Build the keyword reference manual as a single hyperlinked PDF (markdown → doxygen → LaTeX). |
| `generate_keywords` | Generate `Keywords.txt` and `Keywords.Vim` (vim/gedit syntax). |
| `cache` | Inspect and manage `~/.cache/trustify/`. |
| `install-completion` | Install shell tab-completion (needs the `[optional]` extra). |

`trustify <command> --help` always prints the full set of options.

### `trustify generate_schema`

```
trustify generate_schema [--out DIR]
```

Walks `--trust-root` and each `--projects` baltik for `// XD`
declarations, builds a `TRAD2Content`, then code-generates the parser
and pydantic model. Outputs five files in
`~/.cache/trustify/<hash>/` (or `--out DIR` if given):

- `TRAD2_trustify` — canonical flat-text schema.
- `source_locations.json` — `block_name → {project, path, line}` map.
- `trustify_gen_pyd_<digest>.py` — pydantic data model. The digest is
  the same content hash that names the cache directory, so multiple
  schemas loaded in one Python process keep distinct `sys.modules`
  entries.
- `trustify_gen.py` — parser module.
- `provenance.json` — generation timestamp + version + step timings.

The hash is computed from the trustify version + canonical TRAD2 text
+ `source_locations.json`. Re-running with the same inputs is a
no-op (cache hit). Two clones of the same project set at different
absolute paths share one cache entry (no absolute paths feed the
hash, and the digested pyd filename is content-derived too).

Prints the cache directory to stdout. Exits 0 on success, non-zero on
schema-extraction failure.

**Options**

- `--out DIR`, `-o DIR` — Write into `DIR` instead of the cache. When
  set, the auto-trim post-step is **skipped** (the user's directory is
  considered user-controlled territory).

### `trustify check`

```
trustify check [--no-skip] DATA_FILE
```

Parses `DATA_FILE` against the schema, then writes it back and asserts
the round-trip is byte-identical (modulo trailing content past the
`end`/`fin` keyword). Prints one line:

```
PASSED: /path/to/foo.data
FAILED: /path/to/foo.data — TrustifyParseError: ...
SKIPPED: /path/to/foo.data — dataset opted out via 'TRUSTIFY NOT' marker — sed template, not a real dataset
```

A dataset carrying a `# TRUSTIFY NOT #` comment is reported as
SKIPPED unless `--no-skip` is given. Use it to opt out files that
aren't valid TRUST input (templates, intentional counter-examples,
etc.). Two forms are accepted:

- `# TRUSTIFY NOT #` — bare sentinel.
- `# TRUSTIFY NOT: <justification text> #` — inline form; the
  justification is whitespace-collapsed and appended to every CLI
  line that mentions the opt-out (SKIPPED summary, live progress,
  `--summary-out` file, and the `--no-skip` obsolete-marker warning).

A malformed sentinel — `TRUSTIFY NOT` followed by tag text but no
colon, e.g. `# TRUSTIFY NOT toto #` — is reported as **FAILED** with
a clear "ill-formed TRUSTIFY NOT tag at line N" message rather than
silently ignored. A typo of the colon would otherwise let an invalid
dataset slip through the gate.

**Exit codes**: `0` on PASSED or SKIPPED, `1` on FAILED.

**Options**

- `--no-skip` — Run the check even on files carrying the `TRUSTIFY
  NOT` marker. If the parse passes anyway, the trailing line warns
  that the marker is now obsolete and can be removed.

### `trustify batch-check`

```
trustify batch-check [--exclude PATTERN]... [--follow-symlinks] [--no-skip] [--quiet] [--summary-out PATH] [--jobs N] DATA_PATH [DATA_PATH...]
```

Runs `check` on every `.data` file under each `DATA_PATH`. Each path
may be a file or a directory; directories are recursed.

By default emits one progress line per file:

```
[1/1986] PASSED /abs/path/to/file.data
[2/1986] FAILED /abs/path/to/other.data: TrustifyParseError: ...
[3/1986] SKIPPED /.../template.data: dataset opted out via 'TRUSTIFY NOT' marker — sed template, not a real dataset
...
Summary: 1981/1986 test(s) passed successfully (5 skipped, 0 failed)
  [SKIPPED] /.../template.data: dataset opted out via 'TRUSTIFY NOT' marker — sed template, not a real dataset
  [FAILED]  /.../broken.data: ...
```

Each `[SKIPPED]` line carries the trailing justification when the
opt-out tag uses the `# TRUSTIFY NOT: <text> #` inline form; the
bare `# TRUSTIFY NOT #` form prints just the static message. The
same suffix appears in the `--summary-out` file and on the
`[obsolete] ...` warning lines emitted under `--no-skip`.

**Exit codes**: `0` when zero failures, `1` otherwise. SKIPPED files
do not contribute to the failure count.

**Options**

- `--exclude PATTERN` — fnmatch glob to skip while recursing
  directories. Repeatable. A pattern without `/` matches basenames
  anywhere; with `/` it's anchored to the listed root; a trailing `/`
  restricts to directories. Files passed explicitly on the command
  line bypass `--exclude`.
- `--no-skip` — disable the `TRUSTIFY NOT` opt-out marker.
- `--quiet`, `-q` — suppress per-file progress lines. The final
  summary still prints.
- `--summary-out PATH` — also write the final summary (Summary
  line + SKIPPED/FAILED list + any obsolete-marker warnings) to
  `PATH`.
- `--jobs N`, `-j N` — parallelism (default `1`). `N>=2` farms the
  per-file checks out to a fork-based `multiprocessing.Pool`; workers
  inherit the already-loaded schema via fork+CoW (no re-import cost).
  Results stream back in submission order so the `[i/N]` lines stay
  correctly numbered. Ctrl-C aborts the pool cleanly.

### `trustify batch-check-projects`

```
trustify batch-check-projects [--only ROLE[,ROLE...]] [--subdir REL_PATH] [--exclude PATTERN]... [--follow-symlinks] [--no-skip] [--quiet] [--summary-out PATH] [--jobs N]
```

Convenience wrapper over `batch-check` for "check the datasets shipped
by the projects in my workspace". It resolves the same project set that
`trustify projects` dumps, appends `--subdir` to each project root, and
runs `batch-check` over the resulting directories. Equivalent to:

```
trustify projects --only=dependency | sed 's:$:/tests/Reference:' | xargs trustify batch-check
```

The directories selected are listed up front before the run begins.

**Options** (in addition to all `batch-check` options above):

- `--only ROLE[,ROLE...]` — keep only projects with one of these roles
  (`trust_root`, `primary`, `dependency`), exactly like
  `trustify projects --only`. With no `--only`, **every** resolved
  project is checked. Use `--only dependency` for just the transitive
  dependencies. (Role filtering is `--only` only — the `batch-check`
  `--exclude` here is still the file-glob filter, not a role filter.)
- `--subdir REL_PATH` — directory appended to each selected project root
  to locate datasets. Default `tests/Reference`; e.g.
  `--subdir tests/Validation`.

Projects that lack the `--subdir` directory are skipped with a one-line
note to stderr. Needs the workspace to resolve (run inside a baltik or
pass `--projects` / `--trust-root`); errors otherwise.

### `trustify format`

```
trustify format [--out PATH | --in-place] DATA_FILE
```

Reformats `DATA_FILE` using the canonical rules: 4-space indentation
proportional to brace depth, consecutive blank lines collapsed to one,
blank lines immediately after `{` or before `}` removed. **No schema
is needed** — this is pure text manipulation that runs on TRUST
syntax, baltik datasets, or templated `.data` files alike (the
formatter passes `${key}` / `$key` placeholders through unchanged).

By default prints the formatted text to stdout. Mutually-exclusive
output flags:

- `--out PATH`, `-o PATH` — write to `PATH`.
- `--in-place`, `-i` — rewrite the input file.

The formatter is intentionally **parameter-free** — every TRUST
dataset is reformatted the same way regardless of editor settings,
keeping files byte-identical across the ecosystem.

### `trustify batch-format`

```
trustify batch-format [--apply] [--exclude PATTERN]... [--follow-symlinks] DATA_PATH [DATA_PATH...]
```

Same path expansion + `--exclude` semantics as `batch-check`.

**Dry-run by default** — prints the unified diff for every file
that would change and exits 1 if there is at least one such file
(exit 0 when every input is already canonical). `--apply` switches
to in-place rewrite — every file is atomically replaced
(temp-file + `os.replace`, so a crash mid-write leaves the source
untouched). `TRUSTIFY_FORCE_FORMAT` set to any non-empty value in
the env has the same effect as `--apply` — useful from Makefile
targets that want a single env knob.

`--follow-symlinks` opts into following symbolic links during
directory recursion (both dir and file symlinks). Default skips
them — a symlinked file at a leaf emits a one-line skip warning
to stderr so out-of-scope links don't sneak in. Explicit symlink
args on the command line are always processed.

Mirrors `modernize`'s pattern: destructive operations are opt-in,
never the default. No schema is needed; no stdout / `--out` flags
(this is the fleet operation — `format` is the single-file
scalpel).

### `trustify modernize`

```
trustify modernize [--apply]
```

Rewrites legacy XD-tag forms in the **current top-level project's**
`<project>/src/**/*.{cpp,xd}` tree. Three rules apply per file, in
order:

1. **Brace flag**: numeric `-3/-2/-1/0/1` on `XD` block headers →
   `BRACE` / `NO_BRACE` / `INHERITS_BRACE`.
2. **Opt flag**: numeric `0/1` on `XD attr` lines → `REQ` / `OPT`.
3. **Long-line split**: any XD opener (incl. `XD_ADD_P`,
   `XD_ADD_DICO`, and the `2`/`3`-prefixed variants) exceeding 120
   chars is split into a base line + `XD_CONT` continuation(s) via
   a greedy O(words) word-packer. The opener tag stays on the first
   line; only the description tail moves into continuations. Lines
   whose structured prefix (C++ prefix + `// <opener> <structured
   tokens>`) alone exceeds 120 chars are reported as unsplittable
   and left as-is.

Idempotent: existing `XD_CONT` blocks fold first, then re-split
canonically. Running modernize on its own output is a no-op.

Dry-run by default — unified diff on stdout, one-line summary on
stderr (`summary: N files would change — A brace flags, B opt flags,
C splits, U unsplittable`). Exits 1 when changes are pending,
matching `ruff format --check` semantics so CI can gate on it.

**Scope.** Modernize rewrites sources in place, so it deliberately
opts OUT of the four-signal resolution other commands use. Rules
(first match wins):

- `--projects A [B ...]` — those exact projects, **no** transitive
  `[dependencies]` expansion, **no** `--trust-root` overlay. An
  end-of-run warning prints on stderr (dry-run and `--apply` both)
  to flag the multi-project scope explicitly.
- Otherwise the **auto-detected baltik** — `project.cfg` walk-up
  from CWD or `$project_directory` env. Modernize that baltik only;
  `--trust-root` / `$TRUST_ROOT` are ignored.
- Otherwise, when only `--trust-root` / `$TRUST_ROOT` resolves (no
  baltik, no `$project_directory`) — modernize the TRUST tree itself.

**Failure modes** (exit 2):

- `.trustify.json` contributing `workspace.projects` /
  `workspace.trust_root` — modernize refuses. `.trustify.json` is
  for user-environment management; driving in-place source rewrites
  from it is the wrong tool. Get inside a baltik or pass
  `--projects` explicitly.
- Nothing resolves at all — error spells out the four options.

When CWD lives under `$TRUST_ROOT`, the implicit baltik signals
(auto-detected `project.cfg`, `$project_directory`) are suppressed so
a stale env from a previous `source <baltik>/env.sh` cannot smuggle
a baltik rewrite into a plain modernize-from-trust.

**Options**

- `--apply` — write the rewrite in place. Summary on stderr, no diff
  on stdout, exits 0.

**Cross-link**: `generate_schema` separately emits a WARNING for each
modernizable pattern it sees (long lines, legacy brace flags, legacy
opt flags), each suggesting `trustify modernize` as the fix.

### `trustify init-config`

```
trustify init-config [--force] [DIR]
```

Resolves the workspace as usual (CLI flags > `.trustify.json` >
baltik auto-detection > env vars, with transitive `[dependencies]`
expanded), then writes the resolved `workspace.projects` +
`workspace.trust_root` (plus a placeholder `lsp` block) into a
fresh `.trustify.json` in `DIR` (default: current directory).

The result is a config frozen from the live env — think
`pip freeze > requirements.txt` for trustify workspaces. Future
invocations from inside that directory will use the captured
project set without any env setup.

**Options**

- `DIR` — target directory (default: cwd). The file is always named
  `.trustify.json`.
- `--force` — overwrite an existing `.trustify.json`. Without it,
  the command exits non-zero if the file already exists.

### `trustify projects`

```
trustify projects [--json] [--only ROLE[,ROLE...] | --exclude ROLE[,ROLE...]]
```

Dumps the resolved project list in **overlay order**: `trust_root`
first if any, then transitive dependencies in post-order, then
explicitly-listed primaries last. One absolute path per line —
designed for shell scripting:

```make
BALTIK_DEPENDENCIES := $(shell trustify projects --only=dependency 2>/dev/null)
```

Each entry carries a **role**:

- `trust_root` — the TRUST source tree.
- `primary` — explicitly listed by the user (CLI flag,
  `.trustify.json`, baltik auto-detection, or `$project_directory`).
- `dependency` — pulled in transitively via `project.cfg
  [dependencies]`.

**Exit codes**: `0` on success (including empty output after
filtering), `2` when nothing resolves at all.

**Options**

- `--json` — emit `[{path, role}, ...]` instead of one-path-per-line.
- `--only ROLE[,ROLE...]` — keep only matching roles. Mutually
  exclusive with `--exclude`.
- `--exclude ROLE[,ROLE...]` — drop matching roles.

### `trustify generate_markdown`

```
trustify generate_markdown --out DIR [--kw-ref-path PATH]
```

Renders the keyword reference manual (one markdown file per top-level
class, plus a master index) into `--out DIR`. The directory must not
exist OR be empty.

`generate_markdown` consumes a schema — pass `--schema DIR` to point
at one already generated, or omit and let it run `generate_schema`
internally. Extras (LaTeX-style `\input` directives, image embeds)
are pulled from each project's `docs/trustify/extras/` directory.

**Options**

- `--out DIR`, `-o DIR` — required output directory.
- `--kw-ref-path PATH` — override the master index location (default:
  `<out>/keyword_reference.md`).

**Exit codes**: `0` on success, `2` on precondition errors (bad
`--out`), `1` on content errors (missing referenced extras, duplicate
extras, broken image copy).

### `trustify generate_pdf`

```
trustify generate_pdf --out PATH.pdf [--from-markdown DIR] [--keep-build] [--doxygen PATH]
```

Builds a single hyperlinked PDF of the keyword reference manual. By
default it runs the whole pipeline end-to-end — schema → markdown →
doxygen → LaTeX → PDF — using throwaway temporary directories, so a
bare `trustify generate_pdf --out manual.pdf` is all you need. Pass
`--from-markdown DIR` to reuse an existing `generate_markdown` output
and skip schema resolution + generation.

Requires `doxygen` and a LaTeX toolchain (`pdflatex`, `make`) on
`PATH`; both are probed up-front so a missing install fails fast,
before the long doxygen pass. The orchestration lives in
`src/trustify/pdf.py`; see @ref Trustify_KeywordsReferencePDF for the
standalone, install-free variant of the same pipeline.

**Options**

- `--out PATH.pdf`, `-o PATH.pdf` — required PDF output path (a file).
  The parent directory must already exist.
- `--from-markdown DIR` — reuse an existing `generate_markdown`
  output instead of generating one; skips schema resolution.
- `--keep-build` — keep the scratch build directory (and the temp
  markdown directory in end-to-end mode) for inspection; their paths
  are printed to stderr.
- `--doxygen PATH` — explicit doxygen binary (overrides
  `$TRUST_DOXYGEN_BINARY` and `PATH`).

**Exit codes**: `0` on success, `2` on precondition errors (missing
`--out` parent, an invalid `--from-markdown` directory), `1` on
toolchain / build failures (doxygen or LaTeX errors).

### `trustify generate_keywords`

```
trustify generate_keywords [--out DIR]
```

Emits two files for downstream syntax-highlighting tooling:

- `Keywords.txt` — one keyword per line, prefixed with `|`. Source of
  truth for vim/gedit syntax generators.
- `Keywords.Vim` — a single `syntax keyword TRUSTLanguageKeywords  ...`
  line, keywords double-space-separated.

The keyword set includes every block name + its synonyms, every
attribute name + its synonyms, and every `chaine(into=[...])` enum
value. Abstract `_base` / `_deriv` classes are dropped, as are tokens
containing `{}#*` (which would break downstream syntax files).

**Options**

- `--out DIR`, `-o DIR` — output directory (default: cwd).

### `trustify cache`

```
trustify cache <action> [<args>]
```

Inspect and manage the schema cache at `~/.cache/trustify/` (or
`$TRUSTIFY_CACHE_DIR/` when that env var is set to a non-empty
path — useful on HPC nodes with a quota'd / read-only `$HOME`,
ephemeral CI runners, and multi-user shared TRUST installs that
want one shared cache). Entries are 16-char hash directories; any
unique **prefix** of a hash works as an ID (git-short-SHA style).

#### `cache list`

```
trustify cache list [--json]
```

Table of entries with creation time, trustify version, project set,
and on-disk size. `--json` emits the same data as machine-readable
JSON.

#### `cache clean`

```
trustify cache clean [--force] [ID...]
```

- With **no IDs**: wipes every entry. Prompts for confirmation;
  `--force` skips the prompt.
- With **explicit IDs**: deletes only those entries, no prompt.

Ambiguous prefixes are rejected with the full candidate list.

#### `cache copy`

```
trustify cache copy ID DEST_DIR
```

Copies one entry's generated files into `DEST_DIR`. The destination
is created if absent; if it exists it must be empty. `__pycache__`
and `*.pyc` are filtered out.

#### `cache path`

```
trustify cache path ID
```

Prints the absolute path of one entry. Useful for
`--schema $(trustify cache path <id>)` plumbing.

#### `cache trim`

```
trustify cache trim
```

Groups entries by project set and deletes all but the newest in each
group. The same logic runs **automatically** after each
`generate_schema` against the default cache root (unless
`cache.auto_trim: false` in `.trustify.json` or
`TRUSTIFY_NO_AUTO_TRIM` set to any non-empty value in the env opts
out — the literal value is not parsed, so `=0` also disables).

### `trustify install-completion`

```
trustify install-completion [--shell {bash,zsh,fish}] [--dest PATH] [--print] [--force]
```

Installs a shell completion script generated via `argcomplete`.
Requires the `[optional]` extra (`pip install trustify[optional]`). In
a TRUST-bundled python, `configurer_env` already attempted this for you
on a best-effort basis (warns rather than fails if your install host
had no internet).

Defaults pick the right per-shell lazy-load destination:
- bash → `$XDG_DATA_HOME/bash-completion/completions/trustify`
- fish → `$XDG_CONFIG_HOME/fish/completions/trustify.fish`
- zsh → `$XDG_DATA_HOME/zsh/site-functions/_trustify` (you must add
  the parent directory to `$fpath` before `compinit`; the command
  prints the exact incantation).

**Options**

- `--shell {bash,zsh,fish}` — target shell. Default: detected from
  `$SHELL` (fallback bash).
- `--dest PATH` — override the destination file path.
- `--print` — dump the script to stdout instead of writing it.
- `--force`, `-f` — overwrite an existing destination file.

---

## Miscellaneous

### Pydantic version

trustify uses [`pydantic`](https://docs.pydantic.dev/latest) v2 for the
generated data model.

### TRUST daily checks

trustify is wired into the daily TRUST regression suite via
`make check` (which runs `trustify batch-check $TRUST_ROOT/tests
--jobs $(nproc)`). The check catches schema/dataset drift — if a
keyword's documentation in the C++ source disagrees with how the
dataset actually uses it, the round-trip fails. Same logic mirrored
in each baltik via the generated `trustify_check` Make target.

### Known limitations

- **Repeated keywords**: datasets repeating the same keyword/option
  within a block are tolerated by TRUST itself but not by trustify.
  Only the last instance is kept on read; on write-back, only that
  single instance reappears. Avoid `dt_max 1.0 dt_max 2.0` style
  patterns.

### Running the tests

```bash
make install-dev    # editable install + ruff
make test           # full unit suite (+ Jupyter smoke unless TRUST_DISABLE_JUPYTER=1)
make lint           # ruff check + ruff format --check
make check          # batch-check every .data under $TRUST_ROOT/tests
```

The hello-world notebook lives at
[`docs/examples/hello_world/HelloWorld.ipynb`](docs/examples/hello_world/HelloWorld.ipynb) —
useful as a starting point for scripting against the generated model.

### Documentation pipeline

To update the readthedocs documentation site, maintainers:
- run `make doc` in TRUST to generate the per-keyword markdown
  files,
- clone `github.com:cea-trust-platform/trust-documentation.git` and
  run its `update_from_TRUST.sh` script to import the freshly
  generated content.
