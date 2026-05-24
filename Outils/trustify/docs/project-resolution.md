@page Trustify_ProjectResolution Project resolution in trustify

This doc describes how trustify decides **which TRUST sources to scan**
(`--trust-root`) and **which baltik projects to overlay on top**
(`--projects`) when you run a command. The same logic applies whether
you call the CLI (`trustify <cmd>`) or the programmatic API
(`from trustify import generate_schema, ...`), but the CLI does some
extra work behind the scenes that the API doesn't — see
[CLI vs programmatic API](#trustify_cli_vs_api).

## The four signals

For each of the two project-list fields (`projects`, `trust_root`),
trustify consults up to **four** sources, in strictly decreasing
priority:

| # | Source | When it kicks in |
|---|---|---|
| 1 | **CLI flag** (`--projects` / `--trust-root`) | always wins when present |
| 2 | **`.trustify.json`** in an ancestor of the current path | nearest one wins; CLI auto-discovers, API on opt-in |
| 3 | **Baltik auto-detection** via `project.cfg` in an ancestor | only for `--projects`; CLI auto-detects, API on opt-in |
| 4 | **Env var** (`$TRUST_ROOT` / `$project_directory`) | always consulted last |

Crucially: **each field is resolved independently**. A `.trustify.json`
that declares only `projects` does NOT prevent `$TRUST_ROOT` from
filling in `trust_root`. This is the "per-field strict-replace" rule —
once a higher-priority source sets a field, no lower-priority source
can change it; but absent fields cascade down to the next layer.

### Worked example

You source the TRUST env (`$TRUST_ROOT` and `$project_directory` are
set), `cd` into your dataset folder where you have a `.trustify.json`
declaring `projects: ["/abs/baltik-A", "/abs/baltik-B"]` but **no**
`trust_root`. You run `trustify check foo.data`:

- `projects` resolution: CLI flag absent → config declares
  `["/abs/baltik-A", "/abs/baltik-B"]` → that's the answer. Baltik
  detection and `$project_directory` are skipped.
- `trust_root` resolution: CLI flag absent → config has no
  `trust_root` → no auto-detection mechanism for `trust_root` →
  `$TRUST_ROOT` env wins.

Effective schema scan: `$TRUST_ROOT` + baltik-A + baltik-B.

## The `.trustify.json` file

The config file is a JSON object with three optional top-level
sections — `workspace` (project resolution), `cache`, and `lsp`. Every
field within each section is optional. Newer keys may be added later;
unknown keys are preserved verbatim under `raw` for the LSP and future
consumers.

```json
{
  "workspace": {
    "projects": ["/abs/path/to/baltik-A", "/abs/path/to/baltik-B"],
    "trust_root": "/abs/path/to/trust",
    "auto_resolve_dependencies": true
  },
  "cache": {
    "auto_trim": true
  },
  "lsp": {
    "schema_dir": "/abs/path/to/pre-generated-schema",
    "enum_dedup_threshold": 4
  }
}
```

- `workspace.projects`, `workspace.trust_root`: as described above.
- `workspace.auto_resolve_dependencies` (default `true`): when set to
  `false`, trustify uses only the projects you list explicitly without
  transitively pulling in each project's `[dependencies]` from
  `project.cfg`. Equivalent to the `--no-auto-deps` CLI flag.
- `cache.auto_trim` (default `true`): when set to `false`,
  `generate_schema` will NOT evict superseded same-project-set cache
  entries after writing a new one. Useful when you want to keep
  intermediate snapshots around for comparison. Equivalent to setting
  `TRUSTIFY_NO_AUTO_TRIM` to any non-empty value in the environment.
- `lsp.*`: forwarded as-is to the trustify-lsp server. The CLI doesn't
  read these.

### Where to put it

Two intended use cases drive two placement patterns:

1. **Power user / multi-baltik / datasets anywhere.** You work with
   `.data` files in a folder that isn't inside any baltik (e.g.
   `~/work/datasets-X/`). Drop a `.trustify.json` there listing the
   permanent baltik installs your datasets use:
   ```json
   {"workspace": {"projects": ["/opt/baltiks/A", "/opt/baltiks/B"], "trust_root": "/opt/trust-1.x"}}
   ```
   Now `trustify check` / `batch-check` / etc. find this config when
   walking up from your `.data` paths and use the listed installs —
   regardless of whatever `$TRUST_ROOT` / `$project_directory` you may
   have inherited from sourcing some other env.

2. **Developer inside a baltik (or TRUST itself).** A `project.cfg`
   in the baltik root already identifies the baltik as the active
   project — you don't need to declare `workspace.projects` in
   `.trustify.json`. You might still want a `.trustify.json` for LSP
   customization (`enum_dedup_threshold`, ...) or to opt out of
   `cache.auto_trim`. Keep the file minimal:
   ```json
   {"lsp": {"enum_dedup_threshold": 6}}
   ```

   **Hard rule:** a `.trustify.json` sitting next to a `project.cfg`
   MUST NOT declare `workspace.projects` or `workspace.trust_root`.
   The project is already auto-detected from the baltik marker;
   redeclaring it here is ambiguous and trustify rejects the
   combination with a `ConfigError` pointing at the offending file.
   Strip those keys, or move the `.trustify.json` outside the baltik
   tree if you intended it to override the auto-detected project.

### Gitignore convention

Workspace `.trustify.json` files are usually personal: they encode
your local install paths, not project policy. The TRUST repo's root
`.gitignore` ignores `.trustify.json` so accidental copies don't get
committed. The repo deliberately tracks **one** `.trustify.json` at
the very root, used to ship the LSP defaults that apply to anyone
working in core TRUST. Baltiks are encouraged to do the same.

## Baltik auto-detection

Walking up from the relevant path, trustify treats the first directory
that contains a `project.cfg` file as the active baltik. This is how
`trustify check $YOUR_BALTIK/datasets/foo.data` knows to add your
baltik to `--projects` without you having to set anything up.

Every directory passed as `--projects` (or listed under `projects` in a
`.trustify.json`, or auto-detected) **must** be a proper baltik —
i.e. it must contain a top-level `project.cfg` with `[description].name`
set. `bin/trust baltik` and `baltik_configure` create this file when
they scaffold a baltik. Pointing `--projects` at a bare directory now
fails fast with a `ConfigError` rather than silently scanning it as
an unnamed project: run `baltik_configure` inside the directory to fix
it, or pass it via `--trust-root` if you actually meant to point at a
TRUST source tree.

The baltik detector and the `$project_directory` env var are
**redundant on purpose**. Detection is the more reliable path (works
even if you didn't source the baltik's env), but the env var is kept
as a backup for users who run trustify from outside any baltik
checkout while still having `$project_directory` set from an earlier
shell session.

### Implicit-resolution summary on stderr

Whenever any source other than a CLI flag contributed to the resolved
workspace (auto-discovered `.trustify.json`, baltik auto-detection,
`$TRUST_ROOT` / `$project_directory` env, or transitive
`[dependencies]` expansion), the CLI prints a one-shot summary to
stderr right after resolution:

```
trustify: resolved workspace:
  trust_root  /abs/path/to/trust       ($TRUST_ROOT env)
  project     /abs/path/to/main_baltik (auto-detected (project.cfg walk-up))
  dependency  /abs/path/to/dep_a       (transitive [dependencies])
  dependency  /abs/path/to/dep_b       (transitive [dependencies])
```

The block is silent when every entry came from a CLI flag, and for the
`projects` subcommand (which is itself a stdout dump of the same data).

### Anchor-inside-`$TRUST_ROOT` warning

When every path arg (or, in the no-path case, the CWD) sits inside the
resolved `$TRUST_ROOT` and an **implicit** baltik signal (auto-detected
`project.cfg`, or the `$project_directory` env-var fallback) supplied
the overlay, trustify keeps the overlay but prints a one-line WARNING
to stderr:

```
trustify: WARNING: dataset anchor is inside $TRUST_ROOT (/abs/path/to/trust) but a
baltik overlay is active (/abs/path/to/baltik, via $project_directory env). Proceeding
with the overlay. If this is a stale environment, pass --trust-root alone (or unset
$project_directory) to check against plain TRUST.
```

This used to be a hard suppression — both implicit signals were
silently dropped so a stale `$project_directory` left over from a
previous `source <baltik>/env.sh` could not overlay the wrong baltik.
That also blocked the legitimate case: a baltik whose sources live
*inside* the TRUST tree (e.g. ICoCo) is a perfectly valid overlay
target when you want to check whether the baltik's XD tags break TRUST
datasets. The call is now left to the user — the overlay applies and
the warning makes an accidental stale env visible.

The explicit signals (`--projects` flag and a
`.trustify.json workspace.projects` list) never trigger the warning —
the user clearly asked for the overlay.

## Transitive dependency resolution

Once the (projects, trust_root) pair is decided, trustify expands each
project's `[dependencies]` section transitively — same model as
`baltik_configure`. A baltik whose `project.cfg` declares:

```ini
[dependencies]
solver_kit : ./vendored_kit
shared     : `pwd`/../shared
```

contributes its own keywords *plus* every keyword from `solver_kit`
and `shared`, recursively. Paths support the same conventions as
baltik scripts:

- `` `pwd` `` substitutes the declaring baltik's directory.
- `$VAR` / `${VAR}` / `~` expand normally.
- Relative paths resolve against the declaring baltik's root.

The resolved list is **post-order**: each dependency appears before
the baltik that declares it, so the inner baltik shadows its deps in
`scanSourceFiles`. Diamond dependencies are deduplicated; each unique
resolved path is scanned exactly once.

Two rigid checks (mirroring `baltik_configure`):

1. Every declared dep path must exist and be a directory.
2. The declared `dep_name` must match the dep's own
   `[description].name`. Mismatch raises `ConfigError` with both names
   in the message.

### Freezing the resolved workspace

`trustify init-config <DIR>` runs the full resolution chain
(CLI flags → `.trustify.json` → baltik auto-detection → env vars → transitive
dependencies) and writes the result into `<DIR>/.trustify.json`. The
intended workflow:

```bash
# Inside your baltik (or with $TRUST_ROOT/$project_directory set), in
# the directory where your datasets live:
trustify init-config .
```

After this, that directory carries a self-contained config: future
trustify invocations from there work without any env setup, against
the exact same project list that was active when you ran init-config.
Same idea as `pip freeze > requirements.txt` for a Python env.

Because resolution validates every path, init-config refuses to write
an "aspirational" config pointing at directories that don't exist —
fix the projects, then re-run.

### Opting out

Two ways to disable transitive expansion (the resolved list is then
exactly what the four-signal step produced, with no further
expansion):

- CLI: `trustify --no-auto-deps <command> ...` (global flag).
- Config: `{"workspace": {"auto_resolve_dependencies": false}}` in
  `.trustify.json`. The CLI flag wins when both are set.

The CLI flag is intended for one-off invocations; the config switch
fits workflows where you've pinned every baltik path explicitly in
`.trustify.json` and don't want trustify second-guessing.

## CLI vs programmatic API {#trustify_cli_vs_api}

| | CLI (`trustify <cmd>`) | Programmatic API (`from trustify import ...`) |
|---|---|---|
| `.trustify.json` discovery | **automatic**, walks up from path args (or CWD if none) | **opt-in** — caller invokes `discover_config(path)` and threads result through |
| Baltik auto-detection (`project.cfg`) | **automatic** | **opt-in** — caller invokes `detect_baltik_root(path)` |
| `$TRUST_ROOT` / `$project_directory` fallback | yes | yes (via `effective_*` helpers, used internally) |
| Multi-path discovery anchor | per path arg, see below | N/A (single call site) |

The library API is **side-effect-free by design**: no implicit
filesystem walks, no surprises. If you want CLI-style auto-discovery
in a script, do it explicitly:

```python
from pathlib import Path
import trustify

cfg = trustify.discover_config(Path("~/datasets/projectX/foo.data").expanduser())
kwargs = cfg.to_api_kwargs() if cfg else {}

# Env-var fallback still works inside generate_schema, so you don't
# need to pass anything if the user has $TRUST_ROOT set.
schema_dir = trustify.generate_schema(**kwargs)
```

The returned `TrustifyConfig` also carries the cache and LSP settings,
useful when the LSP needs to know e.g. `cfg.lsp.enum_dedup_threshold`.

## Multi-path commands

`trustify batch-check`, `trustify batch-format`, and any other command
that takes multiple paths use **per-path discovery**: each path arg is
treated as an independent anchor for the upward walk. If they
disagree — e.g. two paths under different `.trustify.json` files or
under different baltiks — the CLI errors out with a clear message
rather than silently picking one and running with mismatched context.

If you genuinely need to run trustify across two workspaces, do it as
two separate invocations.

## Inspecting the resolved list from scripts

`trustify projects` dumps the full resolved workspace as a script-friendly
list of paths — useful when another tool (Doxygen, a build system, a
CI driver) needs the same dependency chain trustify computed.

```bash
trustify projects                       # one path per line, overlay order
trustify projects --json                # [{path, role}, ...]
trustify projects --only=dependency     # filter by role
trustify projects --exclude=trust_root,primary
```

Roles:

- **`trust_root`** — the TRUST source tree.
- **`primary`** — explicitly listed by the user (CLI flag, config,
  baltik auto-detection, or `$project_directory`).
- **`dependency`** — pulled in transitively from `project.cfg`
  `[dependencies]`.

Errors (exit 2) when nothing can be resolved at all — matches every
other schema-needing command. Filtering down to an empty list (e.g.
`--only=dependency` on a baltik with no deps) is **not** an error; it
prints nothing and exits 0.

The canonical use case is the baltik scaffold's `docs/Makefile`, which
asks trustify for the list of baltiks this one depends on, so
Doxygen's `INPUT` can be extended to include their C++ sources:

```make
BALTIK_DEPENDENCIES ?= $(shell trustify --trust-root $(TRUST_ROOT) -p $(project_directory) projects --only=dependency 2>/dev/null)
```

## Cheat sheet

```text
projects:    CLI flag  >  .trustify.json[projects]  >  project.cfg auto-detect  >  $project_directory
trust_root:  CLI flag  >  .trustify.json[trust_root]                            >  $TRUST_ROOT
each project gets its [dependencies] expanded transitively, unless
  --no-auto-deps  or  .trustify.json {"dependencies": {"auto_resolve": false}}
```

```text
CLI: auto-discovers .trustify.json + project.cfg from path args (or CWD).
     Multi-path = per-path discovery, error on disagreement.
     Auto-resolves transitive [dependencies] (off with --no-auto-deps).
API: env-var fallback only. Opt into discovery via trustify.discover_config(),
     opt into dependency expansion via trustify.expand_dependencies().
```

## `trustify modernize` is an exception

The four-signal resolver above applies to every command that *reads*
sources for schema purposes — `check`, `batch-check`,
`generate_schema`, `generate_markdown`, `generate_keywords`,
`init-config`, `projects`.

`trustify modernize` is the one command that *writes* sources back in
place, so its scope is narrower and bypasses the resolver entirely
(rules live in `cli/__init__.py::_scope_modernize`):

- `--projects A [B ...]` — those exact projects, **no** transitive
  `[dependencies]` expansion, **no** `--trust-root` overlay. An
  end-of-run warning prints on stderr (dry-run and `--apply` both).
- Otherwise the auto-detected baltik (`project.cfg` walk-up from CWD
  or `$project_directory` env) — that one baltik only.
- Otherwise only `--trust-root` / `$TRUST_ROOT` resolves (no baltik,
  no `$project_directory`) — modernize the TRUST tree itself.
- `.trustify.json` contributing `workspace.projects` /
  `workspace.trust_root` — **error**. `.trustify.json` captures a
  user environment (decoupling a dataset folder from whatever env is
  active); driving in-place source rewrites from it is the wrong
  tool. Run modernize from inside a baltik or pass `--projects`
  explicitly.
- Nothing resolves at all — error.

The anchor-inside-`$TRUST_ROOT` warning applies here too: when CWD sits
inside `$TRUST_ROOT` and an implicit baltik signal (auto-detected
`project.cfg` or `$project_directory`) is active, modernize targets
that baltik and prints a WARNING — it no longer suppresses the signal
in favour of the whole TRUST tree. This is what lets you modernize a
baltik whose sources live *inside* TRUST (e.g. ICoCo). The warning is
sterner than the schema resolver's because modernize rewrites `.cpp` /
`.xd` sources in place, so a stale `$project_directory` from a previous
`source <baltik>/env.sh` could rewrite the wrong tree — but, per the
"leave it to the user" policy, it warns rather than silently dropping
the overlay.

Note modernize's scope is decided entirely here, in the CLI:
`api.modernize` / `modernize._resolve_src_dirs` do NOT re-apply the
`$project_directory` / `$TRUST_ROOT` env fallbacks. Passing
`trust_root=None` means "do not include TRUST", even when `$TRUST_ROOT`
is set in the shell — otherwise modernizing a single baltik would fold
in all of TRUST.

## Appendix — `baltik_configure` parity

trustify reads `project.cfg` files via Python's `configparser`, while
`baltik_configure` parses them through a hand-rolled `sed` / `awk` /
`eval` pipeline. The two have meaningfully different default
semantics — trustify configures `configparser` to match `baltik_configure`
on the cases that occur in practice. This appendix is the reference for
anyone debugging a "baltik accepts this, trustify doesn't" (or vice
versa) report.

### What `baltik_configure` actually does

The parsing primitives live in
`bin/baltik/share/baltik/bin/baltik_configuration_parsing`:

```bash
extract_section() {
    cat "$1" |
      sed 's/#.*//' |                  # strip # to EOL on every line
        sed -n '/^ *$/!p' |            # drop empty / whitespace-only lines
          sed 's/ *[:=][ \t]*/:/' |    # normalize " = " / " : " to ":"
            tr '[]' '\n ' |
              sed 's/ *$//' |
                sed -n "/^$2$/,/^$/ p" |  # from header to next empty line
                  grep -v "$2" |
                    sed -n '/^ *$/!p'
}

extract_field() {
    echo "$1" |
      sed -n "/$2/p" |                 # substring match on field name
        cut -d : -f 2                  # only the 2nd `:`-split column
}
```

Dependency-graph traversal lives in
`bin/baltik/share/baltik/bin/baltik_dependencies_management`
(`add_dependencies` recurses, `check_dependencies` errors on
same-name-different-path).

The behaviour that matters for trustify parity:

| `baltik_configure` behaviour | Mechanism |
|---|---|
| Strip inline `#` comments on every line | `sed 's/#.*//'` |
| `=` and `:` are interchangeable separators | `sed 's/ *[:=][ \t]*/:/' ` |
| Section ends at the **first empty line** (not next `[header]`) | `sed -n "/^$2$/,/^$/ p"` |
| `[DEFAULT]` is just a regular (ignored) section | sed treats it like any other header |
| Field-name lookup is **case-sensitive substring** match | `sed -n "/$2/p"` |
| Field value is the second `:`-split column | `cut -d : -f 2` |
| Continuation lines (` indented`) are not interpreted | each line stands alone |
| Dependency value resolved via shell `eval "echo …"` | shell power: `$(...)`, backticks, quoting, concat |
| Dependency path resolved relative to declaring baltik's CWD | `cd ${dependency_orig_path}` before the `cd $declared_path` |
| Two distinct paths with the same `[description].name` is an error | `check_dependencies` → `too_many_paths_error` |
| Declared dep name must match the dep's own `[description].name` | `invalid_dependency_name_error` |
| Unknown `[description]` fields rejected | whitelist: `name\|author\|executable\|kernel\|cpp_flags\|ld_flags` |

### What trustify mirrors

Implemented in `src/trustify/projects.py` (`_read_baltik_config`,
`_resolve_dep_path`, `expand_dependencies`):

- Inline `#` comments stripped via `inline_comment_prefixes=("#",)`.
- `[DEFAULT]` treated as a regular section via
  `default_section=_UNUSED_DEFAULT_SECTION` — no inheritance magic.
- Option keys stay case-sensitive (`cp.optionxform = str`), so
  `Name = Foo` is **not** recognised as `name` (matches baltik).
- Dependency keys keep their declared case (`TrioCFD : ...` is
  preserved verbatim).
- Surrounding `"…"` / `'…'` quotes on dependency paths are stripped —
  matches what `eval "echo $entry"` does naturally.
- Same `[description].name` on two distinct paths raises `ConfigError`
  (mirrors `too_many_paths_error`).
- Declared dep name must equal the dep's own `[description].name`
  (mirrors `invalid_dependency_name_error`).
- `=` vs `:` separators both supported (configparser default).
- Dep path resolution: relative paths resolve against the declaring
  baltik's directory (same as baltik); `` `pwd` ``, `$VAR`, `${VAR}`,
  and `~` are expanded.

### What trustify diverges from intentionally

These corner cases are not chased bug-for-bug; either nobody writes
them in practice, or trustify is strictly better.

- **Empty line inside a section.** baltik truncates the section at the
  first empty line; configparser reads through to the next `[header]`.
  A `[dependencies]` with a blank line in the middle gives trustify
  more entries than baltik. There is no clean configparser knob for
  this and the case has never been seen in real configs.
- **Continuation lines.** A line beginning with whitespace is
  concatenated onto the previous value by configparser. baltik would
  treat it as an unrelated (unmatched) line. Realistic `project.cfg`
  files never use continuation lines for paths.
- **Shell substitution beyond the supported subset.** baltik's
  `eval "echo $entry"` runs arbitrary shell — `$(cmd)`, non-`pwd`
  backticks, command substitution, string concatenation like
  `"foo"/"bar"`. trustify supports only `` `pwd` ``, `$VAR`, `${VAR}`,
  `~`, and surrounding quote stripping. See `_resolve_dep_path`'s
  docstring. No known real-world `project.cfg` uses anything outside
  this subset.
- **Extra `[description]` fields.** baltik rejects field names outside
  `name|author|executable|kernel|cpp_flags|ld_flags`; trustify only
  reads `name` and silently ignores the rest. trustify is more
  lenient — won't break a valid baltik config.
- **`extract_field`'s substring bug.** baltik's `sed -n "/name/p"`
  matches any line containing the substring `name` (`nickname`,
  `username`, …), and `cut -d : -f 2` keeps only the second
  `:`-separated column. trustify uses an exact-key lookup and the full
  value. trustify is stricter and more correct; nobody should be
  relying on the baltik quirk.

### How to extend this in the future

If a real-world `project.cfg` ever surfaces that baltik accepts but
trustify rejects (or vice versa):

1. Reproduce in `tests/test_projects.py` — under `TestProjectNameFor`
   for `[description].name` issues, under `TestExpandDependencies` for
   `[dependencies]` issues. Each existing test cites the baltik shell
   construct it mirrors.
2. The parser entry points are `project_name_for` and
   `_parse_dependencies` (both via `_read_baltik_config`), and
   `_resolve_dep_path` for path semantics.
3. Update the tables above so the next debugging session starts from
   an accurate baseline.
