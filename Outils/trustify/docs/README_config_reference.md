@page Trustify_ConfigReference trustify — .trustify.json reference

`.trustify.json` is the per-workspace config. It's auto-discovered
by the CLI by walking upward from the dataset path (or CWD when no
path is given). The programmatic API does NOT auto-discover —
callers thread it in explicitly via `discover_config()`.

The file is a JSON object with **three optional top-level sections**.
Every field within each section is optional. Unknown keys are
preserved verbatim under `TrustifyConfig.raw` for forward-compatibility.

## Full scaffold

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
    "...": "documented separately when the LSP server lands"
  }
}
```

## `workspace.*` — project resolution

### `workspace.projects` (list of strings)

Absolute (or `$VAR` / `~`-expandable) paths of baltik directories to
overlay on top of `trust_root`. Each entry must contain a
`project.cfg`. Order matters — later entries override earlier ones
on basename collisions in `src/`.

When this key is set, transitive `[dependencies]` from each
project's `project.cfg` are expanded automatically (unless
`auto_resolve_dependencies` is false). You don't need to list
dependencies explicitly.

A `.trustify.json` placed **next to** a `project.cfg` (i.e. inside a
baltik) **cannot** declare `workspace.projects` or
`workspace.trust_root` — that combination is rejected as ambiguous.
Use the baltik's auto-detection instead.

### `workspace.trust_root` (string)

Absolute path of the TRUST source tree. When unset, `$TRUST_ROOT`
from the environment is used.

### `workspace.auto_resolve_dependencies` (bool, default `true`)

When `false`, trustify uses only the projects you list explicitly
without transitively pulling in each project's `[dependencies]`.
Equivalent to the `--no-auto-deps` CLI flag.

## `cache.*` — schema cache

### `cache.auto_trim` (bool, default `true`)

When `true` (the default), `generate_schema` deletes any sibling
cache entries that came from the **exact same project set** every
time it writes a fresh entry to the default cache root
(`~/.cache/trustify/`). Rationale: during active TRUST/baltik
development, the same project set produces a new cache entry on
every source change (line numbers in `source_locations.json`
shift), so without this the cache would balloon.

Set to `false` to keep intermediate snapshots around for
comparison. Equivalent to setting `TRUSTIFY_NO_AUTO_TRIM` to any
non-empty value in the env (the literal value is not parsed —
follows the NO_COLOR convention). `trustify cache trim` is the
manual lever when auto-trim is off.

Auto-trim only fires against the default cache root: passing
`--out DIR` to `generate_schema` skips it.

## `lsp.*` — LSP server config

Forwarded as-is to the [trustify-lsp](../trustify-lsp/) server.
trustify itself parses these keys (so the LSP doesn't need its own
discovery logic) but doesn't consume them. The full set of LSP
keys is documented on the LSP doc side.

## Four-signal precedence

For each of `projects` and `trust_root`, trustify consults up to
**four sources**, in strictly decreasing priority:

| # | Source | Notes |
|---|---|---|
| 1 | CLI flag (`--projects` / `--trust-root`) | always wins when present |
| 2 | `.trustify.json` in an ancestor of the current path | nearest one wins; CLI auto-discovers, API on opt-in via `discover_config` |
| 3 | Baltik auto-detection (`project.cfg` in an ancestor) | only for `projects`; CLI auto-detects, API on opt-in via `detect_baltik_root` |
| 4 | Env var (`$TRUST_ROOT` / `$project_directory`) | always consulted last |

Each field resolves **independently**: a `.trustify.json` that
declares only `projects` does NOT block `$TRUST_ROOT` from filling
in `trust_root`.

**Anchor-inside-`$TRUST_ROOT` warning.** When every path arg (or CWD,
when none is given) sits inside the resolved `$TRUST_ROOT` and an
implicit baltik signal — auto-detection (#3) or the
`$project_directory` env var (#4) — supplied the overlay, trustify
keeps the overlay but prints a WARNING to stderr. This used to be a
hard suppression that silently dropped both implicit signals; it now
only warns, because overlaying a baltik whose sources live inside the
TRUST tree (e.g. ICoCo) is a legitimate way to test whether the
baltik's XD tags break TRUST datasets. `--projects` and
`.trustify.json workspace.projects` are explicit, so they never warn.

Library callers that want NO baltik overlay (and no `$project_directory`
fallback) pass `projects=[]` to `generate_schema` / `load_dataset` /
`check` / `batch_check`. The empty list is the canonical "no overlay"
answer; `None` (the default) keeps the env-var fallback active.

**`trustify modernize` exception.** Modernize rewrites sources in
place and refuses to be steered by `.trustify.json` — if a discovered
config contributes `workspace.projects` or `workspace.trust_root`,
modernize errors out and asks the user to run from inside a baltik
(`project.cfg` ancestry) or pass `--projects` explicitly. The rule
is "in-place source rewrites belong in a developer environment, not
a user-environment config." See @ref Trustify_ProjectResolution for
the full modernize-scope rules.

Full discussion + worked examples in
@ref Trustify_ProjectResolution "the project-resolution reference".

## Freezing the current workspace

```bash
trustify init-config .
```

Writes a `.trustify.json` containing whatever projects + trust_root
trustify just resolved (transitive `[dependencies]` included, paths
absolute). Future trustify invocations from that directory work
without any env setup, against the exact same project list.

## Environment variables (summary)

| Variable | Effect |
|---|---|
| `TRUST_ROOT` | Fallback `trust_root` when not set by CLI/config. |
| `project_directory` | Fallback first entry of `projects` when not set by CLI/config/baltik-detection. |
| `TRUSTIFY_CACHE_DIR` | Redirect the schema cache root away from `~/.cache/trustify/`. Useful on HPC sites with a quota'd / read-only `$HOME`, ephemeral CI runners, and multi-user shared TRUST installs that want a shared cache. `~` is expanded; whitespace-only values fall back to the default. |
| `TRUSTIFY_NO_AUTO_TRIM` | When set to any non-empty value, equivalent to `cache.auto_trim: false`. Follows the NO_COLOR convention — the value is not parsed (so `=0` still disables; empty string is treated as unset). |
| `TRUSTIFY_FORCE_FORMAT` | When set to any non-empty value, `trustify batch-format` writes in place instead of running in its default dry-run mode (same as passing `--apply`). Lets Makefile targets enable the destructive path with one env knob. Same NO_COLOR-style semantic. |
| `TRUSTIFY_MP_START_METHOD` | Override the `multiprocessing` start method `batch-check --jobs N` uses (`fork`, `forkserver`, or `spawn`). Default picks `forkserver` if available, `spawn` otherwise — never `fork` silently, because fork-after-pydantic-import can deadlock on C-extension internal locks. Set this to `fork` to opt back into the older behaviour, or to `spawn` to debug forkserver issues. |
| `TRUSTIFY_DEBUG` | When `=1`, enables verbose debug logging from the parser/generator. Very noisy — for development of trustify itself. |
| `TRUST_DISABLE_JUPYTER` | When `=1`, `make test` skips the Jupyter smoke notebook. |
