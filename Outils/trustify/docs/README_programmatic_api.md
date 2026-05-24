@page Trustify_ProgrammaticAPI trustify — Programmatic API

Everything below is importable directly from `trustify`:

```python
from trustify import (
    # Schema + dataset operations
    generate_schema, load_dataset, check, batch_check, init_config,
    # Formatter
    format_dataset, format_dataset_file,
    # Workspace resolution
    discover_config, detect_baltik_root, expand_dependencies,
    # Types
    TrustifyConfig, ConfigError, TrustifyParseError, SourceRange,
    # Low-level parser primitives
    Dataset_Parser, TRUSTParser, TRUSTStream,
)
```

`trustify.api.*` and `trustify.projects.*` are the same symbols, but
**the only stable import path is the top-level `trustify` package**.
Sub-modules (`trustify.core.*`, `trustify.cli.*`, `trustify.formatter`,
`trustify.doc.*`) are implementation details and can change without
notice.

## Schema and dataset operations

### `generate_schema(projects=None, trust_root=None, out=None) -> Path`

Generate (or reuse-cached) the TRUST schema for the given project
set. Returns the directory containing the five generated files
(`TRAD2_trustify`, `source_locations.json`, `trustify_gen_pyd.py`,
`trustify_gen.py`, `provenance.json`).

- `projects` — list of baltik paths. `None` falls back to
  `$project_directory`. Pass `[]` to explicitly request "no baltik
  overlay" with **no** env-var fallback — useful when you've already
  resolved the workspace yourself and don't want trustify to second-
  guess by reading `$project_directory`.
- `trust_root` — TRUST source tree. `None` falls back to
  `$TRUST_ROOT`.
- `out` — override the cache directory. When set, auto-trim is
  skipped.

Mirrors the `trustify generate_schema` CLI exactly.

### `load_dataset(filename, projects=None, trust_root=None, schema=None) -> Any`

Load a TRUST `.data` file and return it as a pydantic model rooted at
the `Dataset` class of the appropriate schema. Two modes:

- `schema=<directory>` — reuse an already-generated schema.
- `projects=` / `trust_root=` — generate (or reuse-cached) the
  schema, then load. Falls back to env vars when both are `None`.

Returns the populated `Dataset` instance; raises
`TrustifyParseError` on syntax / structural errors.

### `check(data_file, projects=None, trust_root=None, schema=None, no_skip=False) -> CheckResult`

Parse `data_file` against the schema and verify the byte-identical
round-trip. Returns a `CheckResult` (dataclass) with:

- `path: str` — the input path.
- `status: "PASSED" | "FAILED" | "SKIPPED"`.
- `message: str` — error or skip reason; empty when PASSED.
- `obsolete_marker: bool` — True when `no_skip=True` was passed AND
  the dataset carries a `TRUSTIFY NOT` marker AND the parse passed.
  Signals the marker can be removed.

Does NOT raise on parse errors — they're captured in the returned
result. This lets `batch_check` keep iterating after a single
failure.

### `batch_check(data_files, projects=None, trust_root=None, schema=None, no_skip=False, verbose=False, jobs=1) -> BatchCheckResult`

Run `check` on every entry in `data_files`. The schema is generated
(or pulled from cache) ONCE before iteration.

Returns a `BatchCheckResult` with `total`, `passed`, `failed`,
`skipped`, and `results: list[CheckResult]`.

- `verbose=True` — print one `[i/N] STATUS path[: message]` line per
  finished dataset to stdout (default: silent; the CLI flips this
  on).
- `jobs >= 2` — process-level parallelism via a fork-based
  `multiprocessing.Pool`. Results stream back in submission order.

### `init_config(directory='.', projects=None, trust_root=None, overwrite=False) -> Path`

Write a `.trustify.json` template into `directory`. Raises
`FileExistsError` if the file already exists and `overwrite=False`.
Returns the path to the written file.

Emits the `workspace` block (frozen from the passed `projects` /
`trust_root`) plus a placeholder `lsp` block. `cache` keys are not
written — they default to their declared values
(see @ref Trustify_ConfigReference "the .trustify.json reference")
and can be added by hand when needed.

## Formatter

### `format_dataset(text, start_context=None) -> str`

Format TRUST `.data` text. Pure text manipulation — no schema
needed. Pass any `str`: file content, a string built in memory, the
result of a `string.Template` substitution.

`start_context` is an advanced hook for the LSP server (lets it
format a partial document mid-edit). End users can ignore it.

### `format_dataset_file(filename) -> str`

Thin convenience wrapper that reads `filename` and returns the
formatted text.

## Workspace resolution

These are **opt-in** for programmatic callers. The API itself never
walks the filesystem looking for configs or baltik markers — call
these explicitly when you want CLI-style discovery in a script.

### `discover_config(start=None) -> TrustifyConfig | None`

Walk up from `start` (file or directory; defaults to CWD) looking
for a `.trustify.json` and return its parsed `TrustifyConfig`.
Returns `None` when no config is found. Raises `ConfigError` on
malformed content.

```python
from pathlib import Path
import trustify

cfg = trustify.discover_config(Path("~/datasets/foo.data").expanduser())
kwargs = cfg.to_api_kwargs() if cfg else {}
schema_dir = trustify.generate_schema(**kwargs)
```

### `detect_baltik_root(start=None) -> Path | None`

Walk up from `start` looking for a directory containing a
`project.cfg`. Returns the absolute path of that directory, or
`None`.

### `expand_dependencies(projects) -> list[str] | None`

Recursively expand each project's `[dependencies]` from its
`project.cfg`. Returns the full overlay-correct list in post-order
(every dependency appears before the baltik that declares it).
Validates that each declared dep exists and that its
`[description].name` matches what the parent declared.

`None` / empty input returns unchanged.

## Types and exceptions

### `TrustifyConfig`

Frozen dataclass returned by `discover_config` / `TrustifyConfig.from_file`.
Fields:

- `source_path: Path` — the file it was loaded from.
- `workspace: WorkspaceConfig` — `projects`, `trust_root`,
  `auto_resolve_dependencies`.
- `cache: CacheConfig` — `auto_trim`.
- `lsp: LspConfig` — opaque to trustify; consumed by the LSP server.
- `raw: dict` — the original JSON payload, exposed for forward-compat
  with keys not yet modelled here.

Helper: `to_api_kwargs() -> dict` returns
`{"projects": ..., "trust_root": ...}` ready to splat into
`generate_schema(**kwargs)` / `load_dataset(**kwargs)` / etc.

See @ref Trustify_ConfigReference "the .trustify.json reference" for
the on-disk schema.

### `ConfigError`

Raised when a `.trustify.json` is malformed (wrong type for a known
key, conflicting `project.cfg` + `workspace.projects`, missing
dependency, etc.). Subclass of `Exception`.

### `TrustifyParseError`

Raised by the schema parser on syntactic / structural errors in a
`.data` file. Carries `file_name`, `line`, `col`, `end_line`,
`end_col`, `token`, `attr_name`, `kind` for callers that want to
report editor-friendly diagnostics. `check()` catches this and
folds it into the returned `CheckResult.message`.

### `SourceRange`

Frozen 4-tuple `(start_line, start_char, end_line, end_char)`
attached to parsed nodes by `TRUSTStream`. Used by the LSP to map
pydantic model fields back to source positions.

## Low-level parser primitives

These are the entry points the generated schema uses internally.
Most callers don't need them — `load_dataset` and `check` are the
right level of abstraction. Documented here because they're part of
the public surface.

### `Dataset_Parser`

Root parser class. Every generated schema's `Dataset` class
subclasses this. `Dataset_Parser.ReadFromTokens(stream)` is the
entry-point parse call.

### `TRUSTParser`

Tokenizer for TRUST `.data` text. `parser.tokenize(text)` populates
its token lists; pair with `TRUSTStream` for parser-friendly
consumption.

### `TRUSTStream`

Wraps a `TRUSTParser`'s token list with the bookkeeping the schema
parsers expect (current position, source-range tracking, error
context).
