"""trustify.core — internal modules: TRAD2 parsing, pydantic generation,
parser/streaming primitives, runtime hacks, low-level utilities.

These modules are not part of the public API. Use `trustify.api` (or the
top-level `from trustify import generate_schema, check, batch_check, load_dataset`)
for stable programmatic access.

Top-level shim modules (`trustify/<modname>.py`) re-export from this
subpackage for backward compatibility with pre-Phase-2 imports; those
shims will be removed in Phase 3 once external callers (`trustify-lsp`,
`bin/trust`, baltik integration) migrate to the public API.
"""
