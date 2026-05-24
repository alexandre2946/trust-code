# Filter for doxygen stdout - collapses noisy per-item progress lines
# (5000+ "Patching output file N/M", 3000+ "Parsing file X", ...) into
# one informational line each, while preserving high-level phase
# markers, the cache stats summary, and the final "finished..." line.
#
# Usage (in docs/common.mk):
#   doxygen > >( tee build/doxygen.log | awk -f filter_doxygen.awk )
#
# Full output (unfiltered) is still in build/doxygen.log.

function once(label, key) {
    if (!seen[key]) {
        print label
        seen[key] = 1
    }
}

/^Patching output file [0-9]/                  { once("Patching output files (silent) ...", "patch");      next }
/^Running dot for graph [0-9]/                 { once("Running dot for graphs (silent) ...", "dot");       next }
/^Parsing file /                               { once("Parsing files (silent) ...", "parse");              next }
/^Preprocessing /                              { once("Preprocessing files (silent) ...", "preproc");      next }
/^Reading \//                                  { once("Reading files (silent) ...", "read");               next }
/^Generating code for file /                   { once("Generating code for files (silent) ...", "code");   next }
/^Generating docs for /                        { once("Generating docs for symbols (silent) ...", "docs"); next }
/^Generating dependency graph for directory /  { once("Generating dependency graphs (silent) ...", "depgraph"); next }
/^Searching for files in directory /           { once("Searching directories (silent) ...", "search");     next }

# Pass everything else through (phase markers, cache stats, finished, ...).
{ print }
