# TRUST documentation build system
#
# Builds HTML documentation from Doxygen-annotated C++ sources and
# Markdown pages, including the auto-generated keyword reference manual
# and the trustutils Python API reference.
#
# All targets must be invoked from this directory (docs/).
#
# Prerequisites:
#   - python3 with the `trustify` package on PATH (the `trustify_doc` target
#     pip-installs it from $(TRUST_ROOT)/Outils/trustify if missing — both a
#     system Python and the conda env that ./configure sets up work).
#   - graphviz (optional — required only for class/call graphs)
#
# No ./configure step is required: TRUST_ROOT defaults to the repo root,
# and doxygen is unpacked from externalpackages/ so all users get the same
# version.  Set TRUST_DOXYGEN_BINARY=/path/to/doxygen to use a different
# binary (e.g. when the bundled Linux binary is not portable to your OS).
#
# Shared rules (DOXYGEN resolution, install_doxygen, install_trustify,
# clean_doxygen, the doxygen recipe + post-processing) live in common.mk
# and are also pulled in by the baltik docs Makefile template.

# Resolve TRUST_ROOT once (ifndef avoids re-running the shell on every
# $(TRUST_ROOT) reference, which `?=` would).
#
# Derive from THIS Makefile's path (`$(lastword $(MAKEFILE_LIST))` is
# this very file, since the line runs before any `include`), NOT from
# CWD.  CWD-based `cd .. && pwd` works for the documented invocations
# (`cd docs && make`, `make -C docs`) but silently picks the WRONG
# parent for `make -f docs/Makefile` from the repo root — and there's no
# error, just cascading "file not found" later.  Audit M1.
ifndef TRUST_ROOT
TRUST_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
endif
export TRUST_ROOT

# -----------------------------------------------------------------------
# Layout
#
# common.mk derives HTML_DIR, KW_DST_DIR, GENERATED_CONTENT_DIR from
# BUILD_DIR (set below) - we only spell out trust-specific paths.
# -----------------------------------------------------------------------
BUILD_DIR  = $(CURDIR)/build

# -----------------------------------------------------------------------
# trustify (keyword reference)
#
# `trustify generate_markdown` writes kw_*.md, keyword_reference.md, and
# a figures/ subdir directly into the --out directory.  No staging step.
# -----------------------------------------------------------------------
KW_REF_MD  = $(GENERATED_CONTENT_DIR)/keyword_reference.md

# -----------------------------------------------------------------------
# trustutils (validation form API)
#
# trustutils_doc invokes a sub-make in $(TRUSTUTILS_DOC_DIR) that writes
# intermediate Markdown into $(TU_GEN_DIR) (owned by that sub-make, not
# wiped by `make clean_build`).  copy_trustutils brings the final
# pages into $(GENERATED_CONTENT_DIR) so doxygen finds them under
# build/generated_content/ alongside the keyword reference.
# -----------------------------------------------------------------------
TRUSTUTILS_DOC_DIR = $(TRUST_ROOT)/Validation/Outils/documentation
TU_GEN_DIR         = $(TRUSTUTILS_DOC_DIR)/build/md
TU_DST_DIR         = $(GENERATED_CONTENT_DIR)/validation_api
TU_REF_MD          = $(GENERATED_CONTENT_DIR)/validation_form_api.md

# -----------------------------------------------------------------------
# Help preamble - consumed by the `help` target in common.mk
# -----------------------------------------------------------------------
define HELP_PREAMBLE
Usage: make [target]   (run from docs/)

Combo targets:
  all          Regenerate everything and build HTML (default)
  kw           Regenerate keyword reference, then build HTML
  api          Regenerate trustutils API docs, then build HTML
  run_doxygen  Build HTML from existing sources (no regeneration)

For more info: $(CURDIR)/README.md
               $(CURDIR)/How_to_write_doc.md
endef

# `all` is the default goal. We set it explicitly because `include
# common.mk` (below) introduces targets, and Make would otherwise pick
# the first one in the include as the default.
.DEFAULT_GOAL := all

include $(TRUST_ROOT)/docs/common.mk

# -----------------------------------------------------------------------
# Phony declarations
# -----------------------------------------------------------------------
.PHONY: all kw api \
        clean \
        trustify_doc \
        trustutils_doc copy_trustutils

# =======================================================================
# High-level combo targets
#
# Each combo runs in two phases via $(MAKE):
#   1. regenerate the relevant markdown content
#   2. run_doxygen (which scans the now-stable content)
# Splitting into sub-make invocations creates an explicit phase barrier
# that Make's plain prerequisite list (`all: a b`) does NOT provide
# under `make -j` - the two would otherwise be scheduled in parallel.
#
# None of the combos cleans first.  Doxygen overwrites its outputs on
# every run, and our post-processors (flatten_navtree.py, patch_filter.py)
# operate on those freshly-overwritten files, so re-runs are correct AND
# much faster than a clean build (typically ~10x on a populated tree).
# Run `make clean all` when you genuinely need a fresh start - e.g.
# after removing a class/keyword/page whose stale HTML you want to
# evict from the output tree.
# =======================================================================

all:
	@$(MAKE) --no-print-directory trustify_doc copy_trustutils
	@$(MAKE) --no-print-directory run_doxygen

kw:
	@$(MAKE) --no-print-directory trustify_doc
	@$(MAKE) --no-print-directory run_doxygen

api:
	@$(MAKE) --no-print-directory copy_trustutils
	@$(MAKE) --no-print-directory run_doxygen

# =======================================================================
# Low-level targets (clean)
# =======================================================================
#
# Everything we generate lives under $(BUILD_DIR), so `clean` is just an
# alias for the shared `clean_build`.  Finer-grained `clean_html` and
# `clean_kw_ref` come from common.mk for fast iteration.

clean: clean_build                            ## Remove all generated artefacts (alias for clean_build)

# =======================================================================
# Low-level targets (trustify — keyword reference)
# =======================================================================

trustify_doc: install_trustify                ## Generate keyword reference docs via trustify CLI
	@echo "[doc] Generating keyword reference documentation ..."
	@rm -rf $(KW_DST_DIR) $(KW_REF_MD)
	@mkdir -p $(GENERATED_CONTENT_DIR)
	@$(TRUSTIFY) --trust-root $(TRUST_ROOT) generate_markdown \
	    --out $(KW_DST_DIR) \
	    --kw-ref-path $(KW_REF_MD)
	@echo "[doc] Keyword reference docs written to $(KW_DST_DIR)/"

# =======================================================================
# Low-level targets (trustutils — validation form API)
# =======================================================================

trustutils_doc:                               ## Generate trustutils Markdown API docs
	@echo "[doc] Generating trustutils API documentation ..."
	@$(MAKE) -C $(TRUSTUTILS_DOC_DIR) md

copy_trustutils: trustutils_doc               ## Copy trustutils docs into the doxygen tree
	@echo "[doc] Copying trustutils docs into the doxygen source tree ..."
	@test -d "$(TU_GEN_DIR)" || { \
	    echo "ERROR: $(TU_GEN_DIR) does not exist."; \
	    echo "       Run 'make trustutils_doc' first."; \
	    exit 1; }
	@ls $(TU_GEN_DIR)/trustutils_*.md > /dev/null 2>&1 || { \
	    echo "ERROR: no trustutils_*.md found in $(TU_GEN_DIR)" >&2; \
	    echo "       (trustutils_doc produced no module pages - check $(TRUSTUTILS_DOC_DIR))" >&2; \
	    exit 1; }
	@rm -rf $(TU_DST_DIR)
	@mkdir -p $(TU_DST_DIR)
	@cp $(TU_GEN_DIR)/trustutils_*.md      $(TU_DST_DIR)/
	@cp $(TU_GEN_DIR)/validation_form_api.md $(TU_REF_MD)
	@echo "[doc] Trustutils copy complete."

# =======================================================================
# `help` is provided by common.mk (driven by the HELP_PREAMBLE above).
