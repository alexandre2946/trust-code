# TRUST documentation build — shared rules
#
# Included by:
#   $(TRUST_ROOT)/docs/Makefile                            (TRUST itself)
#   $(TRUST_ROOT)/bin/baltik/templates/basic/docs/Makefile (every baltik)
#
# Contract — the includer must set these BEFORE the `include`:
#   TRUST_ROOT      path to the TRUST checkout
#   BUILD_DIR       root for all generated artefacts (typically docs/build/)
#
# A `Doxyfile` must exist in the includer's $(CURDIR).
#
# The includer can then either:
#   * use `run_doxygen` as their build target directly, OR
#   * extend it with extra prereqs, e.g. `run_doxygen: my_extra_prereq`.

# Strict CWD: the includer's Makefile must be invoked from its own
# directory.  Several paths resolve relative to CWD — the includer's
# `BUILD_DIR = $(CURDIR)/build` default, the `Doxyfile` read by the
# doxygen recipe, the Doxyfile's own `OUTPUT_DIRECTORY = build` /
# `@INCLUDE_PATH = .` — so running from elsewhere silently writes to
# (and `make clean` would `rm -rf`) the wrong tree. pre-fix,
# `make -f docs/Makefile clean` from the repo root would have wiped
# $TRUST_ROOT/build, the C++ build directory.
#
# `firstword $(MAKEFILE_LIST)` is the includer's Makefile (the file
# passed to make); its directory is where we require make to be running.
# `realpath` on both sides normalises symlinks so a symlinked checkout
# is not falsely rejected.
_includer_dir := $(realpath $(dir $(firstword $(MAKEFILE_LIST))))
ifneq ($(realpath $(CURDIR)),$(_includer_dir))
$(error This Makefile must be invoked from $(_includer_dir). Use `make -C $(_includer_dir)` or `cd $(_includer_dir) && make`. CWD=$(CURDIR))
endif

# Enforce the contract above: abort with a clear message if a required
# variable is unset or empty.  Make evaluates $(error) at parse time, so
# this fires before any recipe runs.
$(foreach v,TRUST_ROOT BUILD_DIR, \
    $(if $(strip $($v)),, \
        $(error common.mk requires $v to be set before `include common.mk`)))

TRUST_VERSION ?= $(shell cat $(TRUST_ROOT)/VERSION)
TRUST_DOC_DIR ?= $(TRUST_ROOT)/docs
TRUSTIFY_DIR  ?= $(TRUST_ROOT)/Outils/trustify

export BUILD_DIR
export TRUST_VERSION

# Resolve the interpreter to install trustify into + run it from
# (TRUSTIFY_PYTHON / TRUSTIFY_PY_KIND / TRUSTIFY_USER_FLAG / TRUSTIFY).
# Shared with the TRUST root Makefile so both pick the same Python and
# feed the same flags to `make -C trustify install`.
include $(TRUSTIFY_DIR)/python-select.mk

# Use bash so the `run_doxygen` recipe below can rely on process
# substitution (`> >(tee ... | awk ...)`) to tee + filter doxygen's
# stdout without losing its exit code.
SHELL := /bin/bash

# -----------------------------------------------------------------------
# Required tools (parse-time check)
#
# Probe early so a missing tool surfaces immediately with an actionable
# message — otherwise the build dies later inside a recipe with cryptic
# shell output like "make: /bin/bash: No such file or directory", and the
# user has no hint where the dependency came from.
#
# Why each:
#   /bin/bash : SHELL is hard-set above; the run_doxygen recipe needs
#               process substitution, which POSIX /bin/sh does not
#               support.  Pre-installed on Linux and macOS; Alpine /
#               busybox containers need `apk add bash`.
#   curl      : NOT checked here — install_doxygen / install_mathjax
#               probe it themselves and only fail when the corresponding
#               externalpackages/ tarball is also absent.  A user with
#               both tarballs but no curl should still be able to build.
# -----------------------------------------------------------------------
ifeq ($(wildcard /bin/bash),)
$(error /bin/bash not found - required by the run_doxygen recipe for process substitution. See "Prerequisites" in $(TRUST_DOC_DIR)/README.md)
endif

# Generated artefact layout (derived from BUILD_DIR).  The includer may
# override any of these before the `include` to shape the build tree;
# defaults shown below are what TRUST and every baltik use.
#
#   docs/build/                       ($(BUILD_DIR))
#     ├── generated_content/          ($(GENERATED_CONTENT_DIR))
#     │   ├── kw_ref/                 ($(KW_DST_DIR))
#     │   └── ... (other generated md added by the includer)
#     ├── html/                       ($(HTML_DIR))    doxygen output
#     └── Doxyfile.chain              (baltik only, written by baltik.mk)
#
# generated_content/ mirrors docs/content/ in role: both hold markdown
# pages that doxygen scans via INPUT.  The split keeps source pages
# (under content/) cleanly separated from generated ones.
GENERATED_CONTENT_DIR ?= $(BUILD_DIR)/generated_content
KW_DST_DIR            ?= $(GENERATED_CONTENT_DIR)/kw_ref
HTML_DIR              ?= $(BUILD_DIR)/html

# Favicon: no Makefile involvement.  Each project declares its own in its
# own Doxyfile via `HTML_EXTRA_FILES += <path>/favicon.ico` (page chrome,
# like the HTML header and any sidebar logo) — see the project Doxyfiles.

# -----------------------------------------------------------------------
# Doxygen binary resolution (three-way fallback, in order):
#
#   1. TRUST_DOXYGEN_BINARY=/path/to/doxygen
#        → user-provided binary, skip extraction.
#   2. $(DOXY_TARBALL) present (i.e. externalpackages/ was cloned)
#        → extract the bundled tarball into $(TRUST_ROOT)/exec/doxygen/.
#   3. Network available
#        → download the official Linux release from doxygen.nl and
#          extract it into $(TRUST_ROOT)/exec/doxygen/.
#
# Case 3 covers readthedocs / lightweight clones that don't pull
# externalpackages/ (~hundreds of MB).  Requires curl on PATH.
# -----------------------------------------------------------------------
DOXY_VERSION = 1.16.1
DOXY_TARBALL = $(TRUST_ROOT)/externalpackages/doxygen/doxygen-$(DOXY_VERSION).linux.bin.tar.gz
DOXY_URL     = https://www.doxygen.nl/files/doxygen-$(DOXY_VERSION).linux.bin.tar.gz
DOXY_INSTALL = $(TRUST_ROOT)/exec/doxygen/bin/doxygen

ifneq ($(strip $(TRUST_DOXYGEN_BINARY)),)
DOXYGEN := $(TRUST_DOXYGEN_BINARY)
else
DOXYGEN := $(DOXY_INSTALL)
endif

.PHONY: install_trustify install_doxygen install_mathjax clean_doxygen clean_mathjax clean_build clean_html clean_kw_ref run_doxygen help

# -----------------------------------------------------------------------
# Doxygen install (extract bundled tarball on first use, or validate the
# user-provided TRUST_DOXYGEN_BINARY)
# -----------------------------------------------------------------------
ifneq ($(strip $(TRUST_DOXYGEN_BINARY)),)
install_doxygen:
	@test -x "$(DOXYGEN)" || { \
	    echo "ERROR: TRUST_DOXYGEN_BINARY=$(DOXYGEN) is not an executable file" >&2; \
	    exit 1; }
	@"$(DOXYGEN)" --version 2>/dev/null | head -n 1 | grep -qE '^[0-9]+\.[0-9]+' || { \
	    echo "ERROR: TRUST_DOXYGEN_BINARY=$(DOXYGEN) does not look like a doxygen binary" >&2; \
	    echo "       (\`--version\` did not produce a recognisable version line)" >&2; \
	    exit 1; }
	@echo "[doc] Using TRUST_DOXYGEN_BINARY=$(DOXYGEN)"
else
install_doxygen: $(DOXY_INSTALL)
endif

# Source the tarball from externalpackages/ if available, else download
# the official release from doxygen.nl.  Either way, extract into a
# scratch dir (cleaned up via trap, even on failure) and copy out the
# binary.
#
# Post-copy: validate the extracted binary actually runs on THIS host.
# The bundled tarball is Linux x86_64 — on macOS, ARM Linux, or any
# non-glibc system the cp succeeds but `--version` would later die with
# "exec format error" or "cannot execute binary file", and Make would
# halt during run_doxygen with cryptic shell output.  Probing here
# turns that into an actionable error with a pointer to the override.
# Same probe shape as the TRUST_DOXYGEN_BINARY branch above.
$(DOXY_INSTALL):
	@tmp=$$(mktemp -d); \
	 trap 'rm -rf "$$tmp"' EXIT; \
	 if [ -f "$(DOXY_TARBALL)" ]; then \
	     echo "[doc] Unpacking doxygen-$(DOXY_VERSION) from externalpackages ..."; \
	     tar xzf "$(DOXY_TARBALL)" -C "$$tmp" || exit $$?; \
	 else \
	     echo "[doc] $(DOXY_TARBALL) not found - downloading from doxygen.nl ..."; \
	     command -v curl > /dev/null 2>&1 || { \
	         echo "ERROR: curl is required to download doxygen when externalpackages/ is absent." >&2; \
	         echo "       Either install curl, populate externalpackages/, or" >&2; \
	         echo "       set TRUST_DOXYGEN_BINARY=/path/to/doxygen." >&2; \
	         exit 1; }; \
	     curl -fsSL "$(DOXY_URL)" | tar xzf - -C "$$tmp" || { \
	         echo "ERROR: failed to download/extract $(DOXY_URL)" >&2; \
	         echo "       Set TRUST_DOXYGEN_BINARY=/path/to/doxygen to bypass this step." >&2; \
	         exit 1; }; \
	 fi; \
	 mkdir -p $(dir $@) && \
	 cp "$$tmp/doxygen-$(DOXY_VERSION)/bin/doxygen" $@
	@"$@" --version 2>/dev/null | head -n 1 | grep -qE '^[0-9]+\.[0-9]+' || { \
	    rm -f "$@"; \
	    echo "ERROR: extracted $@ does not run on this host." >&2; \
	    echo "       The bundled tarball is a Linux x86_64 binary; on macOS," >&2; \
	    echo "       ARM, or non-glibc systems it cannot execute." >&2; \
	    echo "       Set TRUST_DOXYGEN_BINARY=/path/to/doxygen to use a native build," >&2; \
	    echo "       or see 'Building doxygen from source' in $(TRUST_DOC_DIR)/README.md." >&2; \
	    exit 1; \
	}
	@echo "[doc] Installed $@"

clean_doxygen:                                ## Remove the extracted doxygen binary
	@echo "[doc] Removing $(dir $(DOXY_INSTALL)) ..."
	@rm -rf $(dir $(DOXY_INSTALL))

# -----------------------------------------------------------------------
# MathJax 3 bundle (offline math rendering).
#
# Without a local copy the Doxyfile's MATHJAX_RELPATH would point at
# jsdelivr — meaning the browser must fetch MathJax over the network at
# every page view, and offline readers see raw `\(...\)` placeholders
# instead of formulas.  install_mathjax extracts the bundle into
# $(MATHJAX_CACHE); the run_doxygen recipe then mirrors es5/ into
# $(HTML_DIR)/mathjax/ so the generated tree is self-contained.
#
# Resolution mirrors install_doxygen: prefer the tarball in
# externalpackages/, else download from npm if curl is on PATH.
# On version bump, `make clean_mathjax` to evict the cached copy.
# -----------------------------------------------------------------------
MATHJAX_VERSION_NUM = 3.2.2
MATHJAX_TARBALL     = $(TRUST_ROOT)/externalpackages/mathjax/mathjax-$(MATHJAX_VERSION_NUM).tgz
MATHJAX_URL         = https://registry.npmjs.org/mathjax/-/mathjax-$(MATHJAX_VERSION_NUM).tgz
MATHJAX_CACHE       = $(TRUST_ROOT)/exec/mathjax
MATHJAX_INSTALL     = $(MATHJAX_CACHE)/es5/tex-chtml.js

install_mathjax: $(MATHJAX_INSTALL)

$(MATHJAX_INSTALL):
	@tmp=$$(mktemp -d); \
	 trap 'rm -rf "$$tmp"' EXIT; \
	 if [ -f "$(MATHJAX_TARBALL)" ]; then \
	     echo "[doc] Unpacking mathjax-$(MATHJAX_VERSION_NUM) from externalpackages ..."; \
	     tar xzf "$(MATHJAX_TARBALL)" -C "$$tmp" || exit $$?; \
	 else \
	     echo "[doc] $(MATHJAX_TARBALL) not found - downloading from npm ..."; \
	     command -v curl > /dev/null 2>&1 || { \
	         echo "ERROR: curl is required to download MathJax when externalpackages/ is absent." >&2; \
	         echo "       Either install curl, or drop mathjax-$(MATHJAX_VERSION_NUM).tgz into" >&2; \
	         echo "       $(dir $(MATHJAX_TARBALL))." >&2; \
	         exit 1; }; \
	     curl -fsSL "$(MATHJAX_URL)" | tar xzf - -C "$$tmp" || { \
	         echo "ERROR: failed to download/extract $(MATHJAX_URL)" >&2; \
	         echo "       Drop mathjax-$(MATHJAX_VERSION_NUM).tgz into $(dir $(MATHJAX_TARBALL)) to bypass this step." >&2; \
	         exit 1; }; \
	 fi; \
	 rm -rf $(MATHJAX_CACHE) && \
	 mkdir -p $(MATHJAX_CACHE) && \
	 mv "$$tmp/package/es5" $(MATHJAX_CACHE)/es5
	@echo "[doc] Installed MathJax 3 bundle at $(MATHJAX_CACHE)/"

clean_mathjax:                                ## Remove the extracted MathJax bundle
	@echo "[doc] Removing $(MATHJAX_CACHE)/ ..."
	@rm -rf $(MATHJAX_CACHE)

# -----------------------------------------------------------------------
# Shared clean components.
#
# These wipe individual stages of the pipeline.  `clean_build` is the
# universal nuke (everything lives under $(BUILD_DIR)).  The finer-
# grained targets are kept for fast iteration during development.
# -----------------------------------------------------------------------
clean_build:                                  ## Remove the entire build/ directory
	@echo "[doc] Removing $(BUILD_DIR)/ ..."
	@rm -rf $(BUILD_DIR)

clean_html:                                   ## Remove the HTML output directory
	@echo "[doc] Removing $(HTML_DIR)/ ..."
	@rm -rf $(HTML_DIR)

clean_kw_ref:                                 ## Remove the generated keyword reference
	@echo "[doc] Removing $(KW_DST_DIR)/ ..."
	@rm -rf $(KW_DST_DIR)

# -----------------------------------------------------------------------
# trustify install
# -----------------------------------------------------------------------
install_trustify:                             ## Install trustify into the resolved Python (TRUSTIFY_PYTHON) if its CLI is unavailable
	@# Fast path: if trustify already imports from our interpreter, skip
	@# (covers the conda env where ./configure pre-installed it, and RTD
	@# workspace reuse).  Probe with `-m trustify`, NOT `command -v
	@# trustify`: a stray PATH binary of that name must not masquerade, and
	@# `-m` pins the probe to TRUSTIFY_PYTHON. (3.7)
	@#
	@# No locking: make builds this prerequisite at most once per invocation
	@# (even under -j — a shared prereq is deduped), ./configure pre-installs
	@# trustify into the conda env so the fast path skips the install in the
	@# normal case, and the only residual race (two SEPARATE make processes
	@# pip-installing into one site-packages at once) is rare and was never
	@# covered by the old per-BUILD_DIR lock anyway.  Dropping flock also
	@# removes a hard parse-time dependency that broke `make` (even clean /
	@# help) on macOS / Alpine.
	@#
	@# The opt-in gate, pydantic precondition and offline-first two-attempt
	@# pip live in the trustify Makefile's `install` target (single source of
	@# truth); we just hand it the resolved interpreter + flags, and let its
	@# errors propagate.
	@if $(TRUSTIFY) --version >/dev/null 2>&1; then \
	    echo "[doc] trustify already available via '$(TRUSTIFY)'."; \
	else \
	    echo "[doc] Installing trustify into $(TRUSTIFY_PYTHON) [$(TRUSTIFY_PY_KIND)] ..."; \
	    $(MAKE) -C $(TRUSTIFY_DIR) install PYTHON=$(TRUSTIFY_PYTHON) USER_FLAG=$(TRUSTIFY_USER_FLAG) PY_KIND=$(TRUSTIFY_PY_KIND); \
	fi

# -----------------------------------------------------------------------
# Run doxygen + post-process the generated HTML
#
# After doxygen runs, copy the trustify-generated figures/ subdir into
# $(HTML_DIR)/figures.  The generated keyword markdown references images
# as `figures/<file>` (relative path), and Doxygen does NOT propagate
# markdown-embedded image directories to the HTML output by itself.
# -----------------------------------------------------------------------
run_doxygen: install_doxygen install_mathjax
	@mkdir -p $(BUILD_DIR)
	@echo "[doc] Running doxygen ($(DOXYGEN)) - typically a few minutes."
	@echo "[doc]   full log:      $(BUILD_DIR)/doxygen.log"
	@echo "[doc]   warnings only: $(BUILD_DIR)/doxygen_warnings.log"
	@# Tee full stdout to the log, but on terminal collapse the per-item
	@# progress spam (Patching/Parsing/dot/...) into one line each via
	@# filter_doxygen.awk.  Process substitution preserves doxygen's exit
	@# code; warnings are already routed to a separate file via
	@# WARN_LOGFILE in the Doxyfile.
	@$(DOXYGEN) Doxyfile > >( \
	    tee $(BUILD_DIR)/doxygen.log | awk -f $(TRUST_DOC_DIR)/scripts/filter_doxygen.awk \
	)
	@echo "[doc] HTML output written to $(HTML_DIR)/"
	@if [ -d "$(KW_DST_DIR)/figures" ]; then \
	    rm -rf $(HTML_DIR)/figures; \
	    cp -r $(KW_DST_DIR)/figures $(HTML_DIR)/figures; \
	    echo "[doc] Copied $(KW_DST_DIR)/figures/ into $(HTML_DIR)/figures/."; \
	fi
	@echo "[doc] Bundling MathJax 3 into HTML output ..."
	@rm -rf $(HTML_DIR)/mathjax
	@mkdir -p $(HTML_DIR)/mathjax
	@cp -r $(MATHJAX_CACHE)/es5 $(HTML_DIR)/mathjax/es5
	@echo "[doc] Flattening navtree ..."
	@python3 $(TRUST_DOC_DIR)/scripts/flatten_navtree.py $(HTML_DIR)
	@echo "[doc] Navtree flattened."
	@python3 $(TRUST_DOC_DIR)/scripts/patch_filter.py $(HTML_DIR)
	@echo "[doc] Search bar customized."

# -----------------------------------------------------------------------
# Help
#
# The preamble (Usage line, combo-target summary, "see also" pointers)
# is supplied by the includer via the HELP_PREAMBLE variable (use
# `define ... endef` for multi-line content; this file `export`s it).
# The low-level target list is scraped from $(MAKEFILE_LIST) - any
# `target:    ## description` line in either the includer's Makefile or
# any included `.mk` file shows up automatically.
# -----------------------------------------------------------------------
export HELP_PREAMBLE

help:
	@echo ""
	@printf '%s\n' "$$HELP_PREAMBLE"
	@echo ""
	@echo "Low-level targets:"
	@awk 'BEGIN {FS = ":.*##"} /^[a-zA-Z_-]+:.*##/ { printf "  %-24s %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""
	@echo "Environment overrides:"
	@echo "  TRUST_DOXYGEN_BINARY=/path/to/doxygen"
	@echo "                       Use a custom doxygen binary instead of the bundled"
	@echo "                       $(DOXY_VERSION) Linux binary from externalpackages/."
	@echo "                       See 'Building doxygen from source' in"
	@echo "                       \$$TRUST_ROOT/docs/README.md."
	@echo "  TRUSTIFY_USER_INSTALL=1"
	@echo "                       Allow installing trustify into a bare system Python"
	@echo "                       (no TRUST conda env, no active venv) via 'pip --user'."
	@echo "                       Resolved interpreter: $(TRUSTIFY_PYTHON) [$(TRUSTIFY_PY_KIND)]."
	@echo ""
