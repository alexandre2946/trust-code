# TRUST documentation build - baltik-specific shared rules
#
# Included by every baltik's docs/Makefile (alongside common.mk).
# Trust itself does NOT include this file - the `Doxyfile` rule below
# would regenerate it from the baltik template, which is wrong for the
# trust tree.
#
# Must be included AFTER common.mk (uses install_trustify, KW_DST_DIR).
#
# Contract - the includer must set these BEFORE the include (in
# addition to common.mk's contract):
#   PROJECT_NAME         this baltik's name
#   project_directory    path to the baltik's root (for trustify -p)
#
# Optional:
#   SCRIPT_EXPAND        path to expand_env.py.  Defaults to
#                        $(TRUST_ROOT)/scripts/expand_env.py; override
#                        only if expand_env.py lives elsewhere in your
#                        tree.  The default lives here (not in the
#                        baltik template) so future relocations of
#                        expand_env.py do not break existing baltiks.
#   BALTIK_DEPENDENCIES  override the trustify auto-resolution.
#                        Set on the command line / in env to force a
#                        specific list, e.g.
#                          make BALTIK_DEPENDENCIES="/abs/dep_a /abs/dep_b"
#                        When unset/empty, the doxyfile_chain recipe
#                        resolves dependencies via trustify.

# Enforce the contract above.  TRUST_ROOT is also checked by common.mk
# when it is loaded first (the expected order), but we re-check it here
# in case a baltik drops common.mk or swaps the include order - the
# recipes below dereference $(TRUST_ROOT) directly.
$(foreach v,TRUST_ROOT PROJECT_NAME project_directory, \
    $(if $(strip $($v)),, \
        $(error baltik.mk requires $v to be set before `include baltik.mk`)))

SCRIPT_EXPAND            ?= $(TRUST_ROOT)/scripts/expand_env.py
BALTIK_DOXYFILE_TEMPLATE = $(TRUST_ROOT)/bin/baltik/templates/basic/docs/Doxyfile

.PHONY: trustify_doc doxyfile_chain

# Generate the keyword reference for this baltik (TRUST + this baltik's
# project.cfg [dependencies], transitively).  `trustify generate_markdown`
# already expands [dependencies] internally when we pass this baltik via
# -p, so we only need to point at this baltik's root.
trustify_doc: install_trustify                ## Generate the keyword reference for this baltik (TRUST + deps + own)
	@echo "[doc] Generating keyword reference (TRUST + $(PROJECT_NAME)) ..."
	@rm -rf $(KW_DST_DIR)
	@mkdir -p $(GENERATED_CONTENT_DIR)
	@$(TRUSTIFY) --trust-root $(TRUST_ROOT) -p $(project_directory) \
	    generate_markdown --out $(KW_DST_DIR)
	@echo "[doc] Keyword reference written to $(KW_DST_DIR)/"

# Generate the Doxyfile from the TRUST baltik template only when the
# baltik does not already have a local copy. This lets users customise
# their Doxyfile (e.g. tweak INPUT, GENERATE_LATEX, ...) without it
# being overwritten on every build. Delete docs/Doxyfile to regenerate
# from the current template.
#
# This is a file target with no prereqs: Make runs the recipe iff the
# file is missing. We intentionally do NOT depend on the template so
# user customisations survive template updates.
#
# The recipe is atomic: expand_env.py writes to Doxyfile.tmp first, and
# we only `mv` it to the final name on success.  Reason: expand_env.py
# exits 1 (and Make halts) when a @@VAR@@ token is unresolved, but
# STILL writes the partially-substituted output to its target file.
# Without the tmp+mv dance, the next `make` would see a Doxyfile on
# disk, skip this rule (no prereqs), and feed a broken Doxyfile to
# doxygen.  With the dance, a failed run leaves no final Doxyfile and
# the next `make` retries cleanly from the template.
Doxyfile:                                     ## Generate Doxyfile from the TRUST baltik template (only if missing)
	@echo "[doc] Generating Doxyfile from TRUST baltik template ..."
	@env PROJECT_NAME=$(PROJECT_NAME) $(SCRIPT_EXPAND) \
	    $(BALTIK_DOXYFILE_TEMPLATE) -o Doxyfile.tmp \
	    && mv Doxyfile.tmp Doxyfile \
	    || { rm -f Doxyfile.tmp; exit 1; }

# Always (re)regenerate Doxyfile.chain so the C++ API tab picks up
# every baltik this one depends on without us touching the
# user-customisable Doxyfile.
#
# Dependency resolution is done inside the recipe (not at Make parse
# time) so that install_trustify can run first - the trustify CLI is
# what resolves the [dependencies] section of project.cfg, so it has
# to be on PATH before we call it.  A non-empty BALTIK_DEPENDENCIES
# (set in the environment or on the command line) bypasses the
# trustify call entirely.
#
# Iteration is line-based (`read -r` via here-string), NOT
# whitespace-split via `for path in $$deps`.  trustify emits one path
# per line by design (see `trustify projects --help`), and the
# line-based reader preserves embedded spaces in trustify-resolved
# paths.
#
# Caveat: the BALTIK_DEPENDENCIES env-var override remains nominally
# whitespace-separated (its documented contract).  A path with an
# embedded space passed via the override is still split, because the
# Make-side `$(BALTIK_DEPENDENCIES)` substitution feeds an unquoted
# shell expansion that word-splits on IFS.  Users with paths-with-
# spaces should drop the override and let trustify resolve naturally
# (the trustify-output branch IS space-safe).
doxyfile_chain: install_trustify              ## Regenerate Doxyfile.chain (dependency baltik INPUT extension)
	@echo "[doc] Regenerating $(BUILD_DIR)/Doxyfile.chain ..."
	@mkdir -p $(BUILD_DIR)
	@if [ -n "$(BALTIK_DEPENDENCIES)" ]; then \
	     deps=$$(printf '%s\n' $(BALTIK_DEPENDENCIES)); \
	 else \
	     deps=$$($(TRUSTIFY) --trust-root $(TRUST_ROOT) -p $(project_directory) projects --only=dependency) || { \
	         echo "ERROR: trustify failed to resolve dependencies for $(project_directory)" >&2; \
	         echo "       (Doxyfile.chain not written; the C++ API tab would be missing dependency baltiks.)" >&2; \
	         exit 1; \
	     }; \
	 fi; \
	 { echo "# Auto-generated by docs/Makefile - do not edit by hand."; \
	   echo "# Re-run 'make run_doxygen' (or 'make doxyfile_chain') to refresh."; \
	   while IFS= read -r path; do \
	       [ -z "$$path" ] && continue; \
	       echo "INPUT                 += $$path/src"; \
	   done <<< "$$deps"; \
	 } > $(BUILD_DIR)/Doxyfile.chain
