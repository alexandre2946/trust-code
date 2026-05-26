# Makefile for trustutils documentation
#
# Targets
# -------
#   make md      Generate Doxygen-compatible Markdown (default)
#   make html    Build HTML with Sphinx (requires sphinx-build)
#   make clean   Remove build/

TRUST_ROOT    ?= $(shell cd ../../.. && pwd)
TRUSTUTILS_DIR = $(TRUST_ROOT)/Validation/Outils/trustutils
BUILDDIR       = build

# Sphinx settings (only used by the html / other Sphinx targets)
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR      = source

.PHONY: md html clean help

md:                    ## Generate Doxygen-compatible Markdown into build/md/
	@echo "[doc] Generating Markdown API docs for trustutils ..."
	@mkdir -p $(BUILDDIR)/md
	@python3 gen_doc_md.py $(TRUSTUTILS_DIR) $(BUILDDIR)/md
	@echo "[doc] Markdown written to $(BUILDDIR)/md/"

html:                  ## Build HTML documentation with Sphinx
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:                 ## Remove the build/ directory
	@rm -rf $(BUILDDIR)
	@echo "[doc] build/ removed."

help:                  ## Show this help message
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@awk 'BEGIN {FS = ":.*##"} /^[a-zA-Z_-]+:.*##/ { printf "  %-12s %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""
