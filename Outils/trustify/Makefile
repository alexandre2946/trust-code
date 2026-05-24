
SHELL := /bin/bash

# `project_directory` is the lowest-priority `projects` source
# consulted by `effective_projects(None)` (see projects.py). When set
# in the surrounding shell (e.g. after sourcing env_TrioCFD.sh or
# cd'ing into a baltik), it silently overlays that baltik onto every
# trustify invocation the test suite or batch-check launches —
# causing confusing test failures from unrelated baltik content. The
# CLI flags `--projects` / `--trust-root` override this fallback, but
# the trustify tests rely on the default code path. Drop the variable
# from the make-spawned environment so the gotcha cannot bite here.
# Section 7 / CLAUDE.md "`$project_directory` env leak" gotcha.
unexport project_directory

# --force-reinstall rewrites every file trustify owns; --no-deps stops
# pip from also re-resolving + re-downloading transitive deps (notably
# pydantic), which would clobber whatever the conda env already provides.
PIP_OPT ?= --force-reinstall --no-deps

# Interpreter + flags for `install`.  Callers (docs/common.mk, the TRUST
# root Makefile) pass these from python-select.mk so the install targets
# the right interpreter and respects the system-Python opt-in gate; a
# direct `make install` defaults to python3 in the current environment
# (PY_KIND empty -> gate never fires), matching historical behaviour.
#
# USER_FLAG carries pip's `--user` (or empty).  It must NOT be named
# PIP_USER: make auto-exports command-line variable assignments into the
# recipe environment, and pip reads PIP_<OPTION> env vars — an exported
# PIP_USER=--user is parsed as the boolean `user` option and pip aborts
# with "--user is not a valid value for user option".
PYTHON    ?= python3
USER_FLAG ?=
PY_KIND   ?=

# A non-editable install copies src/ into site-packages, so source
# edits silently test the stale copy and waste a debugging cycle
# before anyone notices. Warn (do NOT fail) at both ends of `test`.
define WARN_IF_NOT_EDITABLE
@if ! pip show trustify 2>/dev/null | grep -q '^Editable project location:'; then \
		echo "" >&2; \
		echo "WARNING: trustify is NOT installed as editable (pip install -e)." >&2; \
		echo "         Edits to src/ will NOT be picked up by these tests." >&2; \
		echo "         Run 'make install-dev' to fix." >&2; \
		echo "" >&2; \
	fi
endef

all: install test check

clean:
	@rm -rfv **/*.egg-info build/ tests/generated/*
	@find . -name __pycache__ -exec rm -rfv {} +

install:
	@# Gate: refuse to silently install into a shared / PEP-668 system
	@# Python.  python-select.mk classifies such interpreters as
	@# PY_KIND=system and leaves USER_FLAG empty unless TRUSTIFY_USER_INSTALL
	@# opted in (then USER_FLAG=--user).  A direct `make install` passes no
	@# PY_KIND, so this never fires for plain dev use.
	@if [ "$(PY_KIND)" = "system" ] && [ -z "$(USER_FLAG)" ]; then \
	    echo "ERROR: refusing to install trustify into a bare system Python:"      >&2; \
	    echo "         $(PYTHON)"                                                   >&2; \
	    echo "       (no TRUST conda env and no active virtualenv detected). Either:" >&2; \
	    echo "         - create + activate a venv, then re-run, or"                 >&2; \
	    echo "         - set TRUSTIFY_USER_INSTALL=1 for a 'pip --user' install."   >&2; \
	    exit 1; \
	fi
	@# venv is an isolated env the user created for this build, so install
	@# into it freely: provision the runtime dep (pydantic) AND the build dep
	@# (setuptools, which a stock `python -m venv` omits on 3.12+) up front.
	@# That makes the pydantic check below pass and lets the --no-build-
	@# isolation attempt succeed offline-after, instead of erroring its way to
	@# the build-isolated fallback.  Nothing outside the venv is touched.
	@# conda / system Pythons are externally managed (their pydantic is conda-
	@# or distro-owned), so we never auto-install into them — the gate above
	@# and the pydantic check below keep --no-deps honest there.
	@if [ "$(PY_KIND)" = "venv" ]; then \
	    echo "[trustify] Provisioning venv with pydantic + setuptools ..."; \
	    $(PYTHON) -m pip install "pydantic>=2,<3" setuptools || { \
	        echo "ERROR: could not install pydantic + setuptools into the venv."       >&2; \
	        echo "         $(PYTHON)"                                                   >&2; \
	        echo "       A fresh venv needs network for this one-off provisioning;"     >&2; \
	        echo "       afterwards it builds trustify offline.  Either run it online:" >&2; \
	        echo "         $(PYTHON) -m pip install \"pydantic>=2,<3\" setuptools"       >&2; \
	        echo "       or build from a TRUST conda env, which ships both (offline)."  >&2; \
	        exit 1; }; \
	fi
	@# pydantic >=2,<3 must be importable: we install --no-deps, so pip will
	@# not pull it in, and the trustify CLI needs pydantic v2 at runtime.
	@# Check the MAJOR version, not just importability — a too-old pydantic
	@# (e.g. Debian's v1 python3-pydantic) imports fine and even satisfies
	@# `trustify --version`, but then breaks deep in dataset parsing with a
	@# v1/v2 API traceback.  (Auto-satisfied for venv by provisioning above.)
	@$(PYTHON) -c "import pydantic,sys; sys.exit(0 if pydantic.VERSION.split('.')[0]=='2' else 1)" >/dev/null 2>&1 || { \
	    echo "ERROR: trustify requires pydantic >=2,<3, not satisfied by:"     >&2; \
	    echo "         $(PYTHON)"                                              >&2; \
	    found=$$($(PYTHON) -c "import pydantic; print(pydantic.VERSION)" 2>/dev/null); \
	    if [ -n "$$found" ]; then \
	        echo "       found pydantic $$found (trustify needs major version 2)." >&2; \
	    else \
	        echo "       pydantic is not installed." >&2; \
	    fi; \
	    echo "       Install it: $(PYTHON) -m pip install $(USER_FLAG) \"pydantic>=2,<3\"" >&2; \
	    exit 1; }
	@echo "Installing trustify into $(PYTHON) [$(PIP_OPT) $(USER_FLAG)]"
	@# Choose the build mode up front from whether setuptools is present,
	@# rather than trying --no-build-isolation and falling back — the
	@# fallback prints a confusing "No module named 'setuptools'" even when
	@# the retry succeeds, and stacks two error dumps when both fail.
	@#   setuptools present -> --no-build-isolation (builds offline; true for
	@#                         conda, provisioned venvs, most system Pythons).
	@#   setuptools absent   -> build isolation, so pip provisions setuptools
	@#                         in a temp build env (needs network; one clean
	@#                         "could not find setuptools" error if offline).
	@if $(PYTHON) -c "import setuptools" >/dev/null 2>&1; then \
	    $(PYTHON) -m pip install $(PIP_OPT) $(USER_FLAG) --no-build-isolation --no-cache-dir .; \
	else \
	    $(PYTHON) -m pip install $(PIP_OPT) $(USER_FLAG) .; \
	fi

install-dev:
	@echo "Installing trustify (editable, with [dev] extras)"
	@pip install -e ".[dev]"

test:
	$(WARN_IF_NOT_EDITABLE)
	@echo ""
	@echo "======= Starting unit tests of trustify ========="
	@echo ""
	@python -m unittest discover -s tests -t . -v
ifeq ($(TRUST_DISABLE_JUPYTER),0)
	@jupyter-nbconvert --to html --no-input --output $(TRUST_TMP)/trustify_ex.html --execute docs/examples/hello_world/HelloWorld.ipynb
endif
	$(WARN_IF_NOT_EDITABLE)

test_debug:
	@TRUSTIFY_DEBUG=1 python -m unittest discover -s tests -t . -v

lint:
	@ruff check src/ tests/
	@ruff format --check src/ tests/

# Coverage report. Wires the unittest suite under `coverage run` and
# instruments subprocess-spawning tests (CLI / LANG=C / concurrent
# imports) via PYTHONPATH→sitecustomize so coverage tracks the child
# processes too. `parallel = true` in pyproject.toml means every
# process writes its own data file (.coverage.<host>.<pid>.<rand>);
# `coverage combine` merges them.
COVERAGE_HTML_DIR ?= $(CURDIR)/htmlcov
coverage:
	@rm -rf .coverage .coverage.* $(COVERAGE_HTML_DIR)/
	@COVERAGE_PROCESS_START=$(CURDIR)/pyproject.toml \
	    PYTHONPATH=$(CURDIR)/tests/_coverage_sitecustomize:$$PYTHONPATH \
	    coverage run --rcfile=$(CURDIR)/pyproject.toml -m unittest discover -s tests -t .
	@coverage combine --rcfile=$(CURDIR)/pyproject.toml
	@coverage report --rcfile=$(CURDIR)/pyproject.toml
	@coverage html --rcfile=$(CURDIR)/pyproject.toml -d $(COVERAGE_HTML_DIR)
	@echo "HTML report: $(COVERAGE_HTML_DIR)/index.html"


TRUSTIFY_CHECK_LOG ?= $(CURDIR)/check.log
# Cap at 16 (audit 6.7): `nproc` reports the full physical node on HPC
# login machines (often 128+ cores), but a user allocation typically
# gets only a handful — the implicit fan-out DoS'd the login node.
# Override on the make command line: `make check TRUSTIFY_CHECK_JOBS=N`.
TRUSTIFY_CHECK_JOBS ?= $(shell n=$$(nproc 2>/dev/null || echo 1); [ "$$n" -gt 16 ] && n=16; echo "$$n")

# Probe trustify availability. Fails (exit 1) with a clear message if
# the CLI is missing or broken — used as a make-level dependency by
# every target that needs trustify, so the rest of the recipe never
# runs against a missing binary. Audit 3.10.
check_trustify_install:
	@if ! command -v trustify >/dev/null 2>&1 || ! trustify --version >/dev/null 2>&1; then \
		echo "trustify not available — re-run ./configure or pip install $$TRUST_ROOT/Outils/trustify/" >&2; \
		exit 1; \
	fi

check: check_trustify_install
	@if [ -z "$$TRUST_ROOT" ]; then echo "TRUST_ROOT is not set" >&2; exit 1; fi
	@echo "Running trustify batch-check on $$TRUST_ROOT/tests (--jobs $(TRUSTIFY_CHECK_JOBS)). See $(TRUSTIFY_CHECK_LOG) in cas of failure"
	@rm -f $(TRUSTIFY_CHECK_LOG)
	@trustify batch-check $$TRUST_ROOT/tests --summary-out $(TRUSTIFY_CHECK_LOG) --jobs $(TRUSTIFY_CHECK_JOBS)

help:
	@echo "The following make targets are available"
	@echo "   - all: alias for test"
	@echo "   - clean: cleans build artifacts, generated test artifacts, __pycache__"
	@echo "   - install: pip install trustify (PIP_OPT=--force-reinstall --no-deps by default)"
	@echo "   - install-dev: pip install -e trustify (editable, with [dev] extras: ruff)"
	@echo "   - test: launch the unittests of the trustify package"
	@echo "   - test_debug: launch the unittests in debug mode (very verbose!)"
	@echo "   - lint: run ruff check + ruff format --check on src/ and tests/"
	@echo "   - coverage: run the unittest suite under coverage.py, write htmlcov/"
	@echo "   - check: run trustify batch-check on all .data files in \$$TRUST_ROOT/tests"
	@echo "   - help: displays this message"
	@echo ""
	@echo "For production install / dataset checking, see the 'trustify' console script:"
	@echo "   trustify -p <project_dir> generate_schema"
	@echo "   trustify check <file.data>"
	@echo "   trustify batch-check <file_or_dir> [<file_or_dir> ...]"

.PHONY: help clean test test_debug lint coverage check check_trustify_install all install install-dev
