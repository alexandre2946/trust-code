# Shared Python-interpreter selection for installing + running trustify.
#
# Included by docs/common.mk (TRUST docs + every baltik) and the TRUST
# root Makefile so they pick the SAME interpreter and feed the SAME flags
# to the trustify Makefile's `install` target.  Variables ONLY: this file
# defines no targets, so including it never changes the includer's
# default goal.
#
# Exposes (all prefixed TRUSTIFY_, plus TRUSTIFY itself):
#   TRUSTIFY_PYTHON    interpreter to install into / run trustify from
#   TRUSTIFY_PY_KIND   conda | venv | system   (drives the install gate)
#   TRUSTIFY_USER_FLAG '--user' when a system-Python install is opted in
#                      via TRUSTIFY_USER_INSTALL, else empty
#   TRUSTIFY           '$(TRUSTIFY_PYTHON) -m trustify' (the invocation)
#
# The opt-in gate, the pydantic precondition and the offline-first
# two-attempt pip all live in the trustify Makefile's `install` target,
# which consumes PYTHON / USER_FLAG / PY_KIND.  Keeping policy (here) and
# mechanism (there) split lets configurer_env stay on its own direct-pip
# path while docs + the root Makefile share both halves.
#
# (USER_FLAG, not PIP_USER: make exports command-line variable assignments,
# and pip reads PIP_<OPTION> env vars — an exported PIP_USER=--user would be
# parsed as the boolean `user` option and abort pip.)

# Capture this file's own path on its FIRST line (before any later
# include shifts MAKEFILE_LIST) to derive a TRUST_ROOT default: this file
# lives at $(TRUST_ROOT)/Outils/trustify/python-select.mk.  Callers that
# already know TRUST_ROOT (docs/common.mk, the root Makefile) set it
# first, so the `?=` is a no-op for them.
_TRUSTIFY_SELECT_MK := $(lastword $(MAKEFILE_LIST))
TRUST_ROOT ?= $(realpath $(dir $(_TRUSTIFY_SELECT_MK))../..)

# Resolution order:
#   1. conda  - bundled env at exec/python, IFF a `conda` binary sits
#               next to `python` (configurer_env's own discriminator,
#               cf. `[ -f exec/python/bin/conda ]`).  Ships pydantic +
#               usually trustify; its bundled pip builds offline.
#   2. venv   - an activated virtualenv ($VIRTUAL_ENV).
#   3. system - a bare interpreter.  `./configure -without-conda` only
#               SYMLINKS exec/python/bin/python -> the system python3 (no
#               bundled pip, inherits the caller's env), and "no configure
#               at all" leaves just python3 on PATH.  Both write into a
#               shared / possibly PEP-668 site, so the install is gated
#               behind TRUSTIFY_USER_INSTALL and done with `pip --user`.
#               Prefer the configure-blessed exec/python symlink over PATH
#               when it exists.
_TRUSTIFY_PY    := $(wildcard $(TRUST_ROOT)/exec/python/bin/python)
_TRUSTIFY_CONDA := $(wildcard $(TRUST_ROOT)/exec/python/bin/conda)
ifneq ($(_TRUSTIFY_CONDA),)
TRUSTIFY_PYTHON  := $(_TRUSTIFY_PY)
TRUSTIFY_PY_KIND := conda
else ifneq ($(strip $(VIRTUAL_ENV)),)
TRUSTIFY_PYTHON  := $(VIRTUAL_ENV)/bin/python
TRUSTIFY_PY_KIND := venv
else ifneq ($(_TRUSTIFY_PY),)
TRUSTIFY_PYTHON  := $(_TRUSTIFY_PY)
TRUSTIFY_PY_KIND := system
else
TRUSTIFY_PYTHON  := python3
TRUSTIFY_PY_KIND := system
endif

# --user only for an opted-in system Python; conda/venv own their env.
# When PY_KIND=system and this stays empty, `install` errors with guidance.
TRUSTIFY_USER_FLAG :=
ifeq ($(TRUSTIFY_PY_KIND),system)
ifneq ($(strip $(TRUSTIFY_USER_INSTALL)),)
TRUSTIFY_USER_FLAG := --user
endif
endif

# Always invoke trustify as `<python> -m trustify` so the chosen
# interpreter is honoured regardless of where its console script landed
# (next to python for conda/venv, in the user bin dir for --user).
TRUSTIFY := $(TRUSTIFY_PYTHON) -m trustify
