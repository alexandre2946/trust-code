SHELL=/bin/bash

.PHONY: all tools optim semi_optim opt debug prof gcov semi_opt profiling \
	custom unit opt_avx clean validation check check_optim check_semi_opt \
	check_debug ctest_debug ctest_optim ctest_semi_opt ctest_semi_optim \
	ctest_rerun_failed_debug ctest_rerun_failed_optim \
	ctest_rerun_failed_semi_opt ctest_rerun_failed_semi_optim check_perf \
	install_trustify trustify_check docs _check_env

all: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_opt numeric_kernel_opt standard_kernel_opt opt micro_kernel_debug numeric_kernel_debug standard_kernel_debug debug" monodir)
	
tools: _check_env
	(source ./env_TRUST.sh && compile tools)
optim:	 opt
semi_optim:	 semi_opt
opt: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_opt numeric_kernel_opt standard_kernel_opt opt" monodir)
debug: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_debug numeric_kernel_debug standard_kernel_debug debug" monodir)
prof: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_prof numeric_kernel_prof standard_kernel_prof prof" monodir)
gcov: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_gcov numeric_kernel_gcov standard_kernel_gcov gcov" monodir)
semi_opt: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_semi_opt numeric_kernel_semi_opt standard_kernel_semi_opt semi_opt" monodir)
profiling: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_profiling numeric_kernel_profiling standard_kernel_profiling profiling" monodir)
custom: _check_env
	(source ./env_TRUST.sh && env cibles="micro_kernel_custom numeric_kernel_custom standard_kernel_custom custom" monodir)
unit: _check_env
	(source ./env_TRUST.sh && env cibles="numeric_kernel_debug" monodir)
opt_avx:
	@echo opt_avx target not available since TRUST v1.9.2
clean: _check_env
	(source ./env_TRUST.sh && compile clean)
	(source ./env_TRUST.sh && make -C docs clean)
validation: _check_env
	(source ./env_TRUST.sh && cd Validation/Rapports_automatiques && ./lance_tout)
check: _check_env
	(source ./env_TRUST.sh && echo 0 | lance_test ${exec_opt} ${TRUST_ROOT}/build )
check_optim:	check
check_semi_opt: _check_env
	(source ./env_TRUST.sh && echo 0 | lance_test ${exec_semi_opt} ${TRUST_ROOT}/build )
check_debug: _check_env
	(source ./env_TRUST.sh && echo 0 | lance_test ${exec_debug} ${TRUST_ROOT}/build )
ctest_debug: _check_env
	(source ./env_TRUST.sh && lance_test -ctest ${exec_debug})
ctest_optim: _check_env
	(source ./env_TRUST.sh && lance_test -ctest ${exec_opt})
ctest_semi_opt: _check_env
	(source ./env_TRUST.sh && lance_test -ctest ${exec_semi_opt})
ctest_semi_optim:	 ctest_semi_opt
ctest_rerun_failed_debug: _check_env
	(source ./env_TRUST.sh && env CTEST_OPTS="--rerun-failed" lance_test -ctest ${exec_debug})
ctest_rerun_failed_optim: _check_env
	(source ./env_TRUST.sh && env CTEST_OPTS="--rerun-failed" lance_test -ctest ${exec_opt})
ctest_rerun_failed_semi_opt: _check_env
	(source ./env_TRUST.sh && env CTEST_OPTS="--rerun-failed" lance_test -ctest ${exec_semi_opt})
ctest_rerun_failed_semi_optim:	 ctest_rerun_failed_semi_opt
check_perf: _check_env
	(source ./env_TRUST.sh && compare_TU)


# Resolve the Python to install trustify into (conda exec/python > active
# venv > system Python, with a TRUSTIFY_USER_INSTALL opt-in for the last)
# and share it with the trustify Makefile's `install` target.  Same
# selection docs/common.mk uses; defines variables only (no targets), so
# the default goal (`all`) is unchanged.
TRUST_ROOT := $(CURDIR)
include Outils/trustify/python-select.mk

# Refresh the trustify install in the resolved Python (the bundled conda
# env when configured with conda) without re-running configure. Useful
# after editing Outils/trustify/src/ — alternatively, use `make -C
# Outils/trustify install-dev` for an editable install that picks up edits
# without a re-install. Off-conda (./configure -without-conda, or no
# configure) this targets a venv / system Python and honours
# TRUSTIFY_USER_INSTALL=1 (see Outils/trustify/python-select.mk).
install_trustify: _check_env
	$(MAKE) -C Outils/trustify install \
	    PYTHON=$(TRUSTIFY_PYTHON) USER_FLAG=$(TRUSTIFY_USER_FLAG) PY_KIND=$(TRUSTIFY_PY_KIND)

# equivalent to running trustify batch-check tests -j0
trustify_check: _check_env
	$(MAKE) -C Outils/trustify check

docs:
	$(MAKE) -C docs all



_check_env:
	@test -f ./env_TRUST.sh || { echo >&2 "Error: env_TRUST.sh not found. Please run ./configure first."; exit 1; }



