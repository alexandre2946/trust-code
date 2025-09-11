# Kokkos
lib=$(TRUST_LIB)/src/LIBKOKKOS/$(TRUST_ARCH)$(OPT)/lib64/libkokkoscore.a

all: $(lib) 
$(lib): install_kokkos.sh install_kokkos-kernels.sh install_arborx.sh makefile
	@make clean
	./install_kokkos.sh && ./install_kokkos-kernels.sh && ./install_arborx.sh

clean:
	rm -rf $(TRUST_ROOT)/build/kokkos ${TRUST_ROOT}/build/kokkos-kernels $(TRUST_LIB)/src/LIBKOKKOS
