trg=$(TRUST_ROOT)/exec/lata_tools/bin/compare_lata
tools_src=src
trust_kern=$(TRUST_ROOT)/src/Kernel

# All the sources for lata tools ... TODO - is there a better way using CMake computed dpys?
SRCS=$(wildcard $(tools_src)/trust_commun/*.cpp $(tools_src)/trust_commun/*.h $(tools_src)/tools/*.i $(tools_src)/tools/*.h $(tools_src)/tools/*.cpp)
trust_ut=$(trust_kern)/Utilitaires
EXTRA_SRCS=$(trust_ut)/Static_Int_Lists.cpp $(trust_ut)/Nom.cpp $(trust_ut)/Noms.cpp $(trust_ut)/Motcle.cpp $(trust_ut)/AbstractIO.cpp $(trust_ut)/Entree.cpp $(trust_ut)/EFichier.cpp $(trust_ut)/Entree_Fichier_base.cpp $(trust_ut)/Device.cpp $(trust_ut)/DeviceMemory.cpp $(trust_ut)/Objet_U.cpp
trust_ge=$(trust_kern)/Geometrie
EXTRA_SRCS+=$(trust_ge)/Octree_Double.cpp $(trust_ge)/Octree_Int.cpp $(trust_ge)/Connectivite_som_elem.cpp 
trust_mt=$(trust_kern)/Math
EXTRA_SRCS+=$(trust_mt)/TRUSTArray.cpp $(trust_mt)/TRUSTVect_tools.cpp $(trust_mt)/Array_base.cpp $(trust_mt)/ArrOfBit.cpp
# Add links to .h .hpp .tpp useful for plugin lata
EXTRA_SRCS+=$(trust_ut)/AbstractIO.h $(trust_ut)/arch.h $(trust_ut)/DecBaseInst.h $(trust_ut)/Declare_Base.h $(trust_ut)/Declare_Inst.h $(trust_ut)/Device.h $(trust_ut)/DeviceMemory.h $(trust_ut)/EFichier.h $(trust_ut)/Entree_Fichier_base.h $(trust_ut)/Entree.h $(trust_ut)/EntreeSortie.h $(trust_ut)/generique.h $(trust_ut)/Journal.h $(trust_ut)/kokkos++.h $(trust_ut)/Motcle.h $(trust_ut)/Nom.h $(trust_ut)/Noms.h $(trust_ut)/Objet_U.h $(trust_ut)/Process.h $(trust_ut)/Separateur.h $(trust_ut)/Sortie.h $(trust_ut)/Static_Int_Lists.h $(trust_ut)/TRUST_Vector.h $(trust_ut)/TRUST_Version.h $(trust_ut)/View_Types.h $(trust_ut)/TRUST_Ref.h

EXTRA_SRCS+=$(trust_mt)/Array_base.h $(trust_mt)/ArrOfBit.h $(trust_mt)/Double.h $(trust_mt)/MD_Vector_base.h $(trust_mt)/MD_Vector.h $(trust_mt)/MD_Vector_mono.h $(trust_mt)/MD_Vector_seq.h $(trust_mt)/MD_Vector_std_tools.tpp $(trust_mt)/MD_Vector_tools.h $(trust_mt)/span.hpp $(trust_mt)/TRUSTArray_device.tpp $(trust_mt)/TRUSTArray.h $(trust_mt)/TRUSTArray_kokkos.tpp $(trust_mt)/TRUSTArrays.h $(trust_mt)/TRUSTArray_tools.tpp $(trust_mt)/TRUSTArray.tpp $(trust_mt)/TRUSTList.tpp $(trust_mt)/TRUSTTab.h $(trust_mt)/TRUSTTab_parts.tpp $(trust_mt)/TRUSTTabs_forward.h $(trust_mt)/TRUSTTabs.h $(trust_mt)/TRUSTTab_tools.tpp $(trust_mt)/TRUSTTab.tpp $(trust_mt)/TRUSTTrav.h $(trust_mt)/TRUSTTravPool.h $(trust_mt)/TRUSTVect.h $(trust_mt)/TRUSTVect_tools.tpp $(trust_mt)/TRUSTVect.tpp $(trust_mt)/TVAlloc.h
EXTRA_SRCS+=$(trust_ge)/Connectivite_som_elem.h $(trust_ge)/Octree_Double.h $(trust_ge)/Octree_Double.tpp $(trust_ge)/Octree_Int.h $(trust_ut)/Schema_Comm_Vecteurs.h
all: $(trg)

# Retrieve TRUST_DEFINES variable:
#OPT=_opt
include $(TRUST_ENV)/make.$(TRUST_ARCH_CC)$(OPT)
#include $(TRUST_ENV)/make.$(TRUST_ARCH_CC)
 
$(trg): Makefile $(tools_src)/CMakeLists.txt $(tools_src)/tools/CMakeLists.txt $(tools_src)/scripts/CMakeLists.txt $(tools_src)/build.sh Makefile $(TRUST_MEDCOUPLING_ROOT)/$(TRUST_ARCH)_opt/lib/libmedcoupling.a $(SRCS) $(EXTRA_SRCS) $(trust_ut) $(trust_ge) $(trust_mt)
	(cd $(tools_src);env TRUST_DEFINES="$(Defines)" EXTRA_SRCS="$(EXTRA_SRCS)" ./build.sh)
	(cd $(TRUST_ROOT)/exec;ln -sf lata_tools/bin/lata_analyzer .)
	(cd $(TRUST_ROOT)/exec;ln -sf lata_tools/bin/compare_lata .)
	(cd $(TRUST_ROOT)/exec;ln -sf lata_tools/bin/compare_lata compare_cgns)
	(cd $(TRUST_ROOT)/exec;ln -sf $(TRUST_ROOT)/bin/KSH/compare_lata_gz compare_cgns_gz)
	(cd $(TRUST_ROOT)/exec/lata_tools/;ln -sf ../../Outils/lata_tools/src/tools/env.sh .)

clean:
	rm -rf $(TRUST_ROOT)/build/Outils/lata_tools $(TRUST_ROOT)/exec/lata_tools $(TRUST_ROOT)/exec/lata_analyzer $(TRUST_ROOT)/exec/compare_cgns_gz $(TRUST_ROOT)/exec/compare_cgns $(TRUST_ROOT)/exec/compare_lata

write_extra_src:   # Used by LATA Plugin build in Outils/VisIt/plugin/lata
	@echo $(EXTRA_SRCS) > $(OUTFILE)

.PHONY: clean write_extra_src
