#!/bin/bash

# Some clean due to the rename of lata2dx to lata_tools
if [ -d $TRUST_ROOT/exec/lata2dx_install ]
then
   rm -rf $TRUST_ROOT/exec/lata2dx_install $TRUST_ROOT/exec/lata2dx_reactor
   unlink $TRUST_ROOT/exec/lata2dx
   unlink $TRUST_ROOT/exec/compare_lata
fi

ORG=`pwd`
# Something changed -> we delete all
# Sometime TRUST_CC changes and produces recursive cmake so the last condition IS important
Build=$TRUST_ROOT/build/Outils/lata_tools
if [ build.sh -nt $Build ] || [ CMakeLists.txt -nt $Build ] || [ "`grep $TRUST_CC $Build/CMakeCache.txt`" = "" ]
then
   rm -r -f $Build
fi
mkdir -p $Build
cd $Build
# Better detection of SWIG on Ubuntu 16
SWIG_EXECUTABLE=`type -p swig`
# PL inutile de specifier swig (dans PATH normalement):
build_mode="Release"
if [ "x$OPT" == "x" ]; then
   echo "*****************************************"
   echo "lata_tools will be built in DEBUG mode!!"
   echo "*****************************************"
   build_mode="Debug"
fi
if [ "$TRUST_CC_BASE_EXTP" != "" ] && [ "`basename $TRUST_CC_BASE`" != crayCC ] # Sur adastra GNU ici ne marche pas... tout comme medcoupling
then
   export OMPI_CXX=$TRUST_CC_BASE_EXTP
   export OMPI_CC=$TRUST_cc_BASE_EXTP
   export MPICH_CXX=$OMPI_CXX
   export MPICH_CC=$OMPI_CC
fi
[ "`basename $TRUST_CC_BASE`" != crayCC ] && export CXXFLAGS="-Wno-inconsistent-missing-override -Wno-deprecated-declarations" # to suppress warnings from medcoupling

if [[ $(uname -s) == "Darwin" ]]
then
   PYTHON_VERSION=`python -c "import sys; print(str(sys.version_info.major)+'.'+str(sys.version_info.minor))"`
   export CXXFLAGS="$CXXFLAGS -I${TRUST_ROOT}/exec/python/include/python${PYTHON_VERSION}"
   export LDFLAGS="$LDFLAGS -Wl,-flat_namespace -Wl,-w"
fi

# CMake option + implicit use of TRUST_DEFINES and EXTRA_SRCS env var:
cmake_opt="-DCMAKE_BUILD_TYPE=$build_mode -DCMAKE_INSTALL_PREFIX=$TRUST_ROOT/exec/lata_tools -DCMAKE_C_COMPILER=$TRUST_cc "
cmake_opt="$cmake_opt -DCMAKE_CXX_COMPILER=$TRUST_CC -DSWIG_EXECUTABLE=$SWIG_EXECUTABLE "
cmake $ORG $cmake_opt || exit -1

make -j  $TRUST_NB_PROCS  install || exit -1

# We test everything
$TRUST_ROOT/Outils/lata_tools/src/tests/test_tools.sh

cd ..
rm -rf $Build
echo ""
