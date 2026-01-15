#!/bin/bash
[ "$TRUST_STDCPP" = c++14 ] && exit 0
[ "$TRUST_STDCPP" = c++17 ] && exit 0
archive=$TRUST_ROOT/externalpackages/kokkos/arborx-2.0.1.tar.gz # C++ 20 

build_dir=$TRUST_ROOT/build/arborx
KOKKOS_ROOT_DIR=$TRUST_ROOT/lib/src/LIBKOKKOS
# Log file of the process:
log_file=$TRUST_ROOT/arborx_compile.log

if [ ! -f $KOKKOS_ROOT_DIR/include/ArborX/ArborX_Config.hpp ]
then
    echo "# Installing `basename $archive` ..."
    # Sub shell to get back to correct dir:
    (  
      mkdir -p $build_dir;cd $build_dir
      tar xzf $archive || exit -1
      src_dir=$build_dir/`ls $build_dir | grep -i arborx`

      BUILD_TYPES="Release `[ "$build_debug" = "1" ] && echo Debug`"
      BUILD_TYPES="Release Debug"
      for CMAKE_BUILD_TYPE in $BUILD_TYPES
      do
        rm -rf BUILD;mkdir -p BUILD;cd BUILD
        KOKKOS_INSTALL_DIR=$KOKKOS_ROOT_DIR/$TRUST_ARCH`[ $CMAKE_BUILD_TYPE = Release ] && echo _opt`
        CMAKE_OPT="-DCMAKE_CXX_EXTENSIONS=Off"
        if [ "$TRUST_USE_CUDA" = 1 ]
        then
           CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_COMPILER=$KOKKOS_INSTALL_DIR/bin/nvcc_wrapper"
        fi
        if [ "$TRUST_USE_ROCM" = 1 ]
        then
           CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_COMPILER=hipcc"
        fi   
        CMAKE_OPT=$CMAKE_OPT" -DCMAKE_CXX_FLAGS=-fPIC"
        CMAKE_OPT=$CMAKE_OPT" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DCMAKE_INSTALL_PREFIX=$KOKKOS_INSTALL_DIR -DCMAKE_INSTALL_LIBDIR=lib64"
        CMAKE_OPT=$CMAKE_OPT" -DKokkos_ROOT=$KOKKOS_INSTALL_DIR"
        CMAKE_OPT=$CMAKE_OPT" -DKokkos_DIR=$KOKKOS_INSTALL_DIR/lib64/cmake/Kokkos" # Needed on salacia
        CMAKE_OPT=$CMAKE_OPT" -DARBORX_ENABLE_MPI=ON" # -DARBORX_ENABLE_GPU_AWARE_MPI=ON
        #CMAKE_OPT=$CMAKE_OPT" -DARBORX_ENABLE_EXAMPLES=ON" # Need Boost ?
        echo CMAKE_OPT=$CMAKE_OPT | tee $log_file
        cmake $src_dir $CMAKE_OPT 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when configuring Arborx (CMake) - look at $log_file" && exit -1

        # Build
        make -j$TRUST_NB_PHYSICAL_CORES install 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when compiling Arborx - look at $log_file" && exit -1
        echo "Arborx $CMAKE_BUILD_TYPE installed under $KOKKOS_ROOT_DIR"
        cd ..
      done
      # Clean build:
      rm -rf $build_dir $log_file
    )
else
    echo "# Arborx: already installed. Doing nothing."
fi

