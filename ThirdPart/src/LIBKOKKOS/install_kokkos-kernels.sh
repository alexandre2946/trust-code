#!/bin/bash
# Kokkos-kernels uniquement sur GPU:
[ "$TRUST_USE_GPU" != 1 ] && exit 0

# Kokkos-kernels:
archive=$TRUST_ROOT/externalpackages/kokkos/kokkos-kernels-4.7.00.tar.gz
build_dir=$TRUST_ROOT/build/kokkos-kernels
KOKKOS_ROOT_DIR=$TRUST_ROOT/lib/src/LIBKOKKOS
# Log file of the process:
log_file=$TRUST_ROOT/kokkos-kernels_compile.log

if [ ! -f $KOKKOS_ROOT_DIR/lib64/libkokkoskernels.a ]
then
    echo "# Installing `basename $archive` ..."
    # Sub shell to get back to correct dir:
    (  
      mkdir -p $build_dir;cd $build_dir
      tar xzf $archive || exit -1
      src_dir=$build_dir/`ls $build_dir | grep kokkos-kernels`
      # Set this flag to 1 to have Kokkos compiled/linked in Debug mode for $exec_debug or when developping on GPU:
      if [ $HOST = $TRUST_HOST_ADMIN ] || [ "$TRUST_USE_KOKKOS_SIMD" = 1 ]
      then
         build_debug=1
      else
         build_debug=$TRUST_ENABLE_KOKKOS_DEBUG
      fi
      build_debug=0 # Cause issue when building unit
      BUILD_TYPES="Release `[ "$build_debug" = "1" ] && echo Debug`"
      for CMAKE_BUILD_TYPE in $BUILD_TYPES
      do
        rm -rf BUILD;mkdir -p BUILD;cd BUILD
        echo "Building under `pwd` ..."
        CMAKE_OPT=""
        CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_FLAGS=-fPIC"
	if [ "$TRUST_USE_CUDA" = 1 ]
	then
	   CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_COMPILER=$TRUST_CC_BASE_EXTP"
           # Shit on NVHPC single CUDA bundle: If $CUDA_ROOT/bin/nvcc or $CUDA_ROOT/compilers/bin/nvcc: KK cmake fails...
           [ "$CUDA_VERSION" = "" ] && CMAKE_OPT="-DCMAKE_CUDA_COMPILER=`find $CUDA_ROOT -name nvcc`"
           # What a shit to find cublas/cusolver/cusparse !
           cublas=`find $CUDA_ROOT/../math_libs -name cublas.h`
           rep=`dirname $cublas`
           rep=`dirname $rep`
	   CMAKE_OPT=$CMAKE_OPT" -DCUDAToolkit_ROOT=$rep -DKokkosKernels_CUBLAS_ROOT=$rep -DKokkosKernels_CUSPARSE_ROOT=$rep -DKokkosKernels_CUSOLVER_ROOT=$rep"
        elif [ "$TRUST_USE_ROCM" = 1 ]
        then
	   CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_COMPILER=$TRUST_CC_BASE"
	fi
        CMAKE_INSTALL_PREFIX=$KOKKOS_ROOT_DIR/$TRUST_ARCH`[ $CMAKE_BUILD_TYPE = Release ] && echo _opt`
	export Kokkos_DIR=$CMAKE_INSTALL_PREFIX/lib64/cmake
        CMAKE_OPT=$CMAKE_OPT" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX -DCMAKE_INSTALL_LIBDIR=lib64"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_ENABLE_TESTS=OFF"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_DOUBLE=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_COMPLEX_DOUBLE=ON" # OFF ? 
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_ORDINAL_INT=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_ORDINAL_INT64_T=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_OFFSET_INT=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_OFFSET_SIZE_T=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_LAYOUTLEFT=ON"
        #CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_ADD_DEFAULT_ETI=ON"
        # Configure
	echo $src_dir $CMAKE_OPT
	# cmake --debug-find $src_dir $CMAKE_OPT 1>log 2>&1
        cmake $src_dir $CMAKE_OPT 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when configuring Kokkos (CMake) - look at $log_file" && exit -1

        # Build
        make -j$TRUST_NB_PHYSICAL_CORES install 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when compiling Kokkos-kernels - look at $log_file" && exit -1
        echo "Kokkos-kernels $CMAKE_BUILD_TYPE installed under $CMAKE_INSTALL_PREFIX"
        cd ..
      done
      # Clean build:
      rm -rf $build_dir $log_file
    )
else
    echo "# Kokkos-kernels: already installed. Doing nothing."
fi

