#!/bin/bash
# Scaling a mesh on several GPU
[ "$TRUST_ROOT" = "" ] && echo "TRUST_ROOT empty." && exit

# HOST:   
HOST=${HOST%.intra.cea.fr} && [ "$HOST" = portable ] && HOST=is246827
# ARCH:
GPU_ARCH=""
[ "$TRUST_USE_CUDA" = 1 ] && GPU_ARCH=_cc$TRUST_CUDA_CC
[ "$TRUST_USE_ROCM" = 1 ] && GPU_ARCH=_$ROCM_ARCH

ROOT=`pwd`
log=$ROOT/`basename $ROOT`_SCALING.$HOST$GPU_ARCH

echo -e "Config     [MTET]  [MDOF] TimeStep[s] Solver[s] [its] [ms/it] Kernels[s] RAM[GB] DRAM[GB] [MDOF/s] Solver Kernels Conv Diff Grad  Div" | tee $log
versions=cpu && [ "$TRUST_USE_GPU" = 1 ] && versions=gpu
for version in $versions
do
   # Nombre de GPUs testes:
   if [ $version = cpu ]
   then
      gpus="0"
   elif [ "$TRUST_WITHOUT_HOST" = 0 ]
   then
      gpus="1 2 3 4" && [ "$ROCM_ARCH" = gfx90a ] && gpus="1 2 3 4 5 6 7 8"
   elif [ "`hostname`" = petra ]
   then
      gpus="1 2"
   else
      gpus="1"
   fi   

   # Target problem sizes in MDOF (millions of degrees of freedom).
   # MDOF formula: 1.2*40*(Nx-1)*(Ny-1)*(Nz-1)/1e6
   # Nx, Ny, Nz are scaled uniformly from the base mesh in the .data file.
   
   # Read base Nombre_de_Noeuds from reference data file.
   # Prefer /* Nombre_de_Noeuds X Y Z */ commented template (scaling hint);
   # fall back to the last active Nombre_de_Noeuds line.
   if [ -f $ROOT/OpenMP_Iterateur.data ]
   then
      ref_data=$ROOT/OpenMP_Iterateur_BENCH_PETSc.data
      [ $version = gpu ] && [ "$TRUST_USE_CUDA" = 1 ] && ref_data=$ROOT/OpenMP_Iterateur_BENCH_AmgX.data
   else
      ref_data=$ROOT/`basename $ROOT`.data
   fi
   read Nx0 Ny0 Nz0 <<< $(awk '
      /\/\* Nombre_de_Noeuds [0-9]/          { cx=$3; cy=$4; cz=$5 }
      /^[[:space:]]*Nombre_de_Noeuds [0-9]/  { ax=$2; ay=$3; az=$4 }
      END { if (cx!="") print cx,cy,cz; else print ax,ay,az }
   ' $ref_data)
   [ -z "$Nx0" ] && echo "Error: Nombre_de_Noeuds not found in $ref_data" && exit 1
   echo "# Base mesh from $ref_data: Nombre_de_Noeuds $Nx0 $Ny0 $Nz0"

   for gpu in $gpus
   do
      mpis=$TRUST_NB_PHYSICAL_CORES && [ $version = gpu ] && mpis=$gpu
      mdof_target=0.1
      inc_mdof=0.1
      for mpi in $mpis
      do
         while [ 1 ]
         do
            # Compute new Nx Ny Nz matching the target MDOF via uniform scaling:
            #   alpha = cbrt( target_MDOF*1e6 / (1.2*40*(Nx0-1)*(Ny0-1)*(Nz0-1)) )
            #   Ni = 1 + round( alpha * (Ni0-1) )
            read Nx Ny Nz <<< $(awk -v Nx0=$Nx0 -v Ny0=$Ny0 -v Nz0=$Nz0 -v target=$mdof_target \
               'BEGIN {
                  base  = 1.2 * 40 * (Nx0-1) * (Ny0-1) * (Nz0-1)
                  alpha = (target * 1e6 / base) ^ (1.0/3.0)
                  Nx = int(0.5 + 1 + alpha * (Nx0-1)); if (Nx < 2) Nx = 2
                  Ny = int(0.5 + 1 + alpha * (Ny0-1)); if (Ny < 2) Ny = 2
                  Nz = int(0.5 + 1 + alpha * (Nz0-1)); if (Nz < 2) Nz = 2
                  print Nx, Ny, Nz
               }')
            jdd=$mpi"_"$gpu"_"${Nx}x${Ny}x${Nz}
            mkdir -p $ROOT/scaling/$jdd && cd $ROOT/scaling/$jdd
            # Run ?
            run=1 && [ -f $jdd.out_err ] && run=0
            [ $run = 1 ] && echo "$jdd (target ${mdof_target} MDOF) ..."
            # Creation data
            if [ -f $ROOT/OpenMP_Iterateur.data ]
            then
               [ $version = gpu ] && [ "$TRUST_USE_CUDA" = 1 ] && cp $ROOT/OpenMP_Iterateur_BENCH_AmgX.data $jdd.data
               [ $version = cpu ] && cp $ROOT/OpenMP_Iterateur_BENCH_PETSc.data $jdd.data
            else
               cp $ROOT/`basename $ROOT`.data $jdd.data
            fi
            # Replace /* Nombre_de_Noeuds X Y Z */ comment template with new active line
            sed -i "s|/\* Nombre_de_Noeuds [0-9]\+ [0-9]\+ [0-9]\+ \*/|Nombre_de_Noeuds $Nx $Ny $Nz|" $jdd.data
            # Fallback: replace active Nombre_de_Noeuds matching base values (no comment template)
            sed -i "s/Nombre_de_Noeuds[[:space:]]\+$Nx0[[:space:]]\+$Ny0[[:space:]]\+$Nz0/Nombre_de_Noeuds $Nx $Ny $Nz/" $jdd.data
	    # Important change to rtol for scaling !
	    sed -i "1,$ s? seuil ? rtol ?" $jdd.data
	    sed -i "1,$ s? atol ? rtol ?" $jdd.data
	    # Important change lml en cgns
	    sed -i "1,$ s? lml? cgns?" $jdd.data
	    # Suppress post-processing:
	    sed -i '/Postraitement/d' $jdd.data
            # Decoupage
            [ $run = 1 ] && [ $mpi != 1 ] && (make_PAR.data $jdd $mpi 1>/dev/null 2>&1;cp PAR_$jdd.data $jdd.data)
            # Calcul
            [ $run = 1 ] && (trust $jdd $mpi -ksp_view -journal=0 1>$jdd.out_err 2>&1 || (rm -f *.TU;echo "Error:See "`pwd`/$jdd.out_err))
            [ "`grep 'Arret des process' $jdd.out_err`" = "" ] && break
            # Analyse
	    i=0 && [ "$HOST" = adastra ] && i=1
	    hram=`awk -v i=$i '/RAM taken/ {if ($(13+i)>RAM) RAM=$(13+i)} END {print 0.1*int(0.01*RAM)}' $jdd.out_err`
	    dram=`awk -v i=$i '/RAM allocated on a GPU/ {if ($(1+i)>RAM) RAM=$(1+i)} END {print RAM}' $jdd.out_err`
            row=`awk '/Order of the PETSc matrix/ {print $7;exit}' $jdd.out_err`
            faces=`awk '/Total number of faces/ {print $NF;exit}' $jdd.out_err`
            elems=`awk '/Total number of elements/ {printf($NF);exit}' $jdd.out_err`
            # No better to use dof=row
            dof=$row
            #dof=`echo 1*$faces | bc -l` # En VDF
            #dof=`echo 3*$faces | bc -l` # En VEF
            dt=`TU.sh $jdd.TU -dt`
            its=`TU.sh $jdd.TU -its`            
            awk -v mpi=$mpi -v gpu=$gpu -v elems=$elems -v row=$row -v dof=$dof -v hram=$hram -v dram=$dram -v dt=$dt -v its=$its '\
	    BEGIN {config=mpi"MPI"(gpu==0?"":"+"gpu"GPU");mdof=dof/1e6;mtet=elems/1e6} \
            /Linear solver/       {ts=$6;b=dt-ts;ls=mdof/ts} \
            /Convection operator/ { conv=mdof/$4*$8 } \
            /Diffusion operator/  { diff=mdof/$4*$8 } \
            /Gradient operator/   { grad=mdof/$4*$8 } \
            /Divergence operator/ { dive=mdof/$4*$8 } \
            /Kernels:/            { ks=mdof/$3 } \
            END {printf("%s %7.3f %7.3f %11.3f %9.3f %5d %7.1f %10.3f %7.1f %8.1f %8.1f %6.1f %6.1f %4d %4d %4d %4d\n", \
	             config, mtet, mdof, dt, ts, its, 1000*ts/its,    b,   hram, dram,  mdof/dt,  ls,  ks, conv, diff, grad, dive)}' $jdd.TU | tee -a $log
	    # Clean
	    rm -f *.sauv *.xyz *.*lata*	*.cgns* *.face *.son *.lml    
            mdof_target=`echo $mdof_target+$inc_mdof | bc -l`
         done
      done
   done    
done
echo "$log created."
python3 ./plot_scaling.py
display JEL_bous_SCALING.png
