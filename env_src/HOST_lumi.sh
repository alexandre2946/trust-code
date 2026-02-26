#!/bin/bash
##################################
# Variables for configure script #
##################################
define_modules_config()
{
   env=$TRUST_ROOT/env/machine.env
   # qstat inexistente sur les dernieres machines du CCRT/TGCC
   echo "Command qstat created on $HOST"
   cp $TRUST_ROOT/bin/KSH/qstat_wrapper $TRUST_ROOT/bin/KSH/qstat
   # Initialisation de l environnement module $MODULE_PATH
   echo "source /etc/profile" >> $env
   #
   # Load modules
   if [ "$TRUST_USE_ROCM" = 1 ]
   then
      # ROCM_ARCH par defaut:
      if [ "$ROCM_ARCH" = "" ]
      then
         export ROCM_ARCH=gfx90a # MI250X
      else
         echo "$ROCM_ARCH not supported on lumi!"
      fi
      module="PrgEnv-gnu/8.5.0 craype-accel-amd-$ROCM_ARCH LUMI/24.03 partition/G buildtools/24.03 CrayEnv rocm/6.4.4 gnuplot/5.4.10-cpeGNU-24.03"
   else
      echo "Not configured." && exit -1 
   fi
   module=$module" cray-python/3.11.7" # Pour -without-conda
   #module=$module" swig" # Pour -without-conda
   #
   echo "# Module $module detected and loaded on $HOST."
   echo "module purge 1>/dev/null 2>&1" >> $env
   echo "module load $module 1>/dev/null || exit -1" >> $env
   echo "export TRUST_USE_SACCT=1" >> $env # Energy data 
   #echo "PATH=\$CRAY_MPICH_PREFIX/bin:\$PATH"  >> $env # Pour trouver mpicxx
   . $env
}

##############################
# Variables for trust script #
##############################
define_soumission_batch()
{
   # Documentation lumi: https://docs.lumi-supercomputer.eu/
   soumission=1
   if [ "$prod" = 1 ] || [ $NB_PROCS -gt 256 ]
   then  
      [ $NB_PROCS -le 32 ] && queue=small-g    && qos="" && cpu=4320 # 3 days 
      [ $NB_PROCS -gt 32 ] && queue=standard-g && qos="" && cpu=2880 # 2 days
   else
      [ $NB_PROCS -le 256 ] && queue=dev-g      && qos="" && cpu=30 # 30 mn
      [ $NB_PROCS -le 128 ] && queue=dev-g      && qos="" && cpu=60 # 1h
      [ $NB_PROCS -le 64 ]  && queue=dev-g      && qos="" && cpu=120 # h2
   fi
   project=project_465002428
   if [ "$gpu" = 1 ]
   then
      if [ "$ROCM_ARCH" = gfx90a ] # Partition MI250X (BW: 1600 GB/s)
      then
          #constraint=MI250
          ntasks=64 # Node cores
          cpus_per_task=7 # 56 cores and not 64 cores available on each node
          gpu_per_node=8
      fi
      noeuds=`echo "1+($NB_PROCS-1)/$gpu_per_node" | bc`
      ntasks_per_node=$gpu_per_node && [ $NB_PROCS -lt $gpu_per_node ] && ntasks_per_node=$NB_PROCS
      srun_options="-c $cpus_per_task --gpus-per-node=$ntasks_per_node --ntasks-per-node=$ntasks_per_node $TRUST_ROOT/env_src/lumi_select_gpu.sh" # -cpu-bind=map_cpu:49,57,17,25,1,9,33,41
   else
      echo "ToDo"
   fi
   node=1 # --exclusive
   ram=480g # RAM per node (512 - 32)
   USE_MPIRUN=1
   mpirun="srun $srun_options"
   sub=SLURM
}

