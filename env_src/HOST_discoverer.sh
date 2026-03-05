#!/bin/bash
##################################
# Variables for configure script #
##################################
define_modules_config()
{
   env=$TRUST_ROOT/env/machine.env
   # Initialisation de l environnement module $MODULE_PATH si pas disponible:
   module -v 1>/dev/null 2>&1 || echo $echo "source /etc/profile" >> $env
   #
   # Load modules
   if [ "$TRUST_USE_CUDA" = 1 ]
   then
      TRUST_CUDA_CC=90
      module="nvidia/hpcsdk/nvhpc/25.1 cmake/3/3.31.4"
   else
      echo "Not supported yet"
      exit 1
   fi
   #
   echo "# Module $module detected and loaded on $HOST."
   echo "module purge 1>/dev/null" >> $env
   echo "module load $module 1>/dev/null || exit -1" >> $env
   . $env
   # Creation wrapper qstat -> squeue
   echo "#!/bin/bash
squeue" > $TRUST_ROOT/bin/qstat
   chmod +x $TRUST_ROOT/bin/qstat
}

##############################
# Variables for trust script #
##############################
define_soumission_batch()
{
   soumission=1
   # https://docs.discoverer.bg/resource_overview.html
   if [ "$gpu" = 1 ]
   then
      ntasks=112 # number of cores max
      if [ $NB_PROCS -le 8 ]
      then
         gpus_per_node=$NB_PROCS
      else
         gpus_per_node=8
	 node=1 # --exclusive
      fi	 
      cpus_per_task=14
      # 32 GPUs max
      qos=ehpc-dev-2026d03-190 # && cpu=4320 && [ "$prod" != 1 ] && [ $NB_PROCS -le 32 ] && qos=acc_debug && cpu=120 
      project="ehpc-dev-2026d03-190"
   else
      echo "Not supported yet." && exit 1
   fi
   mpirun="srun -n \$SLURM_NTASKS"
   sub=SLURM
}

