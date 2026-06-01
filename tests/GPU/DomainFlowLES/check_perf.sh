#!/bin/bash
check()
{
   [ "$2" = -nsys ] && exit
   #grep "AMG solver" $1.out_err 2>/dev/null
   if [ ! -f $1.TU ] || [ "`grep 'Stopping processes.' $1.out_err`" = "" ]
   then
      echo "Performance is KO for $1 on $2 : case does not run!"
      exit -1
   fi
   TU=$1.TU
   TU_REF=$1.TU.$2`[ "$np" != "" ] && echo x$np`
   if [ ! -f $TU_REF ]
   then
      mv -f $TU $TU_REF && [ "$TRUST_SCM" = 1 ] && git add $TU_REF
      echo "Creating new reference $TU_REF"
      exit
   fi 
   ref=`TU.sh $TU_REF -dt`
   new=`TU.sh $TU     -dt`
   echo $ref $new | awk '// {if (2*($2-$1)/($1+$2)>0.05) {exit 1}}' # On verifie qu'on ne depasse pas +5% de la performance
   err=$?
   ecart=`echo $ref $new | awk '// {printf("%2.1f\n",200*($2-$1)/($1+$2))}'`
   [ "$UPDATE_REFERENCE" = 1 ] && cp -f $TU $TU_REF && echo "Forced update of $TU_REF"
   if [ $err = 1 ]
   then
      #sdiff -w 200 $TU_REF $TU
      echo "Performance is KO ($ecart%) for $1 on $2 !"
   else
      echo "Performance is OK ($ecart%) $new s < $ref s (reference) for $1 on $2"
      if [ `echo "$ecart<-0.99" | bc -l` = 1 ]
      then
         echo "$TU_REF is updated !" && cp $TU $TU_REF
      else
         # Non regression faiblement testee sur le nombre d'iterations des solveurs
         its_ref=`TU.sh $TU_REF -its`
         its_new=`TU.sh $TU     -its`
         [ $its_ref != $its_new ] && echo "Solver convergence is different ($its_ref != $its_new) ! Possible regression..." && exit -1
      fi
   fi
}
run()
{
   gpu=$1 && [ $gpu = -nsys ] && nsys=$gpu
   np=$2
   [ "$3" != "" ] && jdd=$3
   if [ "$np" = "" ] || [ "$np" = 1 ]
   then
      np=""
   else  
      #echo "Partition for $np MPI..." 
      make_PAR.data $jdd $np 1>/dev/null 2>&1
      jdd=PAR"_"$jdd
   fi
   rm -f $jdd.TU
   # Try to mitigate variablity by setting exclusive mode on GPU (firefox, slack, edge, chrome, use device !)
   [ "$np" = "" ] && [ "$TRUST_WITHOUT_HOST" = 1 ] && [ "`hostname`" = is157091 ] && set_EXCLUSIVE_PROCESS=`sudo ls 2>/dev/null`
   [ "$set_EXCLUSIVE_PROCESS" != "" ] && sudo nvidia-smi -c EXCLUSIVE_PROCESS 1>/dev/null
   trust $nsys $jdd $np 1>$jdd.out_err 2>&1
   [ "$set_EXCLUSIVE_PROCESS" != "" ] && sudo nvidia-smi -c DEFAULT 1>/dev/null         
   check $jdd $gpu
}
# Jeu de donnees
jdd=`pwd`
jdd=`basename $jdd`
cp $jdd.data $jdd"_BENCH".data
jdd=$jdd"_BENCH"
for pattern in Nombre_de_Noeuds raffiner_isotrope
do
   sed -i "/"$pattern"/s/\/\*//g" $jdd.data
   sed -i "/"$pattern"/s/\*\///g" $jdd.data
done
# Replace lml format by lata for faster IO
sed -i "1,$ s?format lml?format lata?g" $jdd.data
# Surcharge CLI: check_perf.sh [-nsys] [$np] [$jdd]
nsys="" && [ "$1" = -nsys ] && nsys="-nsys"
np=1 && [ "$1" != "" ] && np=$1 && shift
[ "$1" != "" ] && jdd=$1 && shift

# Liste des machines:
if [ "$nsys" != "" ]
then
   run $nsys
else
   # MPI Ranks:
   [ "$TRUST_USE_GPU" != 1 ] && np=`echo $TRUST_NB_PHYSICAL_CORES/4 | bc`
   # HOST:   
   HOST=${HOST%.intra.cea.fr} && [ "$HOST" = portable ] && HOST=is246827
   # ARCH:
   GPU_ARCH=""
   [ "$TRUST_USE_CUDA" = 1 ] && GPU_ARCH=_cc$TRUST_CUDA_CC
   [ "$TRUST_USE_ROCM" = 1 ] && GPU_ARCH=_$ROCM_ARCH
   # Run:
   run $HOST$GPU_ARCH $np
   if [ $HOST = topaze ] || [ $HOST = adastra ] || [ $HOST = jean-zay ] || [ $HOST = lumi ] # Test strong scalability on multi-node GPU
   then
      [ "`grep -i 'nb_parts 8' $jdd.data`" != "" ] && run $HOST$GPU_ARCH 8
   fi
fi
# clean
rm -f *.sauv *.lml *.sqlite *.nsys-rep
