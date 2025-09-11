#!/bin/bash
# Support old and new TU format
TU=$1 && [ ! -f $TU ] && echo "$TU not found!" && exit -1
[ "$2" = "" ] && echo "Specify an option" && echo "TU.sh file [-dt] [-total] [-its]" && exit -1
if [ $2 = "-dt" ]
then
   # Time step
   awk '/Secondes/ && /pas de temps/ {val=$NF} /Average time per time step/ {val=$NF} END {print val}' $TU
elif [ $2 = "-total" ]
then
   # Total time
   awk '/Total time for the whole computation/ {val=$NF} /Temps total/ {val+=$NF} END {print val}' $TU
elif [ $2 = "-its" ]
then
   # Iterations solver
   awk '/Iterations/ && /solveur/ {val=$NF} /Average number of iteration of the linear solver per call/ {val=$NF} END {print val}' $TU
elif [ $2 = "-ndt" ]
then
   # Number of time steps
   awk '/Number of time steps:/ {val=$NF} END {print val}' $TU
elif [ $2 = "-solver" ]
then
   # Solver Ax=b
   awk '/Linear solver/ {val=$6} END {print val}' $TU
else
   echo "Option $2 unknown"
   exit -1
fi
