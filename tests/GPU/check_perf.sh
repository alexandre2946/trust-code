#!/bin/bash
rm -f check_perf.log
for rep in `find */ -name check_perf.sh`
do
	rep=`dirname $rep`
	cd $rep
	echo "======================"
	echo $rep
	./check_perf.sh 2>&1 | tee -a ../check_perf.log
	grep GPU: $rep"_BENCH".TU 2>/dev/null
	cd - 1>/dev/null 2>&1
done
echo "File check_perf.log created"

