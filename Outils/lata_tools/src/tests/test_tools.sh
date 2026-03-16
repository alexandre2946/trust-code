#!/bin/bash

echo ""
echo "-------------------------------"
echo "Check if compare_lata works ..."
echo "-------------------------------"

if [ "$TRUST_CGNS_ROOT" != "" ]; then
	echo ""
	echo "###  with CGNS files ..."
	cp $TRUST_ROOT/Outils/lata_tools/src/tests/upwind.cgns .
	$TRUST_ROOT/exec/lata_tools/bin/compare_lata upwind.cgns upwind.cgns --max_delta || exit -1
	echo ""
fi 

echo "###  with LATA files 32b ..."
tar xzf $TRUST_ROOT/Outils/lata_tools/src/tests/vdf32.tar.gz
$TRUST_ROOT/exec/lata_tools/bin/compare_lata vdf32/vdf.lata vdf32/vdf.lata || exit -1
echo ""

echo "###  with LATA files 64b ..."
tar xzf $TRUST_ROOT/Outils/lata_tools/src/tests/upwind64.tar.gz
$TRUST_ROOT/exec/lata_tools/bin/compare_lata upwind64/upwind.lata upwind64/upwind.lata || exit -1
echo ""

echo "###  with LML files ..."
tar xzf $TRUST_ROOT/Outils/lata_tools/src/tests/upwind.tar.gz
$TRUST_ROOT/exec/lata_tools/bin/compare_lata upwind.lml upwind.lml || exit -1
echo ""

if [ "$TRUST_MEDCOUPLING_ROOT" != "" ]; then
  echo "###  with MED files ..."
  cp $TRUST_ROOT/Outils/lata_tools/src/tests/upwind_0000.med .
  $TRUST_ROOT/exec/lata_tools/bin/compare_lata upwind_0000.med upwind_0000.med || exit -1
fi
echo ""

echo "###  with FORT21 files ..."
cp $TRUST_ROOT/Outils/lata_tools/src/tools/FORT21 .
$TRUST_ROOT/exec/lata_tools/bin/compare_lata FORT21 FORT21 --same_mesh || exit -1
echo ""

echo "--------------------------------"
echo "Check if lata_analyzer works ..."
echo "--------------------------------"

if [ "$TRUST_CGNS_ROOT" != "" ]; then
	echo ""
	echo "###  with CGNS files converted to single_lata ..."
	$TRUST_ROOT/exec/lata_tools/bin/lata_analyzer upwind.cgns writelata_convert=cgns2lata || exit -1
	echo ""
fi

echo "###  with LML files converted to single_lata ..."
$TRUST_ROOT/exec/lata_tools/bin/lata_analyzer upwind.lml writelata_convert=test_upwind || exit -1
echo ""

echo "###  with LATA files 32b converted to single_lata ..."
$TRUST_ROOT/exec/lata_tools/bin/lata_analyzer vdf32/vdf.lata write_singlelata=test_vdf32 || exit -1
echo ""

echo "###  with LATA files 64b converted to single_lata ..."
$TRUST_ROOT/exec/lata_tools/bin/lata_analyzer upwind64/upwind.lata write_singlelata=test_upwind64 || exit -1
echo ""

echo "###  with FORT21 files converted to single_lata ..."
$TRUST_ROOT/exec/lata_tools/bin/lata_analyzer FORT21 write_singlelata=testf21 || exit -1
nc=$(grep CHAMP testf21.lata | wc| awk '{print $1}')
ref=62
[ $nc -ne $ref ] && echo invalid number of CHAMP $nc != $ref&& exit -2

if [ "$TRUST_MEDCOUPLING_ROOT" != "" ]; then

	echo ""
	echo "------------------------------------"
	echo "Check if TRUST_Post_Loader works ..."
	echo "------------------------------------"

	if [ "$TRUST_CGNS_ROOT" != "" ]; then
		echo ""
		echo "###  with CGNS files ..."
		$TRUST_ROOT/exec/lata_tools/bin/test_TRUST_Post_Loader upwind.cgns PRESSION_ELEM_dom_ELEM || exit -1
		echo ""
	fi

	echo "###  with LML files ..."
	$TRUST_ROOT/exec/lata_tools/bin/test_TRUST_Post_Loader upwind.lml PRESSION_ELEM_dom || exit -1
	echo ""

	echo "###  with LATA files 64b ..."
	$TRUST_ROOT/exec/lata_tools/bin/test_TRUST_Post_Loader upwind64/upwind.lata PRESSION_ELEM_dom || exit -1
	echo ""

	if [ "$TRUST_CGNS_ROOT" != "" ]; then
		echo ""
		echo "------------------------------------"
		echo "Check if TRUST_Post_Loader works in python ..."
		echo "------------------------------------"

		echo ""
		echo "###  with CGNS files ..."
		source $TRUST_ROOT/env_for_python.sh
		cp $TRUST_ROOT/Outils/lata_tools/src/tests/script_loader.py .
		python script_loader.py || exit -1
		echo ""
	fi
fi