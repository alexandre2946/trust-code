#!/bin/bash

# MEDCoupling charge par defaut est l'optim, possible de charger $TRUST_MEDCOUPLING_ROOT/env_debug.sh
. $TRUST_MEDCOUPLING_ROOT/env.sh

export PYTHONPATH=$TRUST_ROOT/Outils/ICoCo/ICoCo_src/share/swig/install/lib:$PYTHONPATH
#export LD_LIBRARY_PATH=$TRUST_ROOT/exec/python/lib:$LD_LIBRARY_PATH  # Not needed anymore ?

# if with PDI, we need libparaconf.1.dylib ... not from python but from pdi !
if [[ "$TRUST_PDI_ROOT" != "" && "$(uname)" == "Darwin" ]]; then
	if [[ "`echo $DYLD_LIBRARY_PATH | grep -i LIBPDI/lib`" = ""  ]]; then
		export DYLD_LIBRARY_PATH=$TRUST_PDI_ROOT/lib:$DYLD_LIBRARY_PATH
	fi
fi