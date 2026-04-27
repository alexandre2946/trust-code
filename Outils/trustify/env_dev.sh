#!/bin/bash

# Source environment for trustify in DEVELOPPER MODE:
src_dir=$TRUST_ROOT/Outils/trustify/src

# Check resulting PYTHONPATH is OK to avoid ugly mess up with real installed version!
PYTHONPATH=${src_dir}:${PYTHONPATH} pth=`python -c "import trustify;print(trustify.__file__)"`
echo "Path to active trustify is $pth"

if ! [[ "$pth" == *"trustify/src"* ]]; then
    echo "This does not contain trustify/src!! There is probably sth wrong with your PYTHONPATH ... exiting."
else
    export PYTHONPATH=${TRUSTIFY_ROOT_DIR}:${PYTHONPATH}
    echo " -> OK"
fi

