#!/bin/bash

# Helper functions
. $WM_PROJECT_USER_DIR/etc/helpFunctions.sh

setCdTrap "$PWD"

set -x

fromCase="$1"
toCase="$2"

[[ -z "$fromCase" || -z "$toCase" ]] && exit 1


if [[ ! -d "$toCase" ]]; then
    mkdir $toCase || exit 1
fi

cp -r "$fromCase/0" "$toCase"
cp -r "$fromCase/0.org" "$toCase"
cp -r "$fromCase/constant" "$toCase"
cp -r "$fromCase/system" "$toCase"
cp -r "$fromCase/tools" "$toCase"

cp "$fromCase/Allrun"* "$toCase"
cp "$fromCase/Allclean" "$toCase"

cd "$toCase"
exec "./Allclean"
cd -
