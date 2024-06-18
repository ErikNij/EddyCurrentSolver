#!/bin/bash

width=110

cat "$1" | sed -e "s/ at /&\n    at /g" | sed "s/ at \$//" | sed -e "s/.\{$width\}/&\n    /g" > ${1}.aligned
