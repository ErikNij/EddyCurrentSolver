#!/bin/bash

set -x

dataFile="$(basename ${0%.sh}).dat"

flagMake=0
flagUpdate=0
flagPlot=0

plotTime=20

while [[ $# -gt 0 ]]
do
  [[ $1 == 'make' ]] && flagMake=1

  if [[ $1 == 'update' ]]; then
    flagMake=1
    flagUpdate=1
  fi

  [[ $1 == 'plot' ]] && flagPlot=1

  shift
done

if [[ $flagMake -eq 1 ]]; then
  if [[ $flagUpdate -eq 1 && -e "$dataFile" ]]; then
    lastTime=$(tail -1 "$dataFile" | cut -d " " -f 1)
    getFluidDomainVolume "-time" "$lastTime:" "-file" "${dataFile}.update" || exit 1
    tail -n +2 "${dataFile}.update" > "${dataFile}.new"
    mv "$dataFile" "${dataFile}.old"
    cat "${dataFile}.old" "${dataFile}.new" > "$dataFile"
    rm "${dataFile}.update" "${dataFile}.old" "${dataFile}.new"
  else
    getFluidDomainVolume "-file" "$dataFile"
  fi
fi

if [[ $flagPlot -eq 1 ]]; then
  # gnuplot -e "plot '$dataFile' u 1:2, '$dataFile' u 1:3, '$dataFile' u 1:4; pause($plotTime)"
  gnuplot -e "set term wxt size 800,600; plot '$dataFile' u 1:( (\$2 + \$3 + \$4)/3 ) w lp; pause($plotTime)"
fi
