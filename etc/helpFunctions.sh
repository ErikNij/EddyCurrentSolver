#!/bin/bash

#---------------------------------*- sh -*-------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     |
#   \\  /    A nd           | For copyright notice see file Copyright
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     RunFunctions
#
# Description
#
#------------------------------------------------------------------------------

setCdTrap ()
{
    pwd_safe="$1"

    setCdTrapFct ()
    {
        cd "$pwd_safe"
    }

    trap setCdTrapFct INT TERM EXIT
}


setErrorTrap ()
{
    setErrorTrapFct ()
    {
        exitCode="$?"
        if [[ "$exitCode" -ne 0 ]]; then
            echo 'ERROR: Something went wrong. Please check the log-files!'
        fi

        exit "$exitCode"
    }

    trap setErrorTrapFct INT TERM EXIT
}


#####


runApplication ()
{
    LOG_NAME=
    LOG_REGION=

    local OPTIND
    while getopts "r:l:" OPTFLAG; do
        case "$OPTFLAG" in
            r)
                LOG_REGION="$OPTARG"
                ;;
            l)
                LOG_NAME="$OPTARG"
                ;;
        esac
    done
    shift $((OPTIND-1))

    APP_RUN=$1; shift
    APP_NAME=$(basename $APP_RUN)

    if [ -z $LOG_NAME ] ; then
        if [ -z $LOG_REGION ] ; then
            LOG_NAME=log.$APP_NAME
        else
            LOG_NAME=log.$APP_NAME.$LOG_REGION
        fi
    fi

    if [ -f $LOG_NAME ] ; then
        echo "$APP_NAME already run: remove log file to run"
    else
        printf "Running $APP_NAME"
        if [ ! -z $LOG_REGION ]; then
            printf " for region $LOG_REGION\n"
        else
            printf "\n"
        fi
        $APP_RUN $* > $LOG_NAME 2>&1
        return $?
    fi
}


runParallel ()
{
    APP_PROCS=$1; shift

    LOG_NAME=
    LOG_REGION=

    local OPTIND
    while getopts "r:l:" OPTFLAG; do
        case "$OPTFLAG" in
            r)
                LOG_REGION="$OPTARG"
                ;;
            l)
                LOG_NAME="$OPTARG"
                ;;
        esac
    done
    shift $((OPTIND-1))

    APP_RUN=$1; shift
    APP_NAME=$(basename $APP_RUN)

    if [ -z $LOG_NAME ] ; then
        if [ -z $LOG_REGION ] ; then
            LOG_NAME=log.$APP_NAME
        else
            LOG_NAME=log.$APP_NAME.$LOG_REGION
        fi
    fi

    if [ -f $LOG_NAME ] ; then
        echo "$APP_NAME already run: remove log file to run"
    else
        printf "Running $APP_NAME in parallel using $APP_PROCS processes"
        if [ ! -z $LOG_REGION ]; then
            printf " for region $LOG_REGION\n"
        else
            printf "\n"
        fi
    if [ -z "$WM_MPIRUN_PROG" ]
    then
        mpirunProg=mpirun
    else
        # Allow exceentric systems to override the hardcoded mpirun
        mpirunProg=$WM_MPIRUN_PROG
    fi
        ( $mpirunProg -np $APP_PROCS $APP_RUN -parallel $* < /dev/null > $LOG_NAME 2>&1 )
        return $?
    fi
}


#####


getApplication ()
{
    getDictValueByKey 'application' 'system/controlDict'
}


getDictValueByKey ()
{
    grep "$1" "$2" | sed "s/$1\s\+\(\S*\);/\1/"
}


foamDictSubstituteByMatchAllKeys ()
{
    sed -i "s/\($1\s\+\)\($2\)\(\s*;\)/\1$3\3/g"  "$4"
}


foamDictSetAllKeys ()
{
    sed -i "s/\($1\s\+\)\(\S\+\)\(\s*;\)/\1$2\3/g"  "$3"
}


controlDictAscii ()
{
    foamDictSetAllKeys 'writeFormat' 'ascii' 'system/controlDict'
}


controlDictBinary ()
{
    foamDictSetAllKeys 'writeFormat' 'binary' 'system/controlDict'
}

polyMeshBoundaryFormatAscii ()
{
    baseDir='constant'
    [[ ! -z "$1" ]] && baseDir="${baseDir}/$1"

    foamDictSetAllKeys 'format' 'ascii' "$baseDir/polyMesh/boundary"
}


#####


paraviewTouchOpenFOAM ()
{
    if [[ -z "$1" ]]; then
        touch "$(basename $PWD).blockMesh"
        touch "$(basename $PWD).OpenFOAM"
    else
        touch "$(basename $PWD){$(basename $1)}.OpenFOAM"
    fi
}

paraviewTouchOpenFOAMParallel ()
{
    for p in $(ls -1d 'processor'*); do
        cd "$p" > /dev/null 2>&1
        paraviewTouchOpenFOAM "$1"
        cd - > /dev/null 2>&1
    done
}


paraviewTouchFoam ()
{
    if [[ -z "$1" ]]; then
        touch "$(basename $PWD).blockMesh"
        touch "$(basename $PWD).foam"
    else
        touch "$(basename $PWD){$(basename $1)}.foam"
    fi
}


paraviewTouchFoamParallel ()
{
    for p in $(ls -1d 'processor'*); do
        cd "$p" > /dev/null 2>&1
        paraviewTouchFoam "$1"
        cd - > /dev/null 2>&1
    done
}
