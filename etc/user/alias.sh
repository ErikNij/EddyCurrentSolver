#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #
# --- Script ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

alias uapp='cd $FOAM_USER_APP'
alias usol='cd $FOAM_USER_SOLVERS'
alias uutil='cd $FOAM_USER_UTILITIES'
alias ulib='cd $FOAM_USER_LIB'
alias urun='cd $FOAM_USER_RUN'
alias usrc='cd $FOAM_USER_SRC'
alias utools='cd $FOAM_USER_TOOLS'
alias utut='cd $FOAM_USER_TUTORIALS'

alias uruntut='cd ${FOAM_USER_TUTORIALS}${PWD#$FOAM_USER_RUN/tutorials}'
alias ututrun='cd ${FOAM_USER_RUN}/tutorials${PWD#$FOAM_USER_TUTORIALS}'

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
