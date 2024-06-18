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

env_set 'FOAM_USER_APP'       "$WM_PROJECT_USER_DIR/applications"
env_set 'FOAM_USER_SOLVERS'   "$WM_PROJECT_USER_DIR/applications/solvers"
env_set 'FOAM_USER_UTILITIES' "$WM_PROJECT_USER_DIR/applications/utilities"
env_set 'FOAM_USER_LIB'       "$WM_PROJECT_USER_DIR/lib"
env_set 'FOAM_USER_RUN'       "$WM_PROJECT_USER_DIR/run"
env_set 'FOAM_USER_SRC'       "$WM_PROJECT_USER_DIR/src"
env_set 'FOAM_USER_TOOLS'     "$WM_PROJECT_USER_DIR/tools"
env_set 'FOAM_USER_TUTORIALS' "$WM_PROJECT_USER_DIR/tutorials"

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
