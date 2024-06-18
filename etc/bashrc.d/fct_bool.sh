#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #

sourcePSD "$PSD/fct_error.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

bool_unset ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return $?

    eval "unset -f ${prefix}_${name}" || return $?

    return 0
}


bool_set ()
{
    local name="$1"
    local value="$2"
    local prefix="${3:-bool}"

    error_fatalIfEmptyVar 'name value prefix' "$FUNCNAME" || return $?

    if [ $value -eq 0 ]; then
        eval "${prefix}_${name} () { return 1; }" || return $?

    elif [ $value -eq 1 ]; then
        eval "${prefix}_${name} () { return 0; }" || return $?
    else
        error_fatal \
            "Boolean value may only take '0' or '1' (given: '$value')." \
            "$FUNCNAME" || return $?
    fi

    return 0
}


bool_true ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return $?

    bool_set "$name" 1 "$prefix" || return $?

    return 0
}


bool_false ()
{
    local name="$1"
    local prefix="${2:-bool}"

    error_fatalIfEmptyVar 'name prefix' "$FUNCNAME" || return $?

    bool_set "$name" 0 "$prefix" || return $?

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
