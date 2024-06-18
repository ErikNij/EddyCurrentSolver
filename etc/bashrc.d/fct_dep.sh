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

dep_isHashed ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return $?

    hash $bin 2>/dev/null || return $?

    return 0
}


dep_isHashedVerb ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return $?

    if dep_isHashed $bin; then
        return 0
    else
        error_warning "The binary/application '$bin' is not hashed." \
        "$FUNCNAME" || return $?
    fi

    return -1
}


dep_isHashedOrExit ()
{
    local bin="$1"

    error_fatalIfEmptyVar 'bin' "$FUNCNAME" || return $?

    local msg
        msg="This is considered fatal."
        msg="$msg Please make sure to install '$bin'."

    if dep_isHashedVerb $bin; then
        return 0
    else
        error_fatal "$msg" "$FUNCNAME" || return $?
    fi

    return -1
}


dep_listIsHashed ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return $?

    local err=0

    local b
    for b in $binList; do
        dep_isHashed $b
        [ $? -eq 0 ] || err=$?
    done

    return $err
}


dep_listIsHashedVerb ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return $?

    local err=0

    local b
    for b in $binList; do
        dep_isHashedVerb $b
        [ $? -eq 0 ] || err=$?
    done

    return $err
}


dep_listIsHashedOrExit ()
{
    local binList="$1"

    error_fatalIfEmptyVar 'binList' "$FUNCNAME" || return $?

    local b
    for b in $binList; do
        dep_isHashedOrExit $b || return $?
    done

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
