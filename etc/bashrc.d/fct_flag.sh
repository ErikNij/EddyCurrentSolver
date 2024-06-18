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
sourcePSD "$PSD/fct_bool.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

flag_nameOkOrExit ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return $?

    local nameStripped="$name"

    nameStripped=${name%%*'_'}
    nameStripped=${name%%'_'*}
    nameStripped=${name##*'_'}
    nameStripped=${name##'_'*}

    if [ "$name" = "$nameStripped" ]; then
        return 0
    else
        error_fatal "Flag name $name contains 'underscore(s)'." "$FUNCNAME" \
        || return $?
    fi

    return -1
}


flag_init ()
{
    local flags="$1"

    error_fatalIfEmptyVar 'flags' "$FUNCNAME" || return $?

    eval "flag_flagList='$flags'" || return $?

    local f
    for f in $flags; do
        flag_nameOkOrExit "$f" || return $?
        bool_false "$f" 'flag_isSet' || return $?
    done

    return 0
}

flag_initRead ()
{
    local flags="$1"
    shift

    error_fatalIfEmptyVar 'flags' "$FUNCNAME" || return $?

    flag_init "$flags"

    while [ $# -gt 0 ]; do
        flag_enable "$1" || return $?
        shift
    done

    return 0
}


flag_set ()
{
    local name="$1"
    local value="$2"

    error_fatalIfEmptyVar 'name value' "$FUNCNAME" || return $?

    local f
    for f in $flag_flagList; do

        if [ "$f" = "$name" ]; then
            bool_set "$f" "$value" 'flag_isSet' || return $?
            return 0
        fi
    done

    error_fatal "Flag '$name' is unknown. Defined flags: '$flag_flagList'." \
    "$FUNCNAME" || return $?

    return -1
}


flag_enable ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return $?

    flag_set "$name" 1 || return $?

    return 0
}


flag_disable ()
{
    local name="$1"

    error_fatalIfEmptyVar 'name' "$FUNCNAME" || return $?

    flag_set "$name" 0 || return $?

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
