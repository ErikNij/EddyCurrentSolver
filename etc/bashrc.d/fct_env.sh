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
sourcePSD "$PSD/fct_dep.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

env_set ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return $?

    eval "export $var='$key'" || return $?

    return 0
}


env_unset ()
{
    local var="$1"

    error_fatalIfEmptyVar 'var' "$FUNCNAME" || return $?

    unset "$var" || return $?

    return 0
}


env_unsetRegex ()
{
    local keyRegex="$1"

    error_fatalIfEmptyVar 'keyRegex' "$FUNCNAME" || return $?

    if dep_isHashedOrExit 'sed'; then

        local keyExpand
            keyExpand="$( env | \
                sed -n "/^$keyRegex=.*/p" | \
                sed 's/^\([^=]*\)=\(.*\)/\1/g' )" || return $?

        local v
        for v in $keyExpand; do
            unset "$v" || return $?
        done

    else
        return $?
    fi

    return 0
}


env_remove ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return $?

    if dep_listIsHashedOrExit 'sed awk'; then

        local keyRemove
            keyRemove="$var=\$(echo -n \$$var | \
                awk -v RS=: -v ORS=: '\$0 != \"$key\"' | \
                sed 's/:$//')"

        eval "$keyRemove" || return $?

        local value
            value="$(eval "echo \$$var")" || return $?

        if [ -z "$value" ]; then
            unset "$var" || return $?
        fi

    else
        return $?
    fi

    return 0
}


env_removeAny ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return $?

    if dep_isHashedOrExit 'sed'; then

        local keyRemove
            keyRemove="$var=\$(echo -n \$$var | \
                sed 's/^[^:]*$key[^:]*://g' | \
                sed 's/[^:]*$key[^:]*://g' | \
                sed 's/[^:]*$key[^:]*$//g')"

        eval "$keyRemove" || return $?

        local value
            value="$(eval "echo \$$var")" || return $?

        if [ -z "$value" ]; then
            unset "$var" || return $?
        fi

    else
        return $?
    fi

    return 0
}


env_removeAll ()
{
    local key="$1"

    error_fatalIfEmptyVar 'key' "$FUNCNAME" || return $?

    if dep_listIsHashedOrExit 'grep sed'; then

        local varList
            varList="$(env | grep "$key" | sed 's/^\([^=]*\)=\(.*\)/\1/g')" \
            || return $?

        local v
        for v in $varList; do
            env_removeAny "$v" "$key" || return $?
        done

        local varList
            varList="$(env | grep '=$' | sed 's/=//g')" \
            || return $?

        local v
        for v in $varList; do
            env_unset "$v" || return $?
        done

    else
        return $?
    fi

    return 0
}


env_append ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return $?

    env_remove "$var" "$key" || return $?

    local value
        value="$(eval "echo \$$var")" || return $?

    if [ -z "$value" ]; then
        eval "export $var=\"$key\"" || return $?
    else
        eval "export $var=\"\$$var:$key\"" || return $?
    fi

    return 0
}


env_prepend ()
{
    local var="$1"
    local key="$2"

    error_fatalIfEmptyVar 'var key' "$FUNCNAME" || return $?

    env_remove "$var" "$key" || return $?

    local value
        value="$(eval "echo \$$var")" || return $?

    if [ -z "$value" ]; then
        eval "export $var=\"$key\"" || return $?
    else
        eval "export $var=\"$key:\$$var\"" || return $?
    fi

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
