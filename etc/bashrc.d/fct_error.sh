#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #

sourcePSD "$PSD/colors.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

error_fatal ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo -e >&2 "ERROR: $functionStr$message"

    case $- in
        *i*)
            return 1
        ;;
        *)
            exit 1
        ;;
    esac
}


error_fatalIfEmptyVar ()
{
    local variable="$1"
    local function="$2"

    local v
    for v in $variable; do

        local value="$(eval "echo \$$v")"

        [ ! -z "$value" ] || \
            error_fatal "Variable '$v' is empty." "$function" || return $?

    done

    return 0
}


error_warning ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo -e >&2 "WARNING: $functionStr$message"

    return 1
}


error_info ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [ ! -z "$function" ] && functionStr="[$function] "

    echo -e >&2 "INFO: $functionStr$message"

    return 0
}


error_echo ()
{
    local message="$1"
    local function="$2"

    local functionStr=''
    [[ ! -z "$function" ]] && functionStr="[$function] "

    echo -e >&2 "$functionStr$message"

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
