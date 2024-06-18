#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

sourcePSD "$FOAM_USER_ETC/bashrc.d/fct_env.sh"
sourcePSD "$FOAM_USER_ETC/bashrc.d/fct_link.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

foamUserToolBinLinks ()
{
    local binDir="$WM_PROJECT_USER_DIR/bin"
    local toolsDir="$WM_PROJECT_USER_DIR/tools"

    local progtag='foamUserTool'

    link_rmFindDangling "$binDir" || return $?

    if dep_isHashedOrExit 'find'; then

        local toolList
            toolList="$(find "$toolsDir" -maxdepth 2 -type f)"

        local f
        for f in $toolList; do
#
            if [ -x $f ]; then

                local makeLink=1

                local dest="$f"
                local destBase="$(basename $f)"
                local destBaseNoExt="${destBase%.*}"

                local link="$binDir/${progtag}${destBaseNoExt^}"

                if [ -L "$link" ]; then

                    error_info "Removing old link '$link'." "$FUNCNAME" \
                    || return $?

                    rm "$link" || return $?

                elif [ -e "$link" ]; then

                    local msg="File '$link' already exists and is no link."
                        msg="$msg Skipping link for existing file '$dest'."

                    error_info "$msg" "$FUNCNAME" || return $?

                    makeLink=0

                fi

                if [ $makeLink -eq 1 ]; then

                    echo "Linking '$link' to '$dest'."
                    link_lnRel "$dest" "$link" || return $?

                fi

            else
                error_info "Skipping non-executable file '$dest'." \
                || return $?
            fi
        done
    else
        return $?
    fi

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
