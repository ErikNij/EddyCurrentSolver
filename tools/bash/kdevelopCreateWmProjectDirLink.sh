#!/bin/bash

# Link $WM_PROJECT_DIR in $WM_PROJECT_USER_DIR for use in local .kdev_include_paths files (KDevelop)

linkName="$WM_PROJECT_USER_DIR/.wm_project_dir"
linkTarget="$WM_PROJECT_DIR"
echo "Linking $linkName -> $linkTarget"
[[ -L "$linkName" ]] && rm "$linkName"
ln -s "$linkTarget" "$linkName"

exit 0