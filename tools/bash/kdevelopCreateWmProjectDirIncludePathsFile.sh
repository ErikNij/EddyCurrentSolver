#!/bin/bash

swakversion="0.4.0"

# Created/Updated file .kdev_include_paths for KDevelop

if [ -w "$WM_PROJECT_DIR" ]; then

  find "$WM_PROJECT_DIR/src" -name 'lnInclude' -print \
    > "$WM_PROJECT_DIR/.kdev_include_paths"

  find "$WM_PROJECT_DIR/applications" -name 'lnInclude' -print \
    >> "$WM_PROJECT_DIR/.kdev_include_paths"

  swakdir="$WM_PROJECT_DIR/ThirdParty/rpmBuild/BUILD/swak4Foam-$swakversion"
  echo $swakdir
  if [ -w "$swakdir" ]; then
    find "$swakdir" -name 'lnInclude' -print \
      >> "$WM_PROJECT_DIR/.kdev_include_paths"
  fi

  chmod go+r "$WM_PROJECT_DIR/.kdev_include_paths"
  echo "Created/Updated file .kdev_include_paths in directory $WM_PROJECT_DIR."
  exit 0
else
  echo "You do not have write permissions for \$WM_PROJECT_DIR ($WM_PROJECT_DIR)."
  echo "Try to rerun this script with 'sudo -E'."
  echo "WARNING: Pay attention to the -E (preserve-environment) option!"
  exit 1
fi
