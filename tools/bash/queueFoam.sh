#!/bin/bash
# OpenFOAM run PBS tool script
# September 2013
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Info ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# TODO [High]: Fix for foam-extend:
#              - "decomposePar -latestTime" does not exist

# Tool script name

VAR_PROG='qFoam'
VAR_BASH_DEP='awk basename bc cat grep ln mv rm sed sleep qsub'

# Function description

prog_function ()
{
cat << EOF

Function:
Prepare und submit job with PBS.
Needs a 'run.conf' file in a 'qRun' sub-directory
with the following variables set:
 - myPBS_host
 - myPBS_walltime
 - myPBS_nodes
 - myPBS_savetime
 - myPBS_queue
 - myPBS_jobname
 - myPBS_solver
 - myPBS_env

Requirements:
bash
  Apps: $VAR_BASH_DEP

Author:
Pascal Beckstein (p.beckstein@hzdr.de)

EOF
}

# Usage description

prog_usage ()
{
cat << EOF

Usage:
$VAR_PROG [ARGOPTS]

ARGOPTS:

  -h    Print this help

EOF
}

# --------------------------------------------------------------------------- #
# --- Command line options -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Read options

optshift=0
optstore="$@"
while getopts "h" optval "$@"; do
  case "$optval" in
    "h")
      prog_function
      prog_usage
      exit 0
    ;;
    "?")
      echo "ERROR: Unknown option! (HELP: -h)" && exit 1
    ;;
    ":")
      echo "ERROR: Missing argument for option! (HELP: -h)" && exit 1
    ;;
    *)
      echo "ERROR: Unknown error while processing options! (HELP: -h)" && exit 1
    ;;
  esac
done
cloptions=($@)
clcount=${#cloptions[@]}
trueshift=$(($optshift<$clcount?$optshift:$clcount))
shift $trueshift

# Check basic syntax

for i in $(eval echo "{1..$#}"); do
  case ${!i} in
    -*)
      echo "ERROR: Found options behind arguments! (HELP: -h)" && exit 1
    ;;
    *):;;
  esac
done

# --------------------------------------------------------------------------- #
# --- Base variables -------------------------------------------------------- #
# --------------------------------------------------------------------------- #

VAR_LOCALTOOLS_DIR='scripts'
VAR_SCRATCH_DIR="$VAR_LOCALTOOLS_DIR/$VAR_PROG"

VAR_PBS_DIR="$VAR_PROG"
VAR_PBS_QSUB_LOG="$VAR_PBS_DIR/log.qsub"
VAR_PBS_CONFIG="qFoam.conf"
VAR_PBS_PRE_SCRIPT="qFoam.pre"
VAR_PBS_PRE_LOG="$VAR_PBS_DIR/log.pre"
VAR_PBS_INTERMEDIATE_SCRIPT="qFoam.inter"
VAR_PBS_INTERMEDIATE_LOG="$VAR_PBS_DIR/log.inter"
VAR_PBS_POST_SCRIPT="qFoam.post"
VAR_PBS_POST_LOG="$VAR_PBS_DIR/log.post"
VAR_PBS_VARS="$VAR_PBS_DIR/run.vars_pbs"
VAR_PBS_RUN="$VAR_PBS_DIR/run"
VAR_PBS_JOB="$VAR_PBS_DIR/run.job"
VAR_PBS_JOB_LOG="$VAR_PBS_DIR/log.job"

VAR_PBS_STATE_FINISHED="$VAR_PBS_DIR/state.done"
VAR_PBS_STATE_FAILED="$VAR_PBS_DIR/state.error"
VAR_PBS_STATE_LOCKED="$VAR_PBS_DIR/state.locked"

VAR_PBS_CMD_STOP="$VAR_PBS_DIR/cmd.stop"
VAR_PBS_CMD_DISABLE="$VAR_PBS_DIR/cmd.disabled"

VAR_OPENFOAM_DIR="$PWD"
VAR_OPENFOAM_CASE="${PWD##*/}"
VAR_OPENFOAM_VARS="$VAR_PBS_DIR/run.vars_openfoam"
VAR_OPENFOAM_LOG="$VAR_PBS_DIR/log.solver"
VAR_OPENFOAM_DECOMPOSE_LOG="$VAR_PBS_DIR/log.decompose"

# --------------------------------------------------------------------------- #
# --- Print script head ----------------------------------------------------- #
# --------------------------------------------------------------------------- #

echo '/*---------------------------------------------------------------------------*\'
echo '| =========                 |                                                 |'
echo '| \\      /  F ield         | foam-extend: Open Source CFD                    |'
echo "|  \\\\    /   O peration     | Version:     $WM_PROJECT_VERSION                                |"
echo '|   \\  /    A nd           | Web:         http://www.foam-extend.org         |'
echo '|    \\/     M anipulation  | For copyright notice see file Copyright         |'
echo '\*---------------------------------------------------------------------------*/'
echo 'Build  : $WM_PROJECT_VERSION (tool script)'
echo "Exec   : $(basename $0) $optstore"
echo "Date   : $(date +'%b %d %Y')"
echo "Time   : $(date +%X)"
echo 'Host   : "'$HOSTNAME'"'
echo "PID    : $$"
echo "Case   : $VAR_OPENFOAM_DIR"
echo 'nProcs : 1'
echo '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
echo

# --------------------------------------------------------------------------- #
# --- Plausability ---------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Bash requirements
for dep in $VAR_BASH_DEP; do
  if [[ $(type -P $dep) == '' ]]; then
    echo "ERROR: Bash app $dep missing! (HELP: -h)" && exit 1
  fi
done

# --------------------------------------------------------------------------- #
# --- Preparation ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Check scratch
[[ ! -d "$VAR_LOCALTOOLS_DIR" ]] && mkdir "$VAR_LOCALTOOLS_DIR"
[[ ! -d "$VAR_SCRATCH_DIR" ]] && mkdir "$VAR_SCRATCH_DIR"

# Check PBS folder
[[ ! -d "$VAR_PBS_DIR" ]] && mkdir "$VAR_PBS_DIR"

# Read config file
echo "Reading config file $VAR_PBS_CONFIG"
if [[ -e "$VAR_PBS_CONFIG" ]]; then
  . "$VAR_PBS_CONFIG"
else
  echo "ERROR: Config file $VAR_PBS_CONFIG not found! (HELP: -h)"
  exit 1
fi
echo

# Abort if case is locked
if [[ -e "$VAR_PBS_STATE_LOCKED" ]]; then
  echo "ERROR: Case is locked since file $VAR_PBS_STATE_LOCKED was found! Aborting"
  exit 1
fi

# Abort if case failed
if [[ -e "$VAR_PBS_STATE_FAILED" ]]; then
  echo "ERROR: Case indicates failed state since file $VAR_PBS_STATE_FAILED was found! Aborting"
  exit 1
fi

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Openfoam definitions
echo 'Writing openfoam library file'
cat << EOF > "$VAR_OPENFOAM_VARS"
#!/bin/bash

##########################################################################################
#                                    PARAMETERS                                          #
##########################################################################################

##########################################################################################
#                                     FUNCTIONS                                          #
##########################################################################################

##########################################################################################
#                                     FUNCTIONS                                          #
##########################################################################################

# openFOAM Dictionary tools
dictLookup () { echo \$(cat \$1 | grep "^\s*\$2\s\{1,\}" | grep -v '^//' | sed 's|;||g' | awk '{ print \$2 }'); }
dictChange () { sed -i'.bak' "s|^\(\s*\$2\s\{1,\}\)\(\S\{1,\}\);|\1\$3;|g" \$1; }

# openFOAM log-File tools
logCurTime () { tail -500 '$VAR_OPENFOAM_LOG' | grep '^Time =' | awk '{ print \$3 }' | tail -1; }

EOF
echo

# PBS definitions
echo 'Writing library file'
cat << EOF > "$VAR_PBS_VARS"
#!/bin/bash

##########################################################################################
#                                    PARAMETERS                                          #
##########################################################################################

myPBS_savetime_s=\$(echo \$myPBS_savetime | awk -F: '{ print (\$1 * 3600) + (\$2 * 60) + \$3 }')
myPBS_np=\$(dictLookup 'system/decomposeParDict' 'numberOfSubdomains')
myPBS_np_include=\$(dictLookup 'system/decomposeParDict.include' 'numberOfSubdomains')
[[ -z \$myPBS_np ]] && myPBS_np=\$myPBS_np_include
myPBS_stopAt_orig=\$(dictLookup 'system/controlDict' 'stopAt')
myPBS_runTimeModifiable_orig=\$(dictLookup 'system/controlDict' 'runTimeModifiable')

##########################################################################################
#                                     VARIABLES                                          #
##########################################################################################

if [[ ! -z \$PBS_JOBID ]]; then

  myPBS_id="\$(echo \$PBS_JOBID | cut -f1 -d '.')"

fi

##########################################################################################
#                                     FUNCTIONS                                          #
##########################################################################################

# Date
myDate () { date +"[%Y.%m.%d-%H:%M:%S]"; }

# Logfile writing
writeLog () { echo -n "\$0 "; echo "\$(date +"[%Y.%m.%d-%H:%M:%S]") \$2" >> \$1; }
# writeLog () { printf '%10s ', \$0; echo "\$(date +"[%Y.%m.%d-%H:%M:%S]") \$2" >> \$1; }

# Cleanup trap
trapFct ()
{
  # Remove lock file and set failed status
  [[ -e '$VAR_PBS_STATE_LOCKED' ]] && mv '$VAR_PBS_STATE_LOCKED' '$VAR_PBS_STATE_FAILED'

  # Remove stop command file
  [[ -e '$VAR_PBS_CMD_STOP' ]] && rm '$VAR_PBS_CMD_STOP'

  # Make sure that original stopAt value gets set
  dictChange 'system/controlDict' 'stopAt' "\$myPBS_stopAt_orig"
}

# Processing function
jobProc ()
{
  case \$1 in
    'pre')
      script='$VAR_PBS_PRE_SCRIPT'
      log='$VAR_PBS_PRE_LOG'
    ;;
    'intermediate')
      script='$VAR_PBS_INTERMEDIATE_SCRIPT'
      log='$VAR_PBS_INTERMEDIATE_LOG'
    ;;
    'post')
      script='$VAR_PBS_POST_SCRIPT'
      log='$VAR_PBS_POST_LOG'
    ;;
    *)
      exitcode=1
      return
    ;;
  esac

  if [[ ! -e "\$script" ]]; then
    return
  else

    # Cahnge attribute and inform
    [[ ! -x "\$script" ]] && chmod +x "\$script"
    writeLog '$VAR_PBS_JOB_LOG' "Starting \$1 processing"
    writeLog "\$log" "Processing started"

    # Process script
    "./\$script" >> "\$log" &

    # Wait till done and save error
    jobproc_pid=\$!
    wait \$jobproc_pid
    exitcode=\$?
    writeLog "\$log" "Processing finished"
    writeLog '$VAR_PBS_JOB_LOG' "Job \$1 processing finished with exit code \$exitcode"
    exitcodesum=\$((\$exitcodesum+\$exitcode))
  fi
}

EOF
echo

# Include libraries now
. "$VAR_OPENFOAM_VARS"
. "$VAR_PBS_CONFIG"
. "$VAR_PBS_VARS"

# --------------------------------------------------------------------------- #
# --- Further checks and log file init -------------------------------------- #
# --------------------------------------------------------------------------- #

# Test if "runTimeModifiable" option is activated in "system/controlDict"
if [[ "$myPBS_runTimeModifiable_orig" != 'true' ]] && [[ "$myPBS_runTimeModifiable_orig" != 'on' ]] && [[ "$myPBS_runTimeModifiable_orig" != 'yes' ]]; then
  echo 'ERROR: The Option "runTimeModifiable" is not activated in "system/controlDict"! Aborting'
  exit 1
fi

# Test if "stopAt" option is set to "endTime" in "system/controlDict"
if [[ "$myPBS_runTimeModifiable_orig" != 'true' ]] && [[ "$myPBS_runTimeModifiable_orig" != 'on' ]] && [[ "$myPBS_runTimeModifiable_orig" != 'yes' ]]; then
  echo 'ERROR: The Option "stopAt" is not set to "endTime" in "system/controlDict"! Aborting'
  exit 1
fi

# # Test processor folder count
# if [[ $myPBS_np != $(ls -1d processor* 2> /dev/null | wc -l) ]]; then
#   if [[ $(ls -1d processor* 2> /dev/null | wc -l) > 0 ]]; then
#     echo 'ERROR: The processor folder count from decompositioning does not match "numberOfSubdomains" in "system/decomposeParDict"! Aborting'
#     exit 1
#   fi
# fi

# Touch log-files
touch $VAR_PBS_JOB_LOG
touch $VAR_PBS_QSUB_LOG
touch $VAR_OPENFOAM_LOG

# Link job log file
if [[ ! -h job ]]; then
  [[ -e job ]] && mv job job.bak
  ln -s "$VAR_PBS_JOB_LOG" job
fi
# Link openFOAM log file
if [[ ! -h log ]]; then
  [[ -e log ]] && mv log log.bak
  ln -s "$VAR_OPENFOAM_LOG" log
fi

# --------------------------------------------------------------------------- #
# --- Creae run file -------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Write run file
echo 'Writing run file'
cat << EOF > "$VAR_PBS_RUN"
#!/bin/bash

#PBS -d $PWD
#PBS -l walltime=$myPBS_walltime,$myPBS_nodes
#PBS -N ${myPBS_jobname:-"$(basename $PWD)"}
#PBS -q $myPBS_queue
#PBS -S /bin/bash

#PBS -j oe
#PBS -o '$PWD/$VAR_PBS_QSUB_LOG'

# Init error
exitcode=0
exitcodesum=0

# Load environment
. ~/.bashrc
eval "$myPBS_env"

# Change into working directory
cd $PWD

# Include library files
. '$VAR_OPENFOAM_VARS'
. '$VAR_PBS_CONFIG'
. '$VAR_PBS_VARS'

# Set exit trap
trap trapFct SIGHUP SIGINT SIGQUIT SIGTERM
# trap trapFct SIGTERM

# Write info to logfile
echo '##########################################################################################' >> '$VAR_PBS_JOB_LOG'
writeLog '$VAR_PBS_JOB_LOG' "Starting job \$myPBS_jobname aka \$PBS_JOBID"
writeLog '$VAR_PBS_JOB_LOG' "- Host: \$myPBS_host"
writeLog '$VAR_PBS_JOB_LOG' "- Queue: \$myPBS_queue"
writeLog '$VAR_PBS_JOB_LOG' "- Walltime: \$myPBS_walltime"
writeLog '$VAR_PBS_JOB_LOG' "- Nodes: \$myPBS_nodes"
writeLog '$VAR_PBS_JOB_LOG' "- Savetime: \$myPBS_savetime (\$myPBS_savetime_s s)"
writeLog '$VAR_PBS_JOB_LOG' "- Solver: \$myPBS_solver"
writeLog '$VAR_PBS_JOB_LOG' "- Environment: \$myPBS_env"
echo '------------------------------------------------------------------------------------------' >> '$VAR_PBS_JOB_LOG'


# Abort if case is disabled for PBS
if [[ -e '$VAR_PBS_CMD_DISABLE' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' 'Job submitting is disabled for this case since file $VAR_PBS_CMD_DISABLE was found. Aborting'
  exit 1
fi

# Abort if case is locked
if [[ -e '$VAR_PBS_STATE_LOCKED' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' 'Case is locked since file $VAR_PBS_STATE_LOCKED was found. Aborting'
  exit 1
fi

# Remove old state files
writeLog '$VAR_PBS_JOB_LOG' 'Removing old finished state file'
[[ -e '$VAR_PBS_STATE_FINISHED' ]] && rm '$VAR_PBS_STATE_FINISHED'

# Create lockfile
writeLog '$VAR_PBS_JOB_LOG' 'Locking case now'
echo \$PBS_JOBID > '$VAR_PBS_STATE_LOCKED'

# Start pre processing if script given
jobProc 'pre'
exitcode=\$?
exitcodesum=\$((\$exitcodesum+\$exitcode))

# Start batch job
writeLog '$VAR_PBS_JOB_LOG' 'Starting job'
'./$VAR_PBS_JOB'
exitcode=\$?
exitcodesum=\$((\$exitcodesum+\$exitcode))

# Analyse result and create corresponding state file
writeLog '$VAR_PBS_JOB_LOG' 'Checking simulation status'
writeLog '$VAR_PBS_JOB_LOG' "Job has reported exit code \$exitcode"
if [[ \$exitcode == 0 ]]; then
  writeLog '$VAR_PBS_JOB_LOG' "Job ok"
  mv '$VAR_PBS_STATE_LOCKED' '$VAR_PBS_STATE_FINISHED'
else
  writeLog '$VAR_PBS_JOB_LOG' "Job reported errors"
  mv '$VAR_PBS_STATE_LOCKED' '$VAR_PBS_STATE_FAILED'
fi

# Stop if job failed
if [[ -f '$VAR_PBS_STATE_FAILED' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' 'File $VAR_PBS_STATE_FAILED was found! Job chain will be aborted now'
  exit 1
fi

# Stop if stop file found
if [[ -e '$VAR_PBS_CMD_STOP' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' 'Stop file $VAR_PBS_CMD_STOP found! Job chain will be aborted now'
  rm '$VAR_PBS_CMD_STOP'
  exit 0
fi

# Stop if endTime reached
if [[ -e '$VAR_OPENFOAM_LOG' ]]; then
  simtime="\$(logCurTime)"
  writeLog '$VAR_PBS_JOB_LOG' "Current simulation time is: \$simtime"
  endtime="\$(dictLookup 'system/controlDict' 'endTime')"
  writeLog '$VAR_PBS_JOB_LOG' "Target simulation time is: \$endtime"
  if [[ \$(echo "scale=30; \$simtime >= \$endtime" | bc -l) == 1 ]]; then
    writeLog '$VAR_PBS_JOB_LOG' 'Solver endTime has been reached! Job chain will be stopped'

    # Start post processing if script given
    jobProc 'post'
    exitcode=\$?
    exitcodesum=\$((\$exitcodesum+\$exitcode))

    # Exit
    exit \$exitcodesum
  fi
fi

# Restart if not stopped up to now
if [[ -f '$VAR_PBS_STATE_FINISHED' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' 'Job chain will be continued now'
  ssh $myPBS_host "cd \$PWD; \$(which qsub) $VAR_PBS_RUN"
fi

# This should not happen
exit 1

EOF
chmod +x "$VAR_PBS_RUN"
echo

# --------------------------------------------------------------------------- #
# --- Create Run job file --------------------------------------------------- #
# --------------------------------------------------------------------------- #

echo 'Writing job file'
cat << EOF > "$VAR_PBS_JOB"
#!/bin/bash

# Init error
exitcode=0
exitcodesum=0

# Init runtime
runtime=0
runtime_start="\$(date +%s)"

# Include library file
. '$VAR_OPENFOAM_VARS'
. '$VAR_PBS_CONFIG'
. '$VAR_PBS_VARS'

# Create parallel folders if needed
if [[ ! -d processor0 ]]; then
  writeLog '$VAR_PBS_JOB_LOG' "Decomposing latestTime since no processor0 folder was found"
  decomposePar -latestTime >> $VAR_OPENFOAM_DECOMPOSE_LOG
fi

writeLog '$VAR_PBS_JOB_LOG' "Starting solver \$myPBS_solver"
if [[ \$myPBS_np == 1 ]]; then
  \$myPBS_solver >> '$VAR_OPENFOAM_LOG' 2>&1 &
else
  mpirun -np "\$myPBS_np" \$myPBS_solver -parallel >> '$VAR_OPENFOAM_LOG' 2>&1 &
fi

# Save last background pid
solver_pid=\$!
writeLog '$VAR_PBS_JOB_LOG' "Solver 'mpirun' process id is \$solver_pid on \$HOSTNAME"
solver_stopped='false'

# Wait for solver to finish
writeLog '$VAR_PBS_JOB_LOG' 'Waiting for solver to finish'
while [[ -d "/proc/\$solver_pid" ]]; do

  # Get current time
  runtime="\$(echo "\$(date +%s) - \$runtime_start" | bc -l)"
  simtime="\$(logCurTime)"

  # Print info
  writeLog '$VAR_PBS_JOB_LOG' "Current runtime is \$runtime s"
  if [[ ! -z "\$simtime" ]]; then
    writeLog '$VAR_PBS_JOB_LOG' "Current simtime is \$simtime"
  else
    writeLog '$VAR_PBS_JOB_LOG' "Simulation is beeing prepared"
  fi

  # Stop solver if savetime is reached
  if [[ \$(echo "scale=30; \$runtime < \$myPBS_savetime_s" | bc -l) == 0 ]]; then

    # Info
    writeLog '$VAR_PBS_JOB_LOG' "Savetime limit due"
    writeLog '$VAR_PBS_JOB_LOG' "Stopping solver via controlDict"
    solver_stopped='true'

  fi

  # Stop if stop file found
  if [[ -e '$VAR_PBS_CMD_STOP' ]]; then

    # Info
    writeLog '$VAR_PBS_JOB_LOG' 'Stop file $VAR_PBS_CMD_STOP found.'
    writeLog '$VAR_PBS_JOB_LOG' "Stopping solver via controlDict"
    solver_stopped='true'

  fi

  if [[ \$solver_stopped == true ]]; then

    # Change controlDict file to writeNow
    [[ ! -e system/controlDict.orig ]] && cp system/controlDict system/controlDict.orig
    dictChange 'system/controlDict' 'stopAt' 'nextWrite'

    # Break while loop
    break

  fi

  # Sleep
  writeLog '$VAR_PBS_JOB_LOG' "Sleeping for \$myPBS_sleeptime s"
  sleep \$myPBS_sleeptime

done

# Wait for solver process to finish and get exit status
writeLog '$VAR_PBS_JOB_LOG' "Waiting for solver to finish"
wait \$solver_pid
exitcode=\$?
writeLog '$VAR_PBS_JOB_LOG' "Solver finished with exit code \$exitcode"
exitcodesum=\$((\$exitcodesum+\$exitcode))

# Now restore original stopAt value
if [[ \$solver_stopped == 'true' ]]; then
  writeLog '$VAR_PBS_JOB_LOG' "Restoring original controlDict"
  dictChange 'system/controlDict' 'stopAt' "\$myPBS_stopAt_orig"
fi

# Start intermediate processing
jobProc 'intermediate'
exitcode=\$?
exitcodesum=\$((\$exitcodesum+\$exitcode))

# Exit
writeLog '$VAR_PBS_JOB_LOG' "Jobscript finished with exit code \$exitcodesum"
exit \$exitcodesum

EOF
chmod +x "$VAR_PBS_JOB"
echo

# --------------------------------------------------------------------------- #
# --- Start batch job ------------------------------------------------------- #
# --------------------------------------------------------------------------- #

echo -n 'Submitting job to PBS: '
qsub "$VAR_PBS_RUN"
echo

exit 0

# --------------------------------------------------------------------------- #
# --- Configuration file examples ------------------------------------------- #
# --------------------------------------------------------------------------- #

##########################################################################################
#                                 Example qFoam.conf                                     #
##########################################################################################

# myPBS_host='hydra'
# myPBS_walltime='00:15:00'
# myPBS_nodes='nodes=1:ppn=8'
# myPBS_savetime='00:10:00'
# myPBS_queue='short'
# myPBS_jobname='myJob'
# myPBS_solver='icoFoamBF'
# myPBS_env='fe40'
# myPBS_sleeptime='30'

##########################################################################################
#                                Example qFoam.inter                                     #
##########################################################################################

# # Init error
# exitcode=0
# exitcodesum=0
#
# # Include openfoam library
# . "qFoam/run.vars_openfoam"
# . "qFoam.conf"
# . "qFoam/run.vars_pbs"
#
# if [[ $myPBS_np == 1 ]]; then
#
#   # Link bodyforce
#   OF-linkBF.sh
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
#
# else
#
#   # Link bodyforce
#   OF-linkBF.sh -p
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
#
#   # Reconstruct latest time
#   reconstructPar -latestTime
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
#
# fi
#
# # Exit
# exit $exitcodesum

##########################################################################################
#                                Example qFoam.post                                      #
##########################################################################################

# # Init error
# exitcode=0
# exitcodesum=0
#
# # Include openfoam library
# . "qFoam/run.vars_openfoam"
# . "qFoam.conf"
# . "qFoam/run.vars_pbs"
#
# if [[ $myPBS_np == 1 ]]; then
#
#   continue
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
#
# else
#
#   # Reconstruct latest time
#   reconstructPar -latestTime
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
#
# fi
#
# # Calculate Courant number field
# Co -latestTime
# exitcode=$?
# exitcodesum=$(($exitcodesum+$exitcode))
#
# # Calculate yPlus for RAS calcs
# if [[ $(dictLookup 'constant/turbulenceProperties' 'simulationType') == 'RASModel' ]]; then
#   yPlusRAS -latestTime
#   exitcode=$?
#   exitcodesum=$(($exitcodesum+$exitcode))
# fi
#
# # Exit
# exit $exitcodesum
