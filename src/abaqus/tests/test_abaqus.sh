# Make bash script more like high-level languages.
set -Eeuo pipefail

# Get this script's file name
script=`basename "$0"`

# Parse arguments
if [ "$#" -ne 3 ]; then
    echo "${script} USAGE:"
    echo "./${script} ABAQUS_PROGRAM USER_SUBROUTINE INPUT_FILE"
    echo "  REQUIRED POSITIONAL ARGUMENTS:"
    echo "    ABAQUS_PROGRAM: abaqus executable on PATH or absolute path"
    echo "    USER_SUBROUTINE: user subroutine file path, absolute or relative to script location"
    echo "    INPUT_FILE: input file path, absolute or relative to script location"
    exit 1
fi
abaqus_program=$1
input=$2
user=$3

# Verify abaqus path 
if [ ! "$(command -v ${abaqus_program})" ]; then
    echo "Command '${abaqus_program}' not found."
    exit 2
fi

# Verify expected input file path
if [ ! -f "${input}" ]; then
    echo "Input file '${input}' not found."
    exit 3
else
    # Strip paths and extension for abaqus
    basename="${input##*/}"
    job="${basename%.*}"
    # Copy input file to local directory
    cp ${input} .
fi

# Verify user subroutine file
if [ ! -f "${user}" ]; then
    echo "user file ${user} not found"
    exit 4
fi

# Run abaqus against the input file and redirect interactive output to log file
# Use interactive to avoid sleep and wait statements
${abaqus_program} -job ${job} -user ${user} -interactive -ask_delete OFF >> ${job}.log 2>&1

# Test the output
grep "Hello Abaqus" ${job}.log 
grep "COMPLETED SUCCESSFULLY" ${job}.sta

# Check commands with expected error codes
set +e
error_pattern='ERROR'
error_matches=$(grep ${error_pattern} ${job}.log)
if [ ! -z "${error_matches}" ]; then
    echo "Found errors in ${job}.log"
    echo "${error_matches}"
    exit 5
else
    echo "No '${error_pattern}' error pattern found in ${job}.log"
fi
