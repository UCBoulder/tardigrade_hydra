# USAGE:
#
# ./new_build.sh cmake_build_type

# Make bash script more like high-level languages.
set -Eeuxo pipefail

# Get this scripts file name
script=`basename "$0"`

# Parse arguments
if [ "$#" -ne 1 ]; then
    echo "${script} USAGE:"
    echo "./${script} cmake_build_type"
    echo "    cmake_build_type: string for the CMake config -DCMAKE_BUILD_TYPE=<string> option"
    exit 1
fi
cmake_build_type=$1

# Debugging
whoami
ls -l $HOME/include || true
ls -l $HOME/.local/include || true

# Clean and build repo tests
rm -rf build/
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=${cmake_build_type}
cmake --build . --verbose
