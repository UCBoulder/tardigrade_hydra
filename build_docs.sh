# Make bash script more like high-level languages.
set -Eeuxo pipefail

# Clean and build repo tests
rm -rf build/
mkdir build
cd build
cmake ..
cmake --build docs/sphinx --verbose
