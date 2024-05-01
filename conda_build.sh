#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status.

# Activate the build environment.
# Conda-build does this automatically, but if you need to manually adjust paths or
# set environment variables specific to the build, you can do so here.

# Ensure numpy is installed before attempting to run setup.py.
# This step is typically handled by the dependencies in meta.yaml, but is shown here for clarity.
echo "Ensuring numpy is installed..."
"${PYTHON}" -m pip install "numpy>=1.23.0" --no-cache-dir

# Run the package setup.
# This executes the setup.py script, which should handle any compilation steps
# for Fortran extensions or other package-specific build logic.
echo "Running package setup..."
"${PYTHON}" -m pip install . --no-deps -vv

# If there are additional steps required for your package build, they can be added here.
# For example, copying build artifacts to specific locations, additional logging, etc.

echo "Build script completed successfully."
