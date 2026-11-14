#!/bin/bash
# Builds the stress_resolved plugin (compute stress/resolved + kspace_style ewald2)
# into $BUILD_DIR/stress_resolved.so
#
# usage: bash build.sh /path/to/lammps/src
#    or: LAMMPS_SRC=/path/to/lammps/src bash build.sh
# The other variables below can be overridden the same way, e.g.
#   MPICXX=/opt/openmpi/bin/mpicxx bash build.sh /path/to/lammps/src

set -e

# LAMMPS source tree (the "src" folder). Must be the same version as the LAMMPS
# executable/library that will load the plugin, built with the KSPACE and PLUGIN packages
LAMMPS_SRC=${1:-$LAMMPS_SRC}

# MPI compiler wrappers, the same ones used to build LAMMPS
MPICC=${MPICC:-mpicc}
MPICXX=${MPICXX:-mpicxx}

# Release (-O3) or Debug
BUILD_TYPE=${BUILD_TYPE:-Release}

# integer sizes of the LAMMPS build (its -D LAMMPS_SIZES): smallbig, bigbig or smallsmall
LAMMPS_SIZES=${LAMMPS_SIZES:-smallbig}

# where the build happens and the .so is written
BUILD_DIR=${BUILD_DIR:-build}

HERE=$(cd "$(dirname "$0")" && pwd)

if [ -z "$LAMMPS_SRC" ]; then
  echo "error: give the LAMMPS src folder: bash build.sh /path/to/lammps/src" >&2
  exit 1
fi
if [ ! -f "$LAMMPS_SRC/lammpsplugin.h" ] || [ ! -f "$LAMMPS_SRC/KSPACE/ewald.h" ]; then
  echo "error: $LAMMPS_SRC is not a LAMMPS src folder (need lammpsplugin.h and KSPACE/ewald.h)" >&2
  exit 1
fi

cmake -S "$HERE" -B "$BUILD_DIR" \
      -D CMAKE_C_COMPILER="$MPICC" \
      -D CMAKE_CXX_COMPILER="$MPICXX" \
      -D CMAKE_BUILD_TYPE="$BUILD_TYPE" \
      -D LAMMPS_SIZES="$LAMMPS_SIZES" \
      -D LAMMPS_SOURCE_DIR="$(cd "$LAMMPS_SRC" && pwd)"

cmake --build "$BUILD_DIR"

echo "built $(cd "$BUILD_DIR" && pwd)/stress_resolved.so"
