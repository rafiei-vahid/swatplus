#!/usr/bin/env bash
# Build one engine variant from the wt_race worktree.
#   usage: build_variant.sh <name> "<extra ifx flags>"
# Produces build_<name>/swatplus. Never touches production lib/.
set -eo pipefail
NAME="$1"; EXTRA="${2:-}"
ROOT="/data/SWATGenXApp/codes/swatplus_perf/wt_race"
BUILD="$ROOT/build_$NAME"

export SETVARS_COMPLETED=1
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export PKG_CONFIG_PATH=/data/SWATGenXApp/codes/swatplus_perf/deps/netcdf-ifx/lib/pkgconfig:${PKG_CONFIG_PATH:-}

rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"
cmake -G "Unix Makefiles" \
  -DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/2026.0/bin/ifx \
  -DSWATPLUS_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE="-O3" \
  -DCMAKE_Fortran_FLAGS="$EXTRA" \
  "$ROOT" > cmake_cfg.log 2>&1
grep -E "Found netcdf|FFLAGS|RELEASE" cmake_cfg.log | head -5

# POSITIVE CONTROL ON THE BUILD ITSELF. EXTRA used to be injected via
# -DCMAKE_Fortran_FLAGS_RELEASE, and CMakeLists.txt line 219 does
#   set(CMAKE_Fortran_FLAGS_RELEASE "${frelease}")
# which OVERWRITES the -D cache value outright (line 218, commented out, once appended it). So
# every variant built with extra flags silently got the DEFAULT flags, and the build reported
# success. A variant experiment run that way is uncontrolled: it compares a binary against
# itself and calls the null result a finding.
# CMAKE_Fortran_FLAGS is appended to at line 216, so a -D value there survives -- and this
# check proves it landed rather than assuming it.
if [[ -n "$EXTRA" ]]; then
  probe=$(printf '%s' "$EXTRA" | awk '{print $1}')
  if ! grep -q -- "$probe" cmake_cfg.log; then
    echo "BUILD ABORTED: requested flag '$probe' did not reach the compile line." >&2
    echo "  configured FFLAGS: $(grep -m1 'FFLAGS' cmake_cfg.log)" >&2
    exit 1
  fi
  echo "FLAG CONFIRMED in compile line: $probe"
fi
# Bounded parallelism: this is the production web server, and -j$(nproc) is 32 cores of
# compile competing with live user model builds.
make -j"${JOBS:-6}" > build.log 2>&1
EXE=$(find . -maxdepth 1 -type f -executable -name 'swatplus*' | head -1)
echo "BUILT: $BUILD/$(basename "$EXE")"
# DO NOT invoke the binary for its version. The engine parses NO command-line arguments -- zero
# of 4,768 .f90 files call get_command_argument, and the binary has no argv symbols -- so
# `swatplus --version` does not print and exit. It runs a FULL SIMULATION in the current
# directory, truncating simulation.out and erosion.txt on the way in. This script did exactly
# that (CLAUDE.md hard rule #8). MEASURED, rather than assumed: build_ship48 and build_final
# contain NO simulation.out/erosion.txt/area_calc.out, so on those builds it evidently did not
# get that far -- most likely it died at the dynamic loader, since this script sets
# PKG_CONFIG_PATH but never LD_LIBRARY_PATH for the custom netcdf-ifx, and stderr was
# suppressed by 2>/dev/null. That is luck, not design: fix the loader path and the same line
# starts a simulation. Read the string out of the binary instead.
strings "$BUILD/$(basename "$EXE")" 2>/dev/null | grep -m1 "MODULAR Rev" || true
