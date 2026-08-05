#!/usr/bin/env bash
set -eo pipefail
ROOT=/data/SWATGenXApp/codes/swatplus_perf/wt_race
BUILD=$ROOT/build_tsan
export SETVARS_COMPLETED=1
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export PKG_CONFIG_PATH=/data/SWATGenXApp/codes/swatplus_perf/deps/netcdf-ifx/lib/pkgconfig:${PKG_CONFIG_PATH:-}
rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"
# -g BELONGS IN CMAKE_Fortran_FLAGS, NOT _RELEASE: CMakeLists.txt:219 overwrites _RELEASE, so a
# -g placed there vanishes and every frame in a TSan report becomes "??:?" -- the tool then names
# the ROUTINE but not the LINE or the racing variable, which is most of its value. The first TSan
# run produced 248 real reports that could not be localised for exactly that reason.
# NOTE: comments cannot sit INSIDE the backslash-continued cmake invocation below -- doing that
# terminated the command and CMake then treated the build dir as the source dir.
cmake -G "Unix Makefiles" \
  -DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/2026.0/bin/ifx \
  -DSWATPLUS_OPENMP=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE="-O1 -g" \
  -DCMAKE_Fortran_FLAGS="-fsanitize=thread -fpe3 -g -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" \
  "$ROOT" > cmake_cfg.log 2>&1

# THE SANITISER MUST REACH THE COMPILE LINE, NOT JUST THE LINKER.
# This script passed -fsanitize=thread via CMAKE_Fortran_FLAGS_RELEASE, which CMakeLists.txt:219
# OVERWRITES with ${frelease}. The LINKER flag applied, so the produced binary carried 137 __tsan
# symbols and looked instrumented -- while every compile unit was built with
# "-fpp -free -fpe0 -traceback" and no instrumentation at all. A TSan binary in that state reports
# NO RACES WHETHER OR NOT RACES EXIST: pass-signal == failure-signal, in the one tool whose entire
# job is finding races. Same root cause as build_variant.sh (fixed 2026-08-03).
if ! grep -q -- "-fsanitize=thread" cmake_cfg.log; then
  echo "BUILD ABORTED: -fsanitize=thread did not reach the compile line." >&2
  grep -E "FFLAGS|RELEASE" cmake_cfg.log >&2
  exit 1
fi
echo "FLAG CONFIRMED in compile line: -fsanitize=thread"
make -j"${JOBS:-6}" > build.log 2>&1 || { tail -20 build.log; exit 1; }
find . -maxdepth 1 -type f -executable -name 'swatplus*' | head -1
