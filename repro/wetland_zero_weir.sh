#!/bin/bash
# Reproduction for the zero-height weir.res row in wet_initial.f90 (candidate 1 in CANDIDATES.md).
#
# Deck: upstream refdata/Osu_1hru with two edits, so its one wetland is a plain wetland, not a paddy:
#   wetland.wet   rel  weir      -> wetland     (pvol/evol release table already in res_rel.dtl)
#   hru-data.hru  lu   rice140_lum -> agrl_lum  (no weir_height operations)
# weir.res is left as shipped: one row, weir1, height 0.0.
#
# usage: wetland_zero_weir.sh <refdata dir> <exe> [<exe> ...]
# prints basin_wb_aa wetland columns and wetland_aa storage for each exe
set -e
ref=$1; shift
work=$(mktemp -d)
for exe in "$@"; do
  d=$work/$(basename $exe)
  cp -r $ref/Osu_1hru $d
  sed -i 's/paddy0001    initwet1       paddy        weir /paddy0001    initwet1       paddy     wetland /' $d/wetland.wet
  sed -i '3s/rice140_lum/   agrl_lum/' $d/hru-data.hru
  (cd $d && $exe > run.log 2>&1)
  echo "== $exe"
  python3 - $d <<'EOF'
import sys
d = sys.argv[1]
def last(f):
    L = open(f'{d}/{f}').read().splitlines()
    return dict(zip(L[1].split(), L[-1].split()))
b = last('basin_wb_aa.txt'); w = last('wetland_aa.txt')
print('basin_wb_aa (mm):', ' '.join(f'{k}={b[k]}' for k in ('wateryld', 'surq_gen', 'perc', 'et', 'wet_evap', 'wet_oflo', 'wet_stor')))
print('wetland_aa (m3): ', ' '.join(f'{k}={w[k]}' for k in ('flo_stor', 'evap', 'seep')))
EOF
done
rm -rf $work
