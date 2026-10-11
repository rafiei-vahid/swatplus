#!/bin/bash
# Reproduction for the res_control.f90 assignments dropped in 109179b (candidate 2 in CANDIDATES.md).
#
# Deck: upstream refdata/Osu_1hru plus one reservoir below cha01 whose release table is a single
# nonirr-h06 action (alpha 0.85). Principal volume 1.5 ha-m, emergency 2.5 ha-m.
#
# usage: reservoir_h06.sh <refdata dir> <exe> [<exe> ...]
# prints reservoir_aa inflow, outflow, storage and area for each exe
set -e
ref=$1; shift
work=$(mktemp -d)
for exe in "$@"; do
  d=$work/$(basename $exe)
  cp -r $ref/Osu_1hru $d
  sed -i -E '/^connect/{s/(aquifer\.con +null +null +)null/\1reservoir.con/}' $d/file.cio
  sed -i -E '/^reservoir/{s/^reservoir( +)initial\.res( +)null( +)null/reservoir\1initial.res\2reservoir.res\3hydrology.res/}' $d/file.cio
  python3 - $d <<'EOF'
import re, sys
d = sys.argv[1]
L = open(f'{d}/object.cnt').read().splitlines()
h, v = L[1].split(), L[2].split()
v[h.index('obj')] = '5'; v[h.index('res')] = '1'
L[2] = '  '.join(v); open(f'{d}/object.cnt', 'w').write('\n'.join(L) + '\n')
L = open(f'{d}/chandeg.con').read().splitlines()
L[2] = re.sub(r'0\s*$', '1           res         1           tot       1.00000', L[2].rstrip())
open(f'{d}/chandeg.con', 'w').write('\n'.join(L) + '\n')
t = open(f'{d}/res_rel.dtl').read().splitlines()
t[1] = str(int(t[1].split()[0]) + 1)
t += ['', 'name                     conds      alts      acts  ',
      'h06_nonirr                   1         1         1  ',
      'var                        obj   obj_num           lim_var            lim_op     lim_const      alt1  ',
      'month                     null         0              null                 -       0.00000         >  ',
      'act_typ                    obj   obj_num              name            option         const        const2                fp  outcome           ',
      'release                    res         0               h06        nonirr-h06       0.85000       0.00000              null  y   ']
open(f'{d}/res_rel.dtl', 'w').write('\n'.join(t) + '\n')
EOF
  cat > $d/reservoir.con <<'EOF'
reservoir.con: one reservoir below cha01
      id  name                gis_id          area           lat           lon          elev       res               wst       cst      ovfl      rule   out_tot
       1  res001                   1       1.00000      35.52632     127.32895             0         1    s35610n127290e         0         0         0         0
EOF
  cat > $d/reservoir.res <<'EOF'
reservoir.res: one reservoir below cha01
      id  name                  init           hyd           rel           sed           nut
       1  res001              initwet1        hydres1      h06_nonirr      sedwet1       nutwet1
EOF
  cat > $d/hydrology.res <<'EOF'
hydrology.res: one reservoir below cha01
name            yr_op   mon_op   area_ps   vol_ps   area_es   vol_es      k     evap_co   shp_co1   shp_co2
hydres1             0        0      1.00     1.50      1.50     2.50    0.0      0.60      0.00      0.00
EOF
  (cd $d && $exe > run.log 2>&1)
  echo "== $exe"
  python3 - $d <<'EOF'
import sys
d = sys.argv[1]
L = open(f'{d}/reservoir_aa.txt').read().splitlines()
r = dict(zip(L[1].split(), L[-1].split()))
print('reservoir_aa:', ' '.join(f'{k}={r[k]}' for k in ('flo_in', 'flo_out', 'flo_stor', 'area', 'evap')))
EOF
done
rm -rf $work
