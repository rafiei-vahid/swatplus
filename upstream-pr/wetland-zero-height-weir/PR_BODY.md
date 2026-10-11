wet_initial: ignore a zero-height weir.res row instead of setting every wetland's principal volume to zero

**Problem**

`wet_initial.f90:138` takes the weir branch whenever weir.res has a row:

```fortran
if (db_mx%res_weir > 0.and.wet_ob(iihru)%iweir>0) then
  ...
  wet_ob(iihru)%pvol = hru(iihru)%area_ha * res_weir(iweir)%h * 10.**4  ! m3
```

`wet_ob%iweir` defaults to 1 (`reservoir_module.f90:54`). Only a `weir` management operation changes it (`mgt_sched.f90:519`). So every wetland takes row 1 of weir.res.

The SWAT+ Editor writes weir.res with one row of height 0. An example is the Editor v2.3.3 deck on the netCDF-Support branch (`data/netCDF-climate`): `shape1 1.84 2.60 2.50 0.00000`, with 142 wetland HRUs. With that row, every wetland gets pvol = 0 and pdep from hydrology.wet is never used. The `wetland` release table (vol > pvol: release over 10 days) then drains each wetland every day.

**Change**

The weir branch is taken only when the referenced weir.res row has a height > 0. A zero-height row is ignored, and the wetland keeps pvol and the weir height from pdep (the existing else branch).

The nested if avoids indexing `res_weir(0)` and `res_weir(iweir)` with iweir > db_mx%res_weir.

Rows with a real height, paddies included, take the same path as before.

One file: +8 −1.

**Reproduction**

This uses refdata/Osu_1hru with two edits, so the one wetland is a plain wetland:
- wetland.wet: rel `weir` → `wetland`
- hru-data.hru: lu_mgt `rice140_lum` → `agrl_lum` (no weir_height operations)

weir.res is left as shipped (weir1, height 0.0). Run 2010-2023 with print.prt as shipped (nyskip 3), gfortran 13.3 Release.

basin_wb_aa:

| | main 80049bc | this PR |
|---|---|---|
| wateryld (mm) | 576.674 | 387.935 |
| surq_gen (mm) | 576.076 | 386.369 |
| perc (mm) | 224.722 | 474.658 |
| et (mm) | 373.407 | 312.210 |
| wet_evap (mm) | 21.976 | 232.557 |
| wet_stor (mm) | 0.099 | 80.107 |

wetland_aa:

| | main 80049bc | this PR |
|---|---|---|
| evap (m3) | 2.207e3 | 2.338e4 |
| seep (m3) | 1.704e5 | 7.299e4 |

**Other checks**

- **Ames_sub1 (no wetlands):** all 67 output files identical to main.
- **Osu_1hru as shipped:** the paddy's weir height is set by weir_height operations, but weir1 has height 0, so this deck is affected at initialisation only. pvol at initialisation is now pdep instead of 0, and the initial nutrient mass is concentration × pvol.
  - basin_wb_aa is identical to main.
  - Average-annual aquifer flow goes 365.531 → 365.517 mm.
  - The largest change in any `*_aa.txt` value is channel sedp_out, 4.840 → 4.801 (0.8 %).
- **Release build with `-fcheck=bounds,pointer`:** the reproduction deck gives the same numbers, and shipped Osu_1hru runs to completion.
- The Debug build (`-ffpe-trap=...,underflow`) stops on main itself at `aqu_1d_control.f90:160` on shipped Osu_1hru (#288), with or without this change.

**Not changed here**

The WEIR column of wetland.wet is in the file but is not read, so a wetland cannot pick its own weir.res row. That may be the better long-term fix, but it changes the input contract, so I left it out.

dev has the same wet_initial.f90, and the patch applies there unchanged.
