# Q1281: uninitialised locals in the SWAT+ fork

Branch `fix/q1281-uninit-locals` of rafiei-vahid/swatplus, based on `7eca087` (the deployed source). Compile-only audit: no model or engine was run.

## Summary

- gfortran 13.3.0 with `-O2 -Wmaybe-uninitialized -Wuninitialized` on `7eca087` + `c80ab8ce` reports **60 uninitialised-variable warnings in 43 files** (0 errors).
- **42 are REAL** (set on a branch and read outside it) and **18 are false positives**.
- 41 REAL warnings are fixed in **30 commits, one per routine**: an executable assignment at the top of the routine, never a declaration initializer (which would imply SAVE and break OpenMP thread privacy). Two variables gfortran did not flag but that are read uninitialised in the same statements are fixed in the same commits: `istr1` (actions) and `ob2` (hyd_read_connect).
- 1 REAL site is reported and not fixed: `res_control` `irrig_track_b` is cross-call state, not a scratch local (see Follow-ups). That was Vahid's ruling. Per Lane O, it is not static in the shipped binary either, so there it is undefined on every call (see Follow-ups).
- After the fixes the same build reports **19 warnings**: 19 expected (the 18 false positives + res_control), 0 new, and 0 on a fixed variable.
- 4 REAL sites are in routines that nothing calls at `7eca087` (header_cs, recall_read_cs, recall_read_salt, reg_read_elements). The linker discards them under `--gc-sections`. They are fixed anyway so they stay safe if re-enabled.
- Per Vahid's ruling, every fix restores upstream's own initial value. On the uninitialised path **18 fixed sites are still wrong**: 15 index 0 out of bounds and 3 trap under `-fpe0`. The fix makes them deterministic and loud, not correct. Each is listed under Follow-ups.

## Merged branch: `fix/q1281-q1279-merged` (Lane P, step 1)

- Merge commit `28d25c9` joins `fix/q1281-uninit-locals` (`2de4bec`) and `fix/q1279-nut-psed-fractions` (`c80ab8c`) with `--no-ff`. **The merge was clean**: no conflicts, and no file is changed by both branches. Neither source branch was rewritten or force-pushed.
- The q1279 branch carries 5 Q1072 `plant_parm_read` commits beneath `c80ab8ce`, so they come along with the merge. They are not in `7eca087` (they are in `ship/q1130-q1139-combined`): `173befc`; `3115eb3`; `2e0b97e`; `9c4710e`; `c56f19e`.
- Clean builds (fresh build directory, serial `make`) with gfortran 13.3.0, the repo's CMake Release flags minus `-finit-local-zero`, and `-O2 -Wall -Wuninitialized`. Counts are diagnostic lines that start with `Warning:`:

| Warning class | `7eca087` (deployed) | `fix/q1281-uninit-locals` | `fix/q1279-nut-psed-fractions` | `fix/q1281-q1279-merged` |
|---|---|---|---|---|
| `-Wmaybe-uninitialized` | 62 | 22 | 59 | 19 |
| `-Wuninitialized` | 1 | 0 | 1 | 0 |
| `-Wtabs` | 531 | 531 | 531 | 531 |
| `-Wconversion` | 211 | 211 | 211 | 211 |
| `-Wunused-variable` | 207 | 207 | 207 | 207 |
| `-Wunused-label` | 34 | 34 | 34 | 34 |
| `-Wcharacter-truncation` | 12 | 12 | 12 | 12 |
| `-Wunused-dummy-argument` | 7 | 7 | 7 | 7 |
| `-Wunused-function` | 1 | 1 | 1 | 1 |
| legacy/extension/deleted-feature notices | 53 | 53 | 53 | 53 |
| **total** | **1119** | **1078** | **1116** | **1075** |
| errors | 0 | 0 | 0 | 0 |

The merge removes exactly the uninitialised-variable warnings of both branches: 63 at `7eca087` (the 60 audited + nut_psed's 3) become 19 on the merged branch, which are the 18 false positives and `res_control`. Every other class is unchanged.

## Method

1. Source: `7eca087` with `c80ab8ce` (nut_psed, Q1279) cherry-picked on top. The fix/q1279 branch tip also carries five Q1072 commits (`c56f19e`..`173befc`); they are not included, so this is exactly "7eca087 plus c80ab8ce".
2. Build: the repo's own CMake, Release, `SWATPLUS_NETCDF=ON` (netcdf-fortran 4.5.4) and `SWATPLUS_OPENMP=OFF` (the OpenMP directives are ignored; the locals are the same). In a scratch copy of `CMakeLists.txt` only, the Release flags were set to `-O2 -Wmaybe-uninitialized -Wuninitialized`, and `-finit-local-zero -fno-unsafe-math-optimizations -frounding-math` was removed from the gfortran dialect line. `-finit-local-zero` zero-fills locals, which hides exactly this class, and the production ifx build has no equivalent. The serial build log was parsed for `-Wmaybe-uninitialized` / `-Wuninitialized` diagnostics. Control for the two numeric flags: a build that removes only `-finit-local-zero` reports the identical list (60 warnings, same file, line and variable).
3. **Positive control**: `nut_psed.f90` at bare `7eca087` gives 3 warnings (`fr_orgp` 77, `fr_actmin` 78, `fr_stamin` 81), the Q1279 bug. With `c80ab8ce` it gives 0. So the instrument sees a known REAL site and stops seeing it once fixed.
4. Each warning was classified by reading the routine. Reachability was checked with the linker (`-Wl,--gc-sections -Wl,--print-gc-sections` on the same objects; control: `actions` kept, `reg_read_elements`, whose only call is commented out, removed). Upstream's declaration was read at `cb442f7`, the upstream parent of merge `bcb882c`. The commit that removed each initializer was found by walking `7eca087`'s full first-parent history. The session's clone was shallow at `8fa9083`, which made it look like a root commit; the history was deepened before attribution, and `8fa9083`'s parent `051db48` still has the initializers it removed.

## Verification

- **Compile**: all 30 changed files compile. The full serial rebuild of `7eca087` + fixes + `c80ab8ce` with the same flags reports 19 warnings: 0 new, and 0 on a variable this branch assigns.
- **Neutral where the branch was taken**: for every fixed variable, a script checks that the inserted assignment is the variable's first executable reference in its routine. Every pre-existing assignment therefore still runs after it and before the read it fed, and only the path that reached the read with nothing assigned changes. Result: 0 violations. Control: the same script reports `BAD` when the block is moved below an existing assignment. Two fixes assign something other than upstream's 0: `chg_par = val_cur` (all four arms overwrite it) and `koc = cb_wtr_coef%prmt_21` (the only value koc is ever assigned, so it is the same value on the assigned path).
- **The `-fpe0` claims**: in a 10-line standalone program built with `gfortran -ffpe-trap=invalid,zero,overflow` (the gfortran equivalent of ifx `-fpe0`), `log(0.)`, `0.**(-.385)` and `0./0.` each raise SIGFPE, and the controls with 1.0 run. ifx is not installed in this environment.

## All warnings

Line numbers are gfortran's, in `7eca087` numbering (c80ab8ce touches only nut_psed.f90). Tags: `OOB-0` means the restored value 0 is still an out-of-bounds index on the uninitialised path. `FPE` means it traps under `-fpe0`. `dead@7eca087` means the routine is not called.

| # | Warning | Variable | Class | Fix commit | Input condition that reaches the uninitialised read (REAL) / why it is always set (FALSE POSITIVE) | Initializer removed in | Tags |
|---|---|---|---|---|---|---|---|
| 1 | `actions.f90:184` | `irrop` | REAL | `418f098` | An `irr_demand` action (not the `ponding` variant) whose applications already reached its limit (`num_actions(iac) > int(const2)`) with `pco%mgtout = y`: the mgt.out write reads `irrop_db(irrop)`, assigned only in the other branch. | `8fa9083` |  |
| 2 | `actions.f90:953` | `ipud` | REAL | `418f098` | A `puddle` action whose `option` matches no `puddle.ops` name (also read at 955). | `8fa9083` |  |
| 3 | `actions.f90:1200` | `istr` | REAL | `418f098` | Any `grassww_uninstall` action: its unconditional land-use-change write reads `sdr(istr)`, never assigned in that case (unless an earlier structure action in the same call set it). | `8fa9083` |  |
| 4 | `aqu_read_elements.f90:117` | `ihru` | REAL | `56fc5a7` | The aquifer-region block reads a region row with `nspu = 0`, and the `acu_out` block above did not take its own `nspu == 0` branch. | `7e1a882` | OOB-0 |
| 5 | `calsoft_hyd.f90:112` | `chg_val` | FALSE POSITIVE | none needed | `cal_adj = chg_val` runs only when `isim > 0`; `isim = 1` is set in the same block that assigns `chg_val` unconditionally. | `7e1a882` |  |
| 6 | `calsoft_hyd_bfr_et.f90:105` | `chg_val` | FALSE POSITIVE | none needed | Same `isim` pattern as calsoft_hyd: `isim = 1` (line 50) and `chg_val` (56) are set together. | `7e1a882` |  |
| 7 | `calsoft_hyd_bfr_latq.f90:103` | `chg_val` | FALSE POSITIVE | none needed | Same `isim` pattern: `isim = 1` (49) and `chg_val` (55) are set together. | `7e1a882` |  |
| 8 | `calsoft_hyd_bfr_perc.f90:119` | `chg_val` | FALSE POSITIVE | none needed | Same `isim` pattern: `isim = 1` (50) and `chg_val` (57) are set together. | `7e1a882` |  |
| 9 | `calsoft_hyd_bfr_surq.f90:102` | `chg_val` | FALSE POSITIVE | none needed | Same `isim` pattern: `isim = 1` (48) and `chg_val` (54) are set together. | `7e1a882` |  |
| 10 | `calsoft_sed.f90:145` | `diff` | FALSE POSITIVE | none needed | Read inside the slope block's `ireg`/`ilum` loops, which have the same bounds as the first block's loops, where `diff` is assigned unconditionally (line 45) on every iteration. | `7e1a882` |  |
| 11 | `ch_read_elements.f90:70` | `ires` | REAL | `231d0c3` | Any `in_regs%def_cha` row with `nspu = 0`: `ires` is never assigned in the routine. | `7e1a882` | OOB-0 |
| 12 | `ch_read_parms_cal.f90:33` | `mchp` | REAL | `ba87e47` | `ch_sed_parms.sft` exists but ends before its count line (EOF on the title or the count read). | `7e1a882` |  |
| 13 | `ch_rtday.f90:226` | `addp` | FALSE POSITIVE | none needed | `p = addp` only in the `volrt > maxrt` branch, where `sdti = maxrt < volrt` makes the `Do While` body run at least once and assign `addp`; `p` is set on every other path (lines 105, 111, 175). | `7e1a882` |  |
| 14 | `ch_temp.f90:398` | `ilsu` | REAL | `1943133` | `w_temp(0)%sf_on == 1`, any day but jday 366, a channel with no incoming `ru` object (`ru_count = 0`). Upstream has no initializer either, so this predates the strip. | none (no initializer upstream either) |  |
| 15 | `chg_par.f90:52` | `chg_par (result)` | REAL | `7d80016` | A change type other than absval/abschg/pctchg/relchg (misspelt or differently cased) in calibration.cal, ls_parms.sft or a decision-table cn2 change option. No initializer upstream either. | none (no initializer upstream either) |  |
| 16 | `co2_read.f90:114` | `co2_end` | REAL | `ea10cf9` | `co2_yr.dat` whose years have gaps or are out of order, so the record index starts past the last row (`yrc_start - first_year + 1 > rows`) while `yrc_start` is before the last listed year. | `7e1a882` |  |
| 17 | `command_object.f90:370` | `iob` | REAL | `8ec31df` | A delivery-ratio (`dr`) object with no incoming hydrographs (`rcv_tot = 0`) while constituents are simulated (`cs_db%num_tot > 0`). | `8fa9083` | OOB-0 |
| 18 | `conditions.f90:148` | `ipl` | REAL | `92295d7` | A `phu_plant` condition on an object with an empty plant community (`npl = 0`), with no earlier condition in the same call having set `ipl`. | `8fa9083` | OOB-0 |
| 19 | `conditions.f90:191` | `iob` | REAL | `92295d7` | A `phu_base0`, `precip_cur` or `precip_next` condition whose `ob` is neither `hru` nor `hlt` (reads at 161, 176, 191). | `8fa9083` | OOB-0 |
| 20 | `conditions.f90:220` | `ialt` | REAL | `92295d7` | A `plant_name_gro` condition naming a plant that is not in the community, with no earlier `do ialt` loop in the same call. | `8fa9083` | OOB-0 |
| 21 | `conditions.f90:835` | `targ_val` | REAL | `92295d7` | A `vol_wet` condition whose `lim_var` is neither `pvol` nor `evol` (reads at 831-837). | `8fa9083` |  |
| 22 | `flow_dur_curve.f90:124` | `iprv` | FALSE POSITIVE | none needed | The even-years branch is first reached at `iyr = nbyr/2 + 1 >= 2`, after the previous iteration assigned `iprv`; an odd `nbyr` exits through the other branch. | `7e1a882` |  |
| 23 | `gwflow_output.f90:138` | `in_gw` | REAL | `48416ae` | Any gwflow run with `gwflow.wbgroups` present: `in_gw` is never assigned, so the open/read/close at 138-178 used an undefined unit. | `bcb882c` |  |
| 24 | `gwflow_pump_allo.f90:58` | `ob_id_src` | REAL | `ea4f09b` | Every call (gwflow pumping through water allocation): `ob_id_src` is never assigned. Upstream's `= 0` pinned every call to the self-supply branches. | `8fa9083` |  |
| 25 | `gwflow_pump_allo.f90:401` | `cell_id` | REAL | `ea4f09b` | Single-source branch on an unstructured grid; reachable at 7eca087 only through a garbage nonzero `ob_id_src`, and unreachable again once it is 0. | `8fa9083` |  |
| 26 | `gwflow_read.f90:514` | `active_cell` | REAL | `3e3580c` | The first boundary cell (stat 2) with no active cell within 1.0e6 m, e.g. a grid with no active cells. | `7e1a882` | OOB-0 (downstream) |
| 27 | `header_cs.f90:421` | `ipest` | REAL | `ff2bafa` | Would be reached with `pco%hyd%a = y`, `num_paths > 0` and no pesticide header written (`num_pests = 0`). header_cs is not called at 7eca087. | `8a3bf09` | OOB-0, dead@7eca087 |
| 28 | `hru_dtbl_actions_init.f90:50` | `id` | FALSE POSITIVE | none needed | The two `if`s cover every `iauto`: `irr_trn_iauto` is never assigned anywhere in src/ (it keeps its default 0) and `iauto >= 1`. | `8fa9083` |  |
| 29 | `hru_lte_control.f90:120` | `r2` | REAL | `3634d9f` | An HRU-lte with `sw + exp(xx) <= 0.001` (soil water below 0.001 mm and `xx = wrt1 - wrt2*sw` at or below about -6.9). | `8fa9083` |  |
| 30 | `hru_lte_read.f90:252` | `sd_sl` | REAL | `f7f37d1` | Every HRU-lte whose input `tc` is below 1.e-6: `sd_sl` is never assigned (`ch_sl` is assigned on the line before). | `7e1a882` | FPE |
| 31 | `hyd_read_connect.f90:33` | `ob1` | REAL | `3142515` | An object type with count > 0 whose connect file is missing or ends before its header line: the weather crosswalk loop at 337 reads `ob1` (flagged at its declaration) and `ob2`. | `7e1a882` | OOB-0 |
| 32 | `hyd_read_connect.f90:264` | `nmetals` | FALSE POSITIVE | none needed | Assigned (110) under `if (cs_db%num_tot > 0)` earlier in the same loop iteration; read under the same test (248). | `7e1a882` |  |
| 33 | `hyd_read_connect.f90:269` | `nsalts` | FALSE POSITIVE | none needed | Assigned (117) under `if (cs_db%num_tot > 0)` earlier in the same loop iteration; read under the same test. | `7e1a882` |  |
| 34 | `hyd_read_connect.f90:280` | `ncs` | FALSE POSITIVE | none needed | Assigned (144) under `if (cs_db%num_tot > 0)` earlier in the same loop iteration; read under the same test. | `7e1a882` |  |
| 35 | `ls_read_parms_cal.f90:34` | `mlsp` | REAL | `ad8f5b2` | `ls_parms.sft` (`in_chg%wb_parms_sft`) exists but ends before its count line (reads at 34 and 36). | `7e1a882` |  |
| 36 | `lsu_output.f90:76` | `iob` | REAL | `1658178` | The first LSU(s) have no elements (`lsu_out(ilsu)%num_tot = 0`) and any LSU output is on (every gis_id write through 316). | `7e1a882` | OOB-0 |
| 37 | `mgt_sched.f90:255` | `icom` | FALSE POSITIVE | none needed | `icom` is assigned at line 72 for every operation except `fert`; it is read only in the kill/harvest cases. | `8fa9083` |  |
| 38 | `mgt_sched.f90:510` | `idp` | REAL | `c02f025` | Any `dwm` operation with `pco%mgtout = y`: mgt_sched handles one operation per call, and only plnt/harv/kill assign `idp`. | `8fa9083` |  |
| 39 | `mgt_sched.f90:605` | `ipud` | REAL | `c02f025` | A `pudl` operation whose `op_plant` matches no `puddle.ops` name (also 607). | `8fa9083` |  |
| 40 | `nut_orgnc2.f90:156` | `koc` | REAL | `2d916e4` | Layer-1 microbial C at most .01 while a lower layer has microbial C plus percolate of at least .01 and flow (reads at 155-156). | `8fa9083` |  |
| 41 | `recall_read_cs.f90:200` | `istep` | REAL | `1c80240` | Would be reached (a) with no record for the first simulation year (EOF; gfortran's implied-do is guarded and does not evaluate `istep`), or (b) `typ` outside 0-3 with a matching record. recall_read_cs is not called at 7eca087. | `7e1a882` | OOB-0, dead@7eca087 |
| 42 | `recall_read_salt.f90:201` | `istep` | REAL | `0241473` | Same two paths as recall_read_cs. recall_read_salt is not called at 7eca087. | `7e1a882` | OOB-0, dead@7eca087 |
| 43 | `reg_read_elements.f90:71` | `ireg` | REAL | `90474ff` | Would be reached on every pass that reads `in_regs%def_reg` past its count line. reg_read_elements is not called at 7eca087 (its call in proc_cal.f90:25 is commented out). | `7e1a882` | dead@7eca087 |
| 44 | `res_control.f90:83` | `irrig_track_b` | REAL (report only) | not fixed (follow-up) | Every call: `irrig_track_b` is cross-call state (the last demand tracker) that upstream kept through SAVE. A reset at entry would change release numerics for models using `res_irr_dmd`, so per Vahid's ruling it is reported, not fixed. **Lane O:** in the shipped ifx binary it is not static, so there it is undefined on every call. | `8fa9083` | follow-up |
| 45 | `res_cs.f90:93` | `v_settle` | FALSE POSITIVE | none needed | Assigned for `ics` = 1, 2, 3 on every iteration; the loop always runs 1 to n, so for `ics >= 4` it holds iteration 3's value and is never undefined. | `8fa9083` |  |
| 46 | `res_cs.f90:109` | `k_react` | FALSE POSITIVE | none needed | Same as `v_settle` (assigned for `ics` 1-3, loop runs in order). | `8fa9083` |  |
| 47 | `res_hydro.f90:172` | `ic` | REAL | `2166d18` | A `dyrt1` release whose `file_pointer` is not con1/con2/con3 (reads at 172-182). | `8fa9083` | OOB-0 |
| 48 | `res_hydro.f90:184` | `b_lo` | REAL | `2166d18` | A `dyrt1` release whose condition `lim_op` is outside `= * + - /`. | `8fa9083` |  |
| 49 | `res_hydro.f90:207` | `demand` | REAL | `2166d18` | An `irrig_trn` release whose `file_pointer` is neither `wro` nor `hru`. | `8fa9083` |  |
| 50 | `smp_grass_wway.f90:124` | `sf_depth` | REAL | `a35e05e` | `0 < sf_area <= 1.e-6` m2: the read guard is `sf_area > 0`, the assignment guard `sf_area > 1.e-6`. | `8fa9083` | FPE |
| 51 | `smp_grass_wway.f90:128` | `sf_sed` | REAL | `a35e05e` | Same `0 < sf_area <= 1.e-6` window as `sf_depth`. | `8fa9083` |  |
| 52 | `smp_grass_wway.f90:139` | `rcharea` | REAL | `a35e05e` | Surface runoff with `qp_cms <= 0` or NaN (hru km = 0, or tc_gwat < 0), so the depth iteration never runs (reads at 139-140). | `8fa9083` |  |
| 53 | `soil_nutcarb_init.f90:85` | `psp` | REAL | `34e5941` | `bsn_cc%sol_P_model == 1` and a soil whose first layer(s) have clay <= 0 (reads at 85, 87, 93). | `8fa9083` |  |
| 54 | `swift_output.f90:215` | `idb` | FALSE POSITIVE | none needed | The second channel loop has the same bounds as the first, which assigns `idb` on every iteration. | `7e1a882` |  |
| 55 | `unit_hyd.f90:109` | `q` | REAL | `dff7887` | `bsn_cc%uhyd` other than 0 or 1 in codes.bsn (reads at 109-114). | `7e1a882` | FPE |
| 56 | `urban_parm_read.f90:53` | `imax` | REAL | `2c9a4c7` | `urban.urb` exists but ends before its header line. | `aeebf09` (earlier: removed `7e1a882`, restored `ef4dcc4`) |  |
| 57 | `wet_cs.f90:96` | `v_settle` | FALSE POSITIVE | none needed | Same as res_cs `v_settle`. | `8fa9083` |  |
| 58 | `wet_cs.f90:112` | `k_react` | FALSE POSITIVE | none needed | Same as res_cs `k_react`. | `8fa9083` |  |
| 59 | `wet_initial.f90:152` | `iprop` | REAL | `52689c2` | (a) The HRU's hru-data `surf_stor` name matches no wetland.wet entry; (b) an `impound_on` action on an HRU whose hru-data `surf_stor` is "null" (actions sets `dbs%surf_stor`, but wet_initial crosswalks by the name). Reads at 152, 159, 177. | `8fa9083` | OOB-0 |
| 60 | `wet_irrp.f90:48` | `iob` | REAL | `093af86` | A paddy-irrigated HRU (`paddy_irr > 0`, ponding below `irr_hmin`) that is in no routing unit and has no irrigation source in its management (`irr_isc = 0`). | `8fa9083` | OOB-0 |

### Read in the same statements but not flagged by gfortran

| Variable | Routine | Read at | Fix commit | Input condition |
|---|---|---|---|---|
| `istr1` | actions | 1129-1218 (the same write as `istr` at 1200, plus the structure installs) | `418f098` | `grassww_install` (its assignment at 1181 is commented out) and `grassww_uninstall` always. The other installs (tile, septic, fstrip, user_def) when `file_pointer` matches no database entry. |
| `ob2` | hyd_read_connect | 337 (`do i = ob1, ob2`) | `3142515` | Same as `ob1`. |

These two show that gfortran's warning list is a lower bound, not a census.

## Warnings left after the fixes

| Warning (post-fix numbering) | Variable | Expected? |
|---|---|---|
| `calsoft_hyd.f90:112` | `chg_val` | yes |
| `calsoft_hyd_bfr_et.f90:105` | `chg_val` | yes |
| `calsoft_hyd_bfr_latq.f90:103` | `chg_val` | yes |
| `calsoft_hyd_bfr_perc.f90:119` | `chg_val` | yes |
| `calsoft_hyd_bfr_surq.f90:102` | `chg_val` | yes |
| `calsoft_sed.f90:145` | `diff` | yes |
| `ch_rtday.f90:226` | `addp` | yes |
| `flow_dur_curve.f90:124` | `iprv` | yes |
| `hru_dtbl_actions_init.f90:50` | `id` | yes |
| `hyd_read_connect.f90:269` | `nmetals` | yes |
| `hyd_read_connect.f90:274` | `nsalts` | yes |
| `hyd_read_connect.f90:285` | `ncs` | yes |
| `mgt_sched.f90:260` | `icom` | yes |
| `res_control.f90:83` | `irrig_track_b` | yes |
| `res_cs.f90:93` | `v_settle` | yes |
| `res_cs.f90:109` | `k_react` | yes |
| `swift_output.f90:215` | `idb` | yes |
| `wet_cs.f90:96` | `v_settle` | yes |
| `wet_cs.f90:112` | `k_react` | yes |

## Neutrality check output

```
OK  actions.f90              irrop      inserted@114 first-exec-ref@114
OK  actions.f90              ipud       inserted@115 first-exec-ref@115
OK  actions.f90              istr       inserted@116 first-exec-ref@116
OK  actions.f90              istr1      inserted@117 first-exec-ref@117
OK  aqu_read_elements.f90    ihru       inserted@31 first-exec-ref@31
OK  ch_read_elements.f90     ires       inserted@31 first-exec-ref@31
OK  ch_read_parms_cal.f90    mchp       inserted@17 first-exec-ref@17
OK  ch_temp.f90              ilsu       inserted@114 first-exec-ref@114
OK  chg_par.f90              chg_par    inserted@38 first-exec-ref@38
OK  co2_read.f90             co2_end    inserted@51 first-exec-ref@51
OK  command_object.f90       iob        inserted@72 first-exec-ref@72
OK  conditions.f90           ipl        inserted@66 first-exec-ref@66
OK  conditions.f90           iob        inserted@67 first-exec-ref@67
OK  conditions.f90           ialt       inserted@68 first-exec-ref@68
OK  conditions.f90           targ_val   inserted@69 first-exec-ref@69
OK  gwflow_output.f90        in_gw      inserted@34 first-exec-ref@34
OK  gwflow_pump_allo.f90     ob_id_src  inserted@54 first-exec-ref@54
OK  gwflow_pump_allo.f90     cell_id    inserted@55 first-exec-ref@55
OK  gwflow_read.f90          active_cell inserted@200 first-exec-ref@200
OK  header_cs.f90            ipest      inserted@16 first-exec-ref@16
OK  hru_lte_control.f90      r2         inserted@96 first-exec-ref@96
OK  hru_lte_read.f90         sd_sl      inserted@65 first-exec-ref@65
OK  hyd_read_connect.f90     ob1        inserted@50 first-exec-ref@50
OK  hyd_read_connect.f90     ob2        inserted@51 first-exec-ref@51
OK  ls_read_parms_cal.f90    mlsp       inserted@18 first-exec-ref@18
OK  lsu_output.f90           iob        inserted@23 first-exec-ref@23
OK  mgt_sched.f90            idp        inserted@66 first-exec-ref@66
OK  mgt_sched.f90            ipud       inserted@67 first-exec-ref@67
OK  nut_orgnc2.f90           koc        inserted@60 first-exec-ref@60
OK  recall_read_cs.f90       istep      inserted@38 first-exec-ref@38
OK  recall_read_salt.f90     istep      inserted@38 first-exec-ref@38
OK  reg_read_elements.f90    ireg       inserted@36 first-exec-ref@36
OK  res_hydro.f90            ic         inserted@81 first-exec-ref@81
OK  res_hydro.f90            b_lo       inserted@82 first-exec-ref@82
OK  res_hydro.f90            demand     inserted@83 first-exec-ref@83
OK  smp_grass_wway.f90       rcharea    inserted@71 first-exec-ref@71
OK  smp_grass_wway.f90       sf_depth   inserted@72 first-exec-ref@72
OK  smp_grass_wway.f90       sf_sed     inserted@73 first-exec-ref@73
OK  soil_nutcarb_init.f90    psp        inserted@34 first-exec-ref@34
OK  unit_hyd.f90             q          inserted@71 first-exec-ref@71
OK  urban_parm_read.f90      imax       inserted@18 first-exec-ref@18
OK  wet_initial.f90          iprop      inserted@39 first-exec-ref@39
OK  wet_irrp.f90             iob        inserted@26 first-exec-ref@26
violations: 0
```

## Follow-ups (not done here: each one changes behaviour)

**Index 0 is still out of bounds** on the uninitialised path. The upstream value was restored as ruled; the intended logic is a separate change:

- `aqu_read_elements.f90:117`: `acu_reg(i)%num(0)` is a write below `num(1:sp_ob%hru)`. The loop presumably meant `num(iaqu)`, with an allocation sized by `sp_ob%aqu`.
- `ch_read_elements.f90:70`: `ccu_out(i)%num(0)` is a write below `num(1:sp_ob%hru)`. The loop presumably meant `num(icha)`.
- `command_object.f90:370`: `ob(0)%props` is a read below `ob(1:)`. The dr object's own index is `icmd`, which the line above already uses (`dr(ob(icmd)%props)`).
- `conditions.f90:148`: `pcom(ob_num)%plcur(0)`: with `npl = 0` no index is valid. The condition should give "n" without reading `phuacc`.
- `conditions.f90:191`: `ob(0)%wst` is a read below `ob(1:)`, and the station index read from it then indexes `wst`. A weather condition on another object type needs that object's own station.
- `conditions.f90:220`: `act_hit_tl(0)` is a write below `act_hit_tl(1:)`. The plant-absent case presumably meant to set every alternative to "n" in a `do ialt` loop.
- `gwflow_read.f90:514`: 0 is stored into `gw_bound_near`, and `gwflow_rech` then writes `gw_hyd_ss(0)`. A boundary cell with no active cell in range needs a defined fallback.
- `header_cs.f90:421`: `cs_path_solsor(0)` is a read below `(1:num_paths)`. The implied-do variable is `ipath`, the intended index; fixing it changes header text. Dead code at 7eca087.
- `hyd_read_connect.f90:33`: With `ob1 = ob2 = 0` the crosswalk loop runs once on `ob(0)` and writes `ob(0)%wst` through `search`. A missing or empty connect file should stop with a message.
- `lsu_output.f90:76`: `ob(0)%gis_id` is a read below `ob(1:)`. Moving `iob = sp_ob1%ru + ilsu - 1` out of the element loop is exact for LSUs with elements and also correct for empty ones.
- `recall_read_cs.f90:200`: `hd_cs(0, iyrs)` is below the allocation, which starts at 1. Dead code at 7eca087.
- `recall_read_salt.f90:201`: `hd_salt(0, iyrs)` is below the allocation, which starts at 1. Dead code at 7eca087.
- `res_hydro.f90:172`: `d_tbl%cond(0)` is a read below `cond(1:conds)`. An unknown `dyrt1` `file_pointer` should be rejected when the decision tables are read.
- `wet_initial.f90:152`: `wet_dat(0)` is a read below `wet_dat(1:imax)` when wetland.wet exists. For path (b), `hru(iihru)%dbs%surf_stor`, set by `impound_on`, holds the intended index.
- `wet_irrp.f90:48`: `ob(0)%src_tot` is a read below `ob(1:)`, and it then bounds a loop over `ob(0)%obtyp_out`. With no routing unit the code has already set the `unlim` source, which is presumably what should be used.

**Traps under `-fpe0`** on the uninitialised path (measured as above):

- `hru_lte_read.f90:252`: `0.**(-.385)` raises IEEE divide-by-zero for every HRU-lte with input `tc < 1.e-6` (without traps, `tc = +Inf`). The operand is meant to be `ch_sl`, assigned on the line before.
- `smp_grass_wway.f90:124`: `Log(sf_depth)` with `sf_depth = 0.` raises divide-by-zero when `0 < sf_area <= 1.e-6` (without traps, -Inf, and `surq_remove` is clamped to 95). The assignment guard and the read guard should agree.
- `unit_hyd.f90:109`: With `q = 0`, `sumq = 0` and the normalisation `uh / sumq` is 0/0, which raises invalid (NaN without traps). `bsn_cc%uhyd` should be validated where codes.bsn is read.

**Report only:** `res_control.f90:83` `irrig_track_b` is cross-call state: the last `res_ob(jres)%irrig_track` seen, which upstream kept through SAVE. It is shared across reservoirs even upstream. Every call at 7eca087 reads it before setting it. A reset at entry would count `res_irr_dmd` demand every day once the tracker has advanced, and `daily_demand` feeds `d_mon_past`, which the `res_hydro` release rule reads, so that is a numerics change. The thread-safe fix is a per-reservoir tracker component in `res_ob`, which is a reservoir_module type change. **Lane O (engine), relayed by Lane P on 2026-10-05:** in the shipped ifx binary built from `7eca087` none of the 44 routine/variable pairs is static, `irrig_track_b` included. So "kept across calls" describes upstream, not what we ship: there the value is undefined on every call. Not verified here (ifx is not installed in this session). Not changed; Vahid rules.

**Seen in passing, not uninitialised reads:** `swift_output` writes the last channel's name on every `chan_dr.swf` row (`idb` from the previous loop). `lsu_output` prints the previous LSU's `gis_id` for an LSU with no elements. `calsoft_sed` compares `hru(ihru)` (the module-level `ihru`) where `iihru` is meant. `res_cs`/`wet_cs` reuse constituent 3's rates for any constituent beyond 3. `aqu_read_elements` inquires `def_aqu_reg` but opens `def_aqu`.

## Reproduce

```sh
git checkout --detach 7eca087 && git cherry-pick c80ab8ce     # or: fix/q1281-uninit-locals + c80ab8ce
sed -i 's/ -finit-local-zero -fno-unsafe-math-optimizations -frounding-math//' CMakeLists.txt
sed -i 's/set(frelease        "-O")/set(frelease        "-O2 -Wmaybe-uninitialized -Wuninitialized")/' CMakeLists.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran   # needs netcdf-fortran
make -C build -j1 2>&1 | grep -B4 'uninitialized\]'
```

## Limits

- The instrument is one compiler at one optimisation level. `-Wmaybe-uninitialized` depends on gfortran's data-flow analysis and misses derived-type components, arrays and some paths. `istr1` and `ob2` were found by reading, not by the compiler. Treat the 60 as a lower bound.
- Input conditions come from reading the code. None was reproduced by running a model, which was out of scope.
- Reachability (`dead@`) is static: the linker discards code nothing references. A kept routine is referenced from `main`, which does not prove the particular path runs.
- The production build is ifx with `-fpe0`, which is not installed here. The trap behaviour was measured with gfortran's equivalent flags, and gfortran's EOF guard on the recall implied-do was read from its tree dump, not from ifx output.
