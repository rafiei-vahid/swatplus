# Contribution candidates: fork fixes that swat-model/swatplus does not have

Date: 2026-10-11.

Target: upstream main 80049bc9, the merge of #299 on 2026-10-08. I also checked upstream dev f2f6a0ad, which is 20 commits ahead of main and 5 behind. It carries the #254 sync and the 10/01 water-allocation and plants.plt work.

Fork: rafiei-vahid/swatplus, all 70 branches.

Nothing here has been posted upstream. Each PR needs Vahid's approval of that exact PR first.

## How this was checked

- **Inventory.** I took the non-merge commits on every fork branch that are not reachable from upstream/main: 230 commits, deduplicated by patch. The four 2024 legacy branches (Development, jga_develop, CMakePreset_json, gfortranchanges041924) were left out. Commits that are upstream's own work, reaching the fork through a dev-based branch, are not counted.
- **Applies.** `git apply --check` of the commit's `src/` diff against upstream main. "CONFLICT" mostly means the fork's OpenMP initialiser strip changed the context lines. The bug itself was then judged on upstream's own code.
- **Builds.** gfortran 13.3, the upstream CMake Release flags (`-cpp -fimplicit-none -finit-local-zero -frounding-math -O`). Each touched file was patched onto upstream main and compiled against the modules of a full upstream-main build. A full build was done for every candidate that has numbers below.
- **Evidence.** Runs of upstream main versus upstream main plus the one fix, on upstream's refdata (Ames_sub1 1975-2020, Osu_1hru 2010-2023). Where no refdata deck reaches the code, I used a reproduction deck made from Osu_1hru with a few lines edited. The scripts are in `repro/`. Numbers are average-annual (`*_aa.txt`).
- **Upstream status.**
  - Upstream PR pull refs 1-301, fetched with git. None of the heads that contain 109179b carries any of the fixes below.
  - The PR list and issue list pages, 60 issues.
  - "Superseded" means upstream dev already holds an equivalent change.

### gfortran caveat

Upstream's gfortran build zeroes uninitialised locals (`-finit-local-zero`). Locals with a declaration initialiser are implicitly SAVE. So most "uninitialised local" fixes show no difference on a gfortran Release build. Under ifx (no zero-init, `/fpe0`) or the Windows gfortran `-fcheck=all` build they read garbage, crash, or carry the previous call's value.

### Debug build

The upstream-main Debug build already traps on shipped Osu_1hru: SIGFPE at `aqu_1d_control.f90:160`, 2011-01-04. This is underflow; issue #288 tracks it. So the debug checks below used Release plus `-fcheck=bounds,pointer`.

## Ranking (value to SWAT+ users against diff size)

| # | Candidate | Fork commit(s) | Upstream file(s) | Diff | Applies to main / dev | gfortran | Evidence | Upstream status |
|---|---|---|---|---|---|---|---|---|
| 1 | **Wetland: a zero-height weir.res row sets every wetland's principal volume to 0** | 5008b44 (branch fix/wetland-zero-height-weir-row); alt. 9b3b0b5 | wet_initial.f90 | +8 −1 | clean / clean (file identical on dev) | OK (full build, plus bounds-checked build) | repro: water yield 576.7 → 387.9 mm, wetland storage 0.1 → 80.1 mm; Ames_sub1: all 67 output files identical; Osu_1hru as shipped: basin_wb_aa identical, largest aa change 0.8 % | no issue, no PR |
| 2 | **Reservoir: res_control lost its dom/mon/end_of_mo/daily_inflow assignments (109179b), so h06 release rules never release** | 6eb90b8 | res_control.f90 | +8 (4 code lines) | clean / clean (not fixed on dev) | OK (full build) | repro: h06 reservoir outflow 0 → 2.41e5 m³/yr; storage 1.84e7 → 1.00e5 m³; Ames_sub1 (no reservoirs): all 67 output files identical | no issue, no PR |
| 3 | Junk files every run: `null`, `wtps_wuses.wal`, `paddy_test.csv` | 8083d81, 9d25269 | water_allocation_read.f90 (2 inverted inquires), wallo_wtps_wuses.f90, main.f90.in | ~+13 −3 | 8083d81 clean; 9d25269 conflicts (re-port; small) | OK (8083d81) | all three files appear in every upstream run, refdata included | none found (may fit under #291 "Housekeeping after #254") |
| 4 | Water temperature: ch_temp's result never reaches the outflow hydrograph or the outflow-temp column | 56c8dc4 (1 line); then 60bf61c (output column), fa0e2ba (gage window), the ch_temp mixing part of 0e1db00 | sd_channel_control3.f90 (+1); sd_channel_output.f90; ch_temp.f90 | +1 for the core fix | 56c8dc4 conflicts only on context (1 line, re-typed and built on main) | OK (full build of the 1-line version) | Osu_1hru: channel outflow temperature, average-annual, 19.40 °C (air regression 5+0.75·tave) → 24.42 °C (ch_temp's answer) | no issue, no PR |
| 5 | nut_orgnc2: `c_microb_sed` keeps the previous HRU's value (SAVE) when there is no erosion | 221f41f | nut_orgnc2.f90 | +1 code line | clean / **already fixed on dev** (`nut_orgnc2.f90:102`) | OK (full build) | Ames_sub1 (carbon = 2): all 67 output files identical (the stale path is not hit on this deck) | superseded by dev; no PR needed, arrives with the next dev → main merge |
| 6 | chg_par: an unrecognised change type returns an undefined function result | 7d80016 | chg_par.f90 | +4 | clean | OK | static: `-finit-local-zero` does not cover function results, so garbage with either compiler, then clamped to [absmin, absmax] | none found |
| 7 | ch_temp: gage-record window not guarded, so a gage that starts after, or ends before, the simulation reads out of bounds | fa0e2ba | ch_temp.f90 | +9 | clean | OK | static; not reached by refdata (Osu's record covers the run) | none found |
| 8 | res_hydro dyrt1: hard-coded cond(5/4/3) read past a table with fewer conditions | f0343c0 (+ reset 2166d18) | res_hydro.f90 | +13 (+7) | f0343c0 clean; 2166d18 conflicts | OK | static | none found |
| 9 | HRU-lte Kirpich tc divides by a slope that is never set (sd_sl = 0), so tc = Inf or an FP trap | f7f37d1 + 0eba3c2 | hru_lte_read.f90 | +9 −8 | conflict (re-port) | n/a until ported | static | none found |
| 10 | unit_hyd: uhyd other than 0/1 leaves q unset, so NaN unit hydrograph on sub-daily runs | dff7887 + f5f1c7e | unit_hyd.f90 | +12 | dff7887 clean; f5f1c7e conflicts | OK (dff7887) | static | none found |
| 11 | ch_rtmusk: `inflo` keeps the previous call's value on zero-inflow days (feeds the rating curve) | 001d2b4 (part) | ch_rtmusk.f90, sd_channel_control3.f90 | ~4 code lines | clean | OK | static; ephemeral reaches only | none found |
| 12 | pl_nut_demand: `uapd(ipl)` reset with a stale `ipl` (another HRU's demand at dormancy exit; can index mpc+1) | 6322a00 | pl_nut_demand.f90 | +9 −2 | clean | OK | static; perennials and cold annuals | none found |
| 13 | mgt_sched: `itr` (SAVE) is never reset, so a plnt op with op3 = 0 repeats the last transplant | 3418538 | mgt_sched.f90 | +6 −1 | clean | OK | static; models with transplant.ops | none found |
| 14 | hyd_read_connect: a missing first connect file writes ob(0) and crashes later with no message | 3142515 + 4bd1663 | hyd_read_connect.f90 | +14 | 4bd1663 clean; 3142515 conflicts | OK (4bd1663) | static | none found |
| 15 | nut_psed: fr_orgp/fr_actmin/fr_stamin carry the previous call's values when sedp_attach ≤ 1e-9 | c80ab8c | nut_psed.f90 | +7 | clean | OK | static; tiny steep HRUs with exhausted layer-1 P (ifx crash on the fork) | none found |
| 16 | pl_nut_demand `delg_p` and plant_init `igrow` (SAVE, rare paths) | b295c23, b94dee9 | pl_nut_demand.f90, plant_init.f90 | +1, +4 | clean | OK | static; not reached by refdata | none found |
| 17 | co2_read: a co2_yr.dat with year gaps gives CO2 = 0 ppm (no biomass) | ea10cf9 | co2_read.f90 | +4 | clean | OK | static. The fork's reset keeps the 0, so this needs a real fix (hold the last value) before it goes upstream | none found |
| 18 | Rest of the Q1281 audit (fix/q1281-uninit-locals, fix/q1281-followups) | the other ~40 of 48 commits | many | mostly +4 each | mixed | mixed | static; see notes | none found |

The rank is my judgement:

- **1 and 2:** wrong water balances on ordinary inputs; ~8 lines each; numbers on decks made from upstream refdata.
- **3:** touches every run but is cosmetic.
- **4:** high value but changes reported temperatures in every channel model. It needs observed-temperature evidence before upstream will take it. I'd open it as an issue with the Osu numbers first.
- **5:** already fixed on dev.
- **6-17:** real but narrow; most need ifx, a rare input, or a bad input to show. Except where marked, the patches also apply to dev.

## Notes per candidate

### 1. Wetland zero-height weir row (top-ranked; Task 2 branch built)

Upstream main, `wet_initial.f90:138`:

```fortran
if (db_mx%res_weir > 0.and.wet_ob(iihru)%iweir>0) then
  ...
  wet_ob(iihru)%pvol = hru(iihru)%area_ha * res_weir(iweir)%h * 10.**4
```

`wet_ob%iweir` defaults to 1 (`reservoir_module.f90:54`) and is only changed by a `weir` management op (`mgt_sched.f90:519`). The WEIR column of wetland.wet is not read. So whenever weir.res has any row, every wetland takes row 1's height.

SWAT+ Editor writes weir.res with one row of height 0. Upstream's own branch netCDF-Support has an Editor v2.3.3 deck, `data/netCDF-climate`, with exactly that: weir.res `shape1 1.84 2.60 2.50 0.00000`, and 142 wetland HRUs out of 5153. There every wetland gets pvol = 0. pdep from hydrology.wet is ignored, and the `wetland` release table (vol > pvol: release over 10 days) drains the wetland every day.

The fix takes the weir branch only when the referenced row has a height > 0. Rows with a real height, paddies included, are unchanged.

Reproduction: `repro/wetland_zero_weir.sh` (Osu_1hru with rel `weir` → `wetland` and lu `rice140_lum` → `agrl_lum`; weir.res as shipped).

basin_wb_aa:

| | upstream main | + fix |
|---|---|---|
| wateryld (mm) | 576.674 | 387.935 |
| surq_gen (mm) | 576.076 | 386.369 |
| perc (mm) | 224.722 | 474.658 |
| et (mm) | 373.407 | 312.210 |
| wet_evap (mm) | 21.976 | 232.557 |
| wet_stor (mm) | 0.099 | 80.107 |

wetland_aa:

| | upstream main | + fix |
|---|---|---|
| evap (m³) | 2.207e3 | 2.338e4 |
| seep (m³) | 1.704e5 | 7.299e4 |

Other runs:

- **Ames_sub1 (no wetlands):** all 67 output files identical (header line and simulation.out's revision/time lines excluded).
- **Osu_1hru as shipped (paddy, weir height set by weir_height ops; weir1 has height 0, so only the initial pvol and nutrient mass change):**
  - basin_wb_aa is identical.
  - Aquifer flow goes 365.531 → 365.517 mm.
  - The largest change in any `*_aa.txt` value is channel sedp_out, 4.840 → 4.801 (0.8 %).
- **Release with `-fcheck=bounds,pointer`:** the reproduction deck gives the same numbers, and shipped Osu_1hru runs to completion.

Alternative, 9b3b0b5: default `iweir = 0`. It is broader, because every wetland stops using weir.res unless a weir op sets it, and it changes paddies that rely on weir.res at init. Not proposed.

Left for the maintainers: whether wetland.wet's WEIR column should set `iweir` per wetland. It is in the file but not read.

### 2. Reservoir h06 inflow memory (109179b regression)

109179b (Jeff Arnold, 2026-05-07, "Enhance constituent mass tracking in water allocation") removed these lines from `res_control.f90`:

```fortran
dom          = time%day_mo
mon          = time%mo
end_of_mo    = time%end_mo
daily_inflow = ht1%flo
```

The reads below stayed. The four are now SAVEd locals stuck at 0, which has three effects:

- `I_mon_past` and `d_mon_past` are never updated, because `end_of_mo` is never 1.
- `daily_inflow_array` and `daily_demand_array` gain one element per day for the whole run, for every reservoir, because `dom` is never 1. Each day reallocates and copies the array.
- res_hydro's `nonirr-h06` and `irr-h06` get I_mean = 0, so c_ratio = smax/0 and a release of er·0 = 0. The reservoir never releases.

Upstream dev still lacks the lines.

Reproduction: `repro/reservoir_h06.sh` (Osu_1hru plus one reservoir below cha01; pvol 1.5 ha-m, evol 2.5 ha-m; one `nonirr-h06` action, alpha 0.85).

| | upstream main | + fix |
|---|---|---|
| flo_in (m³/yr) | 2.151e5 | 2.151e5 |
| flo_out (m³/yr) | 0 | 2.412e5 |
| flo_stor (m³) | 1.842e7 | 1.004e5 |
| area (ha) | 232.4 | 4.0 |

Basin water balance is identical (the reservoir is downstream of the HRU).

The fixed run's mean storage is still above evol. That is how h06 behaves with these hand-picked parameters, not part of this bug.

Ames_sub1 (no reservoirs): all 67 output files identical to main.

### 3. Junk files

- `water_allocation_read.f90:38` and `:192` test `.not. i_exist .or. name /= "null"`. With `water_rights null null` (both refdata decks) this opens, and so creates, a file literally named `null`.
- `wallo_wtps_wuses.f90:25` compares two literals (always true) and creates an empty `wtps_wuses.wal`.
- `main.f90.in:143` always creates `paddy_test.csv`.

All three are in the Ames_sub1 run directory after an upstream-main run. 8083d81 fixes the POD instance. 9d25269 (the other three) does not apply to main and needs a small re-port.

### 4. Water temperature

`sd_channel_control3.f90:337` copies `ob(icmd)%hd(1) = ht2` before `call ch_temp` (`:357`). So the temperature ch_temp computes (`ht2%temp = tw_final`, `ch_temp.f90:430`) never reaches `hd(1)`. Downstream channels and the outflow-temperature column carry ch_watqual4's `5 + 0.75·tave`. `ch_out_d%temp` is overwritten from `hd(1)` at `:371`/`:433`.

The core fix is one line after `call ch_temp`: `ob(icmd)%hd(1)%temp = ht2%temp`. On Osu_1hru the average-annual outflow temperature goes 19.40 → 24.42 °C.

Follow-ups:

- 60bf61c: the `water_temp` column prints the module scalar `wtemp`, so every channel shows the last-routed channel's value.
- fa0e2ba: the gage window.
- 0e1db00's ch_temp part: the tributary mixing loop keeps only the last inflow, reads `hin_d(0)` out of bounds on a first chandeg inflow, and uses `ht1 = ob(iob)%hd(1)`, the channel's own outflow.
- ch_temp reads only `w_temp(0)`, so temperature.cha rows are never used. That is a separate upstream bug, not fixed on the fork.

The fork's validated heat-budget work (Q380 impoundment heat budget 8685ac3, Q403 relaxation 6ab4865, Q462 0e1db00) is a feature on top of this. It needs its own evidence package (observed temperatures) and is not a small PR.

### 5. c_microb_sed

`nut_orgnc2.f90:55` `real :: c_microb_sed = 0.` is SAVEd. It is set only at `:135` (microb C > 0.01, flow and erosion > 0). On other days `:140` subtracts, and `:148` exports, the previous call's value, which can come from another HRU. carbon = 2 decks reach it; Ames_sub1 is one.

Upstream dev already resets it at the top of the routine (`c_microb_sed = 0.   !! must be reset each call`, `nut_orgnc2.f90:102`), so nothing to send. Ames_sub1 with the fork's one-liner on main: all 67 output files identical, so this deck never reaches the stale case.

### Q1281 audit (rows 6, 8-10, 14, 18)

Of the 48 commits on fix/q1281-uninit-locals and fix/q1281-followups, judged on upstream's code:

- **Already fixed upstream (4):** ch_temp `ilsu`, nut_orgnc2 `koc`, res_hydro `demand`, conditions `plant_name_gro` (8d804d39).
- **Code upstream never calls (6):** header_cs, recall_read_cs/salt (×2 each), reg_read_elements.
- **Exist only because of the fork's initialiser strip (5):** ch_read_parms_cal, gwflow_output_init, gwflow_pump_allo, ls_read_lsparms_cal, urban_parm_read.
- **The rest:** real or stale-SAVE on upstream.

The follow-ups RESULT.md says upstream main still has the original code at all 18 sites. That is wrong for `plant_name_gro`.

Five follow-up guards test for a 0 that only the first branch's reset provides: wet_irrp, hyd_read_connect, smp_grass_wway, wet_initial, lsu_output. Each must go upstream together with its reset.

Worth sending, in this order: chg_par (6), res_hydro dyrt1 (8), hru_lte Kirpich (9), unit_hyd (10), hyd_read_connect (14).

5234386 (weather conditions on non-HRU objects) applies cleanly but does not compile on main. It writes to `act_hit_tl`, a fork name that upstream's conditions.f90 does not have, so it needs a re-port.

Two more upstream bugs the audit found that nobody has filed:

- `grassww_install` never sets the waterway parameters: its `istr1` assignment is commented out.
- The dwm line of mgt.out reads `plcur(j)`/`tot(j)` with `j` the HRU number.

## Already upstream, already proposed, or not for upstream

### Already open as upstream PRs by Vahid

- **#292 plant_parm_read** (fix/q1072-plant-parm-read, bf27af5): parse each plants.plt row explicitly and stop on a bad row. Upstream dev now has be2c1b2 (Taci Ugraskan, 2026-10-01, "Check every plants.plt line and stop on a layout mismatch"), which does the same job in a different way. #292 is likely superseded when dev reaches main; worth a comment pointing to be2c1b2, or closing it.
- **#294 varinit enratio** (68d5324). dev's varinit.f90 is identical to main's, so it still applies.
- **#295 ch_watqual4 gra** (9688e50). dev changed the low-flow guard in the same file, so it may need a rebase.
- **#297 res_sediment nsed** (6132439), for issue #296. dev changed res_sediment.f90 (denormal guards), so it needs a rebase. It already conflicts with main.
- **#301 proc_bsn area_calc.out** (525aced), for issue #298. Same intent as the fork's 67757ff (RECL 2**24); #301's version is the cleaner one.

### Already in upstream

- **Q1083 et/esoil double count** (7eca087): a backport of upstream 2b5236d.
- **mgt_operatn same-date loop** (fix/mgt-operatn-samedate-loop): PR #236, merged.
- **varinit hhsedy row** (perf/varinit-hhsedy-row): PR #220, merged.
- **Q870 co2_annual%yrs** (944681e): upstream already has `integer :: yrs = 0`.
- **ef4dcc4, 03877e0:** fork-only or masked on upstream.

### Closed upstream earlier, unmerged

- #213 NetCDF output.
- #214 print_filter for channel_sd.
- #219 hru_read O(1) lookup.
- #241 and #243 OpenMP race fixes (issue #242, results depend on HRU and channel order, is still open).

### Fork-only features, not small PRs

- OpenMP HRU and wavefront parallelism, and its SAVE-strip, threadprivate and deterministic-reduction commits.
- NetCDF output (Q75, Q430, Q859, Q890).
- MODFLOW 6 coupling: Q1040, Q1130, Q1137, Q1144, and the 06/25 mf6 feats.
- PFAS module, including the daily sediment-PFAS fixes e707183 and 06b6b37.
- The impoundment heat budget, Q380/Q403.
- The 9/27 sync repairs (res_control is candidate 2; the area_calc RECL is #301; cdfout sediment-budget writers are fork-only).
- Build and README changes.

## Task 2: branch for candidate 1

- **Branch:** `claude/upstream-wetland-zero-height-weir`, from upstream main 80049bc9.
- **Contents:** one commit, Vahid's 5008b44 cherry-picked unchanged. The author is Vahid Rafiei. The committer is this session's git identity (Claude <noreply@anthropic.com>), so recommit it under your own name before opening the PR if that matters.
- **Diff:** `src/wet_initial.f90`, +8 −1. The tree is identical to the build the numbers above came from.
- **PR text:** `upstream-pr/wetland-zero-height-weir/PR_BODY.md`, on this branch.
- **Reproduction:** `repro/wetland_zero_weir.sh`. It stays on this branch so the upstream branch carries only the fix.

Run it as:

```
repro/wetland_zero_weir.sh refdata <main exe> <fixed exe>
```
