      subroutine res_temp (jres)
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Impoundment heat budget: warm (or cool) water passing through a reservoir /
!!    pond toward its equilibrium temperature, and hand the result to the outflow
!!    hydrograph so downstream channels see it.
!!
!!    Q380 (WT-M1c). Before this routine, res_*.f90 carried NO temperature state at
!!    all (grep '%temp' over all 29 files returned 0), so a pond passed inflow
!!    temperature through unchanged and the modelled July step across an impoundment
!!    was 0.00 degC by construction.
!!
!!    ~ ~ ~ WHY THE CLOSED-FORM EXPONENTIAL, NOT AN EULER STEP ~ ~ ~
!!    For a well-mixed pond under Edinger equilibrium exchange the outflow satisfies
!!
!!        T_out = T_eq + (T_in - T_eq) * exp(-tau_res / tau_thermal)
!!
!!    with residence time tau_res = A*h/Q and thermal response tau_thermal =
!!    rho*cp*h/K. The MEAN DEPTH h APPEARS IN BOTH AND CANCELS:
!!
!!        tau_res / tau_thermal = A*K / (Q*rho*cp)
!!
!!    That matters twice. First, this model has no bathymetry -- vol_ps/area_ps is
!!    exactly 10.0 for 206 of 265 ponds, a placeholder, not a survey -- so a scheme
!!    that depended on depth would be calibrating a constant. Second, at the real
!!    geometry of the White Cloud impoundment (MDNR 0192/2012: 50 acres, max 8 ft;
!!    20.2 ha, mean depth <= ~1.5 m) tau_res is about 0.90 d at the July mean flow of
!!    3.128 m3/s -- BELOW the daily timestep. An explicit Euler update is unstable in
!!    that regime; the closed form above is exact for a constant T_eq over the step
!!    and is unconditionally stable for any tau_res, which is why it is used here.
!!
!!    ~ ~ ~ WHAT THIS DOES NOT DO ~ ~ ~
!!    Single, fully mixed layer. It cannot reproduce a stratified pond: MDNR measured
!!    72.3 degF surface / 64.9 degF bottom at White Cloud on 9 Aug 1996. Two layers
!!    would need a mixing-depth rule, i.e. a second free parameter calibrated on three
!!    July means, on placeholder bathymetry -- rejected in the Q380 design note and by
!!    review. The drawdown-year signal comes from the AREA term above (draining cuts A,
!!    which cuts the approach fraction), not from a withdrawal-depth rule.
!!
!!    Design note: documents/ops/q380-impoundment-heat-budget-design-2026-09-02.md

      use reservoir_module
      use hydrograph_module
      use climate_module
      use water_body_module
      use channel_data_module

      implicit none

      integer, intent (in) :: jres        !none      |reservoir number
      real, parameter :: rho_cp = 4.186e6 !J/m3/degC |water density * specific heat
      real :: t_in                        !degC      |inflow (mixed) temperature
      real :: q_out                       !m3/s      |outflow rate
      real :: area_m2                     !m2        |pond surface area
      real :: t_ref                       !degC      |reference temp for the exchange coefficients
      real :: tdx_cal, tdx, t_md, ke_beta, f_wind, k_e
      real :: h_sr, e_s, e_a, cloud, e_atm, h_atm, t_equil
      real :: ssff                        !none      |shade factor
      real :: x_ratio                     !none      |A*K/(Q*rho*cp) -- depth-invariant
      real :: f_approach                  !none      |fraction of the way to equilibrium

      !! nothing to do without a surface to exchange across
      area_m2 = res_wat_d(jres)%area_ha * 10000.
      if (area_m2 < 1.e-6) return

      !! Reference temperature is the POND's, not the inflow's. res(jres) is a
      !! hyd_output and already carries %temp, and res_control.f90:76 has already done
      !! `res(jres) = res(jres) + ht1`, whose hydout_add mixes temperature FLOW-WEIGHTED
      !! (hydrograph_module.f90:1374). So by this point res(jres)%temp is the mixture of
      !! yesterday's storage and today's inflow -- the pond's thermal memory, which using
      !! ht1%temp alone would discard. No new state variable is needed for this.
      !! Cold start: %temp defaults to 0, so fall back to the inflow until the pond has a
      !! temperature of its own. With nyskip warm-up and tau_res ~1 d this flushes long
      !! before any printed year.
      t_in = res(jres)%temp
      if (t_in < 0.001) t_in = ht1%temp
      q_out = ht2%flo / 86400.            !! ht2%flo is m3/day

      !! Edinger et al. (1974) equilibrium temperature and surface exchange
      !! coefficient. This duplicates the block in ch_temp.f90 ON PURPOSE: sharing it
      !! would mean editing ch_temp, and the ship gate asserts channel results are
      !! byte-identical to production. A refactor that unifies the two belongs AFTER
      !! both paths are gated, not in the change that introduces the second one.
      t_ref = t_in
      tdx_cal = w%dewpt * w_temp(0)%hex_coef1
      tdx = amin1(tdx_cal, w%tave)
      t_md = (t_ref + tdx) / 2.
      ke_beta = 0.35 + 0.015 * t_md + 0.0012 * (t_md**2)
      f_wind = 9.2 + 0.46 * (w%windsp**2)
      k_e = 4.48 + 0.05 * t_ref + (ke_beta + 0.47) * f_wind

      !! open water: no channel shading. sf_on drives a per-landscape-unit shade table
      !! for channels; a pond surface is not shaded by its banks in the same way, so the
      !! basin default is used and the table is deliberately not consulted here.
      ssff = w_temp(0)%ssff
      h_sr = 0.97 * (w%solrad * 11.57) * (1. - ssff)

      e_s = 6.1275 * Exp(17.62 * w%tave / (237.3 + w%tave))
      e_a = e_s * w%rhum
      if (w%solradmx < 1.e-4) then
        cloud = 0.
      else
        cloud = 0.9 * (w%solrad / w%solradmx) + 0.1
      end if
      e_atm = 0.74 + 0.0065 * e_a * (1. + 0.17 * cloud**2)
      h_atm = 0.96 * e_atm * 5.67e-8 * (w%tave + 273.15)**4

      t_equil = tdx + h_sr / k_e + (h_atm - 305.5 - 4.48 * tdx) / k_e

      !! depth-invariant ratio; guard the no-outflow case, where the pond has all day
      !! to equilibrate and the exponential saturates
      if (q_out > 1.e-9) then
        x_ratio = (area_m2 * k_e) / (q_out * rho_cp)
        if (x_ratio > 50.) then
          f_approach = 1.
        else
          f_approach = 1. - Exp(-x_ratio)
        end if
      else
        f_approach = 1.
      end if

      res(jres)%temp = t_in + (t_equil - t_in) * f_approach

      !! prevent sub-freezing water, as ch_temp does for channels
      if (res(jres)%temp < 0.001) res(jres)%temp = 0.001

      !! the released water carries it downstream -- without this the budget would be
      !! invisible to every channel below the dam, which is the whole point
      ht2%temp = res(jres)%temp

      return
      end subroutine res_temp
