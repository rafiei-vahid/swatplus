      subroutine pfas_sed

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Land-phase sediment-bound PFAS erosion loading for the current HRU
!!    (ihru) in modern (free-form) SWAT+.  Faithful port of the SWAT2012
!!    pfasy.f (Vahid Rafiei) onto the SWAT+ soil/erosion containers, mirroring
!!    the structure of pest_pesty.f90 / pest_enrsb.f90.
!!
!!    The sediment-sorbed PFAS concentration uses the Freundlich solid-phase
!!    isotherm evaluated at the top-layer aqueous equilibrium conc cw solved
!!    in pfas_lch:
!!
!!      conc [kg PFAS / kg soil] = kf * mw * 1.e-9 * cw**nf
!!
!!    where kf is the Freundlich coefficient ((nmol/kg)/(nM)^n), mw the
!!    molecular weight (kg/mol), the 1.e-9 converts nmol -> mol (kg via mw),
!!    and cw the top-layer aqueous conc (nM) left in pfas_soil_hru(j)%ly(1)%cw.
!!
!!    The sediment-bound load is then
!!
!!      pfas_sed = 1000 * sedyld * conc * er / area_ha   [kg PFAS / ha]
!!
!!    enriched by the per-HRU PFAS enrichment ratio (enr) when provided, else
!!    the CREAMS day enrichment ratio (enratio).  The loaded mass is removed
!!    from the top-layer pool (mass-balance bounded to the available mass).
!!
!!    Scope is SURFACE WATER ONLY (HRU sediment yield; no subbasin routing).
!!
!!    ~ ~ ~ REENTRANCY ~ ~ ~
!!    j = ihru is fixed on entry; every read/write touches only HRU j (its
!!    top-layer pool pfas_soil_hru(j)%ly(1) and output slice hpfasb_d(j)).
!!    All scalars are local.  Safe inside the OpenMP parallel land phase.
!!
!!    ~ ~ ~ OUTGOING (per HRU, per PFAS, kg/ha) ~ ~ ~
!!    hpfasb_d(j)%sed(k)  - PFAS lost sorbed to eroded sediment
!!    pfas_soil_hru(j)%ly(1)%sol_pfas(k) - decremented top-layer pool
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use hru_module, only : hru, sedyld, ihru, enratio
      use pfas_module
      use pfas_output_module

      implicit none

      integer :: j = 0          !! none      |HRU number
      integer :: k = 0          !! none      |sequential PFAS counter
      integer :: kk = 0         !! none      |PFAS number from pfas.dat (crosswalk)
      real :: xx = 0.           !! kg/ha     |PFAS mass available in top soil layer
      real :: conc = 0.         !! kg PFAS/kg soil |sediment-phase PFAS concentration
      real :: er = 0.           !! none      |enrichment ratio applied

      j = ihru

      if (npfas == 0) return
      if (.not. allocated(pfas_flag)) return
      if (pfas_flag(j) == 0) return

      do k = 1, npfas
        kk = pfas_num(k)
        if (kk <= 0) cycle

        hpfasb_d(j)%sed(k) = 0.

        xx = pfas_soil_hru(j)%ly(1)%sol_pfas(k)
        if (xx < 1.e-20) cycle

        !! Freundlich solid-phase conc at the top-layer aqueous equilibrium
        !!   kf [(nmol/kg)/(nM)^n] * mw [kg/mol] * 1.e-9 [mol/nmol]
        !!                                      * cw [nM] ** nf  ->  kg/kg
        conc = pfas_soil_hru(j)%ly(1)%kf(k) * pfasdb(kk)%mw * 1.e-9
        conc = conc                                                       &
     &       * (pfas_soil_hru(j)%ly(1)%cw(k) ** pfas_soil_hru(j)%ly(1)%nf(k))

        !! enrichment ratio: per-HRU PFAS value if set, else day's CREAMS ratio
        if (pfas_soil_hru(j)%ly(1)%enr(k) > 0.) then
          er = pfas_soil_hru(j)%ly(1)%enr(k)
        else
          er = enratio
        end if

        !! sediment-bound loading (kg PFAS/ha)
        !!   1000 * sedyld [t] * conc [kg/kg] * er / area_ha [ha]
        hpfasb_d(j)%sed(k) = 1000. * sedyld(j) * conc * er                &
     &                       / hru(j)%area_ha

        if (hpfasb_d(j)%sed(k) < 0.) hpfasb_d(j)%sed(k) = 0.
        if (hpfasb_d(j)%sed(k) > xx) hpfasb_d(j)%sed(k) = xx

        !! decrement the top-layer pool
        pfas_soil_hru(j)%ly(1)%sol_pfas(k) = xx - hpfasb_d(j)%sed(k)
        hpfasb_d(j)%sed_a(k) = hpfasb_d(j)%sed_a(k) + hpfasb_d(j)%sed(k)

      end do

      return
      end subroutine pfas_sed
