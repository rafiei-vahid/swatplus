      subroutine pfas_lch

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Land-phase PFAS leaching / runoff / lateral-flow routine for the
!!    current HRU (ihru) in modern (free-form) SWAT+.  Faithful port of the
!!    SWAT2012 pfaslch.f (Vahid Rafiei) onto the SWAT+ soil/water-flux and
!!    constituent containers.
!!
!!    For each soil layer and each simulated PFAS this routine:
!!      * builds the air-water interfacial area A_aw via pfas_awi,
!!      * solves the three-phase soil equilibrium aqueous concentration cw
!!        via pfas_partition (validated double-precision root-find),
!!      * caps cw at the compound's aqueous solubility,
!!      * converts cw to a water-phase mass concentration co (kg/mm-ha),
!!      * removes dissolved PFAS into the surface-runoff (top layer, scaled
!!        by percop), lateral-flow, and percolation/leach legs,
!!      * decrements the layer pool and cascades leached mass to the layer
!!        below (or emits it as profile leaching at the bottom layer).
!!
!!    Scope is SURFACE WATER ONLY: no tile / aquifer PFAS leg is computed.
!!
!!    ~ ~ ~ REENTRANCY ~ ~ ~
!!    j = ihru is fixed on entry; every read and every write touches only
!!    HRU j (its soil pools pfas_soil_hru(j) and its output slice
!!    hpfasb_d(j)).  All scalars are local.  pfas_partition / pfas_awi are
!!    pure-style.  Safe to call from inside the OpenMP parallel land phase
!!    (one HRU per thread).
!!
!!    ~ ~ ~ OUTGOING (per HRU, per PFAS, kg/ha) ~ ~ ~
!!    hpfasb_d(j)%surq(k)  - PFAS lost in surface runoff
!!    hpfasb_d(j)%latq(k)  - PFAS lost in lateral subsurface flow
!!    hpfasb_d(j)%perc(k)  - PFAS leached from the soil profile
!!    pfas_soil_hru(j)%ly(:)%sol_pfas(:) - updated layer pools
!!    pfas_soil_hru(j)%ly(:)%cw(:)       - last-solved aqueous conc (nM)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      use hru_module, only : surfq, ihru
      use soil_module, only : soil
      use pfas_module
      use pfas_output_module

      implicit none

      integer :: j = 0          !! none      |HRU number
      integer :: ly = 0         !! none      |counter (soil layers)
      integer :: k = 0          !! none      |sequential PFAS counter
      integer :: kk = 0         !! none      |PFAS number from pfas.dat (crosswalk)
      integer :: nly = 0        !! none      |number of soil layers in this HRU
      real :: a_aw = 0.         !! 1/mm      |air-water interfacial area
      real :: cw = 0.           !! nM        |aqueous equilibrium conc (per contaminated site)
      real :: cw_hru = 0.       !! nM        |HRU-scaled aqueous conc (= num_pconta * cw), solubility-capped
      real :: cw_cap = 0.       !! nM        |solubility cap on the aqueous conc
      real :: co = 0.           !! kg/mm-ha  |PFAS water-phase mass concentration
      real :: csurf = 0.        !! kg/mm-ha  |PFAS conc applied to surface runoff / lateral flow
      real :: totmass = 0.      !! kg/ha     |sol_pfas / num_pconta (solver mass term f)
      real :: thick = 0.        !! mm        |thickness of soil layer (solver g)
      real :: vf = 0.           !! mm H2O    |effective flow through the layer (gate)
      real :: qsurf = 0.        !! mm H2O    |surface runoff seen by the layer (top layer only)
      real :: xx = 0.           !! kg/ha     |PFAS mass removed from the layer
      real :: yy = 0.           !! kg/ha     |PFAS mass removed from the layer

      j = ihru

      !! nothing to do if this HRU carries no PFAS
      if (npfas == 0) return
      if (.not. allocated(pfas_flag)) return
      if (pfas_flag(j) == 0) return

      nly = soil(j)%nly

      !! zero the day's leaching / runoff / lateral output for this HRU
      do k = 1, npfas
        hpfasb_d(j)%surq(k) = 0.
        hpfasb_d(j)%latq(k) = 0.
        hpfasb_d(j)%perc(k) = 0.
      end do

      do ly = 1, nly

        thick = soil(j)%phys(ly)%thick

        do k = 1, npfas
          kk = pfas_num(k)
          if (kk <= 0) cycle

          !! reset last-solved aqueous conc for this cell
          pfas_soil_hru(j)%ly(ly)%cw(k) = 0.

          !! surface runoff reaches only the top layer
          if (ly == 1) then
            qsurf = surfq(j)
          else
            qsurf = 0.
          end if

          !! gate: drive partitioning only when there is mass AND outgoing
          !! water AND a defined Freundlich isotherm (matches pfaslch.f)
          vf = max(qsurf, soil(j)%ly(ly)%prk, soil(j)%ly(ly)%flat)

          if (pfas_soil_hru(j)%ly(ly)%sol_pfas(k) >= 1.e-9                 &
     &        .and. vf > 0.                                               &
     &        .and. pfas_soil_hru(j)%ly(ly)%kf(k) > 0.                     &
     &        .and. pfas_soil_hru(j)%ly(ly)%nf(k) > 0.) then

            !! air-water interfacial area for this layer
            a_aw = pfas_awi(soil(j)%phys(ly)%por,                          &
     &                      soil(j)%phys(ly)%st,                           &
     &                      soil(j)%phys(ly)%ul,                           &
     &                      pfas_soil_hru(j)%ly(ly)%sol_d50)

            !! mass term: total layer mass spread over contaminated sites
            totmass = pfas_soil_hru(j)%ly(ly)%sol_pfas(k)                  &
     &              / real(pfas_soil_hru(j)%num_pconta)

            !! solve the three-phase soil equilibrium (per contaminated site)
            cw = pfas_partition(pfasdb(kk)%mw,                             &
     &                          pfas_soil_hru(j)%ly(ly)%kf(k),             &
     &                          pfas_soil_hru(j)%ly(ly)%nf(k),             &
     &                          pfasdb(kk)%lm,                             &
     &                          pfasdb(kk)%kl,                             &
     &                          a_aw,                                      &
     &                          soil(j)%phys(ly)%bd,                       &
     &                          totmass,                                   &
     &                          thick)

            !! rescale to whole-HRU aqueous conc (pfaslch.f: num_pconta * x)
            cw_hru = real(pfas_soil_hru(j)%num_pconta) * cw

            !! enforce maximum aqueous solubility
            !!   cw [nM] * mw [kg/mol] / 1.e3 = mg/L ; cap at pfasdb%sol [mg/L]
            if (pfasdb(kk)%mw > 0.) then
              if (cw_hru * pfasdb(kk)%mw / 1.e3 > pfasdb(kk)%sol) then
                cw_cap = pfasdb(kk)%sol * 1.e3 / pfasdb(kk)%mw
                cw_hru = cw_cap
              end if
            end if

            pfas_soil_hru(j)%ly(ly)%cw(k) = cw_hru

            !! water-phase mass concentration (kg/mm-ha)
            co = cw_hru * pfasdb(kk)%mw / 1.e5

            !! conc carried by surface runoff / lateral flow: top layer is
            !! scaled by the percolation/runoff partition coefficient
            if (ly == 1) then
              csurf = pfasdb(kk)%percop * co
            else
              csurf = co
            end if

            !! ---- PFAS leaching (percolation to the layer below) ----
            xx = co * soil(j)%ly(ly)%prk
            if (xx > pfas_soil_hru(j)%ly(ly)%sol_pfas(k))                  &
     &        xx = pfas_soil_hru(j)%ly(ly)%sol_pfas(k)
            pfas_soil_hru(j)%ly(ly)%sol_pfas(k) =                          &
     &        pfas_soil_hru(j)%ly(ly)%sol_pfas(k) - xx
            if (ly < nly) then
              pfas_soil_hru(j)%ly(ly+1)%sol_pfas(k) =                      &
     &          pfas_soil_hru(j)%ly(ly+1)%sol_pfas(k) + xx
            else
              hpfasb_d(j)%perc(k) = xx
              hpfasb_d(j)%perc_a(k) = hpfasb_d(j)%perc_a(k) + xx
            end if

            !! ---- PFAS lost in surface runoff (top layer only) ----
            if (ly == 1) then
              yy = csurf * surfq(j)
              if (yy > pfas_soil_hru(j)%ly(ly)%sol_pfas(k))               &
     &          yy = pfas_soil_hru(j)%ly(ly)%sol_pfas(k)
              pfas_soil_hru(j)%ly(ly)%sol_pfas(k) =                        &
     &          pfas_soil_hru(j)%ly(ly)%sol_pfas(k) - yy
              hpfasb_d(j)%surq(k) = yy
              hpfasb_d(j)%surq_a(k) = hpfasb_d(j)%surq_a(k) + yy
            end if

            !! ---- PFAS lost in lateral subsurface flow ----
            yy = csurf * soil(j)%ly(ly)%flat
            if (yy > pfas_soil_hru(j)%ly(ly)%sol_pfas(k))                 &
     &        yy = pfas_soil_hru(j)%ly(ly)%sol_pfas(k)
            pfas_soil_hru(j)%ly(ly)%sol_pfas(k) =                          &
     &        pfas_soil_hru(j)%ly(ly)%sol_pfas(k) - yy
            hpfasb_d(j)%latq(k) = hpfasb_d(j)%latq(k) + yy
            hpfasb_d(j)%latq_a(k) = hpfasb_d(j)%latq_a(k) + yy

          end if   ! gate: mass + flow + isotherm

        end do     ! PFAS
      end do        ! soil layers

      return
      end subroutine pfas_lch
