      module pfas_cha_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    In-stream + reservoir PFAS reach state and daily/monthly/yearly
!!    output for the SERIAL channel phase.  Mirrors ch_pesticide_module
!!    exactly (same +, /, // operator idiom) so cha output code can walk
!!    chpfas_*(j)%pfas(:) with the identical pattern used for chpst_*.
!!
!!    Holds the in-stream PFAS database fields that pfas_module%pfasdb does
!!    NOT carry (those are soil/equilibrium params): koc, aq_settle,
!!    aq_resus, ben_bury, ben_act_dep.  These extend the per-PFAS database
!!    with the linear-Koc routing parameters of rtpfas.f / lakeqpfas.f.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      implicit none

      real :: frsol = 0.        !none |fraction of reach PFAS that is soluble
      real :: frsrb = 0.        !none |fraction of reach PFAS that is sorbed

      !! in-stream routing run-cumulative diagnostics (kg) + activity counter
      real, save    :: pfdiag_in   = 0.   !! sum of per-reach-day inflow
      real, save    :: pfdiag_out  = 0.   !! sum of per-reach-day outflow (sol+sor)
      real, save    :: pfdiag_bury = 0.   !! sum of per-reach-day burial (permanent sink)
      integer, save :: pfdiag_active = 0  !! count of reach-days with nonzero PFAS inflow

      !! per-PFAS in-stream routing parameters (extends pfasdb)
      type pfas_cha_db
        character(len=16) :: name = ""    !          |PFAS compound name
        real :: koc       = 0.            !m^3/g     |linear water-sediment partition (Koc)
        real :: aq_settle = 0.            !m/day     |settling velocity of sorbed PFAS
        real :: aq_resus  = 0.            !m/day     |resuspension velocity of bed PFAS
        real :: ben_bury  = 0.            !m/day     |burial velocity in bed sediment
        real :: ben_act_dep = 0.          !m         |active bed-sediment layer depth
      end type pfas_cha_db
      type (pfas_cha_db), dimension(:), allocatable, save :: pfas_chadb

      !! daily reach PFAS balance terms (kg)
      type pfas_cha_processes
        real :: tot_in  = 0.    !kg |total PFAS into reach
        real :: sol_out = 0.    !kg |soluble PFAS out of reach
        real :: sor_out = 0.    !kg |sorbed PFAS out of reach
        real :: settle  = 0.    !kg |PFAS settling to bed sediment
        real :: resus   = 0.    !kg |PFAS resuspended into reach water
        real :: difus   = 0.    !kg |PFAS diffusing between sediment and water
        real :: bury    = 0.    !kg |PFAS buried in bed sediment
        real :: water   = 0.    !kg |PFAS in reach water at end of day
        real :: benthic = 0.    !kg |PFAS in bed sediment at end of day
        real :: conc    = 0.    !ng/L |total PFAS conc of reach outflow (calibration target)
      end type pfas_cha_processes

      type pfas_cha_output
        type (pfas_cha_processes), dimension(:), allocatable :: pfas
      end type pfas_cha_output
      type (pfas_cha_processes) :: ch_pfasbz

      type (pfas_cha_output), dimension(:), allocatable, save :: chpfas_d
      type (pfas_cha_output), dimension(:), allocatable, save :: chpfas_m
      type (pfas_cha_output), dimension(:), allocatable, save :: chpfas_y
      type (pfas_cha_output), dimension(:), allocatable, save :: chpfas_a
      type (pfas_cha_output) :: chpfas, chpfasz

      type pfas_cha_header
          character(len=6)  :: day =      "  jday"
          character(len=6)  :: mo =       "   mon"
          character(len=6)  :: day_mo =   "   day"
          character(len=6)  :: yrc =      "    yr"
          character(len=8)  :: isd =      "   unit "
          character(len=8)  :: id =       " gis_id "
          character(len=16) :: name =     " name           "
          character(len=16) :: pfas =     " pfas"
          character(len=13) :: tot_in =   "tot_in_kg "
          character(len=13) :: sol_out =  "sol_out_kg "
          character(len=14) :: sor_out =  "sor_out_kg "
          character(len=12) :: settle =   "settle_kg "
          character(len=13) :: resus =    "resuspend_kg "
          character(len=12) :: difus =    "diffuse_kg "
          character(len=14) :: bury =     "bury_benth_kg "
          character(len=14) :: water =    "water_stor_kg "
          character(len=12) :: benthic =  "benthic_kg "
      end type pfas_cha_header
      type (pfas_cha_header) :: chpfas_hdr

      interface operator (+)
        module procedure chpfas_add
      end interface
      interface operator (/)
        module procedure chpfas_div
      end interface
      interface operator (//)
        module procedure chpfas_ave
      end interface

      contains

      function chpfas_add(c1, c2) result (c3)
        type (pfas_cha_processes), intent (in) :: c1, c2
        type (pfas_cha_processes) :: c3
        c3%tot_in  = c1%tot_in  + c2%tot_in
        c3%sol_out = c1%sol_out + c2%sol_out
        c3%sor_out = c1%sor_out + c2%sor_out
        c3%settle  = c1%settle  + c2%settle
        c3%resus   = c1%resus   + c2%resus
        c3%difus   = c1%difus   + c2%difus
        c3%bury    = c1%bury    + c2%bury
        c3%water   = c1%water   + c2%water
        c3%benthic = c1%benthic + c2%benthic
      end function chpfas_add

      function chpfas_div(c1, const) result (c2)
        type (pfas_cha_processes), intent (in) :: c1
        real, intent (in) :: const
        type (pfas_cha_processes) :: c2
        c2%tot_in  = c1%tot_in  / const
        c2%sol_out = c1%sol_out / const
        c2%sor_out = c1%sor_out / const
        c2%settle  = c1%settle  / const
        c2%resus   = c1%resus   / const
        c2%difus   = c1%difus   / const
        c2%bury    = c1%bury    / const
        c2%water   = c1%water   / const
        c2%benthic = c1%benthic / const
      end function chpfas_div

      function chpfas_ave(c1, const) result (c2)
        type (pfas_cha_processes), intent (in) :: c1
        real, intent (in) :: const
        type (pfas_cha_processes) :: c2
        c2%tot_in  = const * c1%tot_in
        c2%sol_out = const * c1%sol_out
        c2%sor_out = const * c1%sor_out
        c2%settle  = const * c1%settle
        c2%resus   = const * c1%resus
        c2%difus   = const * c1%difus
        c2%bury    = const * c1%bury
        c2%water   = const * c1%water
        c2%benthic = const * c1%benthic
      end function chpfas_ave

      end module pfas_cha_module
