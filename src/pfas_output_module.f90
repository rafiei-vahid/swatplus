      module pfas_output_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Per-HRU, per-PFAS land-phase loss accumulators for the surface-water-
!!    only PFAS implementation.  These are the SWAT+ homes for the SWAT2012
!!    outputs pfas_surq / lat_pfas / pfassol / pfas_sed of pfaslch.f / pfasy.f.
!!
!!    Storage shape mirrors the SWAT+ pesticide output container
!!    (output_ls_pesticide_module): a per-HRU derived type holding one value
!!    per simulated PFAS, dimensioned by npfas.
!!
!!    ~ ~ ~ REENTRANCY ~ ~ ~
!!    pfas_lch / pfas_sed write ONLY the current HRU slice (hpfasb_d(ihru)).
!!    Under the OpenMP parallel land phase each thread owns a distinct ihru,
!!    so distinct array elements are written and these module arrays are safe.
!!    No element is read/written across HRUs within a day.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      implicit none

      !!    per-HRU PFAS land-phase losses (kg/ha): daily slot (zeroed each day by
      !!    pfas_lch) + run-cumulative slot (never zeroed; each HRU writes only its own
      !!    element, so it is OpenMP-safe in the parallel land phase).
      type pfas_hru_output
        real, dimension(:), allocatable :: surq    !! kg/ha |PFAS lost in surface runoff (top layer), day
        real, dimension(:), allocatable :: latq    !! kg/ha |PFAS lost in lateral subsurface flow, day
        real, dimension(:), allocatable :: perc    !! kg/ha |PFAS leached out the bottom of the profile, day
        real, dimension(:), allocatable :: sed     !! kg/ha |PFAS lost sorbed to eroded sediment, day
        real, dimension(:), allocatable :: surq_a  !! kg/ha |run-cumulative surface-runoff loss
        real, dimension(:), allocatable :: latq_a  !! kg/ha |run-cumulative lateral-flow loss
        real, dimension(:), allocatable :: perc_a  !! kg/ha |run-cumulative leaching loss
        real, dimension(:), allocatable :: sed_a   !! kg/ha |run-cumulative sediment-bound loss
      end type pfas_hru_output
      type (pfas_hru_output), dimension(:), allocatable, save :: hpfasb_d  !! by HRU

      !!    per-HRU initial soil PFAS mass (kg/ha, summed over layers), captured at
      !!    pfas_read -> the reference for the end-of-run mass-balance check.
      real, dimension(:,:), allocatable, save :: pfas_init_hru   !! (hru, pfas) kg/ha

      end module pfas_output_module
