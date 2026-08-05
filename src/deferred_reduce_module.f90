      module deferred_reduce_module

      !! DETERMINISTIC REDUCTION FOR THE BASIN AND REGION CROP-YIELD ACCUMULATORS.
      !!
      !! WHAT WAS WRONG. bsn_crop_yld (basin) and plcal(reg)%lum(lum) (soft-calibration
      !! region) are sums that EVERY HRU contributes to, updated from inside the HRU-parallel
      !! region. Two distinct defects sat on top of that:
      !!
      !!   1. DATA RACE. mgt_sched.f90's two bsn_crop_yld sites and all four plcal sites were
      !!      unguarded read-modify-write on a shared element. Concurrent HRUs lose updates.
      !!
      !!   2. ORDER DEPENDENCE. actions.f90's two bsn_crop_yld sites were made `!$omp atomic`
      !!      (0444d99). Atomic removes the data race but NOT associativity: floating-point
      !!      addition is not associative, so the sum still depends on the order HRUs arrive
      !!      in, which varies with thread count. This is invisible to ThreadSanitizer because
      !!      it is not a race. Measured on the 2026-08-05 acceptance matrix (5 basins x
      !!      {1,2,4,8,16} threads x 3 trials): basin_crop_yld_yr.txt held 12 distinct
      !!      content-hashes and basin_crop_yld_aa.txt 7, while every other output but
      !!      mgt_out.txt was byte-identical.
      !!
      !! THE FIX, which closes both. HRUs never touch the shared accumulators. Each
      !! contribution is appended to the CONTRIBUTING HRU's own buffer; once the parallel
      !! region has closed, dfr_flush replays the buffers in ascending HRU index order and
      !! applies them to the real accumulators. Summation order becomes a property of the
      !! model rather than of the schedule, so the result is identical at every thread count.
      !!
      !! THE DEFERRED PATH IS ALWAYS TAKEN, INCLUDING SERIAL. Keeping a direct path for N=1
      !! would reintroduce exactly the defect being fixed: the acceptance criterion (Vahid) is
      !! ONE equivalence class across all thread counts, and serial is one of them.
      !!
      !! No synchronisation appears below and none is needed: an HRU is processed by exactly
      !! one thread, so only that thread ever appends to that HRU's buffer. That is the whole
      !! point of keying the buffer on the HRU rather than on the thread -- a per-thread
      !! buffer would be just as race-free and just as order-dependent, because the partition
      !! into threads changes with the thread count.

      implicit none

      integer, parameter :: dfr_kind_bsn = 1     !! bsn_crop_yld(i1)  += (area_ha, yield)
      integer, parameter :: dfr_kind_plcal = 2   !! plcal(i1)%lum(i2) += (ha, yield)

      integer, parameter :: dfr_cap0 = 8         !! initial per-HRU capacity; grows by doubling

      type dfr_event
        integer :: kind = 0
        integer :: i1 = 0
        integer :: i2 = 0
        real :: v1 = 0.
        real :: v2 = 0.
      end type dfr_event

      type dfr_hru_buf
        integer :: n = 0
        type (dfr_event), dimension(:), allocatable :: ev
      end type dfr_hru_buf

      type (dfr_hru_buf), dimension(:), allocatable :: dfr_buf
      logical :: dfr_ready = .false.

      contains

      !! Allocate one buffer per HRU. Called once, before the daily loop.
      subroutine dfr_init (nhru)
        integer, intent (in) :: nhru

        if (allocated(dfr_buf)) deallocate (dfr_buf)
        if (nhru < 1) then
          dfr_ready = .false.
          return
        end if
        allocate (dfr_buf(nhru))
        dfr_ready = .true.

      end subroutine dfr_init

      !! Append one contribution to HRU j's own buffer. Safe to call from inside the parallel
      !! region: element j is touched only by the thread that owns HRU j, and growing
      !! dfr_buf(j)%ev reallocates that element's component alone.
      subroutine dfr_add (j, kind, i1, i2, v1, v2)
        integer, intent (in) :: j, kind, i1, i2
        real, intent (in) :: v1, v2
        type (dfr_event), dimension(:), allocatable :: tmp
        integer :: cap

        if (.not. dfr_ready) return
        if (j < 1 .or. j > size(dfr_buf)) return

        if (.not. allocated(dfr_buf(j)%ev)) allocate (dfr_buf(j)%ev(dfr_cap0))
        cap = size(dfr_buf(j)%ev)
        if (dfr_buf(j)%n >= cap) then
          allocate (tmp(2 * cap))
          tmp(1:cap) = dfr_buf(j)%ev(1:cap)
          call move_alloc (tmp, dfr_buf(j)%ev)
        end if

        dfr_buf(j)%n = dfr_buf(j)%n + 1
        dfr_buf(j)%ev(dfr_buf(j)%n) = dfr_event (kind, i1, i2, v1, v2)

      end subroutine dfr_add

      !! Replay every buffered contribution in ascending HRU index order, then reset. Called
      !! once per day from time_control, AFTER the parallel region in command has closed.
      subroutine dfr_flush
        use plant_module, only : bsn_crop_yld
        use calibration_data_module, only : plcal

        integer :: j, k, i1, i2

        if (.not. dfr_ready) return

        do j = 1, size(dfr_buf)
          if (dfr_buf(j)%n == 0) cycle
          do k = 1, dfr_buf(j)%n
            i1 = dfr_buf(j)%ev(k)%i1
            i2 = dfr_buf(j)%ev(k)%i2
            select case (dfr_buf(j)%ev(k)%kind)
            case (dfr_kind_bsn)
              if (allocated(bsn_crop_yld)) then
                if (i1 >= 1 .and. i1 <= size(bsn_crop_yld)) then
                  bsn_crop_yld(i1)%area_ha = bsn_crop_yld(i1)%area_ha + dfr_buf(j)%ev(k)%v1
                  bsn_crop_yld(i1)%yield = bsn_crop_yld(i1)%yield + dfr_buf(j)%ev(k)%v2
                end if
              end if
            case (dfr_kind_plcal)
              if (allocated(plcal)) then
                if (i1 >= 1 .and. i1 <= size(plcal)) then
                  if (allocated(plcal(i1)%lum)) then
                    if (i2 >= 1 .and. i2 <= size(plcal(i1)%lum)) then
                      plcal(i1)%lum(i2)%ha = plcal(i1)%lum(i2)%ha + dfr_buf(j)%ev(k)%v1
                      plcal(i1)%lum(i2)%sim%yield = plcal(i1)%lum(i2)%sim%yield + dfr_buf(j)%ev(k)%v2
                    end if
                  end if
                end if
              end if
            end select
          end do
          dfr_buf(j)%n = 0        !! keep the allocation; only the fill level resets
        end do

      end subroutine dfr_flush

      end module deferred_reduce_module
