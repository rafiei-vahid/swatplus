      subroutine salt_imp_wash (j, fw)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    E2 (road_washoff impervious): washes the fraction fw of HRU j's impervious road-salt store into urbqsalt (kg/ha), which
!!    then takes the surface lag (stor_surfstor) and hru_hyds' routing as surqsalt does; the day's washed mass is kept in
!!    hsaltb_d(j)%salt(:)%impw for output. Writes only HRU j's entries (safe in the HRU-parallel loop); no initialised locals.

      use salt_engine_module, only : rdimp
      use salt_module, only : hsaltb_d
      use hru_module, only : urbqsalt
      use constituent_mass_module, only : cs_db

      implicit none

      integer, intent (in) :: j      !none          |HRU number
      real, intent (in) :: fw        !none          |fraction of the store washed off today
      integer :: isalt               !none          |salt ion counter
      real :: w                      !kg/ha         |salt washed off today

      do isalt = 1, cs_db%num_salts
        w = rdimp(isalt, j) * fw
        rdimp(isalt, j) = rdimp(isalt, j) - w
        hsaltb_d(j)%salt(isalt)%impw = hsaltb_d(j)%salt(isalt)%impw + w
        urbqsalt(j, isalt) = urbqsalt(j, isalt) + w
      end do

      return
      end subroutine salt_imp_wash
