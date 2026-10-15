      subroutine salt_aqu_hold (iaq, m, gw_volume, h)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    E5 (aquifer_salt hold_initial): sets aquifer iaq's mass of salt ion m to C0 x V / 1000 (kg; C0 in g/m3, V in m3) and returns
!!    the mass that setting added (h > 0) or removed (h < 0). Aquifers are routed serially; no initialised locals.

      use salt_engine_module, only : aqu_c0
      use constituent_mass_module, only : cs_aqu

      implicit none

      integer, intent (in) :: iaq        !none     |aquifer number
      integer, intent (in) :: m          !none     |salt ion
      real, intent (in) :: gw_volume     !m3       |groundwater volume today
      real, intent (out) :: h            !kg       |mass set (+ added, - removed)

      h = aqu_c0(m, iaq) * gw_volume / 1000. - cs_aqu(iaq)%salt(m)
      cs_aqu(iaq)%salt(m) = cs_aqu(iaq)%salt(m) + h

      return
      end subroutine salt_aqu_hold
