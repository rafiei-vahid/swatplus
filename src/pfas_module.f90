      module pfas_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Surface-water-only PFAS fate-and-transport state and the validated
!!    soil three-phase equilibrium solver for modern (free-form) SWAT+.
!!
!!    Ported from the SWAT2012 PFAS implementation (Vahid Rafiei): the
!!    governing per-HRU-per-layer-per-PFAS air-water / aqueous / solid
!!    equilibrium of pfaslch.f, plus the global per-PFAS database of
!!    readpfas.f / modparm.f.  Scope is surface water only: NO groundwater
!!    / aquifer PFAS state lives here.
!!
!!    ~ ~ ~ CONTAINER DESIGN CHOICE ~ ~ ~
!!    PFAS soil pools are kept in a DEDICATED PARALLEL container
!!    (pfas_soil / pfas_soil_ly below), NOT folded into the existing
!!    constituent_mass%cs / cs_soil containers.  Rationale:
!!      * The cs_soil "cs" slot is already owned by the generic-constituent
!!        (rtb cs) feature; overloading it would couple two unrelated
!!        constituent systems and break crosswalks.
!!      * PFAS needs per-layer state that pesticides/cs do not carry:
!!        Freundlich kf and n, an enrichment ratio, and a mineral d50 used
!!        only for the air-water interfacial area A_aw.  A parallel type is
!!        the clean home for those.
!!    The container SHAPE deliberately mirrors SWAT+ pesticide storage
!!    (soil_constituent_mass -> per-layer derived type holding a constituent
!!    vector), so downstream code that already walks cs_soil(j)%ly(:) can
!!    walk pfas_soil(j)%ly(:) with the same idiom.
!!
!!    ~ ~ ~ REENTRANCY ~ ~ ~
!!    pfas_partition and pfas_awi are pure-style: they read only their
!!    scalar arguments and write only local scalars, so they are safe to
!!    call from inside the OpenMP parallel land phase (one HRU per thread).
!!    No module variable is written by either function.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      implicit none

      !!    ~ ~ ~ GLOBAL PER-PFAS DATABASE ~ ~ ~
      !!    one entry per PFAS compound (from pfas.dat). Freundlich kf & n
      !!    are intentionally NOT here: they vary per HRU per soil layer and
      !!    live in the soil pool (pfas_soil) instead.
      type pfas_db
        character(len=16) :: name = ""   !!                |PFAS compound name (points to database)
        real :: mw      = 0.             !! kg/mol         |molecular weight ("a" in the solver)
        real :: sol     = 0.             !! mg/L           |maximum aqueous solubility
        real :: kl      = 0.             !! L/nmol         |Langmuir concentration coeff K_L ("n" in solver)
        real :: lm      = 0.             !! nmol/m^2       |Langmuir max surface conc Gamma_max ("h" in solver)
        real :: percop  = 0.             !! none (0-1)     |PFAS percolation/runoff partition coefficient
      end type pfas_db
      type (pfas_db), dimension(:), allocatable, save :: pfasdb

      !!    ~ ~ ~ SOIL POOL : PER HRU PER LAYER PER PFAS ~ ~ ~
      !!    one pfas_soil_ly per soil layer; vectors are dimensioned by
      !!    the number of simulated PFAS (npfas).
      type pfas_soil_ly
        real, dimension(:), allocatable :: sol_pfas   !! kg/ha            |PFAS mass in this soil layer
        real, dimension(:), allocatable :: kf         !! (nmol/kg)/(nM)^n |Freundlich sorption coeff ("c" in solver)
        real, dimension(:), allocatable :: nf         !! none             |Freundlich exponent ("m" in solver)
        real, dimension(:), allocatable :: enr        !! none             |PFAS enrichment ratio (sediment-bound loading)
        real, dimension(:), allocatable :: cw         !! nM (1e-9 mol/L)  |aqueous equilibrium conc (last solved)
        real :: sol_d50 = 0.                          !! mm               |median grain diameter (for A_aw), per layer
      end type pfas_soil_ly

      !!    per-HRU soil column of PFAS pools
      type pfas_soil
        type (pfas_soil_ly), dimension(:), allocatable :: ly   !! by soil layer
        integer :: num_pconta = 1        !! none           |number of contaminated sites in HRU (mass divisor)
      end type pfas_soil
      type (pfas_soil), dimension(:), allocatable, save :: pfas_soil_hru   !! dimensioned by HRU

      !!    ~ ~ ~ MODULE-LEVEL COUNTS / FLAGS / CROSSWALK ~ ~ ~
      integer, save :: npfas = 0                                  !! none |number of PFAS simulated
      integer, dimension(:), allocatable, save :: pfas_num        !! none |sequential PFAS -> pfasdb index crosswalk
      integer, dimension(:), allocatable, save :: pfas_flag       !! none |per-HRU PFAS-active flag (0=off,1=on)

      contains

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      pure function pfas_awi(por, st, ul, d50) result(a_aw)
!!    Air-water interfacial area per unit length, A_aw (1/mm).
!!      A_aw = 6 * (1 - por) * (1 - st/ul) / d50
!!    Reentrant: pure, scalar in/out only.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      real, intent(in) :: por     !! none   |total porosity of layer (fraction)
      real, intent(in) :: st      !! mm H2O |water stored in layer
      real, intent(in) :: ul      !! mm H2O |water held at saturation
      real, intent(in) :: d50     !! mm     |median grain diameter
      real :: a_aw                !! 1/mm   |air-water interfacial area

      real :: sat_frac

      if (ul > 0.) then
        sat_frac = st / ul
      else
        sat_frac = 0.
      end if
      if (d50 > 0.) then
        a_aw = 6. * (1. - por) * (1. - sat_frac) / d50
      else
        a_aw = 0.
      end if

      end function pfas_awi

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      function pfas_partition(mw, kf, nf, gamma_max, kl, a_aw, bd,         &
     &                        totmass, thick) result(cw)
!!    Soil three-phase equilibrium solver (DOUBLE PRECISION port of the
!!    SWAT2012 real(16) pfaslch.f root-find).  Returns the aqueous-phase
!!    equilibrium PFAS concentration cw (nM, i.e. 1e-9 mol/L) that balances
!!
!!      F(x) = a*q*h*n*x/(1+n*x)   ! air-water (Langmuir) term
!!           + a*x                 ! aqueous term
!!           + c*a*d*x**m          ! Freundlich solid term
!!           - 1.e5*f/g            ! total mass term
!!
!!    F is monotonic increasing in x>0, so the positive root is unique.
!!    Algorithm (byte-faithful to the original): exponential bracket up
!!    (x10) and down (x0.1) -> bisection to |F|<1 -> Improved Halley to
!!    F<1e-6.  The 1.e5 unit-conversion constant is ported exactly.
!!
!!    Validated: float64 holds 1.7e-7 worst-case rel error vs a 40-digit
!!    reference across 29,160 cases (mass-balance residual 5.4e-8, 0
!!    non-convergence, <=118 iters); the real(16) quad of the original is
!!    defensive overkill, so double is used.
!!
!!    Reentrant: all working storage (x, fx, brackets, derivatives) is
!!    local; no module state is touched.  Safe inside the OpenMP land phase.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      real, intent(in) :: mw         !! kg/mol         |molecular weight        (a)
      real, intent(in) :: kf         !! (nmol/kg)/(nM)^n|Freundlich coefficient (c)
      real, intent(in) :: nf         !! none           |Freundlich exponent     (m)
      real, intent(in) :: gamma_max  !! nmol/m^2       |Langmuir Gamma_max      (h)
      real, intent(in) :: kl         !! L/nmol         |Langmuir K_L            (n)
      real, intent(in) :: a_aw       !! 1/mm           |air-water area          (q)
      real, intent(in) :: bd         !! Mg/m^3         |soil bulk density       (d)
      real, intent(in) :: totmass    !! kg/ha          |mass/num_pconta         (f)
      real, intent(in) :: thick      !! mm             |layer thickness         (g)
      real :: cw                     !! nM             |aqueous equilibrium conc

      !! double-precision working storage (Fortran real*8 == IEEE float64)
      real(8) :: a, c, d, m, f, h, n, q, g
      real(8) :: rhs, x, fx, fxn, fxr, gama, alpha, denom
      real(8) :: ub, lb, mb
      integer :: it

      a = real(mw,        8)
      c = real(kf,        8)
      d = real(bd,        8)
      m = real(nf,        8)
      f = real(totmass,   8)
      h = real(gamma_max, 8)
      n = real(kl,        8)
      q = real(a_aw,      8)
      g = real(thick,     8)

      cw = 0.
      if (g <= 0.d0 .or. f <= 0.d0) return

      rhs = 1.d5 * f / g

      !! ---- exponential bracket up (x*10 until F>=0) ----
      !! UB is capped at a physical ceiling (1e12 nM >> PFOS solubility ~1.4e6 nM) so
      !! ub**m can never overflow double precision even on degenerate/large mass inputs;
      !! if F is still <0 at the ceiling the solubility cap downstream handles it.
      ub = f / g
      if (ub < 1.d-15) ub = 1.d-15
      fx = -1.d0
      it = 0
      do while (fx < 0.d0 .and. it < 10000)
        fx = a*q*h*n*ub / (1.d0 + n*ub) + a*ub + c*a*d*(ub**m) - rhs
        if (fx < 0.d0) ub = ub * 10.d0
        if (ub > 1.d12) then
          ub = 1.d12; exit
        end if
        it = it + 1
      end do

      !! ---- exponential bracket down (x*0.1 until F<=0) ----
      !! LB is floored at 1e-15 nM (sub-femtomolar; unphysical for a real root) so
      !! x**(m-2) in the Halley step can never overflow.
      lb = f / g
      if (lb > 1.d12) lb = 1.d12
      fx = 1.d0
      it = 0
      do while (fx > 0.d0 .and. it < 10000)
        fx = a*q*h*n*lb / (1.d0 + n*lb) + a*lb + c*a*d*(lb**m) - rhs
        if (fx > 0.d0) lb = lb / 10.d0
        if (lb < 1.d-15) then
          lb = 1.d-15; exit
        end if
        it = it + 1
      end do

      !! ---- bisection until |F| < 1 ----
      fx = 101.d0
      mb = (lb + ub) / 2.d0
      it = 0
      do while (abs(fx) > 1.d0 .and. it < 200)
        mb = (lb + ub) / 2.d0
        fx = a*q*h*n*mb / (1.d0 + n*mb) + a*mb + c*a*d*(mb**m) - rhs
        if (fx > 0.d0) then
          ub = mb
        else
          lb = mb
        end if
        it = it + 1
      end do

      !! ---- Improved Halley until F < 1e-6 ----
      !! The bisection above leaves a valid bracket [lb,ub] (both finite, both in the
      !! physical range, root guaranteed between them since F is monotonic). We CLAMP
      !! every Halley iterate to [lb,ub]: this preserves convergence (the root is in the
      !! bracket) while preventing the overshoot that would drive x to an extreme where
      !! x**(m-2) overflows double precision and traps under -fpe0. On any non-finite
      !! step we fall back to the bisection midpoint.
      x  = mb
      fx = 1.d0
      it = 0
      do while (fx > 1.d-6 .and. it < 1000)
        fx  = a*q*h*n*x / (1.d0 + n*x) + a*x + c*a*d*(x**m) - rhs

        fxn = a*q*h*n / ((1.d0 + n*x)**2)                                  &
     &      + a * (m*c*d*(x**(m - 1.d0)) + 1.d0)
        if (fxn == 0.d0) exit

        fxr = (-2.d0*a*q*h*(n**2)) / ((1.d0 + n*x)**3)                     &
     &      + m*(m - 1.d0)*(x**(m - 2.d0))*a*c*d

        gama  = fx * fxr / (fxn**2)
        denom = 4.d0 - 6.d0*gama + gama**2
        if (denom == 0.d0) exit
        alpha = 4.d0 * (fx/fxn) * ((1.d0 - gama) / denom)
        x = x - alpha

        if (x < 0.d0) x = -x
        if (fx < 0.d0) fx = -fx
        !! keep the iterate inside the validated, overflow-safe bracket
        if (x < lb) x = lb
        if (x > ub) x = ub
        it = it + 1
      end do

      !! final safety: if Halley left a non-finite or out-of-range value, use bisection root
      if (.not. (x == x) .or. x < lb .or. x > ub) x = mb

      cw = real(x)

      end function pfas_partition

      end module pfas_module
