!KROME_DRIVER

subroutine krome_driver(d, e, ge, u, v, w, &
      De, HM, HI, HeI, &
      H2I, HII, HeII, H2II, &
      HeIII,  in, jn, kn, imethod, &
      idual, idim, &
      is, js, ks, ie, je, ke, &
      dt, aye, &
      utem, uxyz, uaye, urho, utim, &
      gamma, fh, dtoh)

  !     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
  !
  !     2014, KROME DEVELOPERS to interface the package with ENZO
  !
  !     PURPOSE:
  !     Solve the multi-species rate and cool equations via KROME.
  !
  !     INPUTS:
  !     in,jn,kn - dimensions of 3D fields
  !
  !     d        - total density field
  !
  !     is,ie    - start and end indices of active region (zero based)
  !     idual    - dual energy formalism flag (0 = off, 1 = on)
  !     idim     - dimensionality (rank) of problem
  !     imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
  !
  !     fh       - Hydrogen mass fraction (typically 0.76)
  !     dtoh     - Deuterium to H mass ratio
  !     dt       - timestep to integrate over
  !     aye      - expansion factor (in code units)
  !
  !     utim     - time units (i.e. code units to CGS conversion factor)
  !     uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
  !     urho     - density units
  !     uxyz     - length units
  !     utem     - temperature(-like) units
  !
  !     OUTPUTS:
  !     update chemical abundances densities (HI, HII, etc) and energy
  !
  !     PARAMETERS:
  !     mh      - H mass in cgs units
  !
  !-----------------------------------------------------------------------
  !     USE KROME
  use krome_main
  use krome_user
  use krome_constants

  implicit none

  real*8,parameter::mh=p_mass !mass_h
  real*8::dom,factor,tgas,tgasold,krome_x(krome_nmols)
  real*8::dt,dt_hydro,idom,edot,krome_tiny
  real*8::d(in,jn,kn),e(in,jn,kn),ge(in,jn,kn)
  real*8::u(in,jn,kn),v(in,jn,kn),w(in,jn,kn)
  real*8::aye,utem,uxyz,uaye,urho,utim,gamma,fh,dtoh
  integer::in,jn,kn,imethod,idual,is,js,ks,ie,je,ke,idim
  integer::i,j,k

  real*8::De(in,jn,kn)
  real*8::HM(in,jn,kn)
  real*8::HI(in,jn,kn)
  real*8::HeI(in,jn,kn)
  real*8::H2I(in,jn,kn)
  real*8::HII(in,jn,kn)
  real*8::HeII(in,jn,kn)
  real*8::H2II(in,jn,kn)
  real*8::HeIII(in,jn,kn)

  !******************************

  !set units
  dom = urho*(aye**3)/mh

  !scaling factor for comoving->proper
  factor = aye**(-3)

  !check minimal value and comoving->proper
  krome_tiny = 1d-30
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        HM(i,j,k) = HM(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        HeIII(i,j,k) = HeIII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        HM(i,j,k) = max(HM(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        HeIII(i,j,k) = max(HeIII(i,j,k), krome_tiny)

      end do
    end do
  end do

  !loop over zones
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1

        !rhogas = #KROME_sum to be removed

        !convert to number densities
        krome_x(krome_idx_E) = De(i,j,k) * dom
        krome_x(krome_idx_Hk) = HM(i,j,k) * dom
        krome_x(krome_idx_H) = HI(i,j,k) * dom
        krome_x(krome_idx_HE) = HeI(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_H2) = H2I(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_Hj) = HII(i,j,k) * dom
        krome_x(krome_idx_HEj) = HeII(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_H2j) = H2II(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_HEjj) = HeIII(i,j,k) * dom * 0.25d0

        call evaluate_tgas(d(i,j,k), e(i,j,k), ge(i,j,k),&
            u(i,j,k), v(i,j,k), w(i,j,k),&
            krome_x(:),imethod,idual,idim,tgas,&
            utem)

        !store old tgas
        tgasold = tgas

        dt_hydro = utim*dt !dt*time_conversion

        !call KROME solver
        call krome(krome_x(:),tgas,dt_hydro)

        idom = 1.d0/dom
        !convert back to code units
        De(i,j,k) = krome_x(krome_idx_E) * idom
        HM(i,j,k) = krome_x(krome_idx_Hk) * idom
        HI(i,j,k) = krome_x(krome_idx_H) * idom
        HeI(i,j,k) = krome_x(krome_idx_HE) * idom * 4d0
        H2I(i,j,k) = krome_x(krome_idx_H2) * idom * 2d0
        HII(i,j,k) = krome_x(krome_idx_Hj) * idom
        HeII(i,j,k) = krome_x(krome_idx_HEj) * idom * 4d0
        H2II(i,j,k) = krome_x(krome_idx_H2j) * idom * 2d0
        HeIII(i,j,k) = krome_x(krome_idx_HEjj) * idom * 4d0

        !evaluate energy from temperature difference
        edot = (tgas - tgasold) * d(i,j,k) &
            / ((gamma - 1.d0) * utem * dt)

        !update internal energy
        e(i,j,k)  = e(i,j,k) + edot / d(i,j,k) * dt
        !when using dual
        if (idual .eq. 1) ge(i,j,k) = ge(i,j,k)+ edot / d(i,j,k) * dt
      end do
    end do
  end do

  !scale comoving<-proper
  factor = aye**3
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        HM(i,j,k) = HM(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        HeIII(i,j,k) = HeIII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        HM(i,j,k) = max(HM(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        HeIII(i,j,k) = max(HeIII(i,j,k), krome_tiny)

      end do
    end do
  end do

end subroutine krome_driver
