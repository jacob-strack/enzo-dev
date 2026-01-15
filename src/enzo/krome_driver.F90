!KROME_DRIVER

subroutine krome_driver(d, e, ge, u, v, w, &
      De, CI, CII,  in, jn, kn, imethod, &
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
  real*8::CI(in,jn,kn)
  real*8::CII(in,jn,kn)

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
        CI(i,j,k) = CI(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)

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
        krome_x(krome_idx_C) = CI(i,j,k) * dom * 0.08333333333333333d0
        krome_x(krome_idx_Cj) = CII(i,j,k) * dom * 0.08333333333333333d0

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
        CI(i,j,k) = krome_x(krome_idx_C) * idom * 12d0
        CII(i,j,k) = krome_x(krome_idx_Cj) * idom * 12d0

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
        CI(i,j,k) = CI(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)

      end do
    end do
  end do

end subroutine krome_driver
