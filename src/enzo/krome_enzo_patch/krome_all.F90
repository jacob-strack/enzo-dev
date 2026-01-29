
!############### MODULE ##############
module krome_commons
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::idx_E=1
  integer,parameter::idx_Hk=2
  integer,parameter::idx_Ck=3
  integer,parameter::idx_Ok=4
  integer,parameter::idx_H=5
  integer,parameter::idx_HE=6
  integer,parameter::idx_H2=7
  integer,parameter::idx_C=8
  integer,parameter::idx_O=9
  integer,parameter::idx_OH=10
  integer,parameter::idx_CO=11
  integer,parameter::idx_CH=12
  integer,parameter::idx_CH2=13
  integer,parameter::idx_C2=14
  integer,parameter::idx_HCO=15
  integer,parameter::idx_H2O=16
  integer,parameter::idx_O2=17
  integer,parameter::idx_CO_total=18
  integer,parameter::idx_H2O_total=19
  integer,parameter::idx_Hj=20
  integer,parameter::idx_HEj=21
  integer,parameter::idx_H2j=22
  integer,parameter::idx_Cj=23
  integer,parameter::idx_Oj=24
  integer,parameter::idx_HOCj=25
  integer,parameter::idx_HCOj=26
  integer,parameter::idx_H3j=27
  integer,parameter::idx_CHj=28
  integer,parameter::idx_CH2j=29
  integer,parameter::idx_COj=30
  integer,parameter::idx_CH3j=31
  integer,parameter::idx_OHj=32
  integer,parameter::idx_H2Oj=33
  integer,parameter::idx_H3Oj=34
  integer,parameter::idx_O2j=35
  integer,parameter::idx_HEjj=36
  integer,parameter::idx_CR=37
  integer,parameter::idx_g=38
  integer,parameter::idx_Tgas=39
  integer,parameter::idx_dummy=40
  integer,parameter::nrea=275
  integer,parameter::nmols=36
  integer,parameter::nspec=40
  integer,parameter::natoms=5
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=0
  integer,parameter::nPhotoRea=0

  !cooling index
  integer,parameter::idx_cool_h2 = 1
  integer,parameter::idx_cool_h2gp = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_cen = 3
  integer,parameter::idx_cool_hd = 4
  integer,parameter::idx_cool_z = 5
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_dh = 6
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_cie = 9
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_co = 14
  integer,parameter::idx_cool_zcie = 15
  integer,parameter::idx_cool_zcienouv = 16
  integer,parameter::idx_cool_zextend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::idx_cool_oh = 19
  integer,parameter::idx_cool_h2o = 20
  integer,parameter::idx_cool_hcn = 21
  integer,parameter::ncools = 21

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_dh = 4
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_photoav = 5
  integer,parameter::idx_heat_av = 5
  integer,parameter::idx_heat_cr = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_zcie = 11
  integer,parameter::nheats = 11

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating), same for chemistry
  integer::krome_thermo_toggle, krome_chemo_toggle
  !$omp threadprivate(krome_thermo_toggle, krome_chemo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_r3(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for frequency bins

  ! Draine dust absorption data loaded from file, via load_kabs
  ! in krome_photo module
  real*8::find_Av_draine_kabs(nPhotoBins)

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_crate,user_Av
  !$omp threadprivate(user_crate,user_Av)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !xsecs from file variables

  ! Gibbs free energy data from file variables

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

  !partition functions common vars

  !verbatim reactions
  character*50::reactionNames(nrea)

  !store evaluate once rate
  real*8::rateEvaluateOnce(nrea)
  !$omp threadprivate(rateEvaluateOnce)

end module krome_commons

!############### MODULE ##############
module krome_constants
  implicit none

  !constants
  real*8,parameter::boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::iboltzmann_eV = 1d0/boltzmann_eV !K / eV
  real*8,parameter::iboltzmann_erg = 1d0/boltzmann_erg !K / erg
  real*8,parameter::planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::planck_J = 6.62606957d-34 !J s
  real*8,parameter::planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::iplanck_eV = 1d0/planck_eV !1 / eV / s
  real*8,parameter::iplanck_J = 1d0/planck_J !1 / J / s
  real*8,parameter::iplanck_erg = 1d0/planck_erg !1 / erg / s
  real*8,parameter::gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::e_mass = 9.10938188d-28 !g
  real*8,parameter::p_mass = 1.67262158d-24 !g
  real*8,parameter::n_mass = 1.674920d-24 !g
  real*8,parameter::ip_mass = 1d0/p_mass !1/g
  real*8,parameter::clight = 2.99792458e10 !cm/s
  real*8,parameter::pi = 3.14159265359d0 !#
  real*8,parameter::eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::km_to_cm = 1d5 !km -> cm
  real*8,parameter::cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass !
  real*8,parameter::pre_kvgas_sqrt = sqrt(8.d0*boltzmann_erg/pi) !
  real*8,parameter::pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter::exp_planck = planck_erg / boltzmann_erg !s*K
  real*8,parameter::stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::N_avogadro = 6.0221d23 !#
  real*8,parameter::Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::hubble = 0.704d0 !dimensionless
  real*8,parameter::Omega0 = 1.0d0 !dimensionless
  real*8,parameter::Omegab = 0.0456d0 !dimensionless
  real*8,parameter::Hubble0 = 1.d2*hubble*km_to_cm*cm_to_Mpc !1/s

end module krome_constants

!############### MODULE ##############
module krome_fit
contains

  !*****************************
  subroutine init_anytab3D(filename,x,y,z,f,xmul,ymul,zmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8::rout(4)
    integer::i,j,k,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(f,1)) then
      print *,"ERROR: in init_anytab3D x size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(f,2)) then
      print *,"ERROR: in init_anytab3D y size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Z input array
    if(size(z).ne.size(f,3)) then
      print *,"ERROR: in init_anytab3D z size differs from f(x,y,z)"
      stop
    end if

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab3D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of grid points"
      print *," per dimension in the format"
      print *,"  XX, YY, ZZ"
      print *,row_string
      stop
    end if

    !loop to read file (3rd dimension of f() is
    ! first in the tables. i.e. tables are z,x,y,
    ! while f() is x,y,z
    do i=1,size(z)
      do j=1,size(x)
        do k=1,size(y)
          read(unit,*,iostat=ios) rout(:)
          y(k) = rout(3)
          f(j,k,i) = rout(4)
        end do
        x(j) = rout(2)
        read(unit,*,iostat=ios) !skip blanks
      end do
      z(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))
    zmul = 1d0/(z(2)-z(1))

  end subroutine init_anytab3D

  !********************************************
  !load 2d tables from filename
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:,:),xmul,ymul
    real*8::rout(3)
    integer::i,j,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
      print *,"ERROR: in init_anytab2D x size differs from z"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
      print *,"ERROR: in init_anytab2D y size differs from z"
      stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of rows and "
      print *," columns in the format"
      print *,"  RR, CC"
      print *,row_string
      stop
    end if

    !loop to read file
    do i=1,size(x)
      do j=1,size(y)
        read(unit,*,iostat=ios) rout(:)
        y(j) = rout(2)
        z(i,j) = rout(3)
      end do
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !********************************************
  !load 1d tables from filename
  subroutine init_anytab1D(filename,x,y,xmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),xmul
    real*8::rout(2)
    integer::i,ios,unit

    !check the size of the X input array
    if(size(x) /= size(y)) then
      print *,"ERROR: in init_anytab1D x size differs from y"
      stop
    end if

    if (krome_mpi_rank <= 1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios /= 0) then
      print *,"ERROR: in init_anytab1D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    ! !check if first line is OK
    ! if(scan(row_string,",")==0) then
    !    print *,"ERROR: file "//filename//" should"
    !    print *," contain the number of rows and "
    !    print *," columns in the format"
    !    print *,"  RR, CC"
    !    print *,row_string
    !    stop
    ! end if

    !loop to read file
    do i=1,size(x)
      read(unit,*,iostat=ios) rout(:)
      y(i) = rout(2)
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios /= 0) exit

    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))

  end subroutine init_anytab1D

  !******************************
  !test 2d fit and save to file
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
      do j=1,size(y)
        xx = x(i)
        yy = y(j)
        zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
        write(unit1,*) xx,yy,zz
        write(unit2,*) x(i),y(j),z(i,j)
      end do
      write(unit1,*)
      write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_anytab2D

  !*****************************
  function fit_anytab3D(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D

  !******************************
  !return 2d fit at xx,yy
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    implicit none
    real*8::fit_anytab2D
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  !return 1d fit at xx
  function fit_anytab1D(x,z,xmul,xx)
    real*8,intent(in)::x(:),z(:),xmul,xx
    real*8::fit_anytab1D,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  function fit_anytab3D_linlinlog(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D_linlinlog,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D_linlog(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D_linlog(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D_linlinlog = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D_linlinlog

  !***************************
  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
      z1 = log10(-z1)
      z2 = log10(-z2)
      fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
      return
    end if

    if(z1>0d0 .and. z2>0d0) then
      z1 = log10(z1)
      z2 = log10(z2)
      fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
      return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
      b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
      ispline = y(1)
      return
    end if

    if(u>=x(n)) then
      ispline = y(n)
      return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
      k = (i+j)/2
      if(u<x(k)) then
        j=k
      else
        i=k
      end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_fit
!This module contains useful routines to get physical
! quantities, like mean molecular weight, mass density,
! mass, jeans length, etc. etc.

!############### MODULE ##############
module krome_getphys
contains

  !*****************************
  !get the mean molecular weight
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = max(sum(n(1:nmols)*m(1:nmols)),1d-40) &
        / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28	!E
    get_mass(2) = 1.6744434563759998d-24	!H-
    get_mass(3) = 2.0077106047315998d-23	!C-
    get_mass(4) = 2.6769171083692d-23	!O-
    get_mass(5) = 1.673532518188d-24	!H
    get_mass(6) = 6.692065036376d-24	!HE
    get_mass(7) = 3.347065036376d-24	!H2
    get_mass(8) = 2.0076195109128d-23	!C
    get_mass(9) = 2.6768260145504d-23	!O
    get_mass(10) = 2.8441792663692003d-23	!OH
    get_mass(11) = 4.6844455254632d-23	!CO
    get_mass(12) = 2.1749727627316d-23	!CH
    get_mass(13) = 2.3423260145503998d-23	!CH2
    get_mass(14) = 4.0152390218256d-23	!C2
    get_mass(15) = 4.8517987772819996d-23	!HCO
    get_mass(16) = 3.011532518188d-23	!H2O
    get_mass(17) = 5.3536520291008d-23	!O2
    get_mass(18) = 0d0	!CO_total
    get_mass(19) = 0d0	!H2O_total
    get_mass(20) = 1.67262158d-24	!H+
    get_mass(21) = 6.691154098188d-24	!HE+
    get_mass(22) = 3.346154098188d-24	!H2+
    get_mass(23) = 2.007528417094d-23	!C+
    get_mass(24) = 2.6767349207316d-23	!O+
    get_mass(25) = 4.8517076834631994d-23	!HOC+
    get_mass(26) = 4.8517076834631994d-23	!HCO+
    get_mass(27) = 5.0196866163760004d-24	!H3+
    get_mass(28) = 2.1748816689128d-23	!CH+
    get_mass(29) = 2.3422349207316d-23	!CH2+
    get_mass(30) = 4.6843544316444d-23	!CO+
    get_mass(31) = 2.5095881725504002d-23	!CH3+
    get_mass(32) = 2.8440881725504d-23	!OH+
    get_mass(33) = 3.0114414243692d-23	!H2O+
    get_mass(34) = 3.178794676188d-23	!H3O+
    get_mass(35) = 5.353560935282d-23	!O2+
    get_mass(36) = 6.690243159999999d-24	!HE++
    get_mass(37) = 0d0	!CR
    get_mass(38) = 0d0	!g
    get_mass(39) = 0d0	!Tgas
    get_mass(40) = 0d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 33132602150543.92	!E
    get_imass_sqrt(2) = 772795806394.0071	!H-
    get_imass_sqrt(3) = 223177004181.41986	!C-
    get_imass_sqrt(4) = 193278051340.87015	!O-
    get_imass_sqrt(5) = 773006102110.9268	!H
    get_imass_sqrt(6) = 386562679981.0883	!HE
    get_imass_sqrt(7) = 546597856701.2171	!H2
    get_imass_sqrt(8) = 223182067345.7445	!C
    get_imass_sqrt(9) = 193281339990.54416	!O
    get_imass_sqrt(10) = 187508740611.18686	!OH
    get_imass_sqrt(11) = 146106959624.0866	!CO
    get_imass_sqrt(12) = 214423849574.04788	!CH
    get_imass_sqrt(13) = 206621889668.37103	!CH2
    get_imass_sqrt(14) = 157813553259.4087	!C2
    get_imass_sqrt(15) = 143565011358.4871	!HCO
    get_imass_sqrt(16) = 182224271009.1322	!H2O
    get_imass_sqrt(17) = 136670546184.13641	!O2
    get_imass_sqrt(18) = 0d0	!CO_total
    get_imass_sqrt(19) = 0d0	!H2O_total
    get_imass_sqrt(20) = 773216569600.4055	!H+
    get_imass_sqrt(21) = 386588992536.287	!HE+
    get_imass_sqrt(22) = 546672253002.7066	!H2+
    get_imass_sqrt(23) = 223187130854.6853	!C+
    get_imass_sqrt(24) = 193284628808.09424	!O+
    get_imass_sqrt(25) = 143566359113.19513	!HOC+
    get_imass_sqrt(26) = 143566359113.19513	!HCO+
    get_imass_sqrt(27) = 446335774600.3221	!H3+
    get_imass_sqrt(28) = 214428340044.27756	!CH+
    get_imass_sqrt(29) = 206625907581.7343	!CH2+
    get_imass_sqrt(30) = 146108380244.195	!CO+
    get_imass_sqrt(31) = 199617572780.53006	!CH3+
    get_imass_sqrt(32) = 187511743462.96832	!OH+
    get_imass_sqrt(33) = 182227027061.27744	!H2O+
    get_imass_sqrt(34) = 177365342346.192	!H3O+
    get_imass_sqrt(35) = 136671708941.89337	!O2+
    get_imass_sqrt(36) = 386615310465.34766	!HE++
    get_imass_sqrt(37) = 0d0	!CR
    get_imass_sqrt(38) = 0d0	!g
    get_imass_sqrt(39) = 0d0	!Tgas
    get_imass_sqrt(40) = 0d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.0977693252662275d+27	!E
    get_imass(2) = 5.9721335838016375d+23	!H-
    get_imass(3) = 4.98079751953935d+22	!C-
    get_imass(4) = 3.7356405130124037d+22	!O-
    get_imass(5) = 5.9753843390072855d+23	!H
    get_imass(6) = 1.4943070555416133d+23	!HE
    get_imass(7) = 2.9876921695036427d+23	!H2
    get_imass(8) = 4.981023518472045d+22	!C
    get_imass(9) = 3.7357676388540333d+22	!O
    get_imass(10) = 3.5159527805593353d+22	!OH
    get_imass(11) = 2.1347243650594476d+22	!CO
    get_imass(12) = 4.5977587266153915d+22	!CH
    get_imass(13) = 4.269260529012849d+22	!CH2
    get_imass(14) = 2.4905117592360225d+22	!C2
    get_imass(15) = 2.0610912486362525d+22	!HCO
    get_imass(16) = 3.320568494480966d+22	!H2O
    get_imass(17) = 1.8678838194270166d+22	!O2
    get_imass(18) = 0d0	!CO_total
    get_imass(19) = 0d0	!H2O_total
    get_imass(20) = 5.978638635046188d+23	!H+
    get_imass(21) = 1.494510491502214d+23	!HE+
    get_imass(22) = 2.988505522030552d+23	!H2+
    get_imass(23) = 4.981249537914642d+22	!C+
    get_imass(24) = 3.7358947733482775d+22	!O+
    get_imass(25) = 2.061129946901891d+22	!HOC+
    get_imass(26) = 2.061129946901891d+22	!HCO+
    get_imass(27) = 1.9921562368806946d+23	!H3+
    get_imass(28) = 4.597951301414432d+22	!CH+
    get_imass(29) = 4.269426568397541d+22	!CH2+
    get_imass(30) = 2.1347658777582276d+22	!CO+
    get_imass(31) = 3.98471753627902d+22	!CH3+
    get_imass(32) = 3.516065393652204d+22	!OH+
    get_imass(33) = 3.320668939159153d+22	!H2O+
    get_imass(34) = 3.1458464665581883d+22	!H3O+
    get_imass(35) = 1.8679156025097617d+22	!O2+
    get_imass(36) = 1.4947139828621716d+23	!HE++
    get_imass(37) = 0d0	!CR
    get_imass(38) = 0d0	!g
    get_imass(39) = 0d0	!Tgas
    get_imass(40) = 0d0	!dummy

  end function get_imass

  !************************
  !species binding energies (surface=BARE), K
  function get_EbindBare()
    use krome_commons
    implicit none
    real*8::get_EbindBare(nspec)

    get_EbindBare(:) = 1d99

    get_EbindBare(idx_H) = 500.0d0
    get_EbindBare(idx_H2) = 300.0d0
    get_EbindBare(idx_O) = 1700.0d0
    get_EbindBare(idx_OH) = 1360.0d0
    get_EbindBare(idx_CO) = 1100.0d0
    get_EbindBare(idx_HCO) = 1100.0d0
    get_EbindBare(idx_H2O) = 4800.0d0
    get_EbindBare(idx_O2) = 1250.0d0

  end function get_EbindBare

  !************************
  !species binding energies (surface=ICE), K
  function get_EbindIce()
    use krome_commons
    implicit none
    real*8::get_EbindIce(nspec)

    get_EbindIce(:) = 1d99

    get_EbindIce(idx_H) = 650.0d0
    get_EbindIce(idx_H2) = 300.0d0
    get_EbindIce(idx_O) = 1700.0d0
    get_EbindIce(idx_OH) = 3500.0d0
    get_EbindIce(idx_CO) = 1300.0d0
    get_EbindIce(idx_HCO) = 3100.0d0
    get_EbindIce(idx_H2O) = 4800.0d0
    get_EbindIce(idx_O2) = 900.0d0

  end function get_EbindIce

  !************************
  function get_kevap70()
    use krome_commons
    implicit none
    real*8::get_kevap70(nspec)

    get_kevap70(idx_E) = 0d0
    get_kevap70(idx_Hk) = 0d0
    get_kevap70(idx_Ck) = 0d0
    get_kevap70(idx_Ok) = 0d0
    get_kevap70(idx_H) = 790490323.1199661
    get_kevap70(idx_HE) = 0d0
    get_kevap70(idx_H2) = 13763786733.050402
    get_kevap70(idx_C) = 0d0
    get_kevap70(idx_O) = 28.369278883298623
    get_kevap70(idx_OH) = 3649.8804308076387
    get_kevap70(idx_CO) = 149751.92964078687
    get_kevap70(idx_CH) = 0d0
    get_kevap70(idx_CH2) = 0d0
    get_kevap70(idx_C2) = 0d0
    get_kevap70(idx_HCO) = 149751.92964078687
    get_kevap70(idx_H2O) = 1.6588493815567146e-18
    get_kevap70(idx_O2) = 17568.77150646201
    get_kevap70(idx_CO_total) = 0d0
    get_kevap70(idx_H2O_total) = 0d0
    get_kevap70(idx_Hj) = 0d0
    get_kevap70(idx_HEj) = 0d0
    get_kevap70(idx_H2j) = 0d0
    get_kevap70(idx_Cj) = 0d0
    get_kevap70(idx_Oj) = 0d0
    get_kevap70(idx_HOCj) = 0d0
    get_kevap70(idx_HCOj) = 0d0
    get_kevap70(idx_H3j) = 0d0
    get_kevap70(idx_CHj) = 0d0
    get_kevap70(idx_CH2j) = 0d0
    get_kevap70(idx_COj) = 0d0
    get_kevap70(idx_CH3j) = 0d0
    get_kevap70(idx_OHj) = 0d0
    get_kevap70(idx_H2Oj) = 0d0
    get_kevap70(idx_H3Oj) = 0d0
    get_kevap70(idx_O2j) = 0d0
    get_kevap70(idx_HEjj) = 0d0
    get_kevap70(idx_CR) = 0d0
    get_kevap70(idx_g) = 0d0
    get_kevap70(idx_Tgas) = 0d0
    get_kevap70(idx_dummy) = 0d0

  end function get_kevap70

  !************************
  !get verbatim reaction names
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    !reaction names are loaded from file
    get_rnames(:) = reactionNames(:)

  end function get_rnames

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "E"
    get_names(2) = "H-"
    get_names(3) = "C-"
    get_names(4) = "O-"
    get_names(5) = "H"
    get_names(6) = "HE"
    get_names(7) = "H2"
    get_names(8) = "C"
    get_names(9) = "O"
    get_names(10) = "OH"
    get_names(11) = "CO"
    get_names(12) = "CH"
    get_names(13) = "CH2"
    get_names(14) = "C2"
    get_names(15) = "HCO"
    get_names(16) = "H2O"
    get_names(17) = "O2"
    get_names(18) = "CO_total"
    get_names(19) = "H2O_total"
    get_names(20) = "H+"
    get_names(21) = "HE+"
    get_names(22) = "H2+"
    get_names(23) = "C+"
    get_names(24) = "O+"
    get_names(25) = "HOC+"
    get_names(26) = "HCO+"
    get_names(27) = "H3+"
    get_names(28) = "CH+"
    get_names(29) = "CH2+"
    get_names(30) = "CO+"
    get_names(31) = "CH3+"
    get_names(32) = "OH+"
    get_names(33) = "H2O+"
    get_names(34) = "H3O+"
    get_names(35) = "O2+"
    get_names(36) = "HE++"
    get_names(37) = "CR"
    get_names(38) = "g"
    get_names(39) = "Tgas"
    get_names(40) = "dummy"

  end function get_names

  !************************
  !get cooling names list (empty element if cooling not present)
  function get_cooling_names()
    use krome_commons
    implicit none
    character*16::get_cooling_names(ncools)

    get_cooling_names(:) = ""

    get_cooling_names(idx_cool_h2) = "H2"
    get_cooling_names(idx_cool_h2gp) = "H2GP"
    get_cooling_names(idx_cool_atomic) = "ATOMIC"
    get_cooling_names(idx_cool_cen) = "CEN"
    get_cooling_names(idx_cool_hd) = "HD"
    get_cooling_names(idx_cool_z) = "Z"
    get_cooling_names(idx_cool_metal) = "METAL"
    get_cooling_names(idx_cool_dh) = "DH"
    get_cooling_names(idx_cool_enthalpic) = "ENTHALPIC"
    get_cooling_names(idx_cool_dust) = "DUST"
    get_cooling_names(idx_cool_compton) = "COMPTON"
    get_cooling_names(idx_cool_cie) = "CIE"
    get_cooling_names(idx_cool_continuum) = "CONTINUUM"
    get_cooling_names(idx_cool_cont) = "CONT"
    get_cooling_names(idx_cool_exp) = "EXP"
    get_cooling_names(idx_cool_expansion) = "EXPANSION"
    get_cooling_names(idx_cool_ff) = "FF"
    get_cooling_names(idx_cool_bss) = "BSS"
    get_cooling_names(idx_cool_custom) = "CUSTOM"
    get_cooling_names(idx_cool_co) = "CO"
    get_cooling_names(idx_cool_zcie) = "ZCIE"
    get_cooling_names(idx_cool_zcienouv) = "ZCIENOUV"
    get_cooling_names(idx_cool_zextend) = "ZEXTEND"
    get_cooling_names(idx_cool_gh) = "GH"
    get_cooling_names(idx_cool_oh) = "OH"
    get_cooling_names(idx_cool_h2o) = "H2O"
    get_cooling_names(idx_cool_hcn) = "HCN"

  end function get_cooling_names

  !************************
  !get heating names list (empty element if heating not present)
  function get_heating_names()
    use krome_commons
    implicit none
    character*16::get_heating_names(nheats)

    get_heating_names(:) = ""

    get_heating_names(idx_heat_chem) = "CHEM"
    get_heating_names(idx_heat_compress) = "COMPRESS"
    get_heating_names(idx_heat_compr) = "COMPR"
    get_heating_names(idx_heat_photo) = "PHOTO"
    get_heating_names(idx_heat_dh) = "DH"
    get_heating_names(idx_heat_enthalpic) = "ENTHALPIC"
    get_heating_names(idx_heat_photoav) = "PHOTOAV"
    get_heating_names(idx_heat_av) = "AV"
    get_heating_names(idx_heat_cr) = "CR"
    get_heating_names(idx_heat_dust) = "DUST"
    get_heating_names(idx_heat_xray) = "XRAY"
    get_heating_names(idx_heat_visc) = "VISC"
    get_heating_names(idx_heat_viscous) = "VISCOUS"
    get_heating_names(idx_heat_custom) = "CUSTOM"
    get_heating_names(idx_heat_zcie) = "ZCIE"

  end function get_heating_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_Hk) + &
        n(idx_H) + &
        n(idx_H2)*2d0 + &
        n(idx_OH) + &
        n(idx_CH) + &
        n(idx_CH2)*2d0 + &
        n(idx_HCO) + &
        n(idx_H2O_total)*2d0 + &
        n(idx_Hj) + &
        n(idx_H2j)*2d0 + &
        n(idx_HOCj) + &
        n(idx_HCOj) + &
        n(idx_H3j)*3d0 + &
        n(idx_CHj) + &
        n(idx_CH2j)*2d0 + &
        n(idx_CH3j)*3d0 + &
        n(idx_OHj) + &
        n(idx_H2Oj)*2d0 + &
        n(idx_H3Oj)*3d0
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 0	!E
    get_zatoms(2) = 1	!H-
    get_zatoms(3) = 6	!C-
    get_zatoms(4) = 8	!O-
    get_zatoms(5) = 1	!H
    get_zatoms(6) = 2	!HE
    get_zatoms(7) = 2	!H2
    get_zatoms(8) = 6	!C
    get_zatoms(9) = 8	!O
    get_zatoms(10) = 9	!OH
    get_zatoms(11) = 14	!CO
    get_zatoms(12) = 7	!CH
    get_zatoms(13) = 8	!CH2
    get_zatoms(14) = 12	!C2
    get_zatoms(15) = 15	!HCO
    get_zatoms(16) = 10	!H2O
    get_zatoms(17) = 16	!O2
    get_zatoms(18) = 14	!CO_total
    get_zatoms(19) = 10	!H2O_total
    get_zatoms(20) = 1	!H+
    get_zatoms(21) = 2	!HE+
    get_zatoms(22) = 2	!H2+
    get_zatoms(23) = 6	!C+
    get_zatoms(24) = 8	!O+
    get_zatoms(25) = 15	!HOC+
    get_zatoms(26) = 15	!HCO+
    get_zatoms(27) = 3	!H3+
    get_zatoms(28) = 7	!CH+
    get_zatoms(29) = 8	!CH2+
    get_zatoms(30) = 14	!CO+
    get_zatoms(31) = 9	!CH3+
    get_zatoms(32) = 9	!OH+
    get_zatoms(33) = 10	!H2O+
    get_zatoms(34) = 11	!H3O+
    get_zatoms(35) = 16	!O2+
    get_zatoms(36) = 2	!HE++
    get_zatoms(37) = 0	!CR
    get_zatoms(38) = 0	!g
    get_zatoms(39) = 0	!Tgas
    get_zatoms(40) = 0	!dummy

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !**************************
  function get_free_fall_time(n)
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),m(nspec)
    real*8::rhogas,get_free_fall_time

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    get_free_fall_time = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time

  !**************************
  function get_free_fall_time_rho(rhogas)
    use krome_constants
    implicit none
    real*8::rhogas,get_free_fall_time_rho

    get_free_fall_time_rho = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time_rho

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    num2col = 1.87d21*(max(ncalc,1d-40)*1d-3)**(2./3.)

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    col2num = 1d3 * (max(ncalc,1d-40)/1.87d21)**1.5

  end function col2num

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  - n(idx_Hk) &
        - n(idx_Ck) &
        - n(idx_Ok) &
        + n(idx_Hj) &
        + n(idx_HEj) &
        + n(idx_H2j) &
        + n(idx_Cj) &
        + n(idx_Oj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_H3j) &
        + n(idx_CHj) &
        + n(idx_CH2j) &
        + n(idx_COj) &
        + n(idx_CH3j) &
        + n(idx_OHj) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + n(idx_O2j) &
        + 2d0*n(idx_HEjj)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = -1.d0 	!E
    get_charges(2) = -1.d0 	!H-
    get_charges(3) = -1.d0 	!C-
    get_charges(4) = -1.d0 	!O-
    get_charges(5) = 0.d0 	!H
    get_charges(6) = 0.d0 	!HE
    get_charges(7) = 0.d0 	!H2
    get_charges(8) = 0.d0 	!C
    get_charges(9) = 0.d0 	!O
    get_charges(10) = 0.d0 	!OH
    get_charges(11) = 0.d0 	!CO
    get_charges(12) = 0.d0 	!CH
    get_charges(13) = 0.d0 	!CH2
    get_charges(14) = 0.d0 	!C2
    get_charges(15) = 0.d0 	!HCO
    get_charges(16) = 0.d0 	!H2O
    get_charges(17) = 0.d0 	!O2
    get_charges(18) = 0.d0 	!CO_total
    get_charges(19) = 0.d0 	!H2O_total
    get_charges(20) = 1.d0 	!H+
    get_charges(21) = 1.d0 	!HE+
    get_charges(22) = 1.d0 	!H2+
    get_charges(23) = 1.d0 	!C+
    get_charges(24) = 1.d0 	!O+
    get_charges(25) = 1.d0 	!HOC+
    get_charges(26) = 1.d0 	!HCO+
    get_charges(27) = 1.d0 	!H3+
    get_charges(28) = 1.d0 	!CH+
    get_charges(29) = 1.d0 	!CH2+
    get_charges(30) = 1.d0 	!CO+
    get_charges(31) = 1.d0 	!CH3+
    get_charges(32) = 1.d0 	!OH+
    get_charges(33) = 1.d0 	!H2O+
    get_charges(34) = 1.d0 	!H3O+
    get_charges(35) = 1.d0 	!O2+
    get_charges(36) = 2.d0 	!HE++
    get_charges(37) = 0.d0 	!CR
    get_charges(38) = 0.d0 	!g
    get_charges(39) = 0.d0 	!Tgas
    get_charges(40) = 0.d0 	!dummy

  end function get_charges

  !*****************************
  ! get metallicity using C as reference
  function get_metallicityC(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityC,zC,nH

    nH = get_Hnuclei(n(:))

    zC = n(idx_Ck) &
        + n(idx_C) &
        + n(idx_CH) &
        + n(idx_CH2) &
        + 2d0*n(idx_C2) &
        + n(idx_HCO) &
        + n(idx_CO_total) &
        + n(idx_Cj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_CHj) &
        + n(idx_CH2j) &
        + n(idx_COj) &
        + n(idx_CH3j)

    zC = max(zC, 0d0)

    get_metallicityC = log10(zC/nH+1d-40) - (-3.5700000000000003)

    phys_metallicity = get_metallicityC

  end function get_metallicityC

  !*****************************
  ! get metallicity using O as reference
  function get_metallicityO(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityO,zO,nH

    nH = get_Hnuclei(n(:))

    zO = n(idx_Ok) &
        + n(idx_O) &
        + n(idx_OH) &
        + n(idx_HCO) &
        + 2d0*n(idx_O2) &
        + n(idx_CO_total) &
        + n(idx_H2O_total) &
        + n(idx_Oj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_COj) &
        + n(idx_OHj) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + 2d0*n(idx_O2j)

    zO = max(zO, 0d0)

    get_metallicityO = log10(zO/nH+1d-40) - (-3.3100000000000005)

    phys_metallicity = get_metallicityO

  end function get_metallicityO

end module krome_getphys
!This module contains the functions and subroutines
! needed to evaluate the adiabatic index.

!############### MODULE ##############
module krome_gadiab
contains

  !#KROME_header

  !**************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part
  function gamma_pop(array_part,dT_part,min_part,Tgasin)
    implicit none
    real*8::array_part(:),dT_part
    real*8::min_part,Tgas,gamma_pop,Tgas2,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !temperature above minimum data point
    inTgas = max(Tgasin,min_part)

    !data index
    idx = (inTgas-min_part)/dT_part+1
    !corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !ln of partition functions (3 points forward)
    logz = log(array_part(idx))
    logz1 = log(array_part(idx+1))
    logz2 = log(array_part(idx+2))

    !derivative for mean energy (2 points forward)
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv1 = (emed2-emed1)/dT_part

    !next point temperature
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of partition functions
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part(idx+3))

    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv2 = (emed2-emed1)/dT_part

    !interpolation for 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns result
    gamma_pop = Cv

  end function gamma_pop

  !*****************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part, for H2 with a ortho/para
  ! ratio of opratio. Needs even and odd partition functions.
  function gamma_pop_H2(array_part_even,array_part_odd,dT_part,&
        min_part,Tgasin,opratio)
    implicit none
    real*8::array_part_even(:),array_part_odd(:),dT_part,zcut(4)
    real*8::min_part,Tgas,opratio,gamma_pop_H2,Tgas2,a,b,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !Tgas above the data limit
    inTgas = max(Tgasin,min_part)

    !exponents for ortho/para ratio
    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    !index in the data for the given Tgas
    idx = (inTgas-min_part)/dT_part+1
    !get the corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !needed for ortho partition function (see Boley+2007)
    zcut(1) = exp(2d0*85.4/Tgas)
    zcut(2) = exp(2d0*85.4/(Tgas+dT_part))
    zcut(3) = exp(2d0*85.4/(Tgas+2d0*dT_part))
    zcut(4) = exp(2d0*85.4/(Tgas+3d0*dT_part))

    !ln of the composite partition function
    logz = log(array_part_even(idx)**b*(3d0*array_part_odd(idx)*zcut(1))**a)
    logz1 = log(array_part_even(idx+1)**b*(3d0*array_part_odd(idx+1)*zcut(2))**a)
    logz2 = log(array_part_even(idx+2)**b*(3d0*array_part_odd(idx+2)*zcut(3))**a)
    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the left point
    Cv1 = (emed2-emed1)/dT_part

    !Tgas of the right point
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of the composite function
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part_even(idx+3)**b*(3d0*array_part_odd(idx+3)*zcut(4))**a)
    !derivative for the mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the right point
    Cv2 = (emed2-emed1)/dT_part

    !interpolation of 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns the result
    gamma_pop_H2 = Cv
  end function gamma_pop_H2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
      j1 = j + 1
      zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
      zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * (zo*exp(2d0*ibTgas))**a) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
      z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

    krome_gamma = 1.66666666667d0

    gamma_index = krome_gamma
  end function gamma_index

end module krome_gadiab
!This module contains functions and subroutines
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.

!############### MODULE ##############
module krome_grfuncs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !**********************
  !get Tdust from tables, K
  function get_table_Tdust(n) result(Tdust)
    use krome_commons
    use krome_fit
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::ntot,Tdust,Tgas

    Tgas = n(idx_Tgas)

    !default, K
    Tdust = 1d0

    !total densitym, cm-3
    ntot = sum(n(1:nmols))

    !zero density returns default
    if(ntot==0d0) return

    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas, &
        log10(ntot), log10(Tgas))

  end function get_table_Tdust

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
        * pre_kvgas_sqrt * ims * sqrTgas &
        * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
        + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
        expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
        + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !******************
  function krate_2bodySi(alpha,Ea,idx1,idx2,n,Tdust) result(krate)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec),Ea,Tdust,alpha
    integer,intent(in)::idx1,idx2
    real*8::krate,amin,amax,pexp,d2g,rho0

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    krate = krate_2body(n(:),idx1,idx2,alpha,amin,amax,pexp,d2g,rho0,Ea,Tdust)

  end function krate_2bodySi

  !********************
  !This routine has been modified to accomodate
  !Semenov framework  for surface chemistry.
  function krate_2body(n,idx1,idx2,alpha,amin,amax,pexp,d2g,rho0, &
        Ea,Tdust) result(krate)
    use krome_commons
    use krome_constants
    use krome_getphys
    implicit none
    integer,intent(in)::idx1,idx2
    real*8,intent(in)::n(nspec),amin,amax,pexp,d2g,rho0,Ea,Tdust,alpha
    real*8::rhog,p3,p4,ndns,krate,mred,fice,fbare,Preac,p
    real*8::iTd23,Ebare(nspec),Eice(nspec),mass(nspec)
    real*8,parameter::app2=(3d-8)**2 !cm^2 (Hocuk+2015)
    real*8,parameter::nu0=1d12 !1/s
    real*8,parameter::hbar=planck_erg/2d0/pi !erg*s
    real*8,parameter::ar=1d-8 !cm

    mass(:) = get_mass()

    !gas density, g/cm3
    rhog = sum(mass(1:nmols)*n(1:nmols))

    !exponentes
    p3 = pexp + 3d0
    p4 = pexp + 4d0

    !number of sites cm-3/mly
    ndns = rhog/(4d0/3d0*rho0*app2)*(amax**p3-amin**p3) &
        / (amax**p4-amin**p4) * p4 / p3

    !reduced mass
    mred = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2))

    !tunneling probability
    Preac = exp(-2d0*ar/hbar*sqrt(2d0*mred*Ea*boltzmann_erg))
    !exponent
    iTd23 = 2d0/3d0/Tdust

    !get Ebind, K
    Ebare(:) = get_Ebind_bare()

    !ice/bare fraction
    fbare = 1d0
    fice = (n(idx_H2O_total)-n(idx_H2O))/ndns
    fbare = 1d0 - fice

    Eice(:) = get_Ebind_ice()

    !compute rate
    krate = fbare*(exp(-Ebare(idx1)*iTd23)+exp(-Ebare(idx2)*iTd23))

    krate = krate + fice*(exp(-Eice(idx1)*iTd23)+exp(-Eice(idx2)*iTd23))

    !rate in cm3/s
    krate = nu0*Preac/ndns*krate

  end function krate_2body

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
      dust_get_inv_phi(i) = 0d0
      dephi = (4d0 * nndust(i) * pi * asize2(i))
      if(dephi.le.0d0) cycle
      dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
      Td100 = Tdust(i) * 1d-2
      dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
          + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !*************************
  function dust_stick(Tgas,Tdust)
    implicit none
    real*8,intent(in)::Tgas,Tdust
    real*8::dust_stick
    real*8::Tg100,Td100

    Tg100 = Tgas * 1d-2
    Td100 = Tdust * 1d-2
    dust_stick = 1d0/(1d0 + 0.4d0*sqrt(Tg100+Td100) &
        + 0.2d0*Tg100 + 0.08d0*Tg100**2)

  end function dust_stick

  !****************************
  !sticking rate (1/s), assuming power-law dust distribution
  ! example rate is
  !  @format:idx,R,P,rate
  !  1,CO,CO_ice,krate_stick(n(:),idx_CO,1d-7,1d-5,-3.5,3d0,1d-2)
  ! n(:): internal status array (number densities, temeperature, etc...)
  ! idx : index of the sticking species, e.g. idx_CO
  ! Tdust: dust temperature (assume same for all bins), K
  ! amin: min grain size, cm
  ! amax: max grain size, cm
  ! pexp: power-law exponent, usually -3.5
  ! rho0: bulk material density, g/cm3, e.g. 3 g/cm3 for silicates
  ! d2g: dust to gass mass ratio, usually 0.01
  function krate_stick(n,idx,Tdust,amin,amax,pexp,rho0,d2g) result(k)
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    real*8,intent(in)::n(nspec),Tdust,amin,amax,pexp,rho0,d2g
    real*8::k,imass(nspec),p4,p3,mass(nspec),rhod
    integer,intent(in)::idx

    !get inverse mass squared
    imass(:) = get_imass_sqrt()
    !get masses
    mass(:) = get_mass()
    !derived exponents
    p3 = pexp + 3.
    p4 = pexp + 4.

    !total dust density, g/cm3
    rhod = sum(n(1:nmols)*mass(1:nmols))*d2g

    !compute rate (1/s) coefficient assuming normalization
    k = pre_kvgas_sqrt*sqrt(n(idx_Tgas)) * imass(idx) &
        * rhod / (4./3.*rho0) * p4 / p3 &
        * (amax**p3-amin**p3) / (amax**p4-amin**p4) &
        * dust_stick(n(idx_Tgas),Tdust)

  end function krate_stick

  !********************************
  !compact version of krate_stick
  function krate_stickSi(n,idx,Tdust) result(k)
    use krome_commons
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,amin,amax,d2g,rho0,pexp

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    k = krate_stick(n(:),idx,Tdust,amin,amax,pexp,rho0,d2g)

  end function krate_stickSi

  !***************************
  !evaporation rate, 1/s
  function krate_evaporation(n,idx,Tdust) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,Ebind(nspec),nu0

    nu0 = 1d12 !1/s
    Ebind(:) = get_EbindBare()

    k = nu0 * exp(-Ebind(idx)/Tdust)

  end function krate_evaporation

  !***************************
  !non-thermal evaporation rate (1/s) following Hollenbach 2009,
  ! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0809.1642
  !Gnot is the habing flux (1.78 is Draine)
  !Av is the visual extinction
  !crflux the ionization flux of cosmic rays, 1/s
  !yield is the efficiency of the photons to desorb the given molecule
  function krate_nonthermal_evaporation(idx, Gnot, Av, crflux, yield) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,parameter::crnot=1.3d-17
    real*8,parameter::Fnot=1d8 !desorbing photons flux, 1/s
    real*8,parameter::ap2=(3d-8)**2 !sites separation squared, cm2
    real*8,intent(in)::Gnot, Av, crflux, yield
    real*8::k,f70,kevap70(nspec)

    f70 = 3.16d-19*crflux/crnot
    kevap70(:) = get_kevap70()

    k = Gnot*Fnot*ap2*yield*exp(-1.8*Av)
    k = k + f70*kevap70(idx)

  end function krate_nonthermal_evaporation

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
      dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

end module krome_grfuncs
!This module mainly contains shielding routine and
! function to initialize radiation background (e.g. Planck).

!############### MODULE ##############
module krome_phfuncs
contains

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
    NHtot  = NHtot + num2col(n(idx_H),n(:))
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))

    shield_dust = exp(-sigma_d*NHtot)

  end function shield_dust

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      planckBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
      a = 2d0/planck_eV**2/clight**2
      planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

  end function krome_fshield

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
      fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
      fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
        0.035d0/sqrt(1d0+x) * &
        exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

end module krome_phfuncs

!############### MODULE ##############
module krome_subs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::invTgas
    real*8::small,nmax
    integer::i
    real*8::Hnuclei  !preproc from coevar
    real*8::ntot !preproc from coevar
    real*8:: T32  !preproc from coevar
    real*8:: Te  !preproc from coevar
    real*8:: invT  !preproc from coevar
    real*8:: lnTe  !preproc from coevar
    real*8:: T  !preproc from coevar
    real*8:: invT32  !preproc from coevar
    real*8::invTe !preproc from coevar
    real*8::logT !preproc from coevar
    real*8::invsqrT !preproc from coevar
    real*8::kl11  !preproc from coevar
    real*8::kh11  !preproc from coevar
    real*8::ncr11  !preproc from coevar
    real*8::ax11 !preproc from coevar
    real*8::a1 !preproc from coevar
    real*8::a2 !preproc from coevar
    real*8::a3 !preproc from coevar
    real*8::a4 !preproc from coevar
    real*8::a5 !preproc from coevar
    real*8::a6 !preproc from coevar
    real*8::a7 !preproc from coevar
    real*8::a8 !preproc from coevar
    real*8::a9 !preproc from coevar
    real*8::a10 !preproc from coevar
    real*8::a11 !preproc from coevar
    real*8::a12 !preproc from coevar
    real*8::asav0 !preproc from coevar
    real*8::asav1 !preproc from coevar
    real*8::asav2 !preproc from coevar
    real*8::asav3 !preproc from coevar
    real*8::asav4 !preproc from coevar
    real*8::asav5 !preproc from coevar
    real*8::asav6 !preproc from coevar
    real*8::asav7 !preproc from coevar
    real*8::bsav0 !preproc from coevar
    real*8::bsav1 !preproc from coevar
    real*8::bsav2 !preproc from coevar
    real*8::bsav3 !preproc from coevar
    real*8::bsav4 !preproc from coevar
    real*8::bsav5 !preproc from coevar
    real*8::bsav6 !preproc from coevar
    real*8::bsav7 !preproc from coevar
    real*8::sumsav !preproc from coevar
    real*8::sumsbv !preproc from coevar
    real*8::kHdiss !preproc from coevar
    real*8::kLdiss !preproc from coevar
    real*8::ncrdiss !preproc from coevar
    real*8::kHdissH2 !preproc from coevar
    real*8::kLdissH2 !preproc from coevar
    real*8::ncrdissH2 !preproc from coevar
    real*8::u1  !preproc from coevar
    real*8::u2  !preproc from coevar
    real*8::u3  !preproc from coevar
    real*8::tabTdust  !preproc from coevar
    real*8:: user_H2self  !preproc from coevar
    real*8:: HnOj  !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 1d-40/(nmax*nmax*nmax)

    invTgas = 1.d0/Tgas !inverse of T (1/K)

    Hnuclei = get_Hnuclei(n(:))
    ntot = sum(n(1:nmols))
    T32 = Tgas/3d2
    Te = Tgas*8.617343d-5
    invT = 1d0/Tgas
    lnTe = log(Te)
    T = Tgas
    invT32 = 1d0/T32
    invTe = 1d0/Te
    logT = log10(Tgas)
    invsqrT = 1d0/sqrt(Tgas)
    kl11 = 1d1**(-27.029d0+3.801d0*logT-29487d0*invT)
    kh11 = 1d1**(-2.729d0-1.75d0*logT-23474d0*invT)
    ncr11 = 1d1**(5.0792d0*(1d0-1.23d-5*(Tgas-2d3)))
    ax11 = 1.d0/(1.d0+(Hnuclei/(ncr11+1d-40)))
    a1 = 1.3500e-09
    a2 = 9.8493e-02
    a3 = 3.2852e-01
    a4 = 5.5610e-01
    a5 = 2.7710e-07
    a6 = 2.1826e+00
    a7 = 6.1910e-03
    a8 = 1.0461e+00
    a9 = 8.9712e-11
    a10 = 3.0424e+00
    a11 = 3.2576e-14
    a12 = 3.7741e+00
    asav0 = -1.9153214d2
    asav1 = 4.0129114d2
    asav2 = -3.7446991d2
    asav3 = 1.9078410d2
    asav4 = -5.7263467d1
    asav5 = 1.0133210d1
    asav6 = -9.8012853d-1
    asav7 = 4.0023414d-2
    bsav0 = -8.8755774d3
    bsav1 = 1.0081246d4
    bsav2 = -4.8606622d3
    bsav3 = 1.2889659d3
    bsav4 = -2.0319575d2
    bsav5 = 1.9057493d1
    bsav6 = -9.8530668d-1
    bsav7 = 2.1675387d-2
    sumsav = asav0+asav1*log10(Tgas)+asav2*(log10(Tgas))**2+asav3*(log10(Tgas))**3+asav4*(log10(Tgas))**4+asav5*(log10(Tgas))**5+asav6*(log10(Tgas))**6+asav7*(log10(Tgas))**7
    sumsbv = bsav0+bsav1*log10(Tgas)+bsav2*(log10(Tgas))**2+bsav3*(log10(Tgas))**3+bsav4*(log10(Tgas))**4+bsav5*(log10(Tgas))**5+bsav6*(log10(Tgas))**6+bsav7*(log10(Tgas))**7
    kHdiss = 3.52d-9*exp(-4.39d4*invT) + 1d-40
    kLdiss = 6.67d-12*sqrt(Tgas)*exp(-(1d0+63590.*invT)) + 1d-40
    ncrdiss = 1d1**(3. - 0.416*log10(Tgas/1d4) - 0.327*log10(Tgas/1d4)**2)
    kHdissH2 = 1.3d-9*exp(-5.33d4*invT) + 1d-40
    kLdissH2 = 5.996d-30*Tgas**4.1881/(1.+6.761d-6*Tgas)**5.6881 * exp(-5.46574d4*invT) + 1d-40
    ncrdissH2 = 1d1**(4.845 - 1.3*log10(Tgas/1d4) + 1.62*log10(Tgas/1d4)**2)
    u1 = 11.26d0*invTe
    u2 = 8.2d0*invTe
    u3 = 13.6*invTe
    tabTdust = get_table_Tdust(n(:))
    user_H2self = fselfH2(1.87d21*(merge(n(idx_H2),n_global(idx_H2),n(idx_H2) > 0.0_8)*1d-3)**(2./3.), 1d5)
    HnOj = fHnOj(user_Av)

    k(:) = small !inizialize coefficients

    !H + E -> H+ + E + E
    k(1) = small + (exp(-32.71396786d0+13.5365560d0&
        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
        *(lnTe**7)-2.03914985d-6*(lnTe**8)))

    !H+ + E -> H
    if(Tgas.GE.2.73d0 .and. Tgas.LE.5.5d3) then
      k(2) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !H+ + E -> H
    if(Tgas.GT.5.5d3 .and. Tgas.LT.1d8) then
      k(3) = small + (exp(-28.61303380689232d0-0.7241125657826851d0&
          *lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0&
          *lnTe**3-0.0003212605213188796d0&
          *lnTe**4-0.00001421502914054107d0&
          *lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7&
          *lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9&
          *lnTe**9))
    end if

    !HE + E -> HE+ + E + E
    k(4) = small + (exp(-44.09864886d0+23.91596563d0&
        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
        *(lnTe**7)-3.64916141d-6*(lnTe**8)))

    !HE+ + E -> HE
    if(Tgas.GE.2.73d0 .and. Tgas.LE.9.28d3) then
      k(5) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !HE+ + E -> HE
    if(Tgas.GT.9.28d3 .and. Tgas.LT.1d8) then
      k(6) = small + (1.54d-9&
          *(1.d0+0.3d0&
          /exp(8.099328789667d0&
          *invTe))&
          /(exp(40.49664394833662d0*invTe)&
          *Te**1.5d0)+3.92d-13&
          /Te**0.6353d0)
    end if

    !HE+ + E -> HE++ + E + E
    k(7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
        *lnTe**7-3.165581065665d-6*lnTe**8))

    !HE+ + H -> HE + H+
    k(8) = small + (1.25d-15*(T32)**.25)

    !HE + H+ -> HE+ + H
    if(Tgas.LT.1d4) then
      k(9) = small + (1.26d-9&
          *Tgas**(-.75)*exp(-1.275d5*invT))
    end if

    !HE + H+ -> HE+ + H
    if(Tgas.GE.1d4) then
      k(10) = small + (4d-37&
          *Tgas**4.74)
    end if

    !H2 + HE -> H + H + HE
    k(11) = small + (kh11**(1.-ax11)&
        *kl11**ax11)

    !H2 + HE+ -> HE + H2+
    k(12) = small + (7.2d-15)

    !H2 + HE+ -> HE + H + H+
    k(13) = small + (3.7d-14*exp(-35d0&
        *invT))

    !H2 + HE+ -> HE+ + H + H
    k(14) = small + (3d-11*sqrt(T32)&
        *exp(-5.2d4*invT))

    !HE++ + E -> HE+
    k(15) = small + (1.891d-10/(sqrt(Tgas&
        /9.37)&
        *(1.+sqrt(Tgas/9.37))**0.2476&
        *(1.+sqrt(Tgas&
        /2.774d6))**1.7524))

    !H + E -> H-
    k(16) = small + (1.4d-18*Tgas**0.928&
        *exp(-Tgas&
        /16200.))

    !H- + H -> H2 + E
    k(17) = small + (a1*(Tgas**a2+a3&
        *Tgas**a4+a5*Tgas**a6)&
        /(1.+a7*Tgas**a8+a9*Tgas**a10+a11&
        *Tgas**a12))

    !H + H+ -> H2+
    if(Tgas.LT.3d1) then
      k(18) = small + (2.10e-20&
          *(Tgas&
          /30.)**(-0.15))
    end if

    !H + H+ -> H2+
    if(Tgas.GE.3d1) then
      k(19) = small + (10**(-18.20-3.194&
          *log10(Tgas)+1.786*(log10(Tgas))**2-0.2072&
          *(log10(Tgas))**3))
    end if

    !H2+ + H -> H2 + H+
    k(20) = small + (6d-10)

    !H2 + H+ -> H2+ + H
    if(Tgas.LT.1d5) then
      k(21) = small + (1d1**sumsav)
    end if

    !H2 + H+ -> H2+ + H
    if(Tgas.GE.1d5) then
      k(22) = small + (1d1**sumsbv)
    end if

    !H2 + E -> H + H + E
    k(23) = small + (4.38d-10*exp(-1.02d5&
        *invT)*Tgas**(0.35))

    !H2 + H -> H + H + H
    k(24) = small + (1d1**(log10(kHdiss)-log10(kHdiss&
        /kLdiss)/(1d0+ntot/ncrdiss)))

    !H- + E -> H + E + E
    k(25) = small + (exp(-18.01849334273d0+2.360852208681d0&
        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
        *lnTe**7-2.631285809207d-6*lnTe**8))

    !H- + H -> H + H + E
    if(Tgas.LE.1.16d3) then
      k(26) = small + (2.56d-9&
          *Te**1.78186)
    end if

    !H- + H -> H + H + E
    if(Tgas.GT.1.16d3) then
      k(27) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
          *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
          *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
          *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
          *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
          *lnTe**9))
    end if

    !H- + H+ -> H + H
    k(28) = small + ((2.96d-6&
        /sqrt(Tgas) - 1.73d-9 + 2.50d-10&
        *sqrt(Tgas)))

    !H- + H+ -> H2+ + E
    k(29) = small + (1.d-8*Tgas**(-0.4))

    !H2+ + E -> H + H
    if(Tgas.LE.6.17d2) then
      k(30) = small + (1.d-8)
    end if

    !H2+ + E -> H + H
    if(Tgas.GT.6.17d2) then
      k(31) = small + (1.32d-6&
          *Tgas**(-0.76))
    end if

    !H2+ + H- -> H + H2
    k(32) = small + (5d-7*sqrt(1d2*invT))

    !H2 + H2 -> H2 + H + H
    k(33) = small + (1d1**(log10(kHdissH2)-log10(kHdissH2&
        /kLdissH2)/(1d0+ntot/ncrdissH2)))

    !H + H + HE -> H2 + HE
    k(34) = small + (6.9d-32&
        *Tgas**(-.4))

    !H + H + H -> H2 + H
    k(35) = small + (6d-32&
        *Tgas**(-.25) + 2d-31*Tgas**(-.5))

    !H2 + H + H -> H2 + H2
    k(36) = small + ((6d-32&
        *Tgas**(-0.25) + 2d-31*Tgas**(-.5)) &
        / 8d0)

    !C+ + E -> C
    if(Tgas.LE.7950d0) then
      k(37) = small + (4.67d-12&
          *(T32)**(-0.6))
    end if

    !C+ + E -> C
    if(Tgas.GT.7950d0 .and. Tgas.LE.21140d0) then
      k(38) = small + (1.23d-17&
          *(T32)**2.49*exp(21845.6d0*invT))
    end if

    !C+ + E -> C
    if(Tgas.GT.21140d0) then
      k(39) = small + (9.62d-8&
          *(T32)**(-1.37)*exp(-115786.2d0*invT))
    end if

    !O+ + E -> O
    if(Tgas.LE.4d2) then
      k(40) = small + (1.30d-10&
          *(Tgas)**(-0.64))
    end if

    !O+ + E -> O
    if(Tgas.GT.4d2) then
      k(41) = small + (1.41d-10&
          *(Tgas)**(-0.66) + 7.4d-4*(Tgas)**(-1.5)*exp(-1.75d5*invT)&
          *(1d0 + 0.062d0*exp(-1.45d5*invT)))
    end if

    !C + E -> C+ + E + E
    k(42) = small + (6.85d-8*u1**0.25&
        *exp(-u1)&
        /(0.193d0+u1))

    !O + E -> O+ + E + E
    k(43) = small + (3.59d-8*u3**0.34&
        *exp(-u3)&
        /(0.073d0+u3))

    !O+ + H -> O + H+
    k(44) = small + (4.99d-11&
        *Tgas**0.405 + 7.54d-10*invT**(0.458))

    !O + H+ -> O+ + H
    k(45) = small + ((1.08d-11&
        *Tgas**0.517 + 4d-10*Tgas**(0.00669))*exp(-2.27d2*invT))

    !O + HE+ -> O+ + HE
    k(46) = small + (4.991d-15*(Tgas&
        *1d-4)**0.3794*exp(-Tgas*8.9206d-7) + 2.78d-15*(Tgas&
        *1d-4)**(-0.2163)*exp(-Tgas*1.2258d-6))

    !C + H+ -> C+ + H
    k(47) = small + (3.9d-16*Tgas**(0.213))

    !C+ + H -> C + H+
    k(48) = small + (6.08d-14*(Tgas&
        *1d-4)**(1.96)*exp(-1.7d5*invT))

    !C + HE+ -> C+ + HE
    if(Tgas.LE.2d2) then
      k(49) = small + (8.58d-17&
          *Tgas**(0.757))
    end if

    !C + HE+ -> C+ + HE
    if(Tgas.GT.2d2 .and. Tgas.LE.2d3) then
      k(50) = small + (3.25d-17&
          *Tgas**(0.968))
    end if

    !C + HE+ -> C+ + HE
    if(Tgas.GT.2d3 .and. Tgas.LT.1d8) then
      k(51) = small + (2.77d-19&
          *Tgas**(1.597))
    end if

    !OH + H -> O + H + H
    k(52) = small + (6d-9*exp(-5.09d4&
        *invT))

    !HOC+ + H2 -> HCO+ + H2
    k(53) = small + (3d-10)

    !HOC+ + CO -> HCO+ + CO
    if(Tgas.LT.1d1) then
      k(54) = small + (1.604d-9)
    end if

    !HOC+ + CO -> HCO+ + CO
    if(Tgas.GE.1d1) then
      k(55) = small + (8.68d-10&
          *(1.+2.42717d-2*sqrt(3e2*invT)+7.1537*invT))
    end if

    !C + H2 -> CH + H
    k(56) = small + (6.64d-10*exp(-11700d0&
        *invT))

    !CH + H -> C + H2
    k(57) = small + (1.31d-10*exp(-8d1*invT))

    !CH + H2 -> CH2 + H
    k(58) = small + (5.46d-10*exp(-1943d0&
        *invT))

    !CH + C -> C2 + H
    k(59) = small + (2.40d-10)

    !CH + O -> CO + H
    k(60) = small + (1.02d-10*exp(-914d0&
        *invT))

    !CH + O -> HCO+ + E
    k(61) = small + (1.9d-11*(T32)**(-2.2)&
        *exp(-165.1d0*invT))

    !CH + O -> OH + C
    k(62) = small + (2.52d-11*exp(-2381d0&
        *invT))

    !CH2 + H -> CH + H2
    k(63) = small + (2.2d-10)

    !CH2 + O -> CO + H + H
    k(64) = small + (2.04d-10*exp(-270d0&
        *invT))

    !CH2 + O -> CO + H2
    k(65) = small + (1.36d-10*exp(-270d0&
        *invT))

    !CH2 + O -> HCO + H
    k(66) = small + (5.01d-11)

    !CH2 + O -> CH + OH
    k(67) = small + (4.98d-10*exp(-6d3&
        *invT))

    !C2 + O -> CO + C
    if(Tgas.LT.3d2) then
      k(68) = small + (2d-12&
          *(T32)**(-.12))
    end if

    !C2 + O -> CO + C
    if(Tgas.GE.3d2) then
      k(69) = small + (2d-12&
          *(T32)**(.757))
    end if

    !O + H2 -> OH + H
    k(70) = small + (1.46e-12*exp(-9650.&
        *invT))

    !OH + H -> O + H2
    if(Tgas.LT.280.d0) then
      k(71) = small + (6.99d-14&
          *T32**2.8*exp(-1950d0*invT))
    end if

    !OH + H -> O + H2
    if(Tgas.GE.280.d0) then
      k(72) = small + (5.45d-17)
    end if

    !H2 + OH -> H2O + H
    k(73) = small + (3.6d-16*T**(1.52)&
        *exp(-1.74d3*invT))

    !C + OH -> H + CO
    if(Tgas.LT.1d1) then
      k(74) = small + (7.051e-11)
    end if

    !C + OH -> H + CO
    if(Tgas.GE.1d1) then
      k(75) = small + (2.25d-11&
          *(T32)**(-.339)*exp(-.108d0*invT))
    end if

    !O + OH -> H + O2
    if(Tgas.GE.150.d0) then
      k(76) = small + (2.4d-11&
          *exp(110d0*invT))
    end if

    !O + OH -> H + O2
    if(Tgas.LT.150.d0) then
      k(77) = small + (4.997d-11)
    end if

    !OH + OH -> H2O + O
    k(78) = small + (1.65d-12*(T32)**1.14&
        *exp(-5d1*invT))

    !H2O + H -> H2 + OH
    k(79) = small + (1.59d-11*(T32)*1.2&
        *exp(-9610.*invT))

    !O2 + H -> OH + O
    k(80) = small + (2.61d-10*1.2*exp(-8156.&
        *invT))

    !O2 + H2 -> OH + OH
    k(81) = small + (3.16d-10*exp(-21890.d0&
        *invT))

    !O2 + C -> CO + O
    if(Tgas.LT.1052d0) then
      k(82) = small + (4.7d-11&
          *T32**(-.34))
    end if

    !O2 + C -> CO + O
    if(Tgas.GE.1052d0) then
      k(83) = small + (2.48d-12&
          *T32**1.54*exp(613d0*invT))
    end if

    !CO + H -> C + OH
    k(84) = small + (1.1d-10*T32**0.5&
        *exp(-77700d0*invT))

    !H2+ + H2 -> H3+ + H
    k(85) = small + (2.24d-9*T32**.042&
        *exp(-Tgas&
        /46600.))

    !H3+ + H -> H2+ + H2
    k(86) = small + (7.7d-9*exp(-17560d0&
        *invT))

    !C + H2+ -> CH+ + H
    k(87) = small + (2.4d-9)

    !C + H3+ -> CH+ + H2
    k(88) = small + ((1.0218d-9 + 7.2733d-11&
        *sqrt(Tgas) + 5.9203d-14*Tgas)&
        /(Tgas**0.1667 + 4.4914d-2&
        *sqrt(Tgas) - 5.9203d-14*Tgas + 2.6397d-6*Tgas**1.5))

    !C + H3+ -> CH2+ + H
    k(89) = small + ((8.5145d-10)&
        /(Tgas**(.1667) + 9.5666d-4&
        *sqrt(Tgas) - 4.404d-5*Tgas + 2.3496d-6 * Tgas**1.5))

    !C+ + H2 -> CH+ + H
    k(90) = small + (1d-10*exp(-4640d0&
        *invT))

    !CH+ + H -> C+ + H2
    k(91) = small + (7.5d-10)

    !CH+ + H2 -> CH2+ + H
    k(92) = small + (1.2d-9)

    !CH+ + O -> CO+ + H
    k(93) = small + (3.5d-10)

    !CH2+ + H -> CH+ + H2
    k(94) = small + (1d-9*exp(-7080d0&
        *invT))

    !CH2+ + H2 -> CH3+ + H
    k(95) = small + (1.6d-9)

    !CH2+ + O -> HCO+ + H
    k(96) = small + (7.5d-10)

    !CH3+ + H -> CH2+ + H2
    k(97) = small + (7.d-10*exp(-10560d0&
        *invT))

    !CH3+ + O -> HOC+ + H2
    k(98) = small + (2.5d-10)

    !CH3+ + O -> HCO+ + H2
    k(99) = small + (2.5d-10)

    !C2 + O+ -> CO+ + C
    k(100) = small + (4.8d-10)

    !O+ + H2 -> H + OH+
    k(101) = small + (1.69d-9)

    !O + H2+ -> H + OH+
    k(102) = small + (1.5d-9)

    !O + H3+ -> H2 + OH+
    k(103) = small + (7.98d-10&
        *T32**(-.156)*exp(-1.41d0*invT))

    !O + H3+ -> H + H2O+
    k(104) = small + (3.42d-10&
        *T32**(-.156)*exp(-1.41d0*invT))

    !OH + H3+ -> H2 + H2O+
    if(Tgas.LT.1d1) then
      k(105) = small + (2.277d-8)
    end if

    !OH + H3+ -> H2 + H2O+
    if(Tgas.GE.1d1) then
      k(106) = small + (1.52d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH + C+ -> H + CO+
    if(Tgas.LT.1d1) then
      k(107) = small + (1.371d-08)
    end if

    !OH + C+ -> H + CO+
    if(Tgas.GE.1d1) then
      k(108) = small + (9.15d-10&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH+ + H2 -> H2O+ + H
    k(109) = small + (1.01d-9)

    !H2O+ + H2 -> H3O+ + H
    k(110) = small + (6.4d-10)

    !H2O + H3+ -> H2 + H3O+
    if(Tgas.LT.1d1) then
      k(111) = small + (2.55d-8)
    end if

    !H2O + H3+ -> H2 + H3O+
    if(Tgas.GE.1d1) then
      k(112) = small + (1.73d-9&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + C+ -> HOC+ + H
    k(113) = small + (1.8d-9)

    !H2O + C+ -> HCO+ + H
    if(Tgas.LT.1d1) then
      k(114) = small + (5.027d-9)
    end if

    !H2O + C+ -> HCO+ + H
    if(Tgas.GE.1d1) then
      k(115) = small + (3.4093e-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + C+ -> H2O+ + C
    k(116) = small + (2.4d-10)

    !H3O+ + C -> HCO+ + H2
    k(117) = small + (1d-11)

    !O2 + C+ -> CO+ + O
    k(118) = small + (3.42d-10)

    !O2 + C+ -> CO + O+
    k(119) = small + (4.53d-10)

    !O2 + CH2+ -> HCO+ + OH
    k(120) = small + (9.1d-10)

    !C + O2+ -> O + CO+
    k(121) = small + (5.2d-11)

    !C + O2+ -> O2 + C+
    k(122) = small + (5.2d-11)

    !CO + H3+ -> H2 + HCO+
    if(Tgas.LT.1d1) then
      k(123) = small + (2.468d-9)
    end if

    !CO + H3+ -> H2 + HCO+
    if(Tgas.GE.1d1) then
      k(124) = small + (1.88055d-9&
          *(1d0 + 0.02427d0*(3d2*invT)**.5 + 1.79558d0*invT))
    end if

    !CO + H3+ -> H2 + HOC+
    if(Tgas.LT.1d1) then
      k(125) = small + (1.421d-10)
    end if

    !CO + H3+ -> H2 + HOC+
    if(Tgas.GE.1d1) then
      k(126) = small + (1.08256d-10&
          *(1d0 + 0.02427d0*(3d2*invT)**.5 + 1.79558d0*invT))
    end if

    !HCO+ + C -> CO + CH+
    k(127) = small + (1.1d-9)

    !HCO+ + H2O -> CO + H3O+
    if(Tgas.LT.1d1) then
      k(128) = small + (7.279e-08)
    end if

    !HCO+ + H2O -> CO + H3O+
    if(Tgas.GE.1d1) then
      k(129) = small + (8.34d-10&
          *(1d0 + 0.5232d0*(3d2*invT)**.5 + 834.165880*invT))
    end if

    !CH + H+ -> CH+ + H
    if(Tgas.LT.1d1) then
      k(130) = small + (3.297d-8)
    end if

    !CH + H+ -> CH+ + H
    if(Tgas.GE.1d1) then
      k(131) = small + (3.54e-09&
          *(0.62d0 + 1.587411d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H2 + CH+
    if(Tgas.LT.1.5d2) then
      k(132) = small + (1.765d-9&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H2 + CH+
    if(Tgas.GE.1.5d2) then
      k(133) = small + (1.765d-9&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.66255d0*invT))
    end if

    !CH2 + H+ -> H + CH2+
    if(Tgas.LT.1.5d2) then
      k(134) = small + (1.765d-9&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H + CH2+
    if(Tgas.GE.1.5d2) then
      k(135) = small + (1.765d-9&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.66255d0*invT))
    end if

    !CH2 + HE+ -> HE + H2 + C+
    if(Tgas.LT.1.5d2) then
      k(136) = small + (9.65d-10&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + HE+ -> HE + H2 + C+
    if(Tgas.GE.1.5d2) then
      k(137) = small + (9.65d-10&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.6625498765d0&
          *invT))
    end if

    !CH2 + HE+ -> HE + H + CH+
    if(Tgas.LT.1.5d2) then
      k(138) = small + (9.65d-10&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + HE+ -> HE + H + CH+
    if(Tgas.GE.1.5d2) then
      k(139) = small + (9.65d-10&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.6625498765d0&
          *invT))
    end if

    !C2 + HE+ -> C+ + C + HE
    k(140) = small + (1.6d-9)

    !OH + H+ -> OH+ + H
    if(Tgas.LT.1d1) then
      k(141) = small + (3.745d-8)
    end if

    !OH + H+ -> OH+ + H
    if(Tgas.GE.1d1) then
      k(142) = small + (2.5d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH + HE+ -> O+ + HE + H
    if(Tgas.LT.1d1) then
      k(143) = small + (2.022d-8)
    end if

    !OH + HE+ -> O+ + HE + H
    if(Tgas.GE.1d1) then
      k(144) = small + (1.35d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !H2O + H+ -> H + H2O+
    if(Tgas.LT.1d1) then
      k(145) = small + (4.202d-8)
    end if

    !H2O + H+ -> H + H2O+
    if(Tgas.GE.1d1) then
      k(146) = small + (2.85d-9&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + OH + H+
    if(Tgas.LT.1d1) then
      k(147) = small + (7.562d-9)
    end if

    !H2O + HE+ -> HE + OH + H+
    if(Tgas.GE.1d1) then
      k(148) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + OH+ + H
    if(Tgas.LT.1d1) then
      k(149) = small + (7.562d-9)
    end if

    !H2O + HE+ -> HE + OH+ + H
    if(Tgas.GE.1d1) then
      k(150) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + H2O+
    if(Tgas.LT.1d1) then
      k(151) = small + (7.56d-9)
    end if

    !H2O + HE+ -> HE + H2O+
    if(Tgas.GE.1d1) then
      k(152) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !O2 + H+ -> O2+ + H
    k(153) = small + (2d-9)

    !O2 + HE+ -> O2+ + HE
    k(154) = small + (3.3d-11)

    !O2 + HE+ -> O+ + HE + O
    k(155) = small + (1.1d-9)

    !CO + HE+ -> C+ + HE + O
    k(156) = small + (1.4d-9&
        *(T32)**(-.5))

    !CO + HE+ -> C + HE + O+
    k(157) = small + (1.4d-16&
        *(T32)**(-.5))

    !CO+ + H -> CO + H+
    k(158) = small + (7.5d-10)

    !C- + H+ -> C + H
    k(159) = small + (2.3d-7*(T32)**(-.5))

    !O- + H+ -> O + H
    k(160) = small + (2.3d-7*(T32)**(-.5))

    !HE+ + H- -> H + HE
    k(161) = small + (2.3d-7*T32**(-.5))

    !H3+ + E -> H2 + H
    k(162) = small + (2.34d-8*T32**(-.52))

    !H3+ + E -> H + H + H
    k(163) = small + (4.36d-8&
        *T32**(-.52))

    !CH+ + E -> C + H
    k(164) = small + (7d-8*T32**(-.5))

    !CH2+ + E -> CH + H
    k(165) = small + (1.6d-7*T32**(-.6))

    !CH2+ + E -> C + H2
    k(166) = small + (7.68d-8*T32**(-.6))

    !CH2+ + E -> C + H + H
    k(167) = small + (4.03d-7&
        *T32**(-.6))

    !CH3+ + E -> CH2 + H
    k(168) = small + (7.75d-8*T32**(-.5))

    !CH3+ + E -> CH + H2
    k(169) = small + (1.95d-7*T32**(-.5))

    !CH3+ + E -> CH + H + H
    k(170) = small + (2d-7*T32**(-.5))

    !OH+ + E -> O + H
    k(171) = small + (6.3d-9*T32**(-.48))

    !H2O+ + E -> O + H2
    k(172) = small + (3.9d-8*T32**(-.5))

    !H2O+ + E -> OH + H
    k(173) = small + (8.6d-8*T32**(-.5))

    !H2O+ + E -> O + H + H
    k(174) = small + (3.05d-7&
        *T32**(-.5))

    !H3O+ + E -> OH + H + H
    k(175) = small + (2.58d-7&
        *T32**(-.5))

    !H3O+ + E -> O + H + H2
    k(176) = small + (5.6d-9&
        *T32**(-.5))

    !H3O+ + E -> H + H2O
    k(177) = small + (1.08d-7*T32**(-.5))

    !H3O+ + E -> OH + H2
    k(178) = small + (6.02d-7*T32**(-.5))

    !O2+ + E -> O + O
    k(179) = small + (1.95d-7*T32**(-.7))

    !CO+ + E -> C + O
    k(180) = small + (2.75d-7*T32**(-.55))

    !HCO+ + E -> CO + H
    k(181) = small + (2.76d-7*T32**(-.64))

    !HCO+ + E -> OH + C
    k(182) = small + (2.4d-8*T32**(-.64))

    !HOC+ + E -> CO + H
    k(183) = small + (1.1d-7*invT32)

    !H- + C -> CH + E
    k(184) = small + (1d-9)

    !H- + O -> OH + E
    k(185) = small + (1d-10)

    !H- + OH -> H2O + E
    k(186) = small + (5d-10)

    !C- + H -> CH + E
    k(187) = small + (1d-13)

    !C- + H2 -> CH2 + E
    k(188) = small + (5d-10)

    !C- + O -> CO + E
    k(189) = small + (5d-10)

    !O- + H -> OH + E
    k(190) = small + (7d-10)

    !O- + H2 -> H2O + E
    k(191) = small + (7d-10)

    !O- + C -> CO + E
    k(192) = small + (5d-10)

    !H2 + H+ -> H + H + H+
    k(193) = small + (3d-11*T32**(.5)&
        *exp(-52000d0*invT))

    !H2 + H+ -> H3+
    k(194) = small + (1d-16)

    !C + E -> C-
    k(195) = small + (2.25d-15)

    !C + H -> CH
    k(196) = small + (1d-17)

    !C + H2 -> CH2
    k(197) = small + (1d-17)

    !C + C -> C2
    k(198) = small + (4.36d-18*T32**.35&
        *exp(-161.3d0*invT))

    !C + O -> CO
    k(199) = small + (3.09d-17*T32**.33&
        *exp(-1629d0*invT))

    !C+ + H -> CH+
    k(200) = small + (4.46d-16*Tgas**(-.5)&
        *exp(-4.93*Tgas**(-.6667)))

    !C+ + H2 -> CH2+
    k(201) = small + (2d-16*T32**(-1.3)&
        *exp(-23d0*invTgas))

    !C+ + O -> CO+
    if(Tgas.LT.3d2) then
      k(202) = small + (2.5d-18)
    end if

    !C+ + O -> CO+
    if(Tgas.GE.3d2) then
      k(203) = small + (3.14d-18&
          *T32**(-.15)*exp(-68d0*invT))
    end if

    !O + E -> O-
    k(204) = small + (1.5d-15)

    !O + H -> OH
    k(205) = small + (9.9d-19*T32**(-.38))

    !O + O -> O2
    k(206) = small + (4.9d-20*T32**(1.58))

    !OH + H -> H2O
    k(207) = small + (5.26d-18*T32**(-5.22)&
        *exp(-9d1*invT))

    !CO -> CO_ice
    k(208) = small + (krate_stickSi(n(:),idx_CO,tabTdust))

    !CO_ice -> CO
    k(209) = small + (krate_evaporation(n(:),idx_CO,tabTdust))

    !H2O -> H2O_ice
    k(210) = small + (krate_stickSi(n(:),idx_H2O,tabTdust))

    !H2O_ice -> H2O
    k(211) = small + (krate_evaporation(n(:),idx_H2O,tabTdust))

    !CO_ice -> CO
    k(212) = small + (krate_nonthermal_evaporation(idx_CO,1.78d0, user_Av, user_crate, 1d-3))

    !H- -> H + E
    k(213) = rateEvaluateOnce(213)

    !H2+ -> H + H+
    k(214) = rateEvaluateOnce(214)

    !H3+ -> H2 + H+
    k(215) = rateEvaluateOnce(215)

    !H3+ -> H2+ + H
    k(216) = rateEvaluateOnce(216)

    !C -> C+ + E
    k(217) = rateEvaluateOnce(217)

    !C- -> C + E
    k(218) = rateEvaluateOnce(218)

    !CH -> C + H
    k(219) = rateEvaluateOnce(219)

    !CH -> CH+ + E
    k(220) = rateEvaluateOnce(220)

    !CH+ -> C + H+
    k(221) = rateEvaluateOnce(221)

    !CH2 -> CH + H
    k(222) = rateEvaluateOnce(222)

    !CH2 -> CH2+ + E
    k(223) = rateEvaluateOnce(223)

    !CH2+ -> CH+ + H
    k(224) = rateEvaluateOnce(224)

    !CH3+ -> CH2+ + H
    k(225) = rateEvaluateOnce(225)

    !CH3+ -> CH+ + H2
    k(226) = rateEvaluateOnce(226)

    !C2 -> C + C
    k(227) = rateEvaluateOnce(227)

    !O- -> O + E
    k(228) = rateEvaluateOnce(228)

    !OH -> O + H
    k(229) = rateEvaluateOnce(229)

    !OH -> OH+ + E
    k(230) = rateEvaluateOnce(230)

    !OH+ -> O + H+
    k(231) = rateEvaluateOnce(231)

    !H2O -> OH + H
    k(232) = rateEvaluateOnce(232)

    !H2O -> H2O+ + E
    k(233) = rateEvaluateOnce(233)

    !O2 -> O2+ + E
    k(234) = rateEvaluateOnce(234)

    !O2 -> O + O
    k(235) = rateEvaluateOnce(235)

    !CO -> C + O
    k(236) = rateEvaluateOnce(236)

    !H2 -> H + H
    k(237) = small + (5.6d-11*exp(-3.74*user_Av)&
        *user_H2self)

    !H2O+ -> H2+ + O
    k(238) = small + (5.d-11*HnOj)

    !H2O+ -> H+ + OH
    k(239) = small + (5.d-11*HnOj)

    !H2O+ -> O+ + H2
    k(240) = small + (5.d-11*HnOj)

    !H2O+ -> OH+ + H
    k(241) = small + (1.5d-10*HnOj)

    !H3O+ -> H+ + H2O
    k(242) = small + (2.5d-11*HnOj)

    !H3O+ -> H2+ + OH
    k(243) = small + (2.5d-11*HnOj)

    !H3O+ -> H2O+ + H
    k(244) = small + (7.5d-12*HnOj)

    !H3O+ -> OH+ + H2
    k(245) = small + (2.5d-11*HnOj)

    !H -> H+ + E
    k(246) = rateEvaluateOnce(246)

    !HE -> HE+ + E
    k(247) = rateEvaluateOnce(247)

    !O -> O+ + E
    k(248) = rateEvaluateOnce(248)

    !CO -> C + O
    k(249) = rateEvaluateOnce(249)

    !CO -> CO+ + E
    k(250) = rateEvaluateOnce(250)

    !C2 -> C + C
    k(251) = rateEvaluateOnce(251)

    !H2 -> H + H
    k(252) = rateEvaluateOnce(252)

    !H2 -> H+ + H-
    k(253) = rateEvaluateOnce(253)

    !H2 -> H2+ + E
    k(254) = rateEvaluateOnce(254)

    !C -> C+ + E
    k(255) = rateEvaluateOnce(255)

    !CH -> C + H
    k(256) = rateEvaluateOnce(256)

    !O2 -> O + O
    k(257) = rateEvaluateOnce(257)

    !O2 -> O2+ + E
    k(258) = rateEvaluateOnce(258)

    !OH -> O + H
    k(259) = rateEvaluateOnce(259)

    !CH2 -> CH2+ + E
    k(260) = rateEvaluateOnce(260)

    !H2O -> OH + H
    k(261) = rateEvaluateOnce(261)

    !HCO -> CO + H
    k(262) = rateEvaluateOnce(262)

    !HCO -> HCO+ + E
    k(263) = rateEvaluateOnce(263)

    !H2 -> H + H+ + E
    k(264) = rateEvaluateOnce(264)

    !C + C -> C2
    if(Tgas.LT.5d3) then
      k(265) = small + (5.99d-33&
          *(Tgas&
          /5d3)**(-1.6)*ntot)
    end if

    !C + C -> C2
    if(Tgas.GE.5d3) then
      k(266) = small + (5.99d-33&
          *(Tgas&
          /5d3)**(-0.64)*exp(5255./Tgas)*ntot)
    end if

    !C + O -> CO
    if(Tgas.LT.2d3) then
      k(267) = small + (6.16d-29&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C + O -> CO
    if(Tgas.GE.2d3) then
      k(268) = small + (2.14d-29&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !C+ + O -> CO+
    if(Tgas.LT.2d3) then
      k(269) = small + (6.16d-27&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C+ + O -> CO+
    if(Tgas.GE.2d3) then
      k(270) = small + (2.14d-27&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !C + O+ -> CO+
    if(Tgas.LT.2d3) then
      k(271) = small + (6.16d-27&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C + O+ -> CO+
    if(Tgas.GE.2d3) then
      k(272) = small + (2.14d-27&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !H + O -> OH
    k(273) = small + (4.33d-32*(T32)**(-1)*ntot)

    !OH + H -> H2O
    k(274) = small + (2.56d-31*(T32)**(-2)*ntot)

    !O + O -> O2
    k(275) = small + (9.2d-34*(T32)**(-1)*ntot)

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if
  end function coe

  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*255::fname,line
    integer::ios,i,nunit

    ! Verbatim reactions filename defaults to `reactions_verbatim.dat`
    fname = "reactions_verbatim.dat"

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: "//trim(fname)//" file not present!"
      stop
    end if

    !load reactions from file
    do i=1,nrea
      read(nunit,'(a)',iostat=ios) line
      if(ios/=0) then
        print *,"ERROR: problem reading "//trim(fname)
        stop
      end if
      reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
    B(:) = ref(:)

    !charge conservation
    x(idx_E) = m(idx_E)*(- 1d0*x(idx_Hk) / m(idx_Hk) &
        - 1d0*x(idx_Ck) / m(idx_Ck) &
        - 1d0*x(idx_Ok) / m(idx_Ok) &
        + 1d0*x(idx_Hj) / m(idx_Hj) &
        + 1d0*x(idx_HEj) / m(idx_HEj) &
        + 1d0*x(idx_H2j) / m(idx_H2j) &
        + 1d0*x(idx_Cj) / m(idx_Cj) &
        + 1d0*x(idx_Oj) / m(idx_Oj) &
        + 1d0*x(idx_HOCj) / m(idx_HOCj) &
        + 1d0*x(idx_HCOj) / m(idx_HCOj) &
        + 1d0*x(idx_H3j) / m(idx_H3j) &
        + 1d0*x(idx_CHj) / m(idx_CHj) &
        + 1d0*x(idx_CH2j) / m(idx_CH2j) &
        + 1d0*x(idx_COj) / m(idx_COj) &
        + 1d0*x(idx_CH3j) / m(idx_CH3j) &
        + 1d0*x(idx_OHj) / m(idx_OHj) &
        + 1d0*x(idx_H2Oj) / m(idx_H2Oj) &
        + 1d0*x(idx_H3Oj) / m(idx_H3Oj) &
        + 1d0*x(idx_O2j) / m(idx_O2j) &
        + 2d0*x(idx_HEjj) / m(idx_HEjj))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(33,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  ! k_rev = k_for * revKc
  ! Note that reaction constant revKc is calculated with
  ! reactants and products from reverse reaction
  function revKc(Tgas,ridx,pidx)
    use krome_constants
    use krome_commons
    implicit none
    real*8::revKc,Tgas,dgibss,stoichiometricChange
    integer::ridx(:),pidx(:),i

    ! when considering forward reaction:
    ! Kc = (P°)**(p+p-r-r) * exp(-dGibss_forward°)
    ! where ° means at standard conditions of
    ! P° = 1 bar = (kb*T/1e6) dyn/cm^2 (cgs)
    ! when considering reverse:
    ! 1/Kc = revKc = (kb*T/1e6)**(p+p-r-r) * exp(-dGibss_reverse°)
    ! kb*T/1e6 is to go from 1 atm pressure to number density cm^-3
    ! When not at standard pressure this does not change:
    ! revKc = P**(p+p-r-r) *exp(-dGibss_reverse° - (p+p-r-r)*ln(P/P°))
    !       = (P°)**(p+p-r-r) * exp(-dGibss_reverse°)

    dgibss = 0.d0 ! Gibbs free energy/(R*T)
    stoichiometricChange = 0d0

    do i=1,size(pidx)
      dgibss = dgibss + revHS(Tgas,pidx(i))
      stoichiometricChange = stoichiometricChange + 1
    end do

    do i=1,size(ridx)
      dgibss = dgibss - revHS(Tgas,ridx(i))
      stoichiometricChange = stoichiometricChange - 1
    end do

    revKc = (boltzmann_erg * Tgas * 1e-6)**(-stoichiometricChange)&
        * exp(-dgibss)

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    use krome_constants
    use krome_fit
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::Tnist,Tnist2,Tnist3,Tnist4,invTnist,invTnist2,lnTnist
    real*8::p1_nasa(40,7), p2_nasa(40,7), Tlim_nasa(40,3), p(7)
    real*8::p1_nist(40,7), p2_nist(40,7), Tlim_nist(40,3)
    integer::idx

    p(:) = 0.d0
    p1_nasa(:,:) = 0.d0
    p2_nasa(:,:) = 0.d0
    Tlim_nasa(:,:) = 0.d0
    p1_nist(:,:) = 0.d0
    p2_nist(:,:) = 0.d0
    Tlim_nist(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)
    ! NIST polynomials are quite differernt
    ! it doesn't like easy stuff...
    Tnist = Tgas * 1.d-3
    Tnist2 = Tnist * Tnist
    Tnist3 = Tnist2 * Tnist
    Tnist4 = Tnist3 * Tnist2
    invTnist = 1d0/Tnist
    invTnist2 = invTnist * invTnist
    lnTnist = log(Tnist)

    p1_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p1_nasa(idx_Ck,:)  = (/2.50025151d0,&
        -1.19774349d-06,&
        2.28919443d-09,&
        -1.98276803d-12,&
        6.44398056d-16,&
        70064.893d0,&
        4.87847086d0/)
    p1_nasa(idx_Ok,:)  = (/2.90805921d0,&
        -0.00169804907d0,&
        2.98069955d-06,&
        -2.43835127d-09,&
        7.61229311d-13,&
        11435.7717d0,&
        2.80339097d0/)
    p1_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1_nasa(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1_nasa(idx_C,:)  = (/2.5542395d0,&
        -0.00032153772d0,&
        7.3379223d-07,&
        -7.3223487d-10,&
        2.6652144d-13,&
        85442.681d0,&
        4.5313085d0/)
    p1_nasa(idx_O,:)  = (/3.1682671d0,&
        -0.00327931884d0,&
        6.64306396d-06,&
        -6.12806624d-09,&
        2.11265971d-12,&
        29122.2592d0,&
        2.05193346d0/)
    p1_nasa(idx_OH,:)  = (/3.99198424d0,&
        -0.00240106655d0,&
        4.61664033d-06,&
        -3.87916306d-09,&
        1.36319502d-12,&
        3368.89836d0,&
        -0.103998477d0/)
    p1_nasa(idx_CO,:)  = (/3.5795335d0,&
        -0.00061035369d0,&
        1.0168143d-06,&
        9.0700586d-10,&
        -9.0442449d-13,&
        -14344.086d0,&
        3.5084093d0/)
    p1_nasa(idx_CH,:)  = (/3.4897583d0,&
        0.0003243216d0,&
        -1.6899751d-06,&
        3.162842d-09,&
        -1.4061803d-12,&
        70660.755d0,&
        2.0842841d0/)
    p1_nasa(idx_CH2,:)  = (/3.84261832d0,&
        -7.36676871d-06,&
        6.16970693d-06,&
        -6.96689962d-09,&
        2.64620979d-12,&
        45863.1528d0,&
        1.2758447d0/)
    p1_nasa(idx_HCO,:)  = (/4.36380907d0,&
        -0.00535204137d0,&
        2.31954508d-05,&
        -2.6610904d-08,&
        1.02711962d-11,&
        25010.8717d0,&
        2.98106307d0/)
    p1_nasa(idx_H2O,:)  = (/4.1986352d0,&
        -0.0020364017d0,&
        6.5203416d-06,&
        -5.4879269d-09,&
        1.771968d-12,&
        -30293.726d0,&
        -0.84900901d0/)
    p1_nasa(idx_O2,:)  = (/3.78245636d0,&
        -0.00299673416d0,&
        9.84730201d-06,&
        -9.68129509d-09,&
        3.24372837d-12,&
        -1063.94356d0,&
        3.65767573d0/)
    p1_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p1_nasa(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p1_nasa(idx_Cj,:)  = (/2.61332254d0,&
        -0.000540148065d0,&
        1.03037233d-06,&
        -8.90092552d-10,&
        2.88500586d-13,&
        216862.274d0,&
        3.8345479d0/)
    p1_nasa(idx_Oj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        187935.284d0,&
        4.39337676d0/)
    p1_nasa(idx_H3j,:)  = (/4.1795698d0,&
        -0.000868875627d0,&
        -1.09017371d-07,&
        4.13349766d-09,&
        -2.37877027d-12,&
        132635.537d0,&
        -5.838001d0/)
    p1_nasa(idx_COj,:)  = (/3.77061642d0,&
        -0.00201773246d0,&
        4.61081738d-06,&
        -2.99175463d-09,&
        6.06065045d-13,&
        149006.795d0,&
        3.38129783d0/)
    p1_nasa(idx_OHj,:)  = (/3.50502572d0,&
        0.000241313747d0,&
        -1.42200948d-06,&
        2.64780232d-09,&
        -1.17038711d-12,&
        155210.676d0,&
        1.97907627d0/)
    p1_nasa(idx_H2Oj,:)  = (/4.02465912d0,&
        -0.00108851414d0,&
        5.13576558d-06,&
        -4.40027838d-09,&
        1.40726746d-12,&
        116895.616d0,&
        0.699968812d0/)
    p1_nasa(idx_H3Oj,:)  = (/3.79295251d0,&
        -0.000910852723d0,&
        1.16363521d-05,&
        -1.21364865d-08,&
        4.26159624d-12,&
        71402.7518d0,&
        1.47156927d0/)
    p1_nasa(idx_O2j,:)  = (/4.61017167d0,&
        -0.00635951952d0,&
        1.42425624d-05,&
        -1.20997923d-08,&
        3.70956878d-12,&
        139742.229d0,&
        -0.201326941d0/)
    p2_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p2_nasa(idx_Ck,:)  = (/2.50001597d0,&
        -1.71721376d-08,&
        6.9283294d-12,&
        -1.20607892d-15,&
        7.60308635d-20,&
        70064.9324d0,&
        4.87955907d0/)
    p2_nasa(idx_Ok,:)  = (/2.54474869d0,&
        -4.66695513d-05,&
        1.84912357d-08,&
        -3.18159223d-12,&
        1.98962956d-16,&
        11504.2089d0,&
        4.52131015d0/)
    p2_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2_nasa(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2_nasa(idx_C,:)  = (/2.605583d0,&
        -0.00019593434d0,&
        1.0673722d-07,&
        -1.642394d-11,&
        8.187058d-16,&
        85411.742d0,&
        4.1923868d0/)
    p2_nasa(idx_O,:)  = (/2.54363697d0,&
        -2.73162486d-05,&
        -4.1902952d-09,&
        4.95481845d-12,&
        -4.79553694d-16,&
        29226.012d0,&
        4.92229457d0/)
    p2_nasa(idx_OH,:)  = (/2.83853033d0,&
        0.00110741289d0,&
        -2.94000209d-07,&
        4.20698729d-11,&
        -2.4228989d-15,&
        3697.80808d0,&
        5.84494652d0/)
    p2_nasa(idx_CO,:)  = (/3.0484859d0,&
        0.0013517281d0,&
        -4.8579405d-07,&
        7.8853644d-11,&
        -4.6980746d-15,&
        -14266.117d0,&
        6.0170977d0/)
    p2_nasa(idx_CH,:)  = (/2.5209369d0,&
        0.0017653639d0,&
        -4.614766d-07,&
        5.9289675d-11,&
        -3.3474501d-15,&
        70994.878d0,&
        7.4051829d0/)
    p2_nasa(idx_CH2,:)  = (/3.11049513d0,&
        0.00373779517d0,&
        -1.37371977d-06,&
        2.23054839d-10,&
        -1.33567178d-14,&
        45971.5953d0,&
        4.62796405d0/)
    p2_nasa(idx_HCO,:)  = (/4.23892214d0,&
        0.0019657617d0,&
        -3.82075171d-07,&
        4.80137647d-11,&
        -3.11176347d-15,&
        24726.1645d0,&
        1.99698242d0/)
    p2_nasa(idx_H2O,:)  = (/2.6770389d0,&
        0.0029731816d0,&
        -7.7376889d-07,&
        9.4433514d-11,&
        -4.2689991d-15,&
        -29885.894d0,&
        6.88255d0/)
    p2_nasa(idx_O2,:)  = (/3.66096065d0,&
        0.000656365811d0,&
        -1.41149627d-07,&
        2.05797935d-11,&
        -1.29913436d-15,&
        -1215.97718d0,&
        3.41536279d0/)
    p2_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p2_nasa(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    p2_nasa(idx_Cj,:)  = (/2.50827618d0,&
        -1.04354146d-05,&
        5.16160809d-09,&
        -1.14187475d-12,&
        9.43539946d-17,&
        216879.645d0,&
        4.3188599d0/)
    p2_nasa(idx_Oj,:)  = (/2.48542028d0,&
        2.56978695d-05,&
        -1.28833378d-08,&
        1.65525487d-12,&
        1.09933344d-16,&
        187940.874d0,&
        4.47425446d0/)
    p2_nasa(idx_H3j,:)  = (/2.01435718d0,&
        0.00415925769d0,&
        -1.42664877d-06,&
        2.22372739d-10,&
        -1.29346518d-14,&
        133230.507d0,&
        5.46168967d0/)
    p2_nasa(idx_COj,:)  = (/2.93062935d0,&
        0.00156033262d0,&
        -6.16246355d-07,&
        1.09957336d-10,&
        -6.66119284d-15,&
        149147.222d0,&
        7.3384673d0/)
    p2_nasa(idx_OHj,:)  = (/2.68358996d0,&
        0.00157006435d0,&
        -5.39972815d-07,&
        9.37643877d-11,&
        -5.70068067d-15,&
        155479.296d0,&
        6.44375894d0/)
    p2_nasa(idx_H2Oj,:)  = (/3.31570445d0,&
        0.00210648746d0,&
        -3.76341515d-07,&
        3.47525972d-11,&
        -1.70335643d-15,&
        117017.475d0,&
        4.03220514d0/)
    p2_nasa(idx_H3Oj,:)  = (/2.49647765d0,&
        0.0057284484d0,&
        -1.83953239d-06,&
        2.73577348d-10,&
        -1.54093917d-14,&
        71624.4227d0,&
        7.45850493d0/)
    p2_nasa(idx_O2j,:)  = (/3.31675922d0,&
        0.00111522244d0,&
        -3.83492556d-07,&
        5.72784687d-11,&
        -2.77648381d-15,&
        139876.823d0,&
        5.44726469d0/)
    Tlim_nasa(idx_Hk,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Ck,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Ok,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_C,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HCO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Cj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_COj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OHj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    ! pick NASA data if present for species
    if (Tlim_nasa(idx,2) /= 0.d0) then
      !select set of NASA polynomials using temperature
      if(Tlim_nasa(idx,1).le.Tgas .and. Tgas.le.Tlim_nasa(idx,2)) then
        p(:) = p1_nasa(idx,:)

      else if(Tlim_nasa(idx,2)<Tgas .and. Tgas.le.Tlim_nasa(idx,3)) then
        p(:) = p2_nasa(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NASA polynomials for enthalpy and enthropy (unitless)
      H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
          p(5)*Tgas4*0.2d0 + p(6)*invT
      S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
          p(5)*Tgas4*0.25d0 + p(7)

      revHS = H - S

      ! else pick NIST data (if present)
    else if (Tlim_nist(idx,2) /= 0.d0) then
      if (Tlim_nist(idx,1) < Tgas .and. Tgas < Tlim_nist(idx,2)) then
        p(:) = p1_nist(idx,:)

      else if (Tlim_nist(idx,2) < Tgas .and. Tgas < Tlim_nist(idx,3)) then
        p(:) = p2_nist(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NIST polynomials for enthalpy and enthropy
      ! H in (kJ/mol)
      H = p(1)*Tnist + p(2)*0.5d0*Tnist2 + p(3)*Tnist3/3.d0 + p(4)*Tnist4*0.25d0&
          - p(5)*invTnist + p(6)
      !  Unitsless
      H = H / (Rgas_kJ * Tgas)

      ! S in (J/mol*K)
      S = p(1)*lnTnist + p(2)*Tnist + p(3)*Tnist2*0.5d0 + p(4)*Tnist3/3.d0&
          - p(5)*invTnist2*0.5d0 + p(7)
      !  Unitless. Note: do not use Tnist
      S = S / Rgas_J

      revHS = H - S

      ! return zero is no data exists
    else
      print *, "No thermochemical data of species index", idx
      revHS = 0.d0

    end if

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_r3(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(1:275) = (/5,20,20,6,21,21,21,21,6,6,7,7,7,7,36,5,2,5&
        ,5,22,7,7,7,7,2,2,2,2,2,22,22,22,7,5,5,7,23,23,23,24,24,8,9&
        ,24,9,9,8,23,8,8,8,10,25,25,25,8,12,12,12,12,12,12,13,13,13&
        ,13,13,14,14,9,10,10,7,8,8,9,9,10,16,17,17,17,17,11,22,27,8,8&
        ,8,23,28,28,28,29,29,29,31,31,31,14,24,9,9,9,10,10,10,10,32&
        ,33,16,16,16,16,16,16,34,17,17,17,8,8,11,11,11,11,26,26,26,12&
        ,12,13,13,13,13,13,13,13,13,14,10,10,10,10,16,16,16,16,16,16&
        ,16,16,17,17,17,11,11,30,3,4,21,27,27,28,29,29,29,31,31,31,32&
        ,33,33,33,34,34,34,34,35,30,26,26,25,2,2,2,3,3,3,4,4,4,7,7,8&
        ,8,8,8,8,23,23,23,23,9,9,9,10,40,40,40,40,40,2,22,27,27,8,3&
        ,12,12,28,13,13,29,31,31,14,4,10,10,32,16,16,17,17,11,7,33,33&
        ,33,33,34,34,34,34,5,6,9,11,11,14,7,7,7,8,12,17,17,10,13,16&
        ,15,15,7,8,8,8,8,23,23,8,8,5,10,9/)
    arr_r2(1:275) = (/1,1&
        ,1,1,1,1,1,5,20,20,6,21,21,21,1,1,5,20,20,5,20,20,1,5,1,5,5&
        ,20,20,1,1,2,7,5,5,5,1,1,1,1,1,1,1,5,20,21,20,5,21,21,21,5,7&
        ,11,11,7,5,7,8,9,9,9,5,9,9,9,9,9,9,7,5,5,10,10,10,10,10,10,5&
        ,5,7,8,8,5,7,5,22,27,27,7,5,7,9,5,7,9,5,9,9,24,7,22,27,27,27&
        ,27,23,23,7,7,27,27,23,23,23,23,8,23,23,29,35,35,27,27,27,27&
        ,8,16,16,20,20,20,20,20,20,21,21,21,21,21,20,20,21,21,20,20&
        ,21,21,21,21,21,21,20,21,21,21,21,5,20,20,2,1,1,1,1,1,1,1,1,1&
        ,1,1,1,1,1,1,1,1,1,1,1,1,1,8,9,10,5,7,9,5,7,8,20,20,1,5,7,8,9&
        ,5,7,9,9,1,5,9,5,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,8,8,9,9,9,9,24,24,9,5,9/)
    arr_r3(1:275) = (/40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,6,5,5,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40&
        ,40/)
    arr_p1(1:275) = (/20,5,5,21,6,6,36,6,21,21,5,6,6,21&
        ,21,2,7,22,22,7,22,22,5,5,5,5,5,5,22,5,5,5,7,7,7,7,8,8,8,9,9&
        ,23,24,9,24,24,23,8,23,23,23,9,26,26,26,12,8,13,14,11,26,10&
        ,12,11,11,15,12,11,11,10,9,9,16,5,5,5,5,16,7,10,10,11,11,8,27&
        ,22,28,28,29,28,23,29,30,28,31,26,29,25,26,30,5,5,7,5,7,7,5,5&
        ,33,34,7,7,25,26,26,33,26,30,11,26,9,17,7,7,7,7,11,11,11,28&
        ,28,7,7,5,5,6,6,6,6,23,32,32,24,24,5,5,6,6,6,6,6,6,35,35,24&
        ,23,8,11,8,9,5,7,5,8,12,8,8,13,12,12,9,9,10,9,10,9,5,10,9,8&
        ,11,10,11,12,10,16,12,13,11,10,16,11,5,27,3,12,13,14,11,28,29&
        ,30,30,4,10,17,16,40,40,40,40,40,5,5,7,22,23,8,8,28,8,12,29&
        ,28,29,28,8,9,9,32,9,10,33,35,9,8,5,22,20,24,32,20,22,33,32&
        ,20,21,24,8,30,8,5,20,22,23,8,9,35,9,29,10,11,26,5,14,14,11&
        ,11,30,30,30,30,10,16,17/)
    arr_p2(1:275) = (/1,40,40,1,40&
        ,40,1,20,5,5,5,22,5,5,40,40,1,40,40,20,5,5,5,5,1,5,5,5,1,5,5&
        ,7,5,6,5,7,40,40,40,40,40,1,1,20,5,6,5,20,6,6,6,5,7,11,11,5,7&
        ,5,5,5,1,8,7,5,7,5,10,8,8,5,7,7,5,11,11,17,17,9,10,9,10,9,9&
        ,10,5,7,5,7,5,5,7,5,5,7,5,5,7,7,7,8,32,32,32,33,33,33,30,30,5&
        ,5,34,34,5,5,5,8,7,9,24,10,30,23,26,26,25,25,28,34,34,5,5,28&
        ,28,29,29,7,7,5,5,8,5,5,6,6,33,33,10,10,32,32,33,33,5,6,6,6,6&
        ,20,5,5,6,5,5,5,5,7,5,5,7,5,5,7,5,5,5,5,16,7,9,9,5,8,5,1,1,1&
        ,1,1,1,1,1,1,5,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,1,20,20,5,1,1,5,1,20,5,1,5,5,7,8,1,5,1,20,5,1,1&
        ,9,9,5,9,10,7,5,16,10,5,7,1,1,1,9,1,8,5,2,1,1,5,9,1,5,1,5,5,1&
        ,20,40,40,40,40,40,40,40,40,40,40,40/)
    arr_p3(1:275) = (/1&
        ,40,40,1,40,40,1,40,40,40,6,40,20,5,40,40,40,40,40,40,40,40,1&
        ,5,1,1,1,40,40,40,40,40,5,40,40,40,40,40,40,40,40,1,1,40,40&
        ,40,40,40,40,40,40,5,40,40,40,40,40,40,40,40,40,40,40,5,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,23,23,28,28,6,40,40,5,5,40,40,20&
        ,20,5,5,40,40,40,40,9,9,24,40,40,40,40,40,5,40,40,40,5,40,40&
        ,5,40,40,40,5,5,7,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,20,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,1,40,40,40,40,40,40&
        ,40,40,40,40,40/)

  end subroutine load_arrays

  ! ************************************
  ! solves linear least squares
  subroutine llsq(n, x, y, a, b)

    !****************************************************
    !
    !! LLSQ solves a linear least squares problem matching a line to data.
    !
    !  Discussion:
    !
    !    A formula for a line of the form Y = A * X + B is sought, which
    !    will minimize the root-mean-square error to N data points
    !    ( X(I), Y(I) );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    In: N, the number of data values.
    !
    !    In: X(N), Y(N), the coordinates of the data points.
    !
    !    Out: A, B, the slope and Y-intercept of the
    !    least-squares approximant to the data.
    !
    implicit none
    integer,intent(in)::n
    real*8,intent(out)::a, b
    real*8,intent(in)::x(n), y(n)
    real*8::bot, top, xbar, ybar

    ! special case
    if(n == 1) then
      a = 0d0
      b = y(1)
      return
    end if

    ! average X and Y
    xbar = sum(x) / n
    ybar = sum(y) / n

    ! compute beta
    top = dot_product(x(:) - xbar, y(:) - ybar)
    bot = dot_product(x(:) - xbar, x(:) - xbar)

    ! if top is zero a is zero
    if(top==0d0) then
      a = 0d0
    else
      a = top / bot
    end if

    b = ybar - a * xbar

  end subroutine llsq

end module krome_subs

!############### MODULE ##############
module krome_stars

end module krome_stars

!############### MODULE ##############
module krome_dust
contains

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_fit
    implicit none

    call init_anytab2D("dust_table_cool.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_cool(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_Tdust.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_H2.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_H2(:,:), dust_mult_ngas, &
        dust_mult_Tgas)

  end subroutine init_dust_tabs

end module krome_dust

!############### MODULE ##############
module krome_photo
contains

end module krome_photo

!############### MODULE ##############
module krome_tabs
contains

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_getphys
    use krome_phfuncs
    use krome_grfuncs
    use krome_constants
    use krome_commons
    use krome_user_commons
    implicit none
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small

    Tgas = max(n(idx_Tgas),phys_Tcmb)
    small = 0d0

    coe_tab(:) = coe(n(:))

  end function coe_tab

  !**************************
  !compute rates that remain constant during solver call
  subroutine makeStoreOnceRates(n)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::small

    small = 0d0
    rateEvaluateOnce(:) = 0d0

    !H- -> H + E
    rateEvaluateOnce(213) = small + (7.1d-7&
        *exp(-.5*user_Av))

    !H2+ -> H + H+
    rateEvaluateOnce(214) = small + (1.1d-9&
        *exp(-1.9*user_Av))

    !H3+ -> H2 + H+
    rateEvaluateOnce(215) = small + (4.9d-13&
        *exp(-1.8*user_Av))

    !H3+ -> H2+ + H
    rateEvaluateOnce(216) = small + (4.9d-13&
        *exp(-2.3*user_Av))

    !C -> C+ + E
    rateEvaluateOnce(217) = small + (3.1d-10&
        *exp(-3.*user_Av))

    !C- -> C + E
    rateEvaluateOnce(218) = small + (2.4d-7&
        *exp(-.9*user_Av))

    !CH -> C + H
    rateEvaluateOnce(219) = small + (8.7d-10&
        *exp(-1.2*user_Av))

    !CH -> CH+ + E
    rateEvaluateOnce(220) = small + (7.7d-10&
        *exp(-2.8*user_Av))

    !CH+ -> C + H+
    rateEvaluateOnce(221) = small + (2.6d-10&
        *exp(-2.5*user_Av))

    !CH2 -> CH + H
    rateEvaluateOnce(222) = small + (7.1d-10&
        *exp(-1.7*user_Av))

    !CH2 -> CH2+ + E
    rateEvaluateOnce(223) = small + (5.9d-10&
        *exp(-2.3*user_Av))

    !CH2+ -> CH+ + H
    rateEvaluateOnce(224) = small + (4.6d-10&
        *exp(-1.7*user_Av))

    !CH3+ -> CH2+ + H
    rateEvaluateOnce(225) = small + (1d-9&
        *exp(-1.7*user_Av))

    !CH3+ -> CH+ + H2
    rateEvaluateOnce(226) = small + (1d-9&
        *exp(-1.7*user_Av))

    !C2 -> C + C
    rateEvaluateOnce(227) = small + (1.5d-10&
        *exp(-2.1*user_Av))

    !O- -> O + E
    rateEvaluateOnce(228) = small + (2.4d-7&
        *exp(-.5*user_Av))

    !OH -> O + H
    rateEvaluateOnce(229) = small + (3.7d-10&
        *exp(-1.7*user_Av))

    !OH -> OH+ + E
    rateEvaluateOnce(230) = small + (1.6d-12&
        *exp(-3.1*user_Av))

    !OH+ -> O + H+
    rateEvaluateOnce(231) = small + (1d-12&
        *exp(-1.8*user_Av))

    !H2O -> OH + H
    rateEvaluateOnce(232) = small + (6d-10&
        *exp(-1.7*user_Av))

    !H2O -> H2O+ + E
    rateEvaluateOnce(233) = small + (3.2d-11&
        *exp(-3.9*user_Av))

    !O2 -> O2+ + E
    rateEvaluateOnce(234) = small + (5.6d-11&
        *exp(-3.7*user_Av))

    !O2 -> O + O
    rateEvaluateOnce(235) = small + (7d-10&
        *exp(-1.8*user_Av))

    !CO -> C + O
    rateEvaluateOnce(236) = small + (2.d-10&
        *exp(-3.53*user_Av))

    !H -> H+ + E
    rateEvaluateOnce(246) = small + (4.6d-1&
        *user_crate)

    !HE -> HE+ + E
    rateEvaluateOnce(247) = small + (5.d-1&
        *user_crate)

    !O -> O+ + E
    rateEvaluateOnce(248) = small + (2.8d0&
        *user_crate)

    !CO -> C + O
    rateEvaluateOnce(249) = small + (5d0&
        *user_crate)

    !CO -> CO+ + E
    rateEvaluateOnce(250) = small + (3d0&
        *user_crate)

    !C2 -> C + C
    rateEvaluateOnce(251) = small + (2.37d2&
        *user_crate)

    !H2 -> H + H
    rateEvaluateOnce(252) = small + (1d-1&
        *user_crate)

    !H2 -> H+ + H-
    rateEvaluateOnce(253) = small + (3d-4&
        *user_crate)

    !H2 -> H2+ + E
    rateEvaluateOnce(254) = small + (9.3d-1&
        *user_crate)

    !C -> C+ + E
    rateEvaluateOnce(255) = small + (1.02d3&
        *user_crate)

    !CH -> C + H
    rateEvaluateOnce(256) = small + (7.3d2&
        *user_crate)

    !O2 -> O + O
    rateEvaluateOnce(257) = small + (7.5d2&
        *user_crate)

    !O2 -> O2+ + E
    rateEvaluateOnce(258) = small + (1.17d2&
        *user_crate)

    !OH -> O + H
    rateEvaluateOnce(259) = small + (5.1d2&
        *user_crate)

    !CH2 -> CH2+ + E
    rateEvaluateOnce(260) = small + (5d2&
        *user_crate)

    !H2O -> OH + H
    rateEvaluateOnce(261) = small + (9.7d2&
        *user_crate)

    !HCO -> CO + H
    rateEvaluateOnce(262) = small + (4.21d2&
        *user_crate)

    !HCO -> HCO+ + E
    rateEvaluateOnce(263) = small + (1.17d3&
        *user_crate)

    !H2 -> H + H+ + E
    rateEvaluateOnce(264) = small + (9.3d-1&
        *user_crate)

  end subroutine makeStoreOnceRates

end module krome_tabs

!############### MODULE ##############
module KROME_coolingGH
end module KROME_coolingGH


!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=0
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
contains

  !*******************
  function cooling(n,inTgas)
    use krome_commons
    implicit none
    real*8::n(:),inTgas,cooling,Tgas

    Tgas = inTgas
    cooling = sum(get_cooling_array(n(:),Tgas))

  end function cooling

  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    real*8::f1,f2,smooth

    f1 = 1d0
    f2 = 1d0

    !returns cooling in erg/cm3/s
    cools(:) = 0d0

    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

    get_cooling_array(:) = cools(:)

  end function get_cooling_array

  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom

    cooling_custom = 0d0

  end function cooling_custom

  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    !use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return

    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0

    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))

    kpla = 1d1**(a0*y + a1) !fit density only

  end function kpla

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)

  end subroutine mylin2

  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if

    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if

    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)

  end subroutine mylin3

  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(33,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool

  !***********************************
  !routine to dump cooling in unit nfile
  subroutine dump_cool(n,Tgas,nfile)
    use krome_commons
    implicit none
    real*8::Tgas,n(:),cools(ncools)
    integer::nfile

    cools(:) = get_cooling_array(n(:),Tgas)
    write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

  end subroutine dump_cool

end module KROME_cooling

!############### MODULE ##############
module KROME_heating
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

    f2 = 1.

    heats(idx_heat_chem) = heatingChem(n(:), Tgas, k(:), nH2dust)

    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array

  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom

    heat_custom = 0d0

  end function heat_custom

  !H2 FORMATION HEATING and other exo/endothermic
  ! processes (including H2 on dust) in erg/cm3/s
  !krome builds the heating/cooling term according
  ! to the chemical network employed
  !*******************************
  function heatingChem(n, Tgas, k, nH2dust)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_subs
    use krome_getphys
    implicit none
    real*8::heatingChem, n(:), Tgas,k(:),nH2dust
    real*8::h2heatfac,HChem,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
    dd = get_Hnuclei(n(:))

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = 1d-40/(nmax*nmax*nmax)

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    n2H = n(idx_H) * n(idx_H)

    !H2 + E -> H + H + E (cooling)
    HChem = HChem + k(23) * (-4.48d0*n(idx_H2)*n(idx_E))
    !H2 + H -> H + H + H (cooling)
    HChem = HChem + k(24) * (-4.48d0*n(idx_H2)*n(idx_H))
    !H2 + H2 -> H2 + H + H (cooling)
    HChem = HChem + k(33) * (-4.48d0*n(idx_H2)*n(idx_H2))

    heatingChem = HChem * eV_to_erg  !erg/cm3/s

  end function heatingChem

end module KROME_heating

!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::dnChem_H2O,dnChem_CO
    real*8::rr
    integer::i,r1,r2,r3,p1,p2,p3

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    dnChem_CO = 0d0  &
        -k(54)*n(idx_HOCj)*n(idx_CO)  &
        +k(54)*n(idx_HOCj)*n(idx_CO)  &
        -k(55)*n(idx_HOCj)*n(idx_CO)  &
        +k(55)*n(idx_HOCj)*n(idx_CO)  &
        +k(60)*n(idx_CH)*n(idx_O)  &
        +k(64)*n(idx_CH2)*n(idx_O)  &
        +k(65)*n(idx_CH2)*n(idx_O)  &
        +k(68)*n(idx_C2)*n(idx_O)  &
        +k(69)*n(idx_C2)*n(idx_O)  &
        +k(74)*n(idx_C)*n(idx_OH)  &
        +k(75)*n(idx_C)*n(idx_OH)  &
        +k(82)*n(idx_O2)*n(idx_C)  &
        +k(83)*n(idx_O2)*n(idx_C)  &
        -k(84)*n(idx_CO)*n(idx_H)  &
        +k(119)*n(idx_O2)*n(idx_Cj)  &
        -k(123)*n(idx_CO)*n(idx_H3j)  &
        -k(124)*n(idx_CO)*n(idx_H3j)  &
        -k(125)*n(idx_CO)*n(idx_H3j)  &
        -k(126)*n(idx_CO)*n(idx_H3j)  &
        +k(127)*n(idx_HCOj)*n(idx_C)  &
        +k(128)*n(idx_HCOj)*n(idx_H2O)  &
        +k(129)*n(idx_HCOj)*n(idx_H2O)  &
        -k(156)*n(idx_CO)*n(idx_HEj)  &
        -k(157)*n(idx_CO)*n(idx_HEj)  &
        +k(158)*n(idx_COj)*n(idx_H)  &
        +k(181)*n(idx_HCOj)*n(idx_E)  &
        +k(183)*n(idx_HOCj)*n(idx_E)  &
        +k(189)*n(idx_Ck)*n(idx_O)  &
        +k(192)*n(idx_Ok)*n(idx_C)  &
        +k(199)*n(idx_C)*n(idx_O)  &
        -k(236)*n(idx_CO)  &
        -k(249)*n(idx_CO)  &
        -k(250)*n(idx_CO)  &
        +k(262)*n(idx_HCO)  &
        +k(267)*n(idx_C)*n(idx_O)  &
        +k(268)*n(idx_C)*n(idx_O)
    dnChem_H2O = 0d0  &
        +k(73)*n(idx_H2)*n(idx_OH)  &
        +k(78)*n(idx_OH)*n(idx_OH)  &
        -k(79)*n(idx_H2O)*n(idx_H)  &
        -k(111)*n(idx_H2O)*n(idx_H3j)  &
        -k(112)*n(idx_H2O)*n(idx_H3j)  &
        -k(113)*n(idx_H2O)*n(idx_Cj)  &
        -k(114)*n(idx_H2O)*n(idx_Cj)  &
        -k(115)*n(idx_H2O)*n(idx_Cj)  &
        -k(116)*n(idx_H2O)*n(idx_Cj)  &
        -k(128)*n(idx_HCOj)*n(idx_H2O)  &
        -k(129)*n(idx_HCOj)*n(idx_H2O)  &
        -k(145)*n(idx_H2O)*n(idx_Hj)  &
        -k(146)*n(idx_H2O)*n(idx_Hj)  &
        -k(147)*n(idx_H2O)*n(idx_HEj)  &
        -k(148)*n(idx_H2O)*n(idx_HEj)  &
        -k(149)*n(idx_H2O)*n(idx_HEj)  &
        -k(150)*n(idx_H2O)*n(idx_HEj)  &
        -k(151)*n(idx_H2O)*n(idx_HEj)  &
        -k(152)*n(idx_H2O)*n(idx_HEj)  &
        +k(177)*n(idx_H3Oj)*n(idx_E)  &
        +k(186)*n(idx_Hk)*n(idx_OH)  &
        +k(191)*n(idx_Ok)*n(idx_H2)  &
        +k(207)*n(idx_OH)*n(idx_H)  &
        -k(232)*n(idx_H2O)  &
        -k(233)*n(idx_H2O)  &
        +k(242)*n(idx_H3Oj)  &
        -k(261)*n(idx_H2O)  &
        +k(274)*n(idx_OH)*n(idx_H)

    !E
    !E
    dn(idx_E) = &
        -k(1)*n(idx_H)*n(idx_E) &
        +2.d0*k(1)*n(idx_H)*n(idx_E) &
        -k(2)*n(idx_Hj)*n(idx_E) &
        -k(3)*n(idx_Hj)*n(idx_E) &
        -k(4)*n(idx_HE)*n(idx_E) &
        +2.d0*k(4)*n(idx_HE)*n(idx_E) &
        -k(5)*n(idx_HEj)*n(idx_E) &
        -k(6)*n(idx_HEj)*n(idx_E) &
        -k(7)*n(idx_HEj)*n(idx_E) &
        +2.d0*k(7)*n(idx_HEj)*n(idx_E) &
        -k(15)*n(idx_HEjj)*n(idx_E) &
        -k(16)*n(idx_H)*n(idx_E) &
        +k(17)*n(idx_Hk)*n(idx_H) &
        -k(23)*n(idx_H2)*n(idx_E) &
        +k(23)*n(idx_H2)*n(idx_E) &
        -k(25)*n(idx_Hk)*n(idx_E) &
        +2.d0*k(25)*n(idx_Hk)*n(idx_E) &
        +k(26)*n(idx_Hk)*n(idx_H) &
        +k(27)*n(idx_Hk)*n(idx_H) &
        +k(29)*n(idx_Hk)*n(idx_Hj) &
        -k(30)*n(idx_H2j)*n(idx_E) &
        -k(31)*n(idx_H2j)*n(idx_E) &
        -k(37)*n(idx_Cj)*n(idx_E) &
        -k(38)*n(idx_Cj)*n(idx_E) &
        -k(39)*n(idx_Cj)*n(idx_E) &
        -k(40)*n(idx_Oj)*n(idx_E) &
        -k(41)*n(idx_Oj)*n(idx_E) &
        -k(42)*n(idx_C)*n(idx_E) &
        +2.d0*k(42)*n(idx_C)*n(idx_E) &
        -k(43)*n(idx_O)*n(idx_E) &
        +2.d0*k(43)*n(idx_O)*n(idx_E) &
        +k(61)*n(idx_CH)*n(idx_O) &
        -k(162)*n(idx_H3j)*n(idx_E) &
        -k(163)*n(idx_H3j)*n(idx_E) &
        -k(164)*n(idx_CHj)*n(idx_E) &
        -k(165)*n(idx_CH2j)*n(idx_E) &
        -k(166)*n(idx_CH2j)*n(idx_E) &
        -k(167)*n(idx_CH2j)*n(idx_E) &
        -k(168)*n(idx_CH3j)*n(idx_E) &
        -k(169)*n(idx_CH3j)*n(idx_E) &
        -k(170)*n(idx_CH3j)*n(idx_E) &
        -k(171)*n(idx_OHj)*n(idx_E) &
        -k(172)*n(idx_H2Oj)*n(idx_E) &
        -k(173)*n(idx_H2Oj)*n(idx_E) &
        -k(174)*n(idx_H2Oj)*n(idx_E) &
        -k(175)*n(idx_H3Oj)*n(idx_E) &
        -k(176)*n(idx_H3Oj)*n(idx_E) &
        -k(177)*n(idx_H3Oj)*n(idx_E) &
        -k(178)*n(idx_H3Oj)*n(idx_E) &
        -k(179)*n(idx_O2j)*n(idx_E) &
        -k(180)*n(idx_COj)*n(idx_E) &
        -k(181)*n(idx_HCOj)*n(idx_E) &
        -k(182)*n(idx_HCOj)*n(idx_E) &
        -k(183)*n(idx_HOCj)*n(idx_E) &
        +k(184)*n(idx_Hk)*n(idx_C) &
        +k(185)*n(idx_Hk)*n(idx_O) &
        +k(186)*n(idx_Hk)*n(idx_OH) &
        +k(187)*n(idx_Ck)*n(idx_H) &
        +k(188)*n(idx_Ck)*n(idx_H2) &
        +k(189)*n(idx_Ck)*n(idx_O) &
        +k(190)*n(idx_Ok)*n(idx_H) &
        +k(191)*n(idx_Ok)*n(idx_H2) &
        +k(192)*n(idx_Ok)*n(idx_C) &
        -k(195)*n(idx_C)*n(idx_E) &
        -k(204)*n(idx_O)*n(idx_E) &
        +k(213)*n(idx_Hk) &
        +k(217)*n(idx_C) &
        +k(218)*n(idx_Ck) &
        +k(220)*n(idx_CH) &
        +k(223)*n(idx_CH2) &
        +k(228)*n(idx_Ok) &
        +k(230)*n(idx_OH) &
        +k(233)*n(idx_H2O) &
        +k(234)*n(idx_O2) &
        +k(246)*n(idx_H) &
        +k(247)*n(idx_HE) &
        +k(248)*n(idx_O) &
        +k(250)*n(idx_CO) &
        +k(254)*n(idx_H2) &
        +k(255)*n(idx_C) &
        +k(258)*n(idx_O2) &
        +k(260)*n(idx_CH2) &
        +k(263)*n(idx_HCO) &
        +k(264)*n(idx_H2)

    !H-
    !H-
    dn(idx_Hk) = &
        +k(16)*n(idx_H)*n(idx_E) &
        -k(17)*n(idx_Hk)*n(idx_H) &
        -k(25)*n(idx_Hk)*n(idx_E) &
        -k(26)*n(idx_Hk)*n(idx_H) &
        -k(27)*n(idx_Hk)*n(idx_H) &
        -k(28)*n(idx_Hk)*n(idx_Hj) &
        -k(29)*n(idx_Hk)*n(idx_Hj) &
        -k(32)*n(idx_H2j)*n(idx_Hk) &
        -k(161)*n(idx_HEj)*n(idx_Hk) &
        -k(184)*n(idx_Hk)*n(idx_C) &
        -k(185)*n(idx_Hk)*n(idx_O) &
        -k(186)*n(idx_Hk)*n(idx_OH) &
        -k(213)*n(idx_Hk) &
        +k(253)*n(idx_H2)

    !C-
    !C-
    dn(idx_Ck) = &
        -k(159)*n(idx_Ck)*n(idx_Hj) &
        -k(187)*n(idx_Ck)*n(idx_H) &
        -k(188)*n(idx_Ck)*n(idx_H2) &
        -k(189)*n(idx_Ck)*n(idx_O) &
        +k(195)*n(idx_C)*n(idx_E) &
        -k(218)*n(idx_Ck)

    !O-
    !O-
    dn(idx_Ok) = &
        -k(160)*n(idx_Ok)*n(idx_Hj) &
        -k(190)*n(idx_Ok)*n(idx_H) &
        -k(191)*n(idx_Ok)*n(idx_H2) &
        -k(192)*n(idx_Ok)*n(idx_C) &
        +k(204)*n(idx_O)*n(idx_E) &
        -k(228)*n(idx_Ok)

    !H
    !H
    dn(idx_H) = &
        -k(1)*n(idx_H)*n(idx_E) &
        +k(2)*n(idx_Hj)*n(idx_E) &
        +k(3)*n(idx_Hj)*n(idx_E) &
        -k(8)*n(idx_HEj)*n(idx_H) &
        +k(9)*n(idx_HE)*n(idx_Hj) &
        +k(10)*n(idx_HE)*n(idx_Hj) &
        +2.d0*k(11)*n(idx_H2)*n(idx_HE) &
        +k(13)*n(idx_H2)*n(idx_HEj) &
        +2.d0*k(14)*n(idx_H2)*n(idx_HEj) &
        -k(16)*n(idx_H)*n(idx_E) &
        -k(17)*n(idx_Hk)*n(idx_H) &
        -k(18)*n(idx_H)*n(idx_Hj) &
        -k(19)*n(idx_H)*n(idx_Hj) &
        -k(20)*n(idx_H2j)*n(idx_H) &
        +k(21)*n(idx_H2)*n(idx_Hj) &
        +k(22)*n(idx_H2)*n(idx_Hj) &
        +2.d0*k(23)*n(idx_H2)*n(idx_E) &
        -k(24)*n(idx_H2)*n(idx_H) &
        +3.d0*k(24)*n(idx_H2)*n(idx_H) &
        +k(25)*n(idx_Hk)*n(idx_E) &
        -k(26)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(26)*n(idx_Hk)*n(idx_H) &
        -k(27)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(27)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(28)*n(idx_Hk)*n(idx_Hj) &
        +2.d0*k(30)*n(idx_H2j)*n(idx_E) &
        +2.d0*k(31)*n(idx_H2j)*n(idx_E) &
        +k(32)*n(idx_H2j)*n(idx_Hk) &
        +2.d0*k(33)*n(idx_H2)*n(idx_H2) &
        -2.d0*k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
        -3.d0*k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
        +k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
        -2.d0*k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
        -k(44)*n(idx_Oj)*n(idx_H) &
        +k(45)*n(idx_O)*n(idx_Hj) &
        +k(47)*n(idx_C)*n(idx_Hj) &
        -k(48)*n(idx_Cj)*n(idx_H) &
        -k(52)*n(idx_OH)*n(idx_H) &
        +2.d0*k(52)*n(idx_OH)*n(idx_H) &
        +k(56)*n(idx_C)*n(idx_H2) &
        -k(57)*n(idx_CH)*n(idx_H) &
        +k(58)*n(idx_CH)*n(idx_H2) &
        +k(59)*n(idx_CH)*n(idx_C) &
        +k(60)*n(idx_CH)*n(idx_O) &
        -k(63)*n(idx_CH2)*n(idx_H) &
        +2.d0*k(64)*n(idx_CH2)*n(idx_O) &
        +k(66)*n(idx_CH2)*n(idx_O) &
        +k(70)*n(idx_O)*n(idx_H2) &
        -k(71)*n(idx_OH)*n(idx_H) &
        -k(72)*n(idx_OH)*n(idx_H) &
        +k(73)*n(idx_H2)*n(idx_OH) &
        +k(74)*n(idx_C)*n(idx_OH) &
        +k(75)*n(idx_C)*n(idx_OH) &
        +k(76)*n(idx_O)*n(idx_OH) &
        +k(77)*n(idx_O)*n(idx_OH) &
        -k(79)*n(idx_H2O)*n(idx_H) &
        -k(80)*n(idx_O2)*n(idx_H) &
        -k(84)*n(idx_CO)*n(idx_H) &
        +k(85)*n(idx_H2j)*n(idx_H2) &
        -k(86)*n(idx_H3j)*n(idx_H) &
        +k(87)*n(idx_C)*n(idx_H2j) &
        +k(89)*n(idx_C)*n(idx_H3j) &
        +k(90)*n(idx_Cj)*n(idx_H2) &
        -k(91)*n(idx_CHj)*n(idx_H) &
        +k(92)*n(idx_CHj)*n(idx_H2) &
        +k(93)*n(idx_CHj)*n(idx_O) &
        -k(94)*n(idx_CH2j)*n(idx_H) &
        +k(95)*n(idx_CH2j)*n(idx_H2) &
        +k(96)*n(idx_CH2j)*n(idx_O) &
        -k(97)*n(idx_CH3j)*n(idx_H) &
        +k(101)*n(idx_Oj)*n(idx_H2) &
        +k(102)*n(idx_O)*n(idx_H2j) &
        +k(104)*n(idx_O)*n(idx_H3j) &
        +k(107)*n(idx_OH)*n(idx_Cj) &
        +k(108)*n(idx_OH)*n(idx_Cj) &
        +k(109)*n(idx_OHj)*n(idx_H2) &
        +k(110)*n(idx_H2Oj)*n(idx_H2) &
        +k(113)*n(idx_H2O)*n(idx_Cj) &
        +k(114)*n(idx_H2O)*n(idx_Cj) &
        +k(115)*n(idx_H2O)*n(idx_Cj) &
        +k(130)*n(idx_CH)*n(idx_Hj) &
        +k(131)*n(idx_CH)*n(idx_Hj) &
        +k(134)*n(idx_CH2)*n(idx_Hj) &
        +k(135)*n(idx_CH2)*n(idx_Hj) &
        +k(138)*n(idx_CH2)*n(idx_HEj) &
        +k(139)*n(idx_CH2)*n(idx_HEj) &
        +k(141)*n(idx_OH)*n(idx_Hj) &
        +k(142)*n(idx_OH)*n(idx_Hj) &
        +k(143)*n(idx_OH)*n(idx_HEj) &
        +k(144)*n(idx_OH)*n(idx_HEj) &
        +k(145)*n(idx_H2O)*n(idx_Hj) &
        +k(146)*n(idx_H2O)*n(idx_Hj) &
        +k(149)*n(idx_H2O)*n(idx_HEj) &
        +k(150)*n(idx_H2O)*n(idx_HEj) &
        +k(153)*n(idx_O2)*n(idx_Hj) &
        -k(158)*n(idx_COj)*n(idx_H) &
        +k(159)*n(idx_Ck)*n(idx_Hj) &
        +k(160)*n(idx_Ok)*n(idx_Hj) &
        +k(161)*n(idx_HEj)*n(idx_Hk) &
        +k(162)*n(idx_H3j)*n(idx_E) &
        +3.d0*k(163)*n(idx_H3j)*n(idx_E) &
        +k(164)*n(idx_CHj)*n(idx_E) &
        +k(165)*n(idx_CH2j)*n(idx_E) &
        +2.d0*k(167)*n(idx_CH2j)*n(idx_E) &
        +k(168)*n(idx_CH3j)*n(idx_E) &
        +2.d0*k(170)*n(idx_CH3j)*n(idx_E) &
        +k(171)*n(idx_OHj)*n(idx_E) &
        +k(173)*n(idx_H2Oj)*n(idx_E) &
        +2.d0*k(174)*n(idx_H2Oj)*n(idx_E) &
        +2.d0*k(175)*n(idx_H3Oj)*n(idx_E) &
        +k(176)*n(idx_H3Oj)*n(idx_E) &
        +k(177)*n(idx_H3Oj)*n(idx_E) &
        +k(181)*n(idx_HCOj)*n(idx_E) &
        +k(183)*n(idx_HOCj)*n(idx_E) &
        -k(187)*n(idx_Ck)*n(idx_H) &
        -k(190)*n(idx_Ok)*n(idx_H) &
        +2.d0*k(193)*n(idx_H2)*n(idx_Hj) &
        -k(196)*n(idx_C)*n(idx_H) &
        -k(200)*n(idx_Cj)*n(idx_H) &
        -k(205)*n(idx_O)*n(idx_H) &
        -k(207)*n(idx_OH)*n(idx_H) &
        +k(213)*n(idx_Hk) &
        +k(214)*n(idx_H2j) &
        +k(216)*n(idx_H3j) &
        +k(219)*n(idx_CH) &
        +k(222)*n(idx_CH2) &
        +k(224)*n(idx_CH2j) &
        +k(225)*n(idx_CH3j) &
        +k(229)*n(idx_OH) &
        +k(232)*n(idx_H2O) &
        +2.d0*k(237)*n(idx_H2) &
        +k(241)*n(idx_H2Oj) &
        +k(244)*n(idx_H3Oj) &
        -k(246)*n(idx_H) &
        +2.d0*k(252)*n(idx_H2) &
        +k(256)*n(idx_CH) &
        +k(259)*n(idx_OH) &
        +k(261)*n(idx_H2O) &
        +k(262)*n(idx_HCO) &
        +k(264)*n(idx_H2) &
        -k(273)*n(idx_H)*n(idx_O) &
        -k(274)*n(idx_OH)*n(idx_H)

    !HE
    !HE
    dn(idx_HE) = &
        -k(4)*n(idx_HE)*n(idx_E) &
        +k(5)*n(idx_HEj)*n(idx_E) &
        +k(6)*n(idx_HEj)*n(idx_E) &
        +k(8)*n(idx_HEj)*n(idx_H) &
        -k(9)*n(idx_HE)*n(idx_Hj) &
        -k(10)*n(idx_HE)*n(idx_Hj) &
        -k(11)*n(idx_H2)*n(idx_HE) &
        +k(11)*n(idx_H2)*n(idx_HE) &
        +k(12)*n(idx_H2)*n(idx_HEj) &
        +k(13)*n(idx_H2)*n(idx_HEj) &
        -k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
        +k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
        +k(46)*n(idx_O)*n(idx_HEj) &
        +k(49)*n(idx_C)*n(idx_HEj) &
        +k(50)*n(idx_C)*n(idx_HEj) &
        +k(51)*n(idx_C)*n(idx_HEj) &
        +k(136)*n(idx_CH2)*n(idx_HEj) &
        +k(137)*n(idx_CH2)*n(idx_HEj) &
        +k(138)*n(idx_CH2)*n(idx_HEj) &
        +k(139)*n(idx_CH2)*n(idx_HEj) &
        +k(140)*n(idx_C2)*n(idx_HEj) &
        +k(143)*n(idx_OH)*n(idx_HEj) &
        +k(144)*n(idx_OH)*n(idx_HEj) &
        +k(147)*n(idx_H2O)*n(idx_HEj) &
        +k(148)*n(idx_H2O)*n(idx_HEj) &
        +k(149)*n(idx_H2O)*n(idx_HEj) &
        +k(150)*n(idx_H2O)*n(idx_HEj) &
        +k(151)*n(idx_H2O)*n(idx_HEj) &
        +k(152)*n(idx_H2O)*n(idx_HEj) &
        +k(154)*n(idx_O2)*n(idx_HEj) &
        +k(155)*n(idx_O2)*n(idx_HEj) &
        +k(156)*n(idx_CO)*n(idx_HEj) &
        +k(157)*n(idx_CO)*n(idx_HEj) &
        +k(161)*n(idx_HEj)*n(idx_Hk) &
        -k(247)*n(idx_HE)

    !H2
    !H2
    dn(idx_H2) = &
        -k(11)*n(idx_H2)*n(idx_HE) &
        -k(12)*n(idx_H2)*n(idx_HEj) &
        -k(13)*n(idx_H2)*n(idx_HEj) &
        -k(14)*n(idx_H2)*n(idx_HEj) &
        +k(17)*n(idx_Hk)*n(idx_H) &
        +k(20)*n(idx_H2j)*n(idx_H) &
        -k(21)*n(idx_H2)*n(idx_Hj) &
        -k(22)*n(idx_H2)*n(idx_Hj) &
        -k(23)*n(idx_H2)*n(idx_E) &
        -k(24)*n(idx_H2)*n(idx_H) &
        +k(32)*n(idx_H2j)*n(idx_Hk) &
        -2.d0*k(33)*n(idx_H2)*n(idx_H2) &
        +k(33)*n(idx_H2)*n(idx_H2) &
        +k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
        +k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
        -k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
        +2.d0*k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
        -k(53)*n(idx_HOCj)*n(idx_H2) &
        +k(53)*n(idx_HOCj)*n(idx_H2) &
        -k(56)*n(idx_C)*n(idx_H2) &
        +k(57)*n(idx_CH)*n(idx_H) &
        -k(58)*n(idx_CH)*n(idx_H2) &
        +k(63)*n(idx_CH2)*n(idx_H) &
        +k(65)*n(idx_CH2)*n(idx_O) &
        -k(70)*n(idx_O)*n(idx_H2) &
        +k(71)*n(idx_OH)*n(idx_H) &
        +k(72)*n(idx_OH)*n(idx_H) &
        -k(73)*n(idx_H2)*n(idx_OH) &
        +k(79)*n(idx_H2O)*n(idx_H) &
        -k(81)*n(idx_O2)*n(idx_H2) &
        -k(85)*n(idx_H2j)*n(idx_H2) &
        +k(86)*n(idx_H3j)*n(idx_H) &
        +k(88)*n(idx_C)*n(idx_H3j) &
        -k(90)*n(idx_Cj)*n(idx_H2) &
        +k(91)*n(idx_CHj)*n(idx_H) &
        -k(92)*n(idx_CHj)*n(idx_H2) &
        +k(94)*n(idx_CH2j)*n(idx_H) &
        -k(95)*n(idx_CH2j)*n(idx_H2) &
        +k(97)*n(idx_CH3j)*n(idx_H) &
        +k(98)*n(idx_CH3j)*n(idx_O) &
        +k(99)*n(idx_CH3j)*n(idx_O) &
        -k(101)*n(idx_Oj)*n(idx_H2) &
        +k(103)*n(idx_O)*n(idx_H3j) &
        +k(105)*n(idx_OH)*n(idx_H3j) &
        +k(106)*n(idx_OH)*n(idx_H3j) &
        -k(109)*n(idx_OHj)*n(idx_H2) &
        -k(110)*n(idx_H2Oj)*n(idx_H2) &
        +k(111)*n(idx_H2O)*n(idx_H3j) &
        +k(112)*n(idx_H2O)*n(idx_H3j) &
        +k(117)*n(idx_H3Oj)*n(idx_C) &
        +k(123)*n(idx_CO)*n(idx_H3j) &
        +k(124)*n(idx_CO)*n(idx_H3j) &
        +k(125)*n(idx_CO)*n(idx_H3j) &
        +k(126)*n(idx_CO)*n(idx_H3j) &
        +k(132)*n(idx_CH2)*n(idx_Hj) &
        +k(133)*n(idx_CH2)*n(idx_Hj) &
        +k(136)*n(idx_CH2)*n(idx_HEj) &
        +k(137)*n(idx_CH2)*n(idx_HEj) &
        +k(162)*n(idx_H3j)*n(idx_E) &
        +k(166)*n(idx_CH2j)*n(idx_E) &
        +k(169)*n(idx_CH3j)*n(idx_E) &
        +k(172)*n(idx_H2Oj)*n(idx_E) &
        +k(176)*n(idx_H3Oj)*n(idx_E) &
        +k(178)*n(idx_H3Oj)*n(idx_E) &
        -k(188)*n(idx_Ck)*n(idx_H2) &
        -k(191)*n(idx_Ok)*n(idx_H2) &
        -k(193)*n(idx_H2)*n(idx_Hj) &
        -k(194)*n(idx_H2)*n(idx_Hj) &
        -k(197)*n(idx_C)*n(idx_H2) &
        -k(201)*n(idx_Cj)*n(idx_H2) &
        +k(215)*n(idx_H3j) &
        +k(226)*n(idx_CH3j) &
        -k(237)*n(idx_H2) &
        +k(240)*n(idx_H2Oj) &
        +k(245)*n(idx_H3Oj) &
        -k(252)*n(idx_H2) &
        -k(253)*n(idx_H2) &
        -k(254)*n(idx_H2) &
        -k(264)*n(idx_H2)

    !C
    !C
    dn(idx_C) = &
        +k(37)*n(idx_Cj)*n(idx_E) &
        +k(38)*n(idx_Cj)*n(idx_E) &
        +k(39)*n(idx_Cj)*n(idx_E) &
        -k(42)*n(idx_C)*n(idx_E) &
        -k(47)*n(idx_C)*n(idx_Hj) &
        +k(48)*n(idx_Cj)*n(idx_H) &
        -k(49)*n(idx_C)*n(idx_HEj) &
        -k(50)*n(idx_C)*n(idx_HEj) &
        -k(51)*n(idx_C)*n(idx_HEj) &
        -k(56)*n(idx_C)*n(idx_H2) &
        +k(57)*n(idx_CH)*n(idx_H) &
        -k(59)*n(idx_CH)*n(idx_C) &
        +k(62)*n(idx_CH)*n(idx_O) &
        +k(68)*n(idx_C2)*n(idx_O) &
        +k(69)*n(idx_C2)*n(idx_O) &
        -k(74)*n(idx_C)*n(idx_OH) &
        -k(75)*n(idx_C)*n(idx_OH) &
        -k(82)*n(idx_O2)*n(idx_C) &
        -k(83)*n(idx_O2)*n(idx_C) &
        +k(84)*n(idx_CO)*n(idx_H) &
        -k(87)*n(idx_C)*n(idx_H2j) &
        -k(88)*n(idx_C)*n(idx_H3j) &
        -k(89)*n(idx_C)*n(idx_H3j) &
        +k(100)*n(idx_C2)*n(idx_Oj) &
        +k(116)*n(idx_H2O)*n(idx_Cj) &
        -k(117)*n(idx_H3Oj)*n(idx_C) &
        -k(121)*n(idx_C)*n(idx_O2j) &
        -k(122)*n(idx_C)*n(idx_O2j) &
        -k(127)*n(idx_HCOj)*n(idx_C) &
        +k(140)*n(idx_C2)*n(idx_HEj) &
        +k(157)*n(idx_CO)*n(idx_HEj) &
        +k(159)*n(idx_Ck)*n(idx_Hj) &
        +k(164)*n(idx_CHj)*n(idx_E) &
        +k(166)*n(idx_CH2j)*n(idx_E) &
        +k(167)*n(idx_CH2j)*n(idx_E) &
        +k(180)*n(idx_COj)*n(idx_E) &
        +k(182)*n(idx_HCOj)*n(idx_E) &
        -k(184)*n(idx_Hk)*n(idx_C) &
        -k(192)*n(idx_Ok)*n(idx_C) &
        -k(195)*n(idx_C)*n(idx_E) &
        -k(196)*n(idx_C)*n(idx_H) &
        -k(197)*n(idx_C)*n(idx_H2) &
        -2.d0*k(198)*n(idx_C)*n(idx_C) &
        -k(199)*n(idx_C)*n(idx_O) &
        -k(217)*n(idx_C) &
        +k(218)*n(idx_Ck) &
        +k(219)*n(idx_CH) &
        +k(221)*n(idx_CHj) &
        +2.d0*k(227)*n(idx_C2) &
        +k(236)*n(idx_CO) &
        +k(249)*n(idx_CO) &
        +2.d0*k(251)*n(idx_C2) &
        -k(255)*n(idx_C) &
        +k(256)*n(idx_CH) &
        -2.d0*k(265)*n(idx_C)*n(idx_C) &
        -2.d0*k(266)*n(idx_C)*n(idx_C) &
        -k(267)*n(idx_C)*n(idx_O) &
        -k(268)*n(idx_C)*n(idx_O) &
        -k(271)*n(idx_C)*n(idx_Oj) &
        -k(272)*n(idx_C)*n(idx_Oj)

    !O
    !O
    dn(idx_O) = &
        +k(40)*n(idx_Oj)*n(idx_E) &
        +k(41)*n(idx_Oj)*n(idx_E) &
        -k(43)*n(idx_O)*n(idx_E) &
        +k(44)*n(idx_Oj)*n(idx_H) &
        -k(45)*n(idx_O)*n(idx_Hj) &
        -k(46)*n(idx_O)*n(idx_HEj) &
        +k(52)*n(idx_OH)*n(idx_H) &
        -k(60)*n(idx_CH)*n(idx_O) &
        -k(61)*n(idx_CH)*n(idx_O) &
        -k(62)*n(idx_CH)*n(idx_O) &
        -k(64)*n(idx_CH2)*n(idx_O) &
        -k(65)*n(idx_CH2)*n(idx_O) &
        -k(66)*n(idx_CH2)*n(idx_O) &
        -k(67)*n(idx_CH2)*n(idx_O) &
        -k(68)*n(idx_C2)*n(idx_O) &
        -k(69)*n(idx_C2)*n(idx_O) &
        -k(70)*n(idx_O)*n(idx_H2) &
        +k(71)*n(idx_OH)*n(idx_H) &
        +k(72)*n(idx_OH)*n(idx_H) &
        -k(76)*n(idx_O)*n(idx_OH) &
        -k(77)*n(idx_O)*n(idx_OH) &
        +k(78)*n(idx_OH)*n(idx_OH) &
        +k(80)*n(idx_O2)*n(idx_H) &
        +k(82)*n(idx_O2)*n(idx_C) &
        +k(83)*n(idx_O2)*n(idx_C) &
        -k(93)*n(idx_CHj)*n(idx_O) &
        -k(96)*n(idx_CH2j)*n(idx_O) &
        -k(98)*n(idx_CH3j)*n(idx_O) &
        -k(99)*n(idx_CH3j)*n(idx_O) &
        -k(102)*n(idx_O)*n(idx_H2j) &
        -k(103)*n(idx_O)*n(idx_H3j) &
        -k(104)*n(idx_O)*n(idx_H3j) &
        +k(118)*n(idx_O2)*n(idx_Cj) &
        +k(121)*n(idx_C)*n(idx_O2j) &
        +k(155)*n(idx_O2)*n(idx_HEj) &
        +k(156)*n(idx_CO)*n(idx_HEj) &
        +k(160)*n(idx_Ok)*n(idx_Hj) &
        +k(171)*n(idx_OHj)*n(idx_E) &
        +k(172)*n(idx_H2Oj)*n(idx_E) &
        +k(174)*n(idx_H2Oj)*n(idx_E) &
        +k(176)*n(idx_H3Oj)*n(idx_E) &
        +2.d0*k(179)*n(idx_O2j)*n(idx_E) &
        +k(180)*n(idx_COj)*n(idx_E) &
        -k(185)*n(idx_Hk)*n(idx_O) &
        -k(189)*n(idx_Ck)*n(idx_O) &
        -k(199)*n(idx_C)*n(idx_O) &
        -k(202)*n(idx_Cj)*n(idx_O) &
        -k(203)*n(idx_Cj)*n(idx_O) &
        -k(204)*n(idx_O)*n(idx_E) &
        -k(205)*n(idx_O)*n(idx_H) &
        -2.d0*k(206)*n(idx_O)*n(idx_O) &
        +k(228)*n(idx_Ok) &
        +k(229)*n(idx_OH) &
        +k(231)*n(idx_OHj) &
        +2.d0*k(235)*n(idx_O2) &
        +k(236)*n(idx_CO) &
        +k(238)*n(idx_H2Oj) &
        -k(248)*n(idx_O) &
        +k(249)*n(idx_CO) &
        +2.d0*k(257)*n(idx_O2) &
        +k(259)*n(idx_OH) &
        -k(267)*n(idx_C)*n(idx_O) &
        -k(268)*n(idx_C)*n(idx_O) &
        -k(269)*n(idx_Cj)*n(idx_O) &
        -k(270)*n(idx_Cj)*n(idx_O) &
        -k(273)*n(idx_H)*n(idx_O) &
        -2.d0*k(275)*n(idx_O)*n(idx_O)

    !OH
    !OH
    dn(idx_OH) = &
        -k(52)*n(idx_OH)*n(idx_H) &
        +k(62)*n(idx_CH)*n(idx_O) &
        +k(67)*n(idx_CH2)*n(idx_O) &
        +k(70)*n(idx_O)*n(idx_H2) &
        -k(71)*n(idx_OH)*n(idx_H) &
        -k(72)*n(idx_OH)*n(idx_H) &
        -k(73)*n(idx_H2)*n(idx_OH) &
        -k(74)*n(idx_C)*n(idx_OH) &
        -k(75)*n(idx_C)*n(idx_OH) &
        -k(76)*n(idx_O)*n(idx_OH) &
        -k(77)*n(idx_O)*n(idx_OH) &
        -2.d0*k(78)*n(idx_OH)*n(idx_OH) &
        +k(79)*n(idx_H2O)*n(idx_H) &
        +k(80)*n(idx_O2)*n(idx_H) &
        +2.d0*k(81)*n(idx_O2)*n(idx_H2) &
        +k(84)*n(idx_CO)*n(idx_H) &
        -k(105)*n(idx_OH)*n(idx_H3j) &
        -k(106)*n(idx_OH)*n(idx_H3j) &
        -k(107)*n(idx_OH)*n(idx_Cj) &
        -k(108)*n(idx_OH)*n(idx_Cj) &
        +k(120)*n(idx_O2)*n(idx_CH2j) &
        -k(141)*n(idx_OH)*n(idx_Hj) &
        -k(142)*n(idx_OH)*n(idx_Hj) &
        -k(143)*n(idx_OH)*n(idx_HEj) &
        -k(144)*n(idx_OH)*n(idx_HEj) &
        +k(147)*n(idx_H2O)*n(idx_HEj) &
        +k(148)*n(idx_H2O)*n(idx_HEj) &
        +k(173)*n(idx_H2Oj)*n(idx_E) &
        +k(175)*n(idx_H3Oj)*n(idx_E) &
        +k(178)*n(idx_H3Oj)*n(idx_E) &
        +k(182)*n(idx_HCOj)*n(idx_E) &
        +k(185)*n(idx_Hk)*n(idx_O) &
        -k(186)*n(idx_Hk)*n(idx_OH) &
        +k(190)*n(idx_Ok)*n(idx_H) &
        +k(205)*n(idx_O)*n(idx_H) &
        -k(207)*n(idx_OH)*n(idx_H) &
        -k(229)*n(idx_OH) &
        -k(230)*n(idx_OH) &
        +k(232)*n(idx_H2O) &
        +k(239)*n(idx_H2Oj) &
        +k(243)*n(idx_H3Oj) &
        -k(259)*n(idx_OH) &
        +k(261)*n(idx_H2O) &
        +k(273)*n(idx_H)*n(idx_O) &
        -k(274)*n(idx_OH)*n(idx_H)

    !CO
    !CO_GAS
    dn(idx_CO) = &
        dnChem_CO &
        -n(idx_CO)*(k(208)+k(212)) &
        +k(212)*n(idx_CO_total)

    !CH
    !CH
    dn(idx_CH) = &
        +k(56)*n(idx_C)*n(idx_H2) &
        -k(57)*n(idx_CH)*n(idx_H) &
        -k(58)*n(idx_CH)*n(idx_H2) &
        -k(59)*n(idx_CH)*n(idx_C) &
        -k(60)*n(idx_CH)*n(idx_O) &
        -k(61)*n(idx_CH)*n(idx_O) &
        -k(62)*n(idx_CH)*n(idx_O) &
        +k(63)*n(idx_CH2)*n(idx_H) &
        +k(67)*n(idx_CH2)*n(idx_O) &
        -k(130)*n(idx_CH)*n(idx_Hj) &
        -k(131)*n(idx_CH)*n(idx_Hj) &
        +k(165)*n(idx_CH2j)*n(idx_E) &
        +k(169)*n(idx_CH3j)*n(idx_E) &
        +k(170)*n(idx_CH3j)*n(idx_E) &
        +k(184)*n(idx_Hk)*n(idx_C) &
        +k(187)*n(idx_Ck)*n(idx_H) &
        +k(196)*n(idx_C)*n(idx_H) &
        -k(219)*n(idx_CH) &
        -k(220)*n(idx_CH) &
        +k(222)*n(idx_CH2) &
        -k(256)*n(idx_CH)

    !CH2
    !CH2
    dn(idx_CH2) = &
        +k(58)*n(idx_CH)*n(idx_H2) &
        -k(63)*n(idx_CH2)*n(idx_H) &
        -k(64)*n(idx_CH2)*n(idx_O) &
        -k(65)*n(idx_CH2)*n(idx_O) &
        -k(66)*n(idx_CH2)*n(idx_O) &
        -k(67)*n(idx_CH2)*n(idx_O) &
        -k(132)*n(idx_CH2)*n(idx_Hj) &
        -k(133)*n(idx_CH2)*n(idx_Hj) &
        -k(134)*n(idx_CH2)*n(idx_Hj) &
        -k(135)*n(idx_CH2)*n(idx_Hj) &
        -k(136)*n(idx_CH2)*n(idx_HEj) &
        -k(137)*n(idx_CH2)*n(idx_HEj) &
        -k(138)*n(idx_CH2)*n(idx_HEj) &
        -k(139)*n(idx_CH2)*n(idx_HEj) &
        +k(168)*n(idx_CH3j)*n(idx_E) &
        +k(188)*n(idx_Ck)*n(idx_H2) &
        +k(197)*n(idx_C)*n(idx_H2) &
        -k(222)*n(idx_CH2) &
        -k(223)*n(idx_CH2) &
        -k(260)*n(idx_CH2)

    !C2
    !C2
    dn(idx_C2) = &
        +k(59)*n(idx_CH)*n(idx_C) &
        -k(68)*n(idx_C2)*n(idx_O) &
        -k(69)*n(idx_C2)*n(idx_O) &
        -k(100)*n(idx_C2)*n(idx_Oj) &
        -k(140)*n(idx_C2)*n(idx_HEj) &
        +k(198)*n(idx_C)*n(idx_C) &
        -k(227)*n(idx_C2) &
        -k(251)*n(idx_C2) &
        +k(265)*n(idx_C)*n(idx_C) &
        +k(266)*n(idx_C)*n(idx_C)

    !HCO
    !HCO
    dn(idx_HCO) = &
        +k(66)*n(idx_CH2)*n(idx_O) &
        -k(262)*n(idx_HCO) &
        -k(263)*n(idx_HCO)

    !H2O
    !H2O_GAS
    dn(idx_H2O) = &
        dnChem_H2O &
        -n(idx_H2O)*(k(210)+k(211)) &
        +k(211)*n(idx_H2O_total)

    !O2
    !O2
    dn(idx_O2) = &
        +k(76)*n(idx_O)*n(idx_OH) &
        +k(77)*n(idx_O)*n(idx_OH) &
        -k(80)*n(idx_O2)*n(idx_H) &
        -k(81)*n(idx_O2)*n(idx_H2) &
        -k(82)*n(idx_O2)*n(idx_C) &
        -k(83)*n(idx_O2)*n(idx_C) &
        -k(118)*n(idx_O2)*n(idx_Cj) &
        -k(119)*n(idx_O2)*n(idx_Cj) &
        -k(120)*n(idx_O2)*n(idx_CH2j) &
        +k(122)*n(idx_C)*n(idx_O2j) &
        -k(153)*n(idx_O2)*n(idx_Hj) &
        -k(154)*n(idx_O2)*n(idx_HEj) &
        -k(155)*n(idx_O2)*n(idx_HEj) &
        +k(206)*n(idx_O)*n(idx_O) &
        -k(234)*n(idx_O2) &
        -k(235)*n(idx_O2) &
        -k(257)*n(idx_O2) &
        -k(258)*n(idx_O2) &
        +k(275)*n(idx_O)*n(idx_O)

    !CO_total
    !CO_TOTAL
    dn(idx_CO_total) = &
        dnChem_CO

    !H2O_total
    !H2O_TOTAL
    dn(idx_H2O_total) = &
        dnChem_H2O

    !H+
    !H+
    dn(idx_Hj) = &
        +k(1)*n(idx_H)*n(idx_E) &
        -k(2)*n(idx_Hj)*n(idx_E) &
        -k(3)*n(idx_Hj)*n(idx_E) &
        +k(8)*n(idx_HEj)*n(idx_H) &
        -k(9)*n(idx_HE)*n(idx_Hj) &
        -k(10)*n(idx_HE)*n(idx_Hj) &
        +k(13)*n(idx_H2)*n(idx_HEj) &
        -k(18)*n(idx_H)*n(idx_Hj) &
        -k(19)*n(idx_H)*n(idx_Hj) &
        +k(20)*n(idx_H2j)*n(idx_H) &
        -k(21)*n(idx_H2)*n(idx_Hj) &
        -k(22)*n(idx_H2)*n(idx_Hj) &
        -k(28)*n(idx_Hk)*n(idx_Hj) &
        -k(29)*n(idx_Hk)*n(idx_Hj) &
        +k(44)*n(idx_Oj)*n(idx_H) &
        -k(45)*n(idx_O)*n(idx_Hj) &
        -k(47)*n(idx_C)*n(idx_Hj) &
        +k(48)*n(idx_Cj)*n(idx_H) &
        -k(130)*n(idx_CH)*n(idx_Hj) &
        -k(131)*n(idx_CH)*n(idx_Hj) &
        -k(132)*n(idx_CH2)*n(idx_Hj) &
        -k(133)*n(idx_CH2)*n(idx_Hj) &
        -k(134)*n(idx_CH2)*n(idx_Hj) &
        -k(135)*n(idx_CH2)*n(idx_Hj) &
        -k(141)*n(idx_OH)*n(idx_Hj) &
        -k(142)*n(idx_OH)*n(idx_Hj) &
        -k(145)*n(idx_H2O)*n(idx_Hj) &
        -k(146)*n(idx_H2O)*n(idx_Hj) &
        +k(147)*n(idx_H2O)*n(idx_HEj) &
        +k(148)*n(idx_H2O)*n(idx_HEj) &
        -k(153)*n(idx_O2)*n(idx_Hj) &
        +k(158)*n(idx_COj)*n(idx_H) &
        -k(159)*n(idx_Ck)*n(idx_Hj) &
        -k(160)*n(idx_Ok)*n(idx_Hj) &
        -k(193)*n(idx_H2)*n(idx_Hj) &
        +k(193)*n(idx_H2)*n(idx_Hj) &
        -k(194)*n(idx_H2)*n(idx_Hj) &
        +k(214)*n(idx_H2j) &
        +k(215)*n(idx_H3j) &
        +k(221)*n(idx_CHj) &
        +k(231)*n(idx_OHj) &
        +k(239)*n(idx_H2Oj) &
        +k(242)*n(idx_H3Oj) &
        +k(246)*n(idx_H) &
        +k(253)*n(idx_H2) &
        +k(264)*n(idx_H2)

    !HE+
    !HE+
    dn(idx_HEj) = &
        +k(4)*n(idx_HE)*n(idx_E) &
        -k(5)*n(idx_HEj)*n(idx_E) &
        -k(6)*n(idx_HEj)*n(idx_E) &
        -k(7)*n(idx_HEj)*n(idx_E) &
        -k(8)*n(idx_HEj)*n(idx_H) &
        +k(9)*n(idx_HE)*n(idx_Hj) &
        +k(10)*n(idx_HE)*n(idx_Hj) &
        -k(12)*n(idx_H2)*n(idx_HEj) &
        -k(13)*n(idx_H2)*n(idx_HEj) &
        -k(14)*n(idx_H2)*n(idx_HEj) &
        +k(14)*n(idx_H2)*n(idx_HEj) &
        +k(15)*n(idx_HEjj)*n(idx_E) &
        -k(46)*n(idx_O)*n(idx_HEj) &
        -k(49)*n(idx_C)*n(idx_HEj) &
        -k(50)*n(idx_C)*n(idx_HEj) &
        -k(51)*n(idx_C)*n(idx_HEj) &
        -k(136)*n(idx_CH2)*n(idx_HEj) &
        -k(137)*n(idx_CH2)*n(idx_HEj) &
        -k(138)*n(idx_CH2)*n(idx_HEj) &
        -k(139)*n(idx_CH2)*n(idx_HEj) &
        -k(140)*n(idx_C2)*n(idx_HEj) &
        -k(143)*n(idx_OH)*n(idx_HEj) &
        -k(144)*n(idx_OH)*n(idx_HEj) &
        -k(147)*n(idx_H2O)*n(idx_HEj) &
        -k(148)*n(idx_H2O)*n(idx_HEj) &
        -k(149)*n(idx_H2O)*n(idx_HEj) &
        -k(150)*n(idx_H2O)*n(idx_HEj) &
        -k(151)*n(idx_H2O)*n(idx_HEj) &
        -k(152)*n(idx_H2O)*n(idx_HEj) &
        -k(154)*n(idx_O2)*n(idx_HEj) &
        -k(155)*n(idx_O2)*n(idx_HEj) &
        -k(156)*n(idx_CO)*n(idx_HEj) &
        -k(157)*n(idx_CO)*n(idx_HEj) &
        -k(161)*n(idx_HEj)*n(idx_Hk) &
        +k(247)*n(idx_HE)

    !H2+
    !H2+
    dn(idx_H2j) = &
        +k(12)*n(idx_H2)*n(idx_HEj) &
        +k(18)*n(idx_H)*n(idx_Hj) &
        +k(19)*n(idx_H)*n(idx_Hj) &
        -k(20)*n(idx_H2j)*n(idx_H) &
        +k(21)*n(idx_H2)*n(idx_Hj) &
        +k(22)*n(idx_H2)*n(idx_Hj) &
        +k(29)*n(idx_Hk)*n(idx_Hj) &
        -k(30)*n(idx_H2j)*n(idx_E) &
        -k(31)*n(idx_H2j)*n(idx_E) &
        -k(32)*n(idx_H2j)*n(idx_Hk) &
        -k(85)*n(idx_H2j)*n(idx_H2) &
        +k(86)*n(idx_H3j)*n(idx_H) &
        -k(87)*n(idx_C)*n(idx_H2j) &
        -k(102)*n(idx_O)*n(idx_H2j) &
        -k(214)*n(idx_H2j) &
        +k(216)*n(idx_H3j) &
        +k(238)*n(idx_H2Oj) &
        +k(243)*n(idx_H3Oj) &
        +k(254)*n(idx_H2)

    !C+
    !C+
    dn(idx_Cj) = &
        -k(37)*n(idx_Cj)*n(idx_E) &
        -k(38)*n(idx_Cj)*n(idx_E) &
        -k(39)*n(idx_Cj)*n(idx_E) &
        +k(42)*n(idx_C)*n(idx_E) &
        +k(47)*n(idx_C)*n(idx_Hj) &
        -k(48)*n(idx_Cj)*n(idx_H) &
        +k(49)*n(idx_C)*n(idx_HEj) &
        +k(50)*n(idx_C)*n(idx_HEj) &
        +k(51)*n(idx_C)*n(idx_HEj) &
        -k(90)*n(idx_Cj)*n(idx_H2) &
        +k(91)*n(idx_CHj)*n(idx_H) &
        -k(107)*n(idx_OH)*n(idx_Cj) &
        -k(108)*n(idx_OH)*n(idx_Cj) &
        -k(113)*n(idx_H2O)*n(idx_Cj) &
        -k(114)*n(idx_H2O)*n(idx_Cj) &
        -k(115)*n(idx_H2O)*n(idx_Cj) &
        -k(116)*n(idx_H2O)*n(idx_Cj) &
        -k(118)*n(idx_O2)*n(idx_Cj) &
        -k(119)*n(idx_O2)*n(idx_Cj) &
        +k(122)*n(idx_C)*n(idx_O2j) &
        +k(136)*n(idx_CH2)*n(idx_HEj) &
        +k(137)*n(idx_CH2)*n(idx_HEj) &
        +k(140)*n(idx_C2)*n(idx_HEj) &
        +k(156)*n(idx_CO)*n(idx_HEj) &
        -k(200)*n(idx_Cj)*n(idx_H) &
        -k(201)*n(idx_Cj)*n(idx_H2) &
        -k(202)*n(idx_Cj)*n(idx_O) &
        -k(203)*n(idx_Cj)*n(idx_O) &
        +k(217)*n(idx_C) &
        +k(255)*n(idx_C) &
        -k(269)*n(idx_Cj)*n(idx_O) &
        -k(270)*n(idx_Cj)*n(idx_O)

    !O+
    !O+
    dn(idx_Oj) = &
        -k(40)*n(idx_Oj)*n(idx_E) &
        -k(41)*n(idx_Oj)*n(idx_E) &
        +k(43)*n(idx_O)*n(idx_E) &
        -k(44)*n(idx_Oj)*n(idx_H) &
        +k(45)*n(idx_O)*n(idx_Hj) &
        +k(46)*n(idx_O)*n(idx_HEj) &
        -k(100)*n(idx_C2)*n(idx_Oj) &
        -k(101)*n(idx_Oj)*n(idx_H2) &
        +k(119)*n(idx_O2)*n(idx_Cj) &
        +k(143)*n(idx_OH)*n(idx_HEj) &
        +k(144)*n(idx_OH)*n(idx_HEj) &
        +k(155)*n(idx_O2)*n(idx_HEj) &
        +k(157)*n(idx_CO)*n(idx_HEj) &
        +k(240)*n(idx_H2Oj) &
        +k(248)*n(idx_O) &
        -k(271)*n(idx_C)*n(idx_Oj) &
        -k(272)*n(idx_C)*n(idx_Oj)

    !HOC+
    !HOC+
    dn(idx_HOCj) = &
        -k(53)*n(idx_HOCj)*n(idx_H2) &
        -k(54)*n(idx_HOCj)*n(idx_CO) &
        -k(55)*n(idx_HOCj)*n(idx_CO) &
        +k(98)*n(idx_CH3j)*n(idx_O) &
        +k(113)*n(idx_H2O)*n(idx_Cj) &
        +k(125)*n(idx_CO)*n(idx_H3j) &
        +k(126)*n(idx_CO)*n(idx_H3j) &
        -k(183)*n(idx_HOCj)*n(idx_E)

    !HCO+
    !HCO+
    dn(idx_HCOj) = &
        +k(53)*n(idx_HOCj)*n(idx_H2) &
        +k(54)*n(idx_HOCj)*n(idx_CO) &
        +k(55)*n(idx_HOCj)*n(idx_CO) &
        +k(61)*n(idx_CH)*n(idx_O) &
        +k(96)*n(idx_CH2j)*n(idx_O) &
        +k(99)*n(idx_CH3j)*n(idx_O) &
        +k(114)*n(idx_H2O)*n(idx_Cj) &
        +k(115)*n(idx_H2O)*n(idx_Cj) &
        +k(117)*n(idx_H3Oj)*n(idx_C) &
        +k(120)*n(idx_O2)*n(idx_CH2j) &
        +k(123)*n(idx_CO)*n(idx_H3j) &
        +k(124)*n(idx_CO)*n(idx_H3j) &
        -k(127)*n(idx_HCOj)*n(idx_C) &
        -k(128)*n(idx_HCOj)*n(idx_H2O) &
        -k(129)*n(idx_HCOj)*n(idx_H2O) &
        -k(181)*n(idx_HCOj)*n(idx_E) &
        -k(182)*n(idx_HCOj)*n(idx_E) &
        +k(263)*n(idx_HCO)

    !H3+
    !H3+
    dn(idx_H3j) = &
        +k(85)*n(idx_H2j)*n(idx_H2) &
        -k(86)*n(idx_H3j)*n(idx_H) &
        -k(88)*n(idx_C)*n(idx_H3j) &
        -k(89)*n(idx_C)*n(idx_H3j) &
        -k(103)*n(idx_O)*n(idx_H3j) &
        -k(104)*n(idx_O)*n(idx_H3j) &
        -k(105)*n(idx_OH)*n(idx_H3j) &
        -k(106)*n(idx_OH)*n(idx_H3j) &
        -k(111)*n(idx_H2O)*n(idx_H3j) &
        -k(112)*n(idx_H2O)*n(idx_H3j) &
        -k(123)*n(idx_CO)*n(idx_H3j) &
        -k(124)*n(idx_CO)*n(idx_H3j) &
        -k(125)*n(idx_CO)*n(idx_H3j) &
        -k(126)*n(idx_CO)*n(idx_H3j) &
        -k(162)*n(idx_H3j)*n(idx_E) &
        -k(163)*n(idx_H3j)*n(idx_E) &
        +k(194)*n(idx_H2)*n(idx_Hj) &
        -k(215)*n(idx_H3j) &
        -k(216)*n(idx_H3j)

    !CH+
    !CH+
    dn(idx_CHj) = &
        +k(87)*n(idx_C)*n(idx_H2j) &
        +k(88)*n(idx_C)*n(idx_H3j) &
        +k(90)*n(idx_Cj)*n(idx_H2) &
        -k(91)*n(idx_CHj)*n(idx_H) &
        -k(92)*n(idx_CHj)*n(idx_H2) &
        -k(93)*n(idx_CHj)*n(idx_O) &
        +k(94)*n(idx_CH2j)*n(idx_H) &
        +k(127)*n(idx_HCOj)*n(idx_C) &
        +k(130)*n(idx_CH)*n(idx_Hj) &
        +k(131)*n(idx_CH)*n(idx_Hj) &
        +k(132)*n(idx_CH2)*n(idx_Hj) &
        +k(133)*n(idx_CH2)*n(idx_Hj) &
        +k(138)*n(idx_CH2)*n(idx_HEj) &
        +k(139)*n(idx_CH2)*n(idx_HEj) &
        -k(164)*n(idx_CHj)*n(idx_E) &
        +k(200)*n(idx_Cj)*n(idx_H) &
        +k(220)*n(idx_CH) &
        -k(221)*n(idx_CHj) &
        +k(224)*n(idx_CH2j) &
        +k(226)*n(idx_CH3j)

    !CH2+
    !CH2+
    dn(idx_CH2j) = &
        +k(89)*n(idx_C)*n(idx_H3j) &
        +k(92)*n(idx_CHj)*n(idx_H2) &
        -k(94)*n(idx_CH2j)*n(idx_H) &
        -k(95)*n(idx_CH2j)*n(idx_H2) &
        -k(96)*n(idx_CH2j)*n(idx_O) &
        +k(97)*n(idx_CH3j)*n(idx_H) &
        -k(120)*n(idx_O2)*n(idx_CH2j) &
        +k(134)*n(idx_CH2)*n(idx_Hj) &
        +k(135)*n(idx_CH2)*n(idx_Hj) &
        -k(165)*n(idx_CH2j)*n(idx_E) &
        -k(166)*n(idx_CH2j)*n(idx_E) &
        -k(167)*n(idx_CH2j)*n(idx_E) &
        +k(201)*n(idx_Cj)*n(idx_H2) &
        +k(223)*n(idx_CH2) &
        -k(224)*n(idx_CH2j) &
        +k(225)*n(idx_CH3j) &
        +k(260)*n(idx_CH2)

    !CO+
    !CO+
    dn(idx_COj) = &
        +k(93)*n(idx_CHj)*n(idx_O) &
        +k(100)*n(idx_C2)*n(idx_Oj) &
        +k(107)*n(idx_OH)*n(idx_Cj) &
        +k(108)*n(idx_OH)*n(idx_Cj) &
        +k(118)*n(idx_O2)*n(idx_Cj) &
        +k(121)*n(idx_C)*n(idx_O2j) &
        -k(158)*n(idx_COj)*n(idx_H) &
        -k(180)*n(idx_COj)*n(idx_E) &
        +k(202)*n(idx_Cj)*n(idx_O) &
        +k(203)*n(idx_Cj)*n(idx_O) &
        +k(250)*n(idx_CO) &
        +k(269)*n(idx_Cj)*n(idx_O) &
        +k(270)*n(idx_Cj)*n(idx_O) &
        +k(271)*n(idx_C)*n(idx_Oj) &
        +k(272)*n(idx_C)*n(idx_Oj)

    !CH3+
    !CH3+
    dn(idx_CH3j) = &
        +k(95)*n(idx_CH2j)*n(idx_H2) &
        -k(97)*n(idx_CH3j)*n(idx_H) &
        -k(98)*n(idx_CH3j)*n(idx_O) &
        -k(99)*n(idx_CH3j)*n(idx_O) &
        -k(168)*n(idx_CH3j)*n(idx_E) &
        -k(169)*n(idx_CH3j)*n(idx_E) &
        -k(170)*n(idx_CH3j)*n(idx_E) &
        -k(225)*n(idx_CH3j) &
        -k(226)*n(idx_CH3j)

    !OH+
    !OH+
    dn(idx_OHj) = &
        +k(101)*n(idx_Oj)*n(idx_H2) &
        +k(102)*n(idx_O)*n(idx_H2j) &
        +k(103)*n(idx_O)*n(idx_H3j) &
        -k(109)*n(idx_OHj)*n(idx_H2) &
        +k(141)*n(idx_OH)*n(idx_Hj) &
        +k(142)*n(idx_OH)*n(idx_Hj) &
        +k(149)*n(idx_H2O)*n(idx_HEj) &
        +k(150)*n(idx_H2O)*n(idx_HEj) &
        -k(171)*n(idx_OHj)*n(idx_E) &
        +k(230)*n(idx_OH) &
        -k(231)*n(idx_OHj) &
        +k(241)*n(idx_H2Oj) &
        +k(245)*n(idx_H3Oj)

    !H2O+
    !H2O+
    dn(idx_H2Oj) = &
        +k(104)*n(idx_O)*n(idx_H3j) &
        +k(105)*n(idx_OH)*n(idx_H3j) &
        +k(106)*n(idx_OH)*n(idx_H3j) &
        +k(109)*n(idx_OHj)*n(idx_H2) &
        -k(110)*n(idx_H2Oj)*n(idx_H2) &
        +k(116)*n(idx_H2O)*n(idx_Cj) &
        +k(145)*n(idx_H2O)*n(idx_Hj) &
        +k(146)*n(idx_H2O)*n(idx_Hj) &
        +k(151)*n(idx_H2O)*n(idx_HEj) &
        +k(152)*n(idx_H2O)*n(idx_HEj) &
        -k(172)*n(idx_H2Oj)*n(idx_E) &
        -k(173)*n(idx_H2Oj)*n(idx_E) &
        -k(174)*n(idx_H2Oj)*n(idx_E) &
        +k(233)*n(idx_H2O) &
        -k(238)*n(idx_H2Oj) &
        -k(239)*n(idx_H2Oj) &
        -k(240)*n(idx_H2Oj) &
        -k(241)*n(idx_H2Oj) &
        +k(244)*n(idx_H3Oj)

    !H3O+
    !H3O+
    dn(idx_H3Oj) = &
        +k(110)*n(idx_H2Oj)*n(idx_H2) &
        +k(111)*n(idx_H2O)*n(idx_H3j) &
        +k(112)*n(idx_H2O)*n(idx_H3j) &
        -k(117)*n(idx_H3Oj)*n(idx_C) &
        +k(128)*n(idx_HCOj)*n(idx_H2O) &
        +k(129)*n(idx_HCOj)*n(idx_H2O) &
        -k(175)*n(idx_H3Oj)*n(idx_E) &
        -k(176)*n(idx_H3Oj)*n(idx_E) &
        -k(177)*n(idx_H3Oj)*n(idx_E) &
        -k(178)*n(idx_H3Oj)*n(idx_E) &
        -k(242)*n(idx_H3Oj) &
        -k(243)*n(idx_H3Oj) &
        -k(244)*n(idx_H3Oj) &
        -k(245)*n(idx_H3Oj)

    !O2+
    !O2+
    dn(idx_O2j) = &
        -k(121)*n(idx_C)*n(idx_O2j) &
        -k(122)*n(idx_C)*n(idx_O2j) &
        +k(153)*n(idx_O2)*n(idx_Hj) &
        +k(154)*n(idx_O2)*n(idx_HEj) &
        -k(179)*n(idx_O2j)*n(idx_E) &
        +k(234)*n(idx_O2) &
        +k(258)*n(idx_O2)

    !HE++
    !HE++
    dn(idx_HEjj) = &
        +k(7)*n(idx_HEj)*n(idx_E) &
        -k(15)*n(idx_HEjj)*n(idx_E)

    !CR

    !CR
    dn(idx_CR) = 0.d0

    !g

    !g
    dn(idx_g) = 0.d0

    !Tgas

    !Tgas
    dn(idx_Tgas) = 0.d0

    !dummy

    !dummy
    dn(idx_dummy) = 0.d0

    krome_gamma = gamma_index(n(:))

    dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
        - cooling(n(:), Tgas)  ) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    krome_gamma = gamma_index(n(:))

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(1)*n(idx_H)  &
          +2.d0*k(1)*n(idx_H)  &
          -k(2)*n(idx_Hj)  &
          -k(3)*n(idx_Hj)  &
          -k(4)*n(idx_HE)  &
          +2.d0*k(4)*n(idx_HE)  &
          -k(5)*n(idx_HEj)  &
          -k(6)*n(idx_HEj)  &
          -k(7)*n(idx_HEj)  &
          +2.d0*k(7)*n(idx_HEj)  &
          -k(15)*n(idx_HEjj)  &
          -k(16)*n(idx_H)  &
          -k(23)*n(idx_H2)  &
          +k(23)*n(idx_H2)  &
          -k(25)*n(idx_Hk)  &
          +2.d0*k(25)*n(idx_Hk)  &
          -k(30)*n(idx_H2j)  &
          -k(31)*n(idx_H2j)  &
          -k(37)*n(idx_Cj)  &
          -k(38)*n(idx_Cj)  &
          -k(39)*n(idx_Cj)  &
          -k(40)*n(idx_Oj)  &
          -k(41)*n(idx_Oj)  &
          -k(42)*n(idx_C)  &
          +2.d0*k(42)*n(idx_C)  &
          -k(43)*n(idx_O)  &
          +2.d0*k(43)*n(idx_O)  &
          -k(162)*n(idx_H3j)  &
          -k(163)*n(idx_H3j)  &
          -k(164)*n(idx_CHj)  &
          -k(165)*n(idx_CH2j)  &
          -k(166)*n(idx_CH2j)  &
          -k(167)*n(idx_CH2j)  &
          -k(168)*n(idx_CH3j)  &
          -k(169)*n(idx_CH3j)  &
          -k(170)*n(idx_CH3j)  &
          -k(171)*n(idx_OHj)  &
          -k(172)*n(idx_H2Oj)  &
          -k(173)*n(idx_H2Oj)  &
          -k(174)*n(idx_H2Oj)  &
          -k(175)*n(idx_H3Oj)  &
          -k(176)*n(idx_H3Oj)  &
          -k(177)*n(idx_H3Oj)  &
          -k(178)*n(idx_H3Oj)  &
          -k(179)*n(idx_O2j)  &
          -k(180)*n(idx_COj)  &
          -k(181)*n(idx_HCOj)  &
          -k(182)*n(idx_HCOj)  &
          -k(183)*n(idx_HOCj)  &
          -k(195)*n(idx_C)  &
          -k(204)*n(idx_O)
      pdj(2) =  &
          +k(16)*n(idx_H)  &
          -k(25)*n(idx_Hk)
      pdj(3) =  &
          +k(195)*n(idx_C)
      pdj(4) =  &
          +k(204)*n(idx_O)
      pdj(5) =  &
          -k(1)*n(idx_H)  &
          +k(2)*n(idx_Hj)  &
          +k(3)*n(idx_Hj)  &
          -k(16)*n(idx_H)  &
          +2.d0*k(23)*n(idx_H2)  &
          +k(25)*n(idx_Hk)  &
          +2.d0*k(30)*n(idx_H2j)  &
          +2.d0*k(31)*n(idx_H2j)  &
          +k(162)*n(idx_H3j)  &
          +3.d0*k(163)*n(idx_H3j)  &
          +k(164)*n(idx_CHj)  &
          +k(165)*n(idx_CH2j)  &
          +2.d0*k(167)*n(idx_CH2j)  &
          +k(168)*n(idx_CH3j)  &
          +2.d0*k(170)*n(idx_CH3j)  &
          +k(171)*n(idx_OHj)  &
          +k(173)*n(idx_H2Oj)  &
          +2.d0*k(174)*n(idx_H2Oj)  &
          +2.d0*k(175)*n(idx_H3Oj)  &
          +k(176)*n(idx_H3Oj)  &
          +k(177)*n(idx_H3Oj)  &
          +k(181)*n(idx_HCOj)  &
          +k(183)*n(idx_HOCj)
      pdj(6) =  &
          -k(4)*n(idx_HE)  &
          +k(5)*n(idx_HEj)  &
          +k(6)*n(idx_HEj)
      pdj(7) =  &
          -k(23)*n(idx_H2)  &
          +k(162)*n(idx_H3j)  &
          +k(166)*n(idx_CH2j)  &
          +k(169)*n(idx_CH3j)  &
          +k(172)*n(idx_H2Oj)  &
          +k(176)*n(idx_H3Oj)  &
          +k(178)*n(idx_H3Oj)
      pdj(8) =  &
          +k(37)*n(idx_Cj)  &
          +k(38)*n(idx_Cj)  &
          +k(39)*n(idx_Cj)  &
          -k(42)*n(idx_C)  &
          +k(164)*n(idx_CHj)  &
          +k(166)*n(idx_CH2j)  &
          +k(167)*n(idx_CH2j)  &
          +k(180)*n(idx_COj)  &
          +k(182)*n(idx_HCOj)  &
          -k(195)*n(idx_C)
      pdj(9) =  &
          +k(40)*n(idx_Oj)  &
          +k(41)*n(idx_Oj)  &
          -k(43)*n(idx_O)  &
          +k(171)*n(idx_OHj)  &
          +k(172)*n(idx_H2Oj)  &
          +k(174)*n(idx_H2Oj)  &
          +k(176)*n(idx_H3Oj)  &
          +2.d0*k(179)*n(idx_O2j)  &
          +k(180)*n(idx_COj)  &
          -k(204)*n(idx_O)
      pdj(10) =  &
          +k(173)*n(idx_H2Oj)  &
          +k(175)*n(idx_H3Oj)  &
          +k(178)*n(idx_H3Oj)  &
          +k(182)*n(idx_HCOj)
      pdj(11) =  &
          +k(181)*n(idx_HCOj)  &
          +k(183)*n(idx_HOCj)
      pdj(12) =  &
          +k(165)*n(idx_CH2j)  &
          +k(169)*n(idx_CH3j)  &
          +k(170)*n(idx_CH3j)
      pdj(13) =  &
          +k(168)*n(idx_CH3j)
      pdj(16) =  &
          +k(177)*n(idx_H3Oj)
      pdj(20) =  &
          +k(1)*n(idx_H)  &
          -k(2)*n(idx_Hj)  &
          -k(3)*n(idx_Hj)
      pdj(21) =  &
          +k(4)*n(idx_HE)  &
          -k(5)*n(idx_HEj)  &
          -k(6)*n(idx_HEj)  &
          -k(7)*n(idx_HEj)  &
          +k(15)*n(idx_HEjj)
      pdj(22) =  &
          -k(30)*n(idx_H2j)  &
          -k(31)*n(idx_H2j)
      pdj(23) =  &
          -k(37)*n(idx_Cj)  &
          -k(38)*n(idx_Cj)  &
          -k(39)*n(idx_Cj)  &
          +k(42)*n(idx_C)
      pdj(24) =  &
          -k(40)*n(idx_Oj)  &
          -k(41)*n(idx_Oj)  &
          +k(43)*n(idx_O)
      pdj(25) =  &
          -k(183)*n(idx_HOCj)
      pdj(26) =  &
          -k(181)*n(idx_HCOj)  &
          -k(182)*n(idx_HCOj)
      pdj(27) =  &
          -k(162)*n(idx_H3j)  &
          -k(163)*n(idx_H3j)
      pdj(28) =  &
          -k(164)*n(idx_CHj)
      pdj(29) =  &
          -k(165)*n(idx_CH2j)  &
          -k(166)*n(idx_CH2j)  &
          -k(167)*n(idx_CH2j)
      pdj(30) =  &
          -k(180)*n(idx_COj)
      pdj(31) =  &
          -k(168)*n(idx_CH3j)  &
          -k(169)*n(idx_CH3j)  &
          -k(170)*n(idx_CH3j)
      pdj(32) =  &
          -k(171)*n(idx_OHj)
      pdj(33) =  &
          -k(172)*n(idx_H2Oj)  &
          -k(173)*n(idx_H2Oj)  &
          -k(174)*n(idx_H2Oj)
      pdj(34) =  &
          -k(175)*n(idx_H3Oj)  &
          -k(176)*n(idx_H3Oj)  &
          -k(177)*n(idx_H3Oj)  &
          -k(178)*n(idx_H3Oj)
      pdj(35) =  &
          -k(179)*n(idx_O2j)
      pdj(36) =  &
          +k(7)*n(idx_HEj)  &
          -k(15)*n(idx_HEjj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(1)*1d-3
      if(dnn>0.d0) then
        nn(1) = n(1) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==2) then
      pdj(1) =  &
          +k(17)*n(idx_H)  &
          -k(25)*n(idx_E)  &
          +2.d0*k(25)*n(idx_E)  &
          +k(26)*n(idx_H)  &
          +k(27)*n(idx_H)  &
          +k(29)*n(idx_Hj)  &
          +k(184)*n(idx_C)  &
          +k(185)*n(idx_O)  &
          +k(186)*n(idx_OH)  &
          +k(213)
      pdj(2) =  &
          -k(17)*n(idx_H)  &
          -k(25)*n(idx_E)  &
          -k(26)*n(idx_H)  &
          -k(27)*n(idx_H)  &
          -k(28)*n(idx_Hj)  &
          -k(29)*n(idx_Hj)  &
          -k(32)*n(idx_H2j)  &
          -k(161)*n(idx_HEj)  &
          -k(184)*n(idx_C)  &
          -k(185)*n(idx_O)  &
          -k(186)*n(idx_OH)  &
          -k(213)
      pdj(5) =  &
          -k(17)*n(idx_H)  &
          +k(25)*n(idx_E)  &
          -k(26)*n(idx_H)  &
          +2.d0*k(26)*n(idx_H)  &
          -k(27)*n(idx_H)  &
          +2.d0*k(27)*n(idx_H)  &
          +2.d0*k(28)*n(idx_Hj)  &
          +k(32)*n(idx_H2j)  &
          +k(161)*n(idx_HEj)  &
          +k(213)
      pdj(6) =  &
          +k(161)*n(idx_HEj)
      pdj(7) =  &
          +k(17)*n(idx_H)  &
          +k(32)*n(idx_H2j)
      pdj(8) =  &
          -k(184)*n(idx_C)
      pdj(9) =  &
          -k(185)*n(idx_O)
      pdj(10) =  &
          +k(185)*n(idx_O)  &
          -k(186)*n(idx_OH)
      pdj(12) =  &
          +k(184)*n(idx_C)
      pdj(16) =  &
          +k(186)*n(idx_OH)
      pdj(20) =  &
          -k(28)*n(idx_Hj)  &
          -k(29)*n(idx_Hj)
      pdj(21) =  &
          -k(161)*n(idx_HEj)
      pdj(22) =  &
          +k(29)*n(idx_Hj)  &
          -k(32)*n(idx_H2j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(2)*1d-3
      if(dnn>0.d0) then
        nn(2) = n(2) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==3) then
      pdj(1) =  &
          +k(187)*n(idx_H)  &
          +k(188)*n(idx_H2)  &
          +k(189)*n(idx_O)  &
          +k(218)
      pdj(3) =  &
          -k(159)*n(idx_Hj)  &
          -k(187)*n(idx_H)  &
          -k(188)*n(idx_H2)  &
          -k(189)*n(idx_O)  &
          -k(218)
      pdj(5) =  &
          +k(159)*n(idx_Hj)  &
          -k(187)*n(idx_H)
      pdj(7) =  &
          -k(188)*n(idx_H2)
      pdj(8) =  &
          +k(159)*n(idx_Hj)  &
          +k(218)
      pdj(9) =  &
          -k(189)*n(idx_O)
      pdj(11) =  &
          +k(189)*n(idx_O)
      pdj(12) =  &
          +k(187)*n(idx_H)
      pdj(13) =  &
          +k(188)*n(idx_H2)
      pdj(20) =  &
          -k(159)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(3)*1d-3
      if(dnn>0.d0) then
        nn(3) = n(3) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==4) then
      pdj(1) =  &
          +k(190)*n(idx_H)  &
          +k(191)*n(idx_H2)  &
          +k(192)*n(idx_C)  &
          +k(228)
      pdj(4) =  &
          -k(160)*n(idx_Hj)  &
          -k(190)*n(idx_H)  &
          -k(191)*n(idx_H2)  &
          -k(192)*n(idx_C)  &
          -k(228)
      pdj(5) =  &
          +k(160)*n(idx_Hj)  &
          -k(190)*n(idx_H)
      pdj(7) =  &
          -k(191)*n(idx_H2)
      pdj(8) =  &
          -k(192)*n(idx_C)
      pdj(9) =  &
          +k(160)*n(idx_Hj)  &
          +k(228)
      pdj(10) =  &
          +k(190)*n(idx_H)
      pdj(11) =  &
          +k(192)*n(idx_C)
      pdj(16) =  &
          +k(191)*n(idx_H2)
      pdj(20) =  &
          -k(160)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(4)*1d-3
      if(dnn>0.d0) then
        nn(4) = n(4) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==5) then
      pdj(1) =  &
          -k(1)*n(idx_E)  &
          +2.d0*k(1)*n(idx_E)  &
          -k(16)*n(idx_E)  &
          +k(17)*n(idx_Hk)  &
          +k(26)*n(idx_Hk)  &
          +k(27)*n(idx_Hk)  &
          +k(187)*n(idx_Ck)  &
          +k(190)*n(idx_Ok)  &
          +k(246)
      pdj(2) =  &
          +k(16)*n(idx_E)  &
          -k(17)*n(idx_Hk)  &
          -k(26)*n(idx_Hk)  &
          -k(27)*n(idx_Hk)
      pdj(3) =  &
          -k(187)*n(idx_Ck)
      pdj(4) =  &
          -k(190)*n(idx_Ok)
      pdj(5) =  &
          -k(1)*n(idx_E)  &
          -k(8)*n(idx_HEj)  &
          -k(16)*n(idx_E)  &
          -k(17)*n(idx_Hk)  &
          -k(18)*n(idx_Hj)  &
          -k(19)*n(idx_Hj)  &
          -k(20)*n(idx_H2j)  &
          -k(24)*n(idx_H2)  &
          +3.d0*k(24)*n(idx_H2)  &
          -k(26)*n(idx_Hk)  &
          +2.d0*k(26)*n(idx_Hk)  &
          -k(27)*n(idx_Hk)  &
          +2.d0*k(27)*n(idx_Hk)  &
          -4.d0*k(34)*n(idx_H)*n(idx_HE)  &
          -9.d0*k(35)*n(idx_H)*n(idx_H)  &
          +3.d0*k(35)*n(idx_H)*n(idx_H)  &
          -4.d0*k(36)*n(idx_H2)*n(idx_H)  &
          -k(44)*n(idx_Oj)  &
          -k(48)*n(idx_Cj)  &
          -k(52)*n(idx_OH)  &
          +2.d0*k(52)*n(idx_OH)  &
          -k(57)*n(idx_CH)  &
          -k(63)*n(idx_CH2)  &
          -k(71)*n(idx_OH)  &
          -k(72)*n(idx_OH)  &
          -k(79)*n(idx_H2O)  &
          -k(80)*n(idx_O2)  &
          -k(84)*n(idx_CO)  &
          -k(86)*n(idx_H3j)  &
          -k(91)*n(idx_CHj)  &
          -k(94)*n(idx_CH2j)  &
          -k(97)*n(idx_CH3j)  &
          -k(158)*n(idx_COj)  &
          -k(187)*n(idx_Ck)  &
          -k(190)*n(idx_Ok)  &
          -k(196)*n(idx_C)  &
          -k(200)*n(idx_Cj)  &
          -k(205)*n(idx_O)  &
          -k(207)*n(idx_OH)  &
          -k(246)  &
          -k(273)*n(idx_O)  &
          -k(274)*n(idx_OH)
      pdj(6) =  &
          +k(8)*n(idx_HEj)  &
          -2.d0*k(34)*n(idx_H)*n(idx_HE)  &
          +2.d0*k(34)*n(idx_H)*n(idx_HE)
      pdj(7) =  &
          +k(17)*n(idx_Hk)  &
          +k(20)*n(idx_H2j)  &
          -k(24)*n(idx_H2)  &
          +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
          +3.d0*k(35)*n(idx_H)*n(idx_H)  &
          -2.d0*k(36)*n(idx_H2)*n(idx_H)  &
          +4.d0*k(36)*n(idx_H2)*n(idx_H)  &
          +k(57)*n(idx_CH)  &
          +k(63)*n(idx_CH2)  &
          +k(71)*n(idx_OH)  &
          +k(72)*n(idx_OH)  &
          +k(79)*n(idx_H2O)  &
          +k(86)*n(idx_H3j)  &
          +k(91)*n(idx_CHj)  &
          +k(94)*n(idx_CH2j)  &
          +k(97)*n(idx_CH3j)
      pdj(8) =  &
          +k(48)*n(idx_Cj)  &
          +k(57)*n(idx_CH)  &
          +k(84)*n(idx_CO)  &
          -k(196)*n(idx_C)
      pdj(9) =  &
          +k(44)*n(idx_Oj)  &
          +k(52)*n(idx_OH)  &
          +k(71)*n(idx_OH)  &
          +k(72)*n(idx_OH)  &
          +k(80)*n(idx_O2)  &
          -k(205)*n(idx_O)  &
          -k(273)*n(idx_O)
      pdj(10) =  &
          -k(52)*n(idx_OH)  &
          -k(71)*n(idx_OH)  &
          -k(72)*n(idx_OH)  &
          +k(79)*n(idx_H2O)  &
          +k(80)*n(idx_O2)  &
          +k(84)*n(idx_CO)  &
          +k(190)*n(idx_Ok)  &
          +k(205)*n(idx_O)  &
          -k(207)*n(idx_OH)  &
          +k(273)*n(idx_O)  &
          -k(274)*n(idx_OH)
      pdj(11) =  &
          -k(84)*n(idx_CO)  &
          +k(158)*n(idx_COj)
      pdj(12) =  &
          -k(57)*n(idx_CH)  &
          +k(63)*n(idx_CH2)  &
          +k(187)*n(idx_Ck)  &
          +k(196)*n(idx_C)
      pdj(13) =  &
          -k(63)*n(idx_CH2)
      pdj(16) =  &
          -k(79)*n(idx_H2O)  &
          +k(207)*n(idx_OH)  &
          +k(274)*n(idx_OH)
      pdj(17) =  &
          -k(80)*n(idx_O2)
      pdj(20) =  &
          +k(1)*n(idx_E)  &
          +k(8)*n(idx_HEj)  &
          -k(18)*n(idx_Hj)  &
          -k(19)*n(idx_Hj)  &
          +k(20)*n(idx_H2j)  &
          +k(44)*n(idx_Oj)  &
          +k(48)*n(idx_Cj)  &
          +k(158)*n(idx_COj)  &
          +k(246)
      pdj(21) =  &
          -k(8)*n(idx_HEj)
      pdj(22) =  &
          +k(18)*n(idx_Hj)  &
          +k(19)*n(idx_Hj)  &
          -k(20)*n(idx_H2j)  &
          +k(86)*n(idx_H3j)
      pdj(23) =  &
          -k(48)*n(idx_Cj)  &
          +k(91)*n(idx_CHj)  &
          -k(200)*n(idx_Cj)
      pdj(24) =  &
          -k(44)*n(idx_Oj)
      pdj(27) =  &
          -k(86)*n(idx_H3j)
      pdj(28) =  &
          -k(91)*n(idx_CHj)  &
          +k(94)*n(idx_CH2j)  &
          +k(200)*n(idx_Cj)
      pdj(29) =  &
          -k(94)*n(idx_CH2j)  &
          +k(97)*n(idx_CH3j)
      pdj(30) =  &
          -k(158)*n(idx_COj)
      pdj(31) =  &
          -k(97)*n(idx_CH3j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(5)*1d-3
      if(dnn>0.d0) then
        nn(5) = n(5) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==6) then
      pdj(1) =  &
          -k(4)*n(idx_E)  &
          +2.d0*k(4)*n(idx_E)  &
          +k(247)
      pdj(5) =  &
          +k(9)*n(idx_Hj)  &
          +k(10)*n(idx_Hj)  &
          +2.d0*k(11)*n(idx_H2)  &
          -2.d0*k(34)*n(idx_H)*n(idx_H)
      pdj(6) =  &
          -k(4)*n(idx_E)  &
          -k(9)*n(idx_Hj)  &
          -k(10)*n(idx_Hj)  &
          -k(11)*n(idx_H2)  &
          +k(11)*n(idx_H2)  &
          -k(34)*n(idx_H)*n(idx_H)  &
          +k(34)*n(idx_H)*n(idx_H)  &
          -k(247)
      pdj(7) =  &
          -k(11)*n(idx_H2)  &
          +k(34)*n(idx_H)*n(idx_H)
      pdj(20) =  &
          -k(9)*n(idx_Hj)  &
          -k(10)*n(idx_Hj)
      pdj(21) =  &
          +k(4)*n(idx_E)  &
          +k(9)*n(idx_Hj)  &
          +k(10)*n(idx_Hj)  &
          +k(247)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(6)*1d-3
      if(dnn>0.d0) then
        nn(6) = n(6) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==7) then
      pdj(1) =  &
          -k(23)*n(idx_E)  &
          +k(23)*n(idx_E)  &
          +k(188)*n(idx_Ck)  &
          +k(191)*n(idx_Ok)  &
          +k(254)  &
          +k(264)
      pdj(2) =  &
          +k(253)
      pdj(3) =  &
          -k(188)*n(idx_Ck)
      pdj(4) =  &
          -k(191)*n(idx_Ok)
      pdj(5) =  &
          +2.d0*k(11)*n(idx_HE)  &
          +k(13)*n(idx_HEj)  &
          +2.d0*k(14)*n(idx_HEj)  &
          +k(21)*n(idx_Hj)  &
          +k(22)*n(idx_Hj)  &
          +2.d0*k(23)*n(idx_E)  &
          -k(24)*n(idx_H)  &
          +3.d0*k(24)*n(idx_H)  &
          +4.d0*k(33)*n(idx_H2)  &
          -2.d0*k(36)*n(idx_H)*n(idx_H)  &
          +k(56)*n(idx_C)  &
          +k(58)*n(idx_CH)  &
          +k(70)*n(idx_O)  &
          +k(73)*n(idx_OH)  &
          +k(85)*n(idx_H2j)  &
          +k(90)*n(idx_Cj)  &
          +k(92)*n(idx_CHj)  &
          +k(95)*n(idx_CH2j)  &
          +k(101)*n(idx_Oj)  &
          +k(109)*n(idx_OHj)  &
          +k(110)*n(idx_H2Oj)  &
          +2.d0*k(193)*n(idx_Hj)  &
          +2.d0*k(237)  &
          +2.d0*k(252)  &
          +k(264)
      pdj(6) =  &
          -k(11)*n(idx_HE)  &
          +k(11)*n(idx_HE)  &
          +k(12)*n(idx_HEj)  &
          +k(13)*n(idx_HEj)
      pdj(7) =  &
          -k(11)*n(idx_HE)  &
          -k(12)*n(idx_HEj)  &
          -k(13)*n(idx_HEj)  &
          -k(14)*n(idx_HEj)  &
          -k(21)*n(idx_Hj)  &
          -k(22)*n(idx_Hj)  &
          -k(23)*n(idx_E)  &
          -k(24)*n(idx_H)  &
          -4.d0*k(33)*n(idx_H2)  &
          +2.d0*k(33)*n(idx_H2)  &
          -k(36)*n(idx_H)*n(idx_H)  &
          +2.d0*k(36)*n(idx_H)*n(idx_H)  &
          -k(53)*n(idx_HOCj)  &
          +k(53)*n(idx_HOCj)  &
          -k(56)*n(idx_C)  &
          -k(58)*n(idx_CH)  &
          -k(70)*n(idx_O)  &
          -k(73)*n(idx_OH)  &
          -k(81)*n(idx_O2)  &
          -k(85)*n(idx_H2j)  &
          -k(90)*n(idx_Cj)  &
          -k(92)*n(idx_CHj)  &
          -k(95)*n(idx_CH2j)  &
          -k(101)*n(idx_Oj)  &
          -k(109)*n(idx_OHj)  &
          -k(110)*n(idx_H2Oj)  &
          -k(188)*n(idx_Ck)  &
          -k(191)*n(idx_Ok)  &
          -k(193)*n(idx_Hj)  &
          -k(194)*n(idx_Hj)  &
          -k(197)*n(idx_C)  &
          -k(201)*n(idx_Cj)  &
          -k(237)  &
          -k(252)  &
          -k(253)  &
          -k(254)  &
          -k(264)
      pdj(8) =  &
          -k(56)*n(idx_C)  &
          -k(197)*n(idx_C)
      pdj(9) =  &
          -k(70)*n(idx_O)
      pdj(10) =  &
          +k(70)*n(idx_O)  &
          -k(73)*n(idx_OH)  &
          +2.d0*k(81)*n(idx_O2)
      pdj(12) =  &
          +k(56)*n(idx_C)  &
          -k(58)*n(idx_CH)
      pdj(13) =  &
          +k(58)*n(idx_CH)  &
          +k(188)*n(idx_Ck)  &
          +k(197)*n(idx_C)
      pdj(16) =  &
          +k(73)*n(idx_OH)  &
          +k(191)*n(idx_Ok)
      pdj(17) =  &
          -k(81)*n(idx_O2)
      pdj(20) =  &
          +k(13)*n(idx_HEj)  &
          -k(21)*n(idx_Hj)  &
          -k(22)*n(idx_Hj)  &
          -k(193)*n(idx_Hj)  &
          +k(193)*n(idx_Hj)  &
          -k(194)*n(idx_Hj)  &
          +k(253)  &
          +k(264)
      pdj(21) =  &
          -k(12)*n(idx_HEj)  &
          -k(13)*n(idx_HEj)  &
          -k(14)*n(idx_HEj)  &
          +k(14)*n(idx_HEj)
      pdj(22) =  &
          +k(12)*n(idx_HEj)  &
          +k(21)*n(idx_Hj)  &
          +k(22)*n(idx_Hj)  &
          -k(85)*n(idx_H2j)  &
          +k(254)
      pdj(23) =  &
          -k(90)*n(idx_Cj)  &
          -k(201)*n(idx_Cj)
      pdj(24) =  &
          -k(101)*n(idx_Oj)
      pdj(25) =  &
          -k(53)*n(idx_HOCj)
      pdj(26) =  &
          +k(53)*n(idx_HOCj)
      pdj(27) =  &
          +k(85)*n(idx_H2j)  &
          +k(194)*n(idx_Hj)
      pdj(28) =  &
          +k(90)*n(idx_Cj)  &
          -k(92)*n(idx_CHj)
      pdj(29) =  &
          +k(92)*n(idx_CHj)  &
          -k(95)*n(idx_CH2j)  &
          +k(201)*n(idx_Cj)
      pdj(31) =  &
          +k(95)*n(idx_CH2j)
      pdj(32) =  &
          +k(101)*n(idx_Oj)  &
          -k(109)*n(idx_OHj)
      pdj(33) =  &
          +k(109)*n(idx_OHj)  &
          -k(110)*n(idx_H2Oj)
      pdj(34) =  &
          +k(110)*n(idx_H2Oj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(7)*1d-3
      if(dnn>0.d0) then
        nn(7) = n(7) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==8) then
      pdj(1) =  &
          -k(42)*n(idx_E)  &
          +2.d0*k(42)*n(idx_E)  &
          +k(184)*n(idx_Hk)  &
          +k(192)*n(idx_Ok)  &
          -k(195)*n(idx_E)  &
          +k(217)  &
          +k(255)
      pdj(2) =  &
          -k(184)*n(idx_Hk)
      pdj(3) =  &
          +k(195)*n(idx_E)
      pdj(4) =  &
          -k(192)*n(idx_Ok)
      pdj(5) =  &
          +k(47)*n(idx_Hj)  &
          +k(56)*n(idx_H2)  &
          +k(59)*n(idx_CH)  &
          +k(74)*n(idx_OH)  &
          +k(75)*n(idx_OH)  &
          +k(87)*n(idx_H2j)  &
          +k(89)*n(idx_H3j)  &
          -k(196)*n(idx_H)
      pdj(6) =  &
          +k(49)*n(idx_HEj)  &
          +k(50)*n(idx_HEj)  &
          +k(51)*n(idx_HEj)
      pdj(7) =  &
          -k(56)*n(idx_H2)  &
          +k(88)*n(idx_H3j)  &
          +k(117)*n(idx_H3Oj)  &
          -k(197)*n(idx_H2)
      pdj(8) =  &
          -k(42)*n(idx_E)  &
          -k(47)*n(idx_Hj)  &
          -k(49)*n(idx_HEj)  &
          -k(50)*n(idx_HEj)  &
          -k(51)*n(idx_HEj)  &
          -k(56)*n(idx_H2)  &
          -k(59)*n(idx_CH)  &
          -k(74)*n(idx_OH)  &
          -k(75)*n(idx_OH)  &
          -k(82)*n(idx_O2)  &
          -k(83)*n(idx_O2)  &
          -k(87)*n(idx_H2j)  &
          -k(88)*n(idx_H3j)  &
          -k(89)*n(idx_H3j)  &
          -k(117)*n(idx_H3Oj)  &
          -k(121)*n(idx_O2j)  &
          -k(122)*n(idx_O2j)  &
          -k(127)*n(idx_HCOj)  &
          -k(184)*n(idx_Hk)  &
          -k(192)*n(idx_Ok)  &
          -k(195)*n(idx_E)  &
          -k(196)*n(idx_H)  &
          -k(197)*n(idx_H2)  &
          -4.d0*k(198)*n(idx_C)  &
          -k(199)*n(idx_O)  &
          -k(217)  &
          -k(255)  &
          -4.d0*k(265)*n(idx_C)  &
          -4.d0*k(266)*n(idx_C)  &
          -k(267)*n(idx_O)  &
          -k(268)*n(idx_O)  &
          -k(271)*n(idx_Oj)  &
          -k(272)*n(idx_Oj)
      pdj(9) =  &
          +k(82)*n(idx_O2)  &
          +k(83)*n(idx_O2)  &
          +k(121)*n(idx_O2j)  &
          -k(199)*n(idx_O)  &
          -k(267)*n(idx_O)  &
          -k(268)*n(idx_O)
      pdj(10) =  &
          -k(74)*n(idx_OH)  &
          -k(75)*n(idx_OH)
      pdj(11) =  &
          +k(74)*n(idx_OH)  &
          +k(75)*n(idx_OH)  &
          +k(82)*n(idx_O2)  &
          +k(83)*n(idx_O2)  &
          +k(127)*n(idx_HCOj)  &
          +k(192)*n(idx_Ok)  &
          +k(199)*n(idx_O)  &
          +k(267)*n(idx_O)  &
          +k(268)*n(idx_O)
      pdj(12) =  &
          +k(56)*n(idx_H2)  &
          -k(59)*n(idx_CH)  &
          +k(184)*n(idx_Hk)  &
          +k(196)*n(idx_H)
      pdj(13) =  &
          +k(197)*n(idx_H2)
      pdj(14) =  &
          +k(59)*n(idx_CH)  &
          +2.d0*k(198)*n(idx_C)  &
          +2.d0*k(265)*n(idx_C)  &
          +2.d0*k(266)*n(idx_C)
      pdj(17) =  &
          -k(82)*n(idx_O2)  &
          -k(83)*n(idx_O2)  &
          +k(122)*n(idx_O2j)
      pdj(20) =  &
          -k(47)*n(idx_Hj)
      pdj(21) =  &
          -k(49)*n(idx_HEj)  &
          -k(50)*n(idx_HEj)  &
          -k(51)*n(idx_HEj)
      pdj(22) =  &
          -k(87)*n(idx_H2j)
      pdj(23) =  &
          +k(42)*n(idx_E)  &
          +k(47)*n(idx_Hj)  &
          +k(49)*n(idx_HEj)  &
          +k(50)*n(idx_HEj)  &
          +k(51)*n(idx_HEj)  &
          +k(122)*n(idx_O2j)  &
          +k(217)  &
          +k(255)
      pdj(24) =  &
          -k(271)*n(idx_Oj)  &
          -k(272)*n(idx_Oj)
      pdj(26) =  &
          +k(117)*n(idx_H3Oj)  &
          -k(127)*n(idx_HCOj)
      pdj(27) =  &
          -k(88)*n(idx_H3j)  &
          -k(89)*n(idx_H3j)
      pdj(28) =  &
          +k(87)*n(idx_H2j)  &
          +k(88)*n(idx_H3j)  &
          +k(127)*n(idx_HCOj)
      pdj(29) =  &
          +k(89)*n(idx_H3j)
      pdj(30) =  &
          +k(121)*n(idx_O2j)  &
          +k(271)*n(idx_Oj)  &
          +k(272)*n(idx_Oj)
      pdj(34) =  &
          -k(117)*n(idx_H3Oj)
      pdj(35) =  &
          -k(121)*n(idx_O2j)  &
          -k(122)*n(idx_O2j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(8)*1d-3
      if(dnn>0.d0) then
        nn(8) = n(8) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==9) then
      pdj(1) =  &
          -k(43)*n(idx_E)  &
          +2.d0*k(43)*n(idx_E)  &
          +k(61)*n(idx_CH)  &
          +k(185)*n(idx_Hk)  &
          +k(189)*n(idx_Ck)  &
          -k(204)*n(idx_E)  &
          +k(248)
      pdj(2) =  &
          -k(185)*n(idx_Hk)
      pdj(3) =  &
          -k(189)*n(idx_Ck)
      pdj(4) =  &
          +k(204)*n(idx_E)
      pdj(5) =  &
          +k(45)*n(idx_Hj)  &
          +k(60)*n(idx_CH)  &
          +2.d0*k(64)*n(idx_CH2)  &
          +k(66)*n(idx_CH2)  &
          +k(70)*n(idx_H2)  &
          +k(76)*n(idx_OH)  &
          +k(77)*n(idx_OH)  &
          +k(93)*n(idx_CHj)  &
          +k(96)*n(idx_CH2j)  &
          +k(102)*n(idx_H2j)  &
          +k(104)*n(idx_H3j)  &
          -k(205)*n(idx_H)  &
          -k(273)*n(idx_H)
      pdj(6) =  &
          +k(46)*n(idx_HEj)
      pdj(7) =  &
          +k(65)*n(idx_CH2)  &
          -k(70)*n(idx_H2)  &
          +k(98)*n(idx_CH3j)  &
          +k(99)*n(idx_CH3j)  &
          +k(103)*n(idx_H3j)
      pdj(8) =  &
          +k(62)*n(idx_CH)  &
          +k(68)*n(idx_C2)  &
          +k(69)*n(idx_C2)  &
          -k(199)*n(idx_C)  &
          -k(267)*n(idx_C)  &
          -k(268)*n(idx_C)
      pdj(9) =  &
          -k(43)*n(idx_E)  &
          -k(45)*n(idx_Hj)  &
          -k(46)*n(idx_HEj)  &
          -k(60)*n(idx_CH)  &
          -k(61)*n(idx_CH)  &
          -k(62)*n(idx_CH)  &
          -k(64)*n(idx_CH2)  &
          -k(65)*n(idx_CH2)  &
          -k(66)*n(idx_CH2)  &
          -k(67)*n(idx_CH2)  &
          -k(68)*n(idx_C2)  &
          -k(69)*n(idx_C2)  &
          -k(70)*n(idx_H2)  &
          -k(76)*n(idx_OH)  &
          -k(77)*n(idx_OH)  &
          -k(93)*n(idx_CHj)  &
          -k(96)*n(idx_CH2j)  &
          -k(98)*n(idx_CH3j)  &
          -k(99)*n(idx_CH3j)  &
          -k(102)*n(idx_H2j)  &
          -k(103)*n(idx_H3j)  &
          -k(104)*n(idx_H3j)  &
          -k(185)*n(idx_Hk)  &
          -k(189)*n(idx_Ck)  &
          -k(199)*n(idx_C)  &
          -k(202)*n(idx_Cj)  &
          -k(203)*n(idx_Cj)  &
          -k(204)*n(idx_E)  &
          -k(205)*n(idx_H)  &
          -4.d0*k(206)*n(idx_O)  &
          -k(248)  &
          -k(267)*n(idx_C)  &
          -k(268)*n(idx_C)  &
          -k(269)*n(idx_Cj)  &
          -k(270)*n(idx_Cj)  &
          -k(273)*n(idx_H)  &
          -4.d0*k(275)*n(idx_O)
      pdj(10) =  &
          +k(62)*n(idx_CH)  &
          +k(67)*n(idx_CH2)  &
          +k(70)*n(idx_H2)  &
          -k(76)*n(idx_OH)  &
          -k(77)*n(idx_OH)  &
          +k(185)*n(idx_Hk)  &
          +k(205)*n(idx_H)  &
          +k(273)*n(idx_H)
      pdj(11) =  &
          +k(60)*n(idx_CH)  &
          +k(64)*n(idx_CH2)  &
          +k(65)*n(idx_CH2)  &
          +k(68)*n(idx_C2)  &
          +k(69)*n(idx_C2)  &
          +k(189)*n(idx_Ck)  &
          +k(199)*n(idx_C)  &
          +k(267)*n(idx_C)  &
          +k(268)*n(idx_C)
      pdj(12) =  &
          -k(60)*n(idx_CH)  &
          -k(61)*n(idx_CH)  &
          -k(62)*n(idx_CH)  &
          +k(67)*n(idx_CH2)
      pdj(13) =  &
          -k(64)*n(idx_CH2)  &
          -k(65)*n(idx_CH2)  &
          -k(66)*n(idx_CH2)  &
          -k(67)*n(idx_CH2)
      pdj(14) =  &
          -k(68)*n(idx_C2)  &
          -k(69)*n(idx_C2)
      pdj(15) =  &
          +k(66)*n(idx_CH2)
      pdj(17) =  &
          +k(76)*n(idx_OH)  &
          +k(77)*n(idx_OH)  &
          +2.d0*k(206)*n(idx_O)  &
          +2.d0*k(275)*n(idx_O)
      pdj(20) =  &
          -k(45)*n(idx_Hj)
      pdj(21) =  &
          -k(46)*n(idx_HEj)
      pdj(22) =  &
          -k(102)*n(idx_H2j)
      pdj(23) =  &
          -k(202)*n(idx_Cj)  &
          -k(203)*n(idx_Cj)  &
          -k(269)*n(idx_Cj)  &
          -k(270)*n(idx_Cj)
      pdj(24) =  &
          +k(43)*n(idx_E)  &
          +k(45)*n(idx_Hj)  &
          +k(46)*n(idx_HEj)  &
          +k(248)
      pdj(25) =  &
          +k(98)*n(idx_CH3j)
      pdj(26) =  &
          +k(61)*n(idx_CH)  &
          +k(96)*n(idx_CH2j)  &
          +k(99)*n(idx_CH3j)
      pdj(27) =  &
          -k(103)*n(idx_H3j)  &
          -k(104)*n(idx_H3j)
      pdj(28) =  &
          -k(93)*n(idx_CHj)
      pdj(29) =  &
          -k(96)*n(idx_CH2j)
      pdj(30) =  &
          +k(93)*n(idx_CHj)  &
          +k(202)*n(idx_Cj)  &
          +k(203)*n(idx_Cj)  &
          +k(269)*n(idx_Cj)  &
          +k(270)*n(idx_Cj)
      pdj(31) =  &
          -k(98)*n(idx_CH3j)  &
          -k(99)*n(idx_CH3j)
      pdj(32) =  &
          +k(102)*n(idx_H2j)  &
          +k(103)*n(idx_H3j)
      pdj(33) =  &
          +k(104)*n(idx_H3j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(9)*1d-3
      if(dnn>0.d0) then
        nn(9) = n(9) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==10) then
      pdj(1) =  &
          +k(186)*n(idx_Hk)  &
          +k(230)
      pdj(2) =  &
          -k(186)*n(idx_Hk)
      pdj(5) =  &
          -k(52)*n(idx_H)  &
          +2.d0*k(52)*n(idx_H)  &
          -k(71)*n(idx_H)  &
          -k(72)*n(idx_H)  &
          +k(73)*n(idx_H2)  &
          +k(74)*n(idx_C)  &
          +k(75)*n(idx_C)  &
          +k(76)*n(idx_O)  &
          +k(77)*n(idx_O)  &
          +k(107)*n(idx_Cj)  &
          +k(108)*n(idx_Cj)  &
          +k(141)*n(idx_Hj)  &
          +k(142)*n(idx_Hj)  &
          +k(143)*n(idx_HEj)  &
          +k(144)*n(idx_HEj)  &
          -k(207)*n(idx_H)  &
          +k(229)  &
          +k(259)  &
          -k(274)*n(idx_H)
      pdj(6) =  &
          +k(143)*n(idx_HEj)  &
          +k(144)*n(idx_HEj)
      pdj(7) =  &
          +k(71)*n(idx_H)  &
          +k(72)*n(idx_H)  &
          -k(73)*n(idx_H2)  &
          +k(105)*n(idx_H3j)  &
          +k(106)*n(idx_H3j)
      pdj(8) =  &
          -k(74)*n(idx_C)  &
          -k(75)*n(idx_C)
      pdj(9) =  &
          +k(52)*n(idx_H)  &
          +k(71)*n(idx_H)  &
          +k(72)*n(idx_H)  &
          -k(76)*n(idx_O)  &
          -k(77)*n(idx_O)  &
          +2.d0*k(78)*n(idx_OH)  &
          +k(229)  &
          +k(259)
      pdj(10) =  &
          -k(52)*n(idx_H)  &
          -k(71)*n(idx_H)  &
          -k(72)*n(idx_H)  &
          -k(73)*n(idx_H2)  &
          -k(74)*n(idx_C)  &
          -k(75)*n(idx_C)  &
          -k(76)*n(idx_O)  &
          -k(77)*n(idx_O)  &
          -4.d0*k(78)*n(idx_OH)  &
          -k(105)*n(idx_H3j)  &
          -k(106)*n(idx_H3j)  &
          -k(107)*n(idx_Cj)  &
          -k(108)*n(idx_Cj)  &
          -k(141)*n(idx_Hj)  &
          -k(142)*n(idx_Hj)  &
          -k(143)*n(idx_HEj)  &
          -k(144)*n(idx_HEj)  &
          -k(186)*n(idx_Hk)  &
          -k(207)*n(idx_H)  &
          -k(229)  &
          -k(230)  &
          -k(259)  &
          -k(274)*n(idx_H)
      pdj(11) =  &
          +k(74)*n(idx_C)  &
          +k(75)*n(idx_C)
      pdj(16) =  &
          +k(73)*n(idx_H2)  &
          +2.d0*k(78)*n(idx_OH)  &
          +k(186)*n(idx_Hk)  &
          +k(207)*n(idx_H)  &
          +k(274)*n(idx_H)
      pdj(17) =  &
          +k(76)*n(idx_O)  &
          +k(77)*n(idx_O)
      pdj(20) =  &
          -k(141)*n(idx_Hj)  &
          -k(142)*n(idx_Hj)
      pdj(21) =  &
          -k(143)*n(idx_HEj)  &
          -k(144)*n(idx_HEj)
      pdj(23) =  &
          -k(107)*n(idx_Cj)  &
          -k(108)*n(idx_Cj)
      pdj(24) =  &
          +k(143)*n(idx_HEj)  &
          +k(144)*n(idx_HEj)
      pdj(27) =  &
          -k(105)*n(idx_H3j)  &
          -k(106)*n(idx_H3j)
      pdj(30) =  &
          +k(107)*n(idx_Cj)  &
          +k(108)*n(idx_Cj)
      pdj(32) =  &
          +k(141)*n(idx_Hj)  &
          +k(142)*n(idx_Hj)  &
          +k(230)
      pdj(33) =  &
          +k(105)*n(idx_H3j)  &
          +k(106)*n(idx_H3j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(10)*1d-3
      if(dnn>0.d0) then
        nn(10) = n(10) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==11) then
      pdj(1) =  &
          +k(250)
      pdj(5) =  &
          -k(84)*n(idx_H)
      pdj(6) =  &
          +k(156)*n(idx_HEj)  &
          +k(157)*n(idx_HEj)
      pdj(7) =  &
          +k(123)*n(idx_H3j)  &
          +k(124)*n(idx_H3j)  &
          +k(125)*n(idx_H3j)  &
          +k(126)*n(idx_H3j)
      pdj(8) =  &
          +k(84)*n(idx_H)  &
          +k(157)*n(idx_HEj)  &
          +k(236)  &
          +k(249)
      pdj(9) =  &
          +k(156)*n(idx_HEj)  &
          +k(236)  &
          +k(249)
      pdj(10) =  &
          +k(84)*n(idx_H)
      pdj(11) =  &
          -k(54)*n(idx_HOCj)  &
          +k(54)*n(idx_HOCj)  &
          -k(55)*n(idx_HOCj)  &
          +k(55)*n(idx_HOCj)  &
          -k(84)*n(idx_H)  &
          -k(123)*n(idx_H3j)  &
          -k(124)*n(idx_H3j)  &
          -k(125)*n(idx_H3j)  &
          -k(126)*n(idx_H3j)  &
          -k(156)*n(idx_HEj)  &
          -k(157)*n(idx_HEj)  &
          -k(236)  &
          -k(249)  &
          -k(250)
      pdj(21) =  &
          -k(156)*n(idx_HEj)  &
          -k(157)*n(idx_HEj)
      pdj(23) =  &
          +k(156)*n(idx_HEj)
      pdj(24) =  &
          +k(157)*n(idx_HEj)
      pdj(25) =  &
          -k(54)*n(idx_HOCj)  &
          -k(55)*n(idx_HOCj)  &
          +k(125)*n(idx_H3j)  &
          +k(126)*n(idx_H3j)
      pdj(26) =  &
          +k(54)*n(idx_HOCj)  &
          +k(55)*n(idx_HOCj)  &
          +k(123)*n(idx_H3j)  &
          +k(124)*n(idx_H3j)
      pdj(27) =  &
          -k(123)*n(idx_H3j)  &
          -k(124)*n(idx_H3j)  &
          -k(125)*n(idx_H3j)  &
          -k(126)*n(idx_H3j)
      pdj(30) =  &
          +k(250)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(11)*1d-3
      if(dnn>0.d0) then
        nn(11) = n(11) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==12) then
      pdj(1) =  &
          +k(61)*n(idx_O)  &
          +k(220)
      pdj(5) =  &
          -k(57)*n(idx_H)  &
          +k(58)*n(idx_H2)  &
          +k(59)*n(idx_C)  &
          +k(60)*n(idx_O)  &
          +k(130)*n(idx_Hj)  &
          +k(131)*n(idx_Hj)  &
          +k(219)  &
          +k(256)
      pdj(7) =  &
          +k(57)*n(idx_H)  &
          -k(58)*n(idx_H2)
      pdj(8) =  &
          +k(57)*n(idx_H)  &
          -k(59)*n(idx_C)  &
          +k(62)*n(idx_O)  &
          +k(219)  &
          +k(256)
      pdj(9) =  &
          -k(60)*n(idx_O)  &
          -k(61)*n(idx_O)  &
          -k(62)*n(idx_O)
      pdj(10) =  &
          +k(62)*n(idx_O)
      pdj(11) =  &
          +k(60)*n(idx_O)
      pdj(12) =  &
          -k(57)*n(idx_H)  &
          -k(58)*n(idx_H2)  &
          -k(59)*n(idx_C)  &
          -k(60)*n(idx_O)  &
          -k(61)*n(idx_O)  &
          -k(62)*n(idx_O)  &
          -k(130)*n(idx_Hj)  &
          -k(131)*n(idx_Hj)  &
          -k(219)  &
          -k(220)  &
          -k(256)
      pdj(13) =  &
          +k(58)*n(idx_H2)
      pdj(14) =  &
          +k(59)*n(idx_C)
      pdj(20) =  &
          -k(130)*n(idx_Hj)  &
          -k(131)*n(idx_Hj)
      pdj(26) =  &
          +k(61)*n(idx_O)
      pdj(28) =  &
          +k(130)*n(idx_Hj)  &
          +k(131)*n(idx_Hj)  &
          +k(220)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(12)*1d-3
      if(dnn>0.d0) then
        nn(12) = n(12) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==13) then
      pdj(1) =  &
          +k(223)  &
          +k(260)
      pdj(5) =  &
          -k(63)*n(idx_H)  &
          +2.d0*k(64)*n(idx_O)  &
          +k(66)*n(idx_O)  &
          +k(134)*n(idx_Hj)  &
          +k(135)*n(idx_Hj)  &
          +k(138)*n(idx_HEj)  &
          +k(139)*n(idx_HEj)  &
          +k(222)
      pdj(6) =  &
          +k(136)*n(idx_HEj)  &
          +k(137)*n(idx_HEj)  &
          +k(138)*n(idx_HEj)  &
          +k(139)*n(idx_HEj)
      pdj(7) =  &
          +k(63)*n(idx_H)  &
          +k(65)*n(idx_O)  &
          +k(132)*n(idx_Hj)  &
          +k(133)*n(idx_Hj)  &
          +k(136)*n(idx_HEj)  &
          +k(137)*n(idx_HEj)
      pdj(9) =  &
          -k(64)*n(idx_O)  &
          -k(65)*n(idx_O)  &
          -k(66)*n(idx_O)  &
          -k(67)*n(idx_O)
      pdj(10) =  &
          +k(67)*n(idx_O)
      pdj(11) =  &
          +k(64)*n(idx_O)  &
          +k(65)*n(idx_O)
      pdj(12) =  &
          +k(63)*n(idx_H)  &
          +k(67)*n(idx_O)  &
          +k(222)
      pdj(13) =  &
          -k(63)*n(idx_H)  &
          -k(64)*n(idx_O)  &
          -k(65)*n(idx_O)  &
          -k(66)*n(idx_O)  &
          -k(67)*n(idx_O)  &
          -k(132)*n(idx_Hj)  &
          -k(133)*n(idx_Hj)  &
          -k(134)*n(idx_Hj)  &
          -k(135)*n(idx_Hj)  &
          -k(136)*n(idx_HEj)  &
          -k(137)*n(idx_HEj)  &
          -k(138)*n(idx_HEj)  &
          -k(139)*n(idx_HEj)  &
          -k(222)  &
          -k(223)  &
          -k(260)
      pdj(15) =  &
          +k(66)*n(idx_O)
      pdj(20) =  &
          -k(132)*n(idx_Hj)  &
          -k(133)*n(idx_Hj)  &
          -k(134)*n(idx_Hj)  &
          -k(135)*n(idx_Hj)
      pdj(21) =  &
          -k(136)*n(idx_HEj)  &
          -k(137)*n(idx_HEj)  &
          -k(138)*n(idx_HEj)  &
          -k(139)*n(idx_HEj)
      pdj(23) =  &
          +k(136)*n(idx_HEj)  &
          +k(137)*n(idx_HEj)
      pdj(28) =  &
          +k(132)*n(idx_Hj)  &
          +k(133)*n(idx_Hj)  &
          +k(138)*n(idx_HEj)  &
          +k(139)*n(idx_HEj)
      pdj(29) =  &
          +k(134)*n(idx_Hj)  &
          +k(135)*n(idx_Hj)  &
          +k(223)  &
          +k(260)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(13)*1d-3
      if(dnn>0.d0) then
        nn(13) = n(13) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==14) then
      pdj(6) =  &
          +k(140)*n(idx_HEj)
      pdj(8) =  &
          +k(68)*n(idx_O)  &
          +k(69)*n(idx_O)  &
          +k(100)*n(idx_Oj)  &
          +k(140)*n(idx_HEj)  &
          +2.d0*k(227)  &
          +2.d0*k(251)
      pdj(9) =  &
          -k(68)*n(idx_O)  &
          -k(69)*n(idx_O)
      pdj(11) =  &
          +k(68)*n(idx_O)  &
          +k(69)*n(idx_O)
      pdj(14) =  &
          -k(68)*n(idx_O)  &
          -k(69)*n(idx_O)  &
          -k(100)*n(idx_Oj)  &
          -k(140)*n(idx_HEj)  &
          -k(227)  &
          -k(251)
      pdj(21) =  &
          -k(140)*n(idx_HEj)
      pdj(23) =  &
          +k(140)*n(idx_HEj)
      pdj(24) =  &
          -k(100)*n(idx_Oj)
      pdj(30) =  &
          +k(100)*n(idx_Oj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(14)*1d-3
      if(dnn>0.d0) then
        nn(14) = n(14) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==15) then
      pdj(1) =  &
          +k(263)
      pdj(5) =  &
          +k(262)
      pdj(11) =  &
          +k(262)
      pdj(15) =  &
          -k(262)  &
          -k(263)
      pdj(26) =  &
          +k(263)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(15)*1d-3
      if(dnn>0.d0) then
        nn(15) = n(15) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==16) then
      pdj(1) =  &
          +k(233)
      pdj(5) =  &
          -k(79)*n(idx_H)  &
          +k(113)*n(idx_Cj)  &
          +k(114)*n(idx_Cj)  &
          +k(115)*n(idx_Cj)  &
          +k(145)*n(idx_Hj)  &
          +k(146)*n(idx_Hj)  &
          +k(149)*n(idx_HEj)  &
          +k(150)*n(idx_HEj)  &
          +k(232)  &
          +k(261)
      pdj(6) =  &
          +k(147)*n(idx_HEj)  &
          +k(148)*n(idx_HEj)  &
          +k(149)*n(idx_HEj)  &
          +k(150)*n(idx_HEj)  &
          +k(151)*n(idx_HEj)  &
          +k(152)*n(idx_HEj)
      pdj(7) =  &
          +k(79)*n(idx_H)  &
          +k(111)*n(idx_H3j)  &
          +k(112)*n(idx_H3j)
      pdj(8) =  &
          +k(116)*n(idx_Cj)
      pdj(10) =  &
          +k(79)*n(idx_H)  &
          +k(147)*n(idx_HEj)  &
          +k(148)*n(idx_HEj)  &
          +k(232)  &
          +k(261)
      pdj(11) =  &
          +k(128)*n(idx_HCOj)  &
          +k(129)*n(idx_HCOj)
      pdj(16) =  &
          -k(79)*n(idx_H)  &
          -k(111)*n(idx_H3j)  &
          -k(112)*n(idx_H3j)  &
          -k(113)*n(idx_Cj)  &
          -k(114)*n(idx_Cj)  &
          -k(115)*n(idx_Cj)  &
          -k(116)*n(idx_Cj)  &
          -k(128)*n(idx_HCOj)  &
          -k(129)*n(idx_HCOj)  &
          -k(145)*n(idx_Hj)  &
          -k(146)*n(idx_Hj)  &
          -k(147)*n(idx_HEj)  &
          -k(148)*n(idx_HEj)  &
          -k(149)*n(idx_HEj)  &
          -k(150)*n(idx_HEj)  &
          -k(151)*n(idx_HEj)  &
          -k(152)*n(idx_HEj)  &
          -k(232)  &
          -k(233)  &
          -k(261)
      pdj(20) =  &
          -k(145)*n(idx_Hj)  &
          -k(146)*n(idx_Hj)  &
          +k(147)*n(idx_HEj)  &
          +k(148)*n(idx_HEj)
      pdj(21) =  &
          -k(147)*n(idx_HEj)  &
          -k(148)*n(idx_HEj)  &
          -k(149)*n(idx_HEj)  &
          -k(150)*n(idx_HEj)  &
          -k(151)*n(idx_HEj)  &
          -k(152)*n(idx_HEj)
      pdj(23) =  &
          -k(113)*n(idx_Cj)  &
          -k(114)*n(idx_Cj)  &
          -k(115)*n(idx_Cj)  &
          -k(116)*n(idx_Cj)
      pdj(25) =  &
          +k(113)*n(idx_Cj)
      pdj(26) =  &
          +k(114)*n(idx_Cj)  &
          +k(115)*n(idx_Cj)  &
          -k(128)*n(idx_HCOj)  &
          -k(129)*n(idx_HCOj)
      pdj(27) =  &
          -k(111)*n(idx_H3j)  &
          -k(112)*n(idx_H3j)
      pdj(32) =  &
          +k(149)*n(idx_HEj)  &
          +k(150)*n(idx_HEj)
      pdj(33) =  &
          +k(116)*n(idx_Cj)  &
          +k(145)*n(idx_Hj)  &
          +k(146)*n(idx_Hj)  &
          +k(151)*n(idx_HEj)  &
          +k(152)*n(idx_HEj)  &
          +k(233)
      pdj(34) =  &
          +k(111)*n(idx_H3j)  &
          +k(112)*n(idx_H3j)  &
          +k(128)*n(idx_HCOj)  &
          +k(129)*n(idx_HCOj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(16)*1d-3
      if(dnn>0.d0) then
        nn(16) = n(16) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==17) then
      pdj(1) =  &
          +k(234)  &
          +k(258)
      pdj(5) =  &
          -k(80)*n(idx_H)  &
          +k(153)*n(idx_Hj)
      pdj(6) =  &
          +k(154)*n(idx_HEj)  &
          +k(155)*n(idx_HEj)
      pdj(7) =  &
          -k(81)*n(idx_H2)
      pdj(8) =  &
          -k(82)*n(idx_C)  &
          -k(83)*n(idx_C)
      pdj(9) =  &
          +k(80)*n(idx_H)  &
          +k(82)*n(idx_C)  &
          +k(83)*n(idx_C)  &
          +k(118)*n(idx_Cj)  &
          +k(155)*n(idx_HEj)  &
          +2.d0*k(235)  &
          +2.d0*k(257)
      pdj(10) =  &
          +k(80)*n(idx_H)  &
          +2.d0*k(81)*n(idx_H2)  &
          +k(120)*n(idx_CH2j)
      pdj(11) =  &
          +k(82)*n(idx_C)  &
          +k(83)*n(idx_C)  &
          +k(119)*n(idx_Cj)
      pdj(17) =  &
          -k(80)*n(idx_H)  &
          -k(81)*n(idx_H2)  &
          -k(82)*n(idx_C)  &
          -k(83)*n(idx_C)  &
          -k(118)*n(idx_Cj)  &
          -k(119)*n(idx_Cj)  &
          -k(120)*n(idx_CH2j)  &
          -k(153)*n(idx_Hj)  &
          -k(154)*n(idx_HEj)  &
          -k(155)*n(idx_HEj)  &
          -k(234)  &
          -k(235)  &
          -k(257)  &
          -k(258)
      pdj(20) =  &
          -k(153)*n(idx_Hj)
      pdj(21) =  &
          -k(154)*n(idx_HEj)  &
          -k(155)*n(idx_HEj)
      pdj(23) =  &
          -k(118)*n(idx_Cj)  &
          -k(119)*n(idx_Cj)
      pdj(24) =  &
          +k(119)*n(idx_Cj)  &
          +k(155)*n(idx_HEj)
      pdj(26) =  &
          +k(120)*n(idx_CH2j)
      pdj(29) =  &
          -k(120)*n(idx_CH2j)
      pdj(30) =  &
          +k(118)*n(idx_Cj)
      pdj(35) =  &
          +k(153)*n(idx_Hj)  &
          +k(154)*n(idx_HEj)  &
          +k(234)  &
          +k(258)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(17)*1d-3
      if(dnn>0.d0) then
        nn(17) = n(17) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==18) then
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(18)*1d-3
      if(dnn>0.d0) then
        nn(18) = n(18) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==19) then
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(19)*1d-3
      if(dnn>0.d0) then
        nn(19) = n(19) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==20) then
      pdj(1) =  &
          -k(2)*n(idx_E)  &
          -k(3)*n(idx_E)  &
          +k(29)*n(idx_Hk)
      pdj(2) =  &
          -k(28)*n(idx_Hk)  &
          -k(29)*n(idx_Hk)
      pdj(3) =  &
          -k(159)*n(idx_Ck)
      pdj(4) =  &
          -k(160)*n(idx_Ok)
      pdj(5) =  &
          +k(2)*n(idx_E)  &
          +k(3)*n(idx_E)  &
          +k(9)*n(idx_HE)  &
          +k(10)*n(idx_HE)  &
          -k(18)*n(idx_H)  &
          -k(19)*n(idx_H)  &
          +k(21)*n(idx_H2)  &
          +k(22)*n(idx_H2)  &
          +2.d0*k(28)*n(idx_Hk)  &
          +k(45)*n(idx_O)  &
          +k(47)*n(idx_C)  &
          +k(130)*n(idx_CH)  &
          +k(131)*n(idx_CH)  &
          +k(134)*n(idx_CH2)  &
          +k(135)*n(idx_CH2)  &
          +k(141)*n(idx_OH)  &
          +k(142)*n(idx_OH)  &
          +k(145)*n(idx_H2O)  &
          +k(146)*n(idx_H2O)  &
          +k(153)*n(idx_O2)  &
          +k(159)*n(idx_Ck)  &
          +k(160)*n(idx_Ok)  &
          +2.d0*k(193)*n(idx_H2)
      pdj(6) =  &
          -k(9)*n(idx_HE)  &
          -k(10)*n(idx_HE)
      pdj(7) =  &
          -k(21)*n(idx_H2)  &
          -k(22)*n(idx_H2)  &
          +k(132)*n(idx_CH2)  &
          +k(133)*n(idx_CH2)  &
          -k(193)*n(idx_H2)  &
          -k(194)*n(idx_H2)
      pdj(8) =  &
          -k(47)*n(idx_C)  &
          +k(159)*n(idx_Ck)
      pdj(9) =  &
          -k(45)*n(idx_O)  &
          +k(160)*n(idx_Ok)
      pdj(10) =  &
          -k(141)*n(idx_OH)  &
          -k(142)*n(idx_OH)
      pdj(12) =  &
          -k(130)*n(idx_CH)  &
          -k(131)*n(idx_CH)
      pdj(13) =  &
          -k(132)*n(idx_CH2)  &
          -k(133)*n(idx_CH2)  &
          -k(134)*n(idx_CH2)  &
          -k(135)*n(idx_CH2)
      pdj(16) =  &
          -k(145)*n(idx_H2O)  &
          -k(146)*n(idx_H2O)
      pdj(17) =  &
          -k(153)*n(idx_O2)
      pdj(20) =  &
          -k(2)*n(idx_E)  &
          -k(3)*n(idx_E)  &
          -k(9)*n(idx_HE)  &
          -k(10)*n(idx_HE)  &
          -k(18)*n(idx_H)  &
          -k(19)*n(idx_H)  &
          -k(21)*n(idx_H2)  &
          -k(22)*n(idx_H2)  &
          -k(28)*n(idx_Hk)  &
          -k(29)*n(idx_Hk)  &
          -k(45)*n(idx_O)  &
          -k(47)*n(idx_C)  &
          -k(130)*n(idx_CH)  &
          -k(131)*n(idx_CH)  &
          -k(132)*n(idx_CH2)  &
          -k(133)*n(idx_CH2)  &
          -k(134)*n(idx_CH2)  &
          -k(135)*n(idx_CH2)  &
          -k(141)*n(idx_OH)  &
          -k(142)*n(idx_OH)  &
          -k(145)*n(idx_H2O)  &
          -k(146)*n(idx_H2O)  &
          -k(153)*n(idx_O2)  &
          -k(159)*n(idx_Ck)  &
          -k(160)*n(idx_Ok)  &
          -k(193)*n(idx_H2)  &
          +k(193)*n(idx_H2)  &
          -k(194)*n(idx_H2)
      pdj(21) =  &
          +k(9)*n(idx_HE)  &
          +k(10)*n(idx_HE)
      pdj(22) =  &
          +k(18)*n(idx_H)  &
          +k(19)*n(idx_H)  &
          +k(21)*n(idx_H2)  &
          +k(22)*n(idx_H2)  &
          +k(29)*n(idx_Hk)
      pdj(23) =  &
          +k(47)*n(idx_C)
      pdj(24) =  &
          +k(45)*n(idx_O)
      pdj(27) =  &
          +k(194)*n(idx_H2)
      pdj(28) =  &
          +k(130)*n(idx_CH)  &
          +k(131)*n(idx_CH)  &
          +k(132)*n(idx_CH2)  &
          +k(133)*n(idx_CH2)
      pdj(29) =  &
          +k(134)*n(idx_CH2)  &
          +k(135)*n(idx_CH2)
      pdj(32) =  &
          +k(141)*n(idx_OH)  &
          +k(142)*n(idx_OH)
      pdj(33) =  &
          +k(145)*n(idx_H2O)  &
          +k(146)*n(idx_H2O)
      pdj(35) =  &
          +k(153)*n(idx_O2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(20)*1d-3
      if(dnn>0.d0) then
        nn(20) = n(20) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==21) then
      pdj(1) =  &
          -k(5)*n(idx_E)  &
          -k(6)*n(idx_E)  &
          -k(7)*n(idx_E)  &
          +2.d0*k(7)*n(idx_E)
      pdj(2) =  &
          -k(161)*n(idx_Hk)
      pdj(5) =  &
          -k(8)*n(idx_H)  &
          +k(13)*n(idx_H2)  &
          +2.d0*k(14)*n(idx_H2)  &
          +k(138)*n(idx_CH2)  &
          +k(139)*n(idx_CH2)  &
          +k(143)*n(idx_OH)  &
          +k(144)*n(idx_OH)  &
          +k(149)*n(idx_H2O)  &
          +k(150)*n(idx_H2O)  &
          +k(161)*n(idx_Hk)
      pdj(6) =  &
          +k(5)*n(idx_E)  &
          +k(6)*n(idx_E)  &
          +k(8)*n(idx_H)  &
          +k(12)*n(idx_H2)  &
          +k(13)*n(idx_H2)  &
          +k(46)*n(idx_O)  &
          +k(49)*n(idx_C)  &
          +k(50)*n(idx_C)  &
          +k(51)*n(idx_C)  &
          +k(136)*n(idx_CH2)  &
          +k(137)*n(idx_CH2)  &
          +k(138)*n(idx_CH2)  &
          +k(139)*n(idx_CH2)  &
          +k(140)*n(idx_C2)  &
          +k(143)*n(idx_OH)  &
          +k(144)*n(idx_OH)  &
          +k(147)*n(idx_H2O)  &
          +k(148)*n(idx_H2O)  &
          +k(149)*n(idx_H2O)  &
          +k(150)*n(idx_H2O)  &
          +k(151)*n(idx_H2O)  &
          +k(152)*n(idx_H2O)  &
          +k(154)*n(idx_O2)  &
          +k(155)*n(idx_O2)  &
          +k(156)*n(idx_CO)  &
          +k(157)*n(idx_CO)  &
          +k(161)*n(idx_Hk)
      pdj(7) =  &
          -k(12)*n(idx_H2)  &
          -k(13)*n(idx_H2)  &
          -k(14)*n(idx_H2)  &
          +k(136)*n(idx_CH2)  &
          +k(137)*n(idx_CH2)
      pdj(8) =  &
          -k(49)*n(idx_C)  &
          -k(50)*n(idx_C)  &
          -k(51)*n(idx_C)  &
          +k(140)*n(idx_C2)  &
          +k(157)*n(idx_CO)
      pdj(9) =  &
          -k(46)*n(idx_O)  &
          +k(155)*n(idx_O2)  &
          +k(156)*n(idx_CO)
      pdj(10) =  &
          -k(143)*n(idx_OH)  &
          -k(144)*n(idx_OH)  &
          +k(147)*n(idx_H2O)  &
          +k(148)*n(idx_H2O)
      pdj(11) =  &
          -k(156)*n(idx_CO)  &
          -k(157)*n(idx_CO)
      pdj(13) =  &
          -k(136)*n(idx_CH2)  &
          -k(137)*n(idx_CH2)  &
          -k(138)*n(idx_CH2)  &
          -k(139)*n(idx_CH2)
      pdj(14) =  &
          -k(140)*n(idx_C2)
      pdj(16) =  &
          -k(147)*n(idx_H2O)  &
          -k(148)*n(idx_H2O)  &
          -k(149)*n(idx_H2O)  &
          -k(150)*n(idx_H2O)  &
          -k(151)*n(idx_H2O)  &
          -k(152)*n(idx_H2O)
      pdj(17) =  &
          -k(154)*n(idx_O2)  &
          -k(155)*n(idx_O2)
      pdj(20) =  &
          +k(8)*n(idx_H)  &
          +k(13)*n(idx_H2)  &
          +k(147)*n(idx_H2O)  &
          +k(148)*n(idx_H2O)
      pdj(21) =  &
          -k(5)*n(idx_E)  &
          -k(6)*n(idx_E)  &
          -k(7)*n(idx_E)  &
          -k(8)*n(idx_H)  &
          -k(12)*n(idx_H2)  &
          -k(13)*n(idx_H2)  &
          -k(14)*n(idx_H2)  &
          +k(14)*n(idx_H2)  &
          -k(46)*n(idx_O)  &
          -k(49)*n(idx_C)  &
          -k(50)*n(idx_C)  &
          -k(51)*n(idx_C)  &
          -k(136)*n(idx_CH2)  &
          -k(137)*n(idx_CH2)  &
          -k(138)*n(idx_CH2)  &
          -k(139)*n(idx_CH2)  &
          -k(140)*n(idx_C2)  &
          -k(143)*n(idx_OH)  &
          -k(144)*n(idx_OH)  &
          -k(147)*n(idx_H2O)  &
          -k(148)*n(idx_H2O)  &
          -k(149)*n(idx_H2O)  &
          -k(150)*n(idx_H2O)  &
          -k(151)*n(idx_H2O)  &
          -k(152)*n(idx_H2O)  &
          -k(154)*n(idx_O2)  &
          -k(155)*n(idx_O2)  &
          -k(156)*n(idx_CO)  &
          -k(157)*n(idx_CO)  &
          -k(161)*n(idx_Hk)
      pdj(22) =  &
          +k(12)*n(idx_H2)
      pdj(23) =  &
          +k(49)*n(idx_C)  &
          +k(50)*n(idx_C)  &
          +k(51)*n(idx_C)  &
          +k(136)*n(idx_CH2)  &
          +k(137)*n(idx_CH2)  &
          +k(140)*n(idx_C2)  &
          +k(156)*n(idx_CO)
      pdj(24) =  &
          +k(46)*n(idx_O)  &
          +k(143)*n(idx_OH)  &
          +k(144)*n(idx_OH)  &
          +k(155)*n(idx_O2)  &
          +k(157)*n(idx_CO)
      pdj(28) =  &
          +k(138)*n(idx_CH2)  &
          +k(139)*n(idx_CH2)
      pdj(32) =  &
          +k(149)*n(idx_H2O)  &
          +k(150)*n(idx_H2O)
      pdj(33) =  &
          +k(151)*n(idx_H2O)  &
          +k(152)*n(idx_H2O)
      pdj(35) =  &
          +k(154)*n(idx_O2)
      pdj(36) =  &
          +k(7)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(21)*1d-3
      if(dnn>0.d0) then
        nn(21) = n(21) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==22) then
      pdj(1) =  &
          -k(30)*n(idx_E)  &
          -k(31)*n(idx_E)
      pdj(2) =  &
          -k(32)*n(idx_Hk)
      pdj(5) =  &
          -k(20)*n(idx_H)  &
          +2.d0*k(30)*n(idx_E)  &
          +2.d0*k(31)*n(idx_E)  &
          +k(32)*n(idx_Hk)  &
          +k(85)*n(idx_H2)  &
          +k(87)*n(idx_C)  &
          +k(102)*n(idx_O)  &
          +k(214)
      pdj(7) =  &
          +k(20)*n(idx_H)  &
          +k(32)*n(idx_Hk)  &
          -k(85)*n(idx_H2)
      pdj(8) =  &
          -k(87)*n(idx_C)
      pdj(9) =  &
          -k(102)*n(idx_O)
      pdj(20) =  &
          +k(20)*n(idx_H)  &
          +k(214)
      pdj(22) =  &
          -k(20)*n(idx_H)  &
          -k(30)*n(idx_E)  &
          -k(31)*n(idx_E)  &
          -k(32)*n(idx_Hk)  &
          -k(85)*n(idx_H2)  &
          -k(87)*n(idx_C)  &
          -k(102)*n(idx_O)  &
          -k(214)
      pdj(27) =  &
          +k(85)*n(idx_H2)
      pdj(28) =  &
          +k(87)*n(idx_C)
      pdj(32) =  &
          +k(102)*n(idx_O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(22)*1d-3
      if(dnn>0.d0) then
        nn(22) = n(22) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==23) then
      pdj(1) =  &
          -k(37)*n(idx_E)  &
          -k(38)*n(idx_E)  &
          -k(39)*n(idx_E)
      pdj(5) =  &
          -k(48)*n(idx_H)  &
          +k(90)*n(idx_H2)  &
          +k(107)*n(idx_OH)  &
          +k(108)*n(idx_OH)  &
          +k(113)*n(idx_H2O)  &
          +k(114)*n(idx_H2O)  &
          +k(115)*n(idx_H2O)  &
          -k(200)*n(idx_H)
      pdj(7) =  &
          -k(90)*n(idx_H2)  &
          -k(201)*n(idx_H2)
      pdj(8) =  &
          +k(37)*n(idx_E)  &
          +k(38)*n(idx_E)  &
          +k(39)*n(idx_E)  &
          +k(48)*n(idx_H)  &
          +k(116)*n(idx_H2O)
      pdj(9) =  &
          +k(118)*n(idx_O2)  &
          -k(202)*n(idx_O)  &
          -k(203)*n(idx_O)  &
          -k(269)*n(idx_O)  &
          -k(270)*n(idx_O)
      pdj(10) =  &
          -k(107)*n(idx_OH)  &
          -k(108)*n(idx_OH)
      pdj(11) =  &
          +k(119)*n(idx_O2)
      pdj(16) =  &
          -k(113)*n(idx_H2O)  &
          -k(114)*n(idx_H2O)  &
          -k(115)*n(idx_H2O)  &
          -k(116)*n(idx_H2O)
      pdj(17) =  &
          -k(118)*n(idx_O2)  &
          -k(119)*n(idx_O2)
      pdj(20) =  &
          +k(48)*n(idx_H)
      pdj(23) =  &
          -k(37)*n(idx_E)  &
          -k(38)*n(idx_E)  &
          -k(39)*n(idx_E)  &
          -k(48)*n(idx_H)  &
          -k(90)*n(idx_H2)  &
          -k(107)*n(idx_OH)  &
          -k(108)*n(idx_OH)  &
          -k(113)*n(idx_H2O)  &
          -k(114)*n(idx_H2O)  &
          -k(115)*n(idx_H2O)  &
          -k(116)*n(idx_H2O)  &
          -k(118)*n(idx_O2)  &
          -k(119)*n(idx_O2)  &
          -k(200)*n(idx_H)  &
          -k(201)*n(idx_H2)  &
          -k(202)*n(idx_O)  &
          -k(203)*n(idx_O)  &
          -k(269)*n(idx_O)  &
          -k(270)*n(idx_O)
      pdj(24) =  &
          +k(119)*n(idx_O2)
      pdj(25) =  &
          +k(113)*n(idx_H2O)
      pdj(26) =  &
          +k(114)*n(idx_H2O)  &
          +k(115)*n(idx_H2O)
      pdj(28) =  &
          +k(90)*n(idx_H2)  &
          +k(200)*n(idx_H)
      pdj(29) =  &
          +k(201)*n(idx_H2)
      pdj(30) =  &
          +k(107)*n(idx_OH)  &
          +k(108)*n(idx_OH)  &
          +k(118)*n(idx_O2)  &
          +k(202)*n(idx_O)  &
          +k(203)*n(idx_O)  &
          +k(269)*n(idx_O)  &
          +k(270)*n(idx_O)
      pdj(33) =  &
          +k(116)*n(idx_H2O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(23)*1d-3
      if(dnn>0.d0) then
        nn(23) = n(23) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==24) then
      pdj(1) =  &
          -k(40)*n(idx_E)  &
          -k(41)*n(idx_E)
      pdj(5) =  &
          -k(44)*n(idx_H)  &
          +k(101)*n(idx_H2)
      pdj(7) =  &
          -k(101)*n(idx_H2)
      pdj(8) =  &
          +k(100)*n(idx_C2)  &
          -k(271)*n(idx_C)  &
          -k(272)*n(idx_C)
      pdj(9) =  &
          +k(40)*n(idx_E)  &
          +k(41)*n(idx_E)  &
          +k(44)*n(idx_H)
      pdj(14) =  &
          -k(100)*n(idx_C2)
      pdj(20) =  &
          +k(44)*n(idx_H)
      pdj(24) =  &
          -k(40)*n(idx_E)  &
          -k(41)*n(idx_E)  &
          -k(44)*n(idx_H)  &
          -k(100)*n(idx_C2)  &
          -k(101)*n(idx_H2)  &
          -k(271)*n(idx_C)  &
          -k(272)*n(idx_C)
      pdj(30) =  &
          +k(100)*n(idx_C2)  &
          +k(271)*n(idx_C)  &
          +k(272)*n(idx_C)
      pdj(32) =  &
          +k(101)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(24)*1d-3
      if(dnn>0.d0) then
        nn(24) = n(24) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==25) then
      pdj(1) =  &
          -k(183)*n(idx_E)
      pdj(5) =  &
          +k(183)*n(idx_E)
      pdj(7) =  &
          -k(53)*n(idx_H2)  &
          +k(53)*n(idx_H2)
      pdj(11) =  &
          -k(54)*n(idx_CO)  &
          +k(54)*n(idx_CO)  &
          -k(55)*n(idx_CO)  &
          +k(55)*n(idx_CO)  &
          +k(183)*n(idx_E)
      pdj(25) =  &
          -k(53)*n(idx_H2)  &
          -k(54)*n(idx_CO)  &
          -k(55)*n(idx_CO)  &
          -k(183)*n(idx_E)
      pdj(26) =  &
          +k(53)*n(idx_H2)  &
          +k(54)*n(idx_CO)  &
          +k(55)*n(idx_CO)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(25)*1d-3
      if(dnn>0.d0) then
        nn(25) = n(25) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==26) then
      pdj(1) =  &
          -k(181)*n(idx_E)  &
          -k(182)*n(idx_E)
      pdj(5) =  &
          +k(181)*n(idx_E)
      pdj(8) =  &
          -k(127)*n(idx_C)  &
          +k(182)*n(idx_E)
      pdj(10) =  &
          +k(182)*n(idx_E)
      pdj(11) =  &
          +k(127)*n(idx_C)  &
          +k(128)*n(idx_H2O)  &
          +k(129)*n(idx_H2O)  &
          +k(181)*n(idx_E)
      pdj(16) =  &
          -k(128)*n(idx_H2O)  &
          -k(129)*n(idx_H2O)
      pdj(26) =  &
          -k(127)*n(idx_C)  &
          -k(128)*n(idx_H2O)  &
          -k(129)*n(idx_H2O)  &
          -k(181)*n(idx_E)  &
          -k(182)*n(idx_E)
      pdj(28) =  &
          +k(127)*n(idx_C)
      pdj(34) =  &
          +k(128)*n(idx_H2O)  &
          +k(129)*n(idx_H2O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(26)*1d-3
      if(dnn>0.d0) then
        nn(26) = n(26) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==27) then
      pdj(1) =  &
          -k(162)*n(idx_E)  &
          -k(163)*n(idx_E)
      pdj(5) =  &
          -k(86)*n(idx_H)  &
          +k(89)*n(idx_C)  &
          +k(104)*n(idx_O)  &
          +k(162)*n(idx_E)  &
          +3.d0*k(163)*n(idx_E)  &
          +k(216)
      pdj(7) =  &
          +k(86)*n(idx_H)  &
          +k(88)*n(idx_C)  &
          +k(103)*n(idx_O)  &
          +k(105)*n(idx_OH)  &
          +k(106)*n(idx_OH)  &
          +k(111)*n(idx_H2O)  &
          +k(112)*n(idx_H2O)  &
          +k(123)*n(idx_CO)  &
          +k(124)*n(idx_CO)  &
          +k(125)*n(idx_CO)  &
          +k(126)*n(idx_CO)  &
          +k(162)*n(idx_E)  &
          +k(215)
      pdj(8) =  &
          -k(88)*n(idx_C)  &
          -k(89)*n(idx_C)
      pdj(9) =  &
          -k(103)*n(idx_O)  &
          -k(104)*n(idx_O)
      pdj(10) =  &
          -k(105)*n(idx_OH)  &
          -k(106)*n(idx_OH)
      pdj(11) =  &
          -k(123)*n(idx_CO)  &
          -k(124)*n(idx_CO)  &
          -k(125)*n(idx_CO)  &
          -k(126)*n(idx_CO)
      pdj(16) =  &
          -k(111)*n(idx_H2O)  &
          -k(112)*n(idx_H2O)
      pdj(20) =  &
          +k(215)
      pdj(22) =  &
          +k(86)*n(idx_H)  &
          +k(216)
      pdj(25) =  &
          +k(125)*n(idx_CO)  &
          +k(126)*n(idx_CO)
      pdj(26) =  &
          +k(123)*n(idx_CO)  &
          +k(124)*n(idx_CO)
      pdj(27) =  &
          -k(86)*n(idx_H)  &
          -k(88)*n(idx_C)  &
          -k(89)*n(idx_C)  &
          -k(103)*n(idx_O)  &
          -k(104)*n(idx_O)  &
          -k(105)*n(idx_OH)  &
          -k(106)*n(idx_OH)  &
          -k(111)*n(idx_H2O)  &
          -k(112)*n(idx_H2O)  &
          -k(123)*n(idx_CO)  &
          -k(124)*n(idx_CO)  &
          -k(125)*n(idx_CO)  &
          -k(126)*n(idx_CO)  &
          -k(162)*n(idx_E)  &
          -k(163)*n(idx_E)  &
          -k(215)  &
          -k(216)
      pdj(28) =  &
          +k(88)*n(idx_C)
      pdj(29) =  &
          +k(89)*n(idx_C)
      pdj(32) =  &
          +k(103)*n(idx_O)
      pdj(33) =  &
          +k(104)*n(idx_O)  &
          +k(105)*n(idx_OH)  &
          +k(106)*n(idx_OH)
      pdj(34) =  &
          +k(111)*n(idx_H2O)  &
          +k(112)*n(idx_H2O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(27)*1d-3
      if(dnn>0.d0) then
        nn(27) = n(27) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==28) then
      pdj(1) =  &
          -k(164)*n(idx_E)
      pdj(5) =  &
          -k(91)*n(idx_H)  &
          +k(92)*n(idx_H2)  &
          +k(93)*n(idx_O)  &
          +k(164)*n(idx_E)
      pdj(7) =  &
          +k(91)*n(idx_H)  &
          -k(92)*n(idx_H2)
      pdj(8) =  &
          +k(164)*n(idx_E)  &
          +k(221)
      pdj(9) =  &
          -k(93)*n(idx_O)
      pdj(20) =  &
          +k(221)
      pdj(23) =  &
          +k(91)*n(idx_H)
      pdj(28) =  &
          -k(91)*n(idx_H)  &
          -k(92)*n(idx_H2)  &
          -k(93)*n(idx_O)  &
          -k(164)*n(idx_E)  &
          -k(221)
      pdj(29) =  &
          +k(92)*n(idx_H2)
      pdj(30) =  &
          +k(93)*n(idx_O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(28)*1d-3
      if(dnn>0.d0) then
        nn(28) = n(28) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==29) then
      pdj(1) =  &
          -k(165)*n(idx_E)  &
          -k(166)*n(idx_E)  &
          -k(167)*n(idx_E)
      pdj(5) =  &
          -k(94)*n(idx_H)  &
          +k(95)*n(idx_H2)  &
          +k(96)*n(idx_O)  &
          +k(165)*n(idx_E)  &
          +2.d0*k(167)*n(idx_E)  &
          +k(224)
      pdj(7) =  &
          +k(94)*n(idx_H)  &
          -k(95)*n(idx_H2)  &
          +k(166)*n(idx_E)
      pdj(8) =  &
          +k(166)*n(idx_E)  &
          +k(167)*n(idx_E)
      pdj(9) =  &
          -k(96)*n(idx_O)
      pdj(10) =  &
          +k(120)*n(idx_O2)
      pdj(12) =  &
          +k(165)*n(idx_E)
      pdj(17) =  &
          -k(120)*n(idx_O2)
      pdj(26) =  &
          +k(96)*n(idx_O)  &
          +k(120)*n(idx_O2)
      pdj(28) =  &
          +k(94)*n(idx_H)  &
          +k(224)
      pdj(29) =  &
          -k(94)*n(idx_H)  &
          -k(95)*n(idx_H2)  &
          -k(96)*n(idx_O)  &
          -k(120)*n(idx_O2)  &
          -k(165)*n(idx_E)  &
          -k(166)*n(idx_E)  &
          -k(167)*n(idx_E)  &
          -k(224)
      pdj(31) =  &
          +k(95)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(29)*1d-3
      if(dnn>0.d0) then
        nn(29) = n(29) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==30) then
      pdj(1) =  &
          -k(180)*n(idx_E)
      pdj(5) =  &
          -k(158)*n(idx_H)
      pdj(8) =  &
          +k(180)*n(idx_E)
      pdj(9) =  &
          +k(180)*n(idx_E)
      pdj(11) =  &
          +k(158)*n(idx_H)
      pdj(20) =  &
          +k(158)*n(idx_H)
      pdj(30) =  &
          -k(158)*n(idx_H)  &
          -k(180)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(30)*1d-3
      if(dnn>0.d0) then
        nn(30) = n(30) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==31) then
      pdj(1) =  &
          -k(168)*n(idx_E)  &
          -k(169)*n(idx_E)  &
          -k(170)*n(idx_E)
      pdj(5) =  &
          -k(97)*n(idx_H)  &
          +k(168)*n(idx_E)  &
          +2.d0*k(170)*n(idx_E)  &
          +k(225)
      pdj(7) =  &
          +k(97)*n(idx_H)  &
          +k(98)*n(idx_O)  &
          +k(99)*n(idx_O)  &
          +k(169)*n(idx_E)  &
          +k(226)
      pdj(9) =  &
          -k(98)*n(idx_O)  &
          -k(99)*n(idx_O)
      pdj(12) =  &
          +k(169)*n(idx_E)  &
          +k(170)*n(idx_E)
      pdj(13) =  &
          +k(168)*n(idx_E)
      pdj(25) =  &
          +k(98)*n(idx_O)
      pdj(26) =  &
          +k(99)*n(idx_O)
      pdj(28) =  &
          +k(226)
      pdj(29) =  &
          +k(97)*n(idx_H)  &
          +k(225)
      pdj(31) =  &
          -k(97)*n(idx_H)  &
          -k(98)*n(idx_O)  &
          -k(99)*n(idx_O)  &
          -k(168)*n(idx_E)  &
          -k(169)*n(idx_E)  &
          -k(170)*n(idx_E)  &
          -k(225)  &
          -k(226)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(31)*1d-3
      if(dnn>0.d0) then
        nn(31) = n(31) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==32) then
      pdj(1) =  &
          -k(171)*n(idx_E)
      pdj(5) =  &
          +k(109)*n(idx_H2)  &
          +k(171)*n(idx_E)
      pdj(7) =  &
          -k(109)*n(idx_H2)
      pdj(9) =  &
          +k(171)*n(idx_E)  &
          +k(231)
      pdj(20) =  &
          +k(231)
      pdj(32) =  &
          -k(109)*n(idx_H2)  &
          -k(171)*n(idx_E)  &
          -k(231)
      pdj(33) =  &
          +k(109)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(32)*1d-3
      if(dnn>0.d0) then
        nn(32) = n(32) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==33) then
      pdj(1) =  &
          -k(172)*n(idx_E)  &
          -k(173)*n(idx_E)  &
          -k(174)*n(idx_E)
      pdj(5) =  &
          +k(110)*n(idx_H2)  &
          +k(173)*n(idx_E)  &
          +2.d0*k(174)*n(idx_E)  &
          +k(241)
      pdj(7) =  &
          -k(110)*n(idx_H2)  &
          +k(172)*n(idx_E)  &
          +k(240)
      pdj(9) =  &
          +k(172)*n(idx_E)  &
          +k(174)*n(idx_E)  &
          +k(238)
      pdj(10) =  &
          +k(173)*n(idx_E)  &
          +k(239)
      pdj(20) =  &
          +k(239)
      pdj(22) =  &
          +k(238)
      pdj(24) =  &
          +k(240)
      pdj(32) =  &
          +k(241)
      pdj(33) =  &
          -k(110)*n(idx_H2)  &
          -k(172)*n(idx_E)  &
          -k(173)*n(idx_E)  &
          -k(174)*n(idx_E)  &
          -k(238)  &
          -k(239)  &
          -k(240)  &
          -k(241)
      pdj(34) =  &
          +k(110)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(33)*1d-3
      if(dnn>0.d0) then
        nn(33) = n(33) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==34) then
      pdj(1) =  &
          -k(175)*n(idx_E)  &
          -k(176)*n(idx_E)  &
          -k(177)*n(idx_E)  &
          -k(178)*n(idx_E)
      pdj(5) =  &
          +2.d0*k(175)*n(idx_E)  &
          +k(176)*n(idx_E)  &
          +k(177)*n(idx_E)  &
          +k(244)
      pdj(7) =  &
          +k(117)*n(idx_C)  &
          +k(176)*n(idx_E)  &
          +k(178)*n(idx_E)  &
          +k(245)
      pdj(8) =  &
          -k(117)*n(idx_C)
      pdj(9) =  &
          +k(176)*n(idx_E)
      pdj(10) =  &
          +k(175)*n(idx_E)  &
          +k(178)*n(idx_E)  &
          +k(243)
      pdj(16) =  &
          +k(177)*n(idx_E)  &
          +k(242)
      pdj(20) =  &
          +k(242)
      pdj(22) =  &
          +k(243)
      pdj(26) =  &
          +k(117)*n(idx_C)
      pdj(32) =  &
          +k(245)
      pdj(33) =  &
          +k(244)
      pdj(34) =  &
          -k(117)*n(idx_C)  &
          -k(175)*n(idx_E)  &
          -k(176)*n(idx_E)  &
          -k(177)*n(idx_E)  &
          -k(178)*n(idx_E)  &
          -k(242)  &
          -k(243)  &
          -k(244)  &
          -k(245)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(34)*1d-3
      if(dnn>0.d0) then
        nn(34) = n(34) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==35) then
      pdj(1) =  &
          -k(179)*n(idx_E)
      pdj(8) =  &
          -k(121)*n(idx_C)  &
          -k(122)*n(idx_C)
      pdj(9) =  &
          +k(121)*n(idx_C)  &
          +2.d0*k(179)*n(idx_E)
      pdj(17) =  &
          +k(122)*n(idx_C)
      pdj(23) =  &
          +k(122)*n(idx_C)
      pdj(30) =  &
          +k(121)*n(idx_C)
      pdj(35) =  &
          -k(121)*n(idx_C)  &
          -k(122)*n(idx_C)  &
          -k(179)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(35)*1d-3
      if(dnn>0.d0) then
        nn(35) = n(35) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==36) then
      pdj(1) =  &
          -k(15)*n(idx_E)
      pdj(21) =  &
          +k(15)*n(idx_E)
      pdj(36) =  &
          -k(15)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(36)*1d-3
      if(dnn>0.d0) then
        nn(36) = n(36) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==37) then
      pdj(39) = 0.d0
    elseif(j==38) then
      pdj(39) = 0.d0
    elseif(j==39) then
      !use fex to compute temperature-dependent Jacobian
      dnn = n(idx_Tgas)*1d-3
      nn(:) = n(:)
      nn(idx_Tgas) = n(idx_Tgas) + dnn
      call fex(neq,tt,nn(:),dn(:))
      do i=1,neq-1
        pdj(i) = dn(i) / dnn
      end do
    elseif(j==40) then
      pdj(39) = 0.d0
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[E_dot]/d[E]
    pd(1,1) =  &
        -k(1)*n(idx_H)  &
        +2.d0*k(1)*n(idx_H)  &
        -k(2)*n(idx_Hj)  &
        -k(3)*n(idx_Hj)  &
        -k(4)*n(idx_HE)  &
        +2.d0*k(4)*n(idx_HE)  &
        -k(5)*n(idx_HEj)  &
        -k(6)*n(idx_HEj)  &
        -k(7)*n(idx_HEj)  &
        +2.d0*k(7)*n(idx_HEj)  &
        -k(15)*n(idx_HEjj)  &
        -k(16)*n(idx_H)  &
        -k(23)*n(idx_H2)  &
        +k(23)*n(idx_H2)  &
        -k(25)*n(idx_Hk)  &
        +2.d0*k(25)*n(idx_Hk)  &
        -k(30)*n(idx_H2j)  &
        -k(31)*n(idx_H2j)  &
        -k(37)*n(idx_Cj)  &
        -k(38)*n(idx_Cj)  &
        -k(39)*n(idx_Cj)  &
        -k(40)*n(idx_Oj)  &
        -k(41)*n(idx_Oj)  &
        -k(42)*n(idx_C)  &
        +2.d0*k(42)*n(idx_C)  &
        -k(43)*n(idx_O)  &
        +2.d0*k(43)*n(idx_O)  &
        -k(162)*n(idx_H3j)  &
        -k(163)*n(idx_H3j)  &
        -k(164)*n(idx_CHj)  &
        -k(165)*n(idx_CH2j)  &
        -k(166)*n(idx_CH2j)  &
        -k(167)*n(idx_CH2j)  &
        -k(168)*n(idx_CH3j)  &
        -k(169)*n(idx_CH3j)  &
        -k(170)*n(idx_CH3j)  &
        -k(171)*n(idx_OHj)  &
        -k(172)*n(idx_H2Oj)  &
        -k(173)*n(idx_H2Oj)  &
        -k(174)*n(idx_H2Oj)  &
        -k(175)*n(idx_H3Oj)  &
        -k(176)*n(idx_H3Oj)  &
        -k(177)*n(idx_H3Oj)  &
        -k(178)*n(idx_H3Oj)  &
        -k(179)*n(idx_O2j)  &
        -k(180)*n(idx_COj)  &
        -k(181)*n(idx_HCOj)  &
        -k(182)*n(idx_HCOj)  &
        -k(183)*n(idx_HOCj)  &
        -k(195)*n(idx_C)  &
        -k(204)*n(idx_O)

    !d[H-_dot]/d[E]
    pd(2,1) =  &
        +k(16)*n(idx_H)  &
        -k(25)*n(idx_Hk)

    !d[C-_dot]/d[E]
    pd(3,1) =  &
        +k(195)*n(idx_C)

    !d[O-_dot]/d[E]
    pd(4,1) =  &
        +k(204)*n(idx_O)

    !d[H_dot]/d[E]
    pd(5,1) =  &
        -k(1)*n(idx_H)  &
        +k(2)*n(idx_Hj)  &
        +k(3)*n(idx_Hj)  &
        -k(16)*n(idx_H)  &
        +2.d0*k(23)*n(idx_H2)  &
        +k(25)*n(idx_Hk)  &
        +2.d0*k(30)*n(idx_H2j)  &
        +2.d0*k(31)*n(idx_H2j)  &
        +k(162)*n(idx_H3j)  &
        +3.d0*k(163)*n(idx_H3j)  &
        +k(164)*n(idx_CHj)  &
        +k(165)*n(idx_CH2j)  &
        +2.d0*k(167)*n(idx_CH2j)  &
        +k(168)*n(idx_CH3j)  &
        +2.d0*k(170)*n(idx_CH3j)  &
        +k(171)*n(idx_OHj)  &
        +k(173)*n(idx_H2Oj)  &
        +2.d0*k(174)*n(idx_H2Oj)  &
        +2.d0*k(175)*n(idx_H3Oj)  &
        +k(176)*n(idx_H3Oj)  &
        +k(177)*n(idx_H3Oj)  &
        +k(181)*n(idx_HCOj)  &
        +k(183)*n(idx_HOCj)

    !d[HE_dot]/d[E]
    pd(6,1) =  &
        -k(4)*n(idx_HE)  &
        +k(5)*n(idx_HEj)  &
        +k(6)*n(idx_HEj)

    !d[H2_dot]/d[E]
    pd(7,1) =  &
        -k(23)*n(idx_H2)  &
        +k(162)*n(idx_H3j)  &
        +k(166)*n(idx_CH2j)  &
        +k(169)*n(idx_CH3j)  &
        +k(172)*n(idx_H2Oj)  &
        +k(176)*n(idx_H3Oj)  &
        +k(178)*n(idx_H3Oj)

    !d[C_dot]/d[E]
    pd(8,1) =  &
        +k(37)*n(idx_Cj)  &
        +k(38)*n(idx_Cj)  &
        +k(39)*n(idx_Cj)  &
        -k(42)*n(idx_C)  &
        +k(164)*n(idx_CHj)  &
        +k(166)*n(idx_CH2j)  &
        +k(167)*n(idx_CH2j)  &
        +k(180)*n(idx_COj)  &
        +k(182)*n(idx_HCOj)  &
        -k(195)*n(idx_C)

    !d[O_dot]/d[E]
    pd(9,1) =  &
        +k(40)*n(idx_Oj)  &
        +k(41)*n(idx_Oj)  &
        -k(43)*n(idx_O)  &
        +k(171)*n(idx_OHj)  &
        +k(172)*n(idx_H2Oj)  &
        +k(174)*n(idx_H2Oj)  &
        +k(176)*n(idx_H3Oj)  &
        +2.d0*k(179)*n(idx_O2j)  &
        +k(180)*n(idx_COj)  &
        -k(204)*n(idx_O)

    !d[OH_dot]/d[E]
    pd(10,1) =  &
        +k(173)*n(idx_H2Oj)  &
        +k(175)*n(idx_H3Oj)  &
        +k(178)*n(idx_H3Oj)  &
        +k(182)*n(idx_HCOj)

    !d[CO_dot]/d[E]
    pd(11,1) =  &
        +k(181)*n(idx_HCOj)  &
        +k(183)*n(idx_HOCj)

    !d[CH_dot]/d[E]
    pd(12,1) =  &
        +k(165)*n(idx_CH2j)  &
        +k(169)*n(idx_CH3j)  &
        +k(170)*n(idx_CH3j)

    !d[CH2_dot]/d[E]
    pd(13,1) =  &
        +k(168)*n(idx_CH3j)

    !d[H2O_dot]/d[E]
    pd(16,1) =  &
        +k(177)*n(idx_H3Oj)

    !d[H+_dot]/d[E]
    pd(20,1) =  &
        +k(1)*n(idx_H)  &
        -k(2)*n(idx_Hj)  &
        -k(3)*n(idx_Hj)

    !d[HE+_dot]/d[E]
    pd(21,1) =  &
        +k(4)*n(idx_HE)  &
        -k(5)*n(idx_HEj)  &
        -k(6)*n(idx_HEj)  &
        -k(7)*n(idx_HEj)  &
        +k(15)*n(idx_HEjj)

    !d[H2+_dot]/d[E]
    pd(22,1) =  &
        -k(30)*n(idx_H2j)  &
        -k(31)*n(idx_H2j)

    !d[C+_dot]/d[E]
    pd(23,1) =  &
        -k(37)*n(idx_Cj)  &
        -k(38)*n(idx_Cj)  &
        -k(39)*n(idx_Cj)  &
        +k(42)*n(idx_C)

    !d[O+_dot]/d[E]
    pd(24,1) =  &
        -k(40)*n(idx_Oj)  &
        -k(41)*n(idx_Oj)  &
        +k(43)*n(idx_O)

    !d[HOC+_dot]/d[E]
    pd(25,1) =  &
        -k(183)*n(idx_HOCj)

    !d[HCO+_dot]/d[E]
    pd(26,1) =  &
        -k(181)*n(idx_HCOj)  &
        -k(182)*n(idx_HCOj)

    !d[H3+_dot]/d[E]
    pd(27,1) =  &
        -k(162)*n(idx_H3j)  &
        -k(163)*n(idx_H3j)

    !d[CH+_dot]/d[E]
    pd(28,1) =  &
        -k(164)*n(idx_CHj)

    !d[CH2+_dot]/d[E]
    pd(29,1) =  &
        -k(165)*n(idx_CH2j)  &
        -k(166)*n(idx_CH2j)  &
        -k(167)*n(idx_CH2j)

    !d[CO+_dot]/d[E]
    pd(30,1) =  &
        -k(180)*n(idx_COj)

    !d[CH3+_dot]/d[E]
    pd(31,1) =  &
        -k(168)*n(idx_CH3j)  &
        -k(169)*n(idx_CH3j)  &
        -k(170)*n(idx_CH3j)

    !d[OH+_dot]/d[E]
    pd(32,1) =  &
        -k(171)*n(idx_OHj)

    !d[H2O+_dot]/d[E]
    pd(33,1) =  &
        -k(172)*n(idx_H2Oj)  &
        -k(173)*n(idx_H2Oj)  &
        -k(174)*n(idx_H2Oj)

    !d[H3O+_dot]/d[E]
    pd(34,1) =  &
        -k(175)*n(idx_H3Oj)  &
        -k(176)*n(idx_H3Oj)  &
        -k(177)*n(idx_H3Oj)  &
        -k(178)*n(idx_H3Oj)

    !d[O2+_dot]/d[E]
    pd(35,1) =  &
        -k(179)*n(idx_O2j)

    !d[HE++_dot]/d[E]
    pd(36,1) =  &
        +k(7)*n(idx_HEj)  &
        -k(15)*n(idx_HEjj)

    !d[Tgas_dot]/d[E]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(1)*1d-3
    if(dnn>0.d0) then
      nn(1) = n(1) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,1) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H-]
    pd(1,2) =  &
        +k(17)*n(idx_H)  &
        -k(25)*n(idx_E)  &
        +2.d0*k(25)*n(idx_E)  &
        +k(26)*n(idx_H)  &
        +k(27)*n(idx_H)  &
        +k(29)*n(idx_Hj)  &
        +k(184)*n(idx_C)  &
        +k(185)*n(idx_O)  &
        +k(186)*n(idx_OH)  &
        +k(213)

    !d[H-_dot]/d[H-]
    pd(2,2) =  &
        -k(17)*n(idx_H)  &
        -k(25)*n(idx_E)  &
        -k(26)*n(idx_H)  &
        -k(27)*n(idx_H)  &
        -k(28)*n(idx_Hj)  &
        -k(29)*n(idx_Hj)  &
        -k(32)*n(idx_H2j)  &
        -k(161)*n(idx_HEj)  &
        -k(184)*n(idx_C)  &
        -k(185)*n(idx_O)  &
        -k(186)*n(idx_OH)  &
        -k(213)

    !d[H_dot]/d[H-]
    pd(5,2) =  &
        -k(17)*n(idx_H)  &
        +k(25)*n(idx_E)  &
        -k(26)*n(idx_H)  &
        +2.d0*k(26)*n(idx_H)  &
        -k(27)*n(idx_H)  &
        +2.d0*k(27)*n(idx_H)  &
        +2.d0*k(28)*n(idx_Hj)  &
        +k(32)*n(idx_H2j)  &
        +k(161)*n(idx_HEj)  &
        +k(213)

    !d[HE_dot]/d[H-]
    pd(6,2) =  &
        +k(161)*n(idx_HEj)

    !d[H2_dot]/d[H-]
    pd(7,2) =  &
        +k(17)*n(idx_H)  &
        +k(32)*n(idx_H2j)

    !d[C_dot]/d[H-]
    pd(8,2) =  &
        -k(184)*n(idx_C)

    !d[O_dot]/d[H-]
    pd(9,2) =  &
        -k(185)*n(idx_O)

    !d[OH_dot]/d[H-]
    pd(10,2) =  &
        +k(185)*n(idx_O)  &
        -k(186)*n(idx_OH)

    !d[CH_dot]/d[H-]
    pd(12,2) =  &
        +k(184)*n(idx_C)

    !d[H2O_dot]/d[H-]
    pd(16,2) =  &
        +k(186)*n(idx_OH)

    !d[H+_dot]/d[H-]
    pd(20,2) =  &
        -k(28)*n(idx_Hj)  &
        -k(29)*n(idx_Hj)

    !d[HE+_dot]/d[H-]
    pd(21,2) =  &
        -k(161)*n(idx_HEj)

    !d[H2+_dot]/d[H-]
    pd(22,2) =  &
        +k(29)*n(idx_Hj)  &
        -k(32)*n(idx_H2j)

    !d[Tgas_dot]/d[H-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(2)*1d-3
    if(dnn>0.d0) then
      nn(2) = n(2) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,2) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C-]
    pd(1,3) =  &
        +k(187)*n(idx_H)  &
        +k(188)*n(idx_H2)  &
        +k(189)*n(idx_O)  &
        +k(218)

    !d[C-_dot]/d[C-]
    pd(3,3) =  &
        -k(159)*n(idx_Hj)  &
        -k(187)*n(idx_H)  &
        -k(188)*n(idx_H2)  &
        -k(189)*n(idx_O)  &
        -k(218)

    !d[H_dot]/d[C-]
    pd(5,3) =  &
        +k(159)*n(idx_Hj)  &
        -k(187)*n(idx_H)

    !d[H2_dot]/d[C-]
    pd(7,3) =  &
        -k(188)*n(idx_H2)

    !d[C_dot]/d[C-]
    pd(8,3) =  &
        +k(159)*n(idx_Hj)  &
        +k(218)

    !d[O_dot]/d[C-]
    pd(9,3) =  &
        -k(189)*n(idx_O)

    !d[CO_dot]/d[C-]
    pd(11,3) =  &
        +k(189)*n(idx_O)

    !d[CH_dot]/d[C-]
    pd(12,3) =  &
        +k(187)*n(idx_H)

    !d[CH2_dot]/d[C-]
    pd(13,3) =  &
        +k(188)*n(idx_H2)

    !d[H+_dot]/d[C-]
    pd(20,3) =  &
        -k(159)*n(idx_Hj)

    !d[Tgas_dot]/d[C-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(3)*1d-3
    if(dnn>0.d0) then
      nn(3) = n(3) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,3) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O-]
    pd(1,4) =  &
        +k(190)*n(idx_H)  &
        +k(191)*n(idx_H2)  &
        +k(192)*n(idx_C)  &
        +k(228)

    !d[O-_dot]/d[O-]
    pd(4,4) =  &
        -k(160)*n(idx_Hj)  &
        -k(190)*n(idx_H)  &
        -k(191)*n(idx_H2)  &
        -k(192)*n(idx_C)  &
        -k(228)

    !d[H_dot]/d[O-]
    pd(5,4) =  &
        +k(160)*n(idx_Hj)  &
        -k(190)*n(idx_H)

    !d[H2_dot]/d[O-]
    pd(7,4) =  &
        -k(191)*n(idx_H2)

    !d[C_dot]/d[O-]
    pd(8,4) =  &
        -k(192)*n(idx_C)

    !d[O_dot]/d[O-]
    pd(9,4) =  &
        +k(160)*n(idx_Hj)  &
        +k(228)

    !d[OH_dot]/d[O-]
    pd(10,4) =  &
        +k(190)*n(idx_H)

    !d[CO_dot]/d[O-]
    pd(11,4) =  &
        +k(192)*n(idx_C)

    !d[H2O_dot]/d[O-]
    pd(16,4) =  &
        +k(191)*n(idx_H2)

    !d[H+_dot]/d[O-]
    pd(20,4) =  &
        -k(160)*n(idx_Hj)

    !d[Tgas_dot]/d[O-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(4)*1d-3
    if(dnn>0.d0) then
      nn(4) = n(4) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,4) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H]
    pd(1,5) =  &
        -k(1)*n(idx_E)  &
        +2.d0*k(1)*n(idx_E)  &
        -k(16)*n(idx_E)  &
        +k(17)*n(idx_Hk)  &
        +k(26)*n(idx_Hk)  &
        +k(27)*n(idx_Hk)  &
        +k(187)*n(idx_Ck)  &
        +k(190)*n(idx_Ok)  &
        +k(246)

    !d[H-_dot]/d[H]
    pd(2,5) =  &
        +k(16)*n(idx_E)  &
        -k(17)*n(idx_Hk)  &
        -k(26)*n(idx_Hk)  &
        -k(27)*n(idx_Hk)

    !d[C-_dot]/d[H]
    pd(3,5) =  &
        -k(187)*n(idx_Ck)

    !d[O-_dot]/d[H]
    pd(4,5) =  &
        -k(190)*n(idx_Ok)

    !d[H_dot]/d[H]
    pd(5,5) =  &
        -k(1)*n(idx_E)  &
        -k(8)*n(idx_HEj)  &
        -k(16)*n(idx_E)  &
        -k(17)*n(idx_Hk)  &
        -k(18)*n(idx_Hj)  &
        -k(19)*n(idx_Hj)  &
        -k(20)*n(idx_H2j)  &
        -k(24)*n(idx_H2)  &
        +3.d0*k(24)*n(idx_H2)  &
        -k(26)*n(idx_Hk)  &
        +2.d0*k(26)*n(idx_Hk)  &
        -k(27)*n(idx_Hk)  &
        +2.d0*k(27)*n(idx_Hk)  &
        -4.d0*k(34)*n(idx_H)*n(idx_HE)  &
        -9.d0*k(35)*n(idx_H)*n(idx_H)  &
        +3.d0*k(35)*n(idx_H)*n(idx_H)  &
        -4.d0*k(36)*n(idx_H2)*n(idx_H)  &
        -k(44)*n(idx_Oj)  &
        -k(48)*n(idx_Cj)  &
        -k(52)*n(idx_OH)  &
        +2.d0*k(52)*n(idx_OH)  &
        -k(57)*n(idx_CH)  &
        -k(63)*n(idx_CH2)  &
        -k(71)*n(idx_OH)  &
        -k(72)*n(idx_OH)  &
        -k(79)*n(idx_H2O)  &
        -k(80)*n(idx_O2)  &
        -k(84)*n(idx_CO)  &
        -k(86)*n(idx_H3j)  &
        -k(91)*n(idx_CHj)  &
        -k(94)*n(idx_CH2j)  &
        -k(97)*n(idx_CH3j)  &
        -k(158)*n(idx_COj)  &
        -k(187)*n(idx_Ck)  &
        -k(190)*n(idx_Ok)  &
        -k(196)*n(idx_C)  &
        -k(200)*n(idx_Cj)  &
        -k(205)*n(idx_O)  &
        -k(207)*n(idx_OH)  &
        -k(246)  &
        -k(273)*n(idx_O)  &
        -k(274)*n(idx_OH)

    !d[HE_dot]/d[H]
    pd(6,5) =  &
        +k(8)*n(idx_HEj)  &
        -2.d0*k(34)*n(idx_H)*n(idx_HE)  &
        +2.d0*k(34)*n(idx_H)*n(idx_HE)

    !d[H2_dot]/d[H]
    pd(7,5) =  &
        +k(17)*n(idx_Hk)  &
        +k(20)*n(idx_H2j)  &
        -k(24)*n(idx_H2)  &
        +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
        +3.d0*k(35)*n(idx_H)*n(idx_H)  &
        -2.d0*k(36)*n(idx_H2)*n(idx_H)  &
        +4.d0*k(36)*n(idx_H2)*n(idx_H)  &
        +k(57)*n(idx_CH)  &
        +k(63)*n(idx_CH2)  &
        +k(71)*n(idx_OH)  &
        +k(72)*n(idx_OH)  &
        +k(79)*n(idx_H2O)  &
        +k(86)*n(idx_H3j)  &
        +k(91)*n(idx_CHj)  &
        +k(94)*n(idx_CH2j)  &
        +k(97)*n(idx_CH3j)

    !d[C_dot]/d[H]
    pd(8,5) =  &
        +k(48)*n(idx_Cj)  &
        +k(57)*n(idx_CH)  &
        +k(84)*n(idx_CO)  &
        -k(196)*n(idx_C)

    !d[O_dot]/d[H]
    pd(9,5) =  &
        +k(44)*n(idx_Oj)  &
        +k(52)*n(idx_OH)  &
        +k(71)*n(idx_OH)  &
        +k(72)*n(idx_OH)  &
        +k(80)*n(idx_O2)  &
        -k(205)*n(idx_O)  &
        -k(273)*n(idx_O)

    !d[OH_dot]/d[H]
    pd(10,5) =  &
        -k(52)*n(idx_OH)  &
        -k(71)*n(idx_OH)  &
        -k(72)*n(idx_OH)  &
        +k(79)*n(idx_H2O)  &
        +k(80)*n(idx_O2)  &
        +k(84)*n(idx_CO)  &
        +k(190)*n(idx_Ok)  &
        +k(205)*n(idx_O)  &
        -k(207)*n(idx_OH)  &
        +k(273)*n(idx_O)  &
        -k(274)*n(idx_OH)

    !d[CO_dot]/d[H]
    pd(11,5) =  &
        -k(84)*n(idx_CO)  &
        +k(158)*n(idx_COj)

    !d[CH_dot]/d[H]
    pd(12,5) =  &
        -k(57)*n(idx_CH)  &
        +k(63)*n(idx_CH2)  &
        +k(187)*n(idx_Ck)  &
        +k(196)*n(idx_C)

    !d[CH2_dot]/d[H]
    pd(13,5) =  &
        -k(63)*n(idx_CH2)

    !d[H2O_dot]/d[H]
    pd(16,5) =  &
        -k(79)*n(idx_H2O)  &
        +k(207)*n(idx_OH)  &
        +k(274)*n(idx_OH)

    !d[O2_dot]/d[H]
    pd(17,5) =  &
        -k(80)*n(idx_O2)

    !d[H+_dot]/d[H]
    pd(20,5) =  &
        +k(1)*n(idx_E)  &
        +k(8)*n(idx_HEj)  &
        -k(18)*n(idx_Hj)  &
        -k(19)*n(idx_Hj)  &
        +k(20)*n(idx_H2j)  &
        +k(44)*n(idx_Oj)  &
        +k(48)*n(idx_Cj)  &
        +k(158)*n(idx_COj)  &
        +k(246)

    !d[HE+_dot]/d[H]
    pd(21,5) =  &
        -k(8)*n(idx_HEj)

    !d[H2+_dot]/d[H]
    pd(22,5) =  &
        +k(18)*n(idx_Hj)  &
        +k(19)*n(idx_Hj)  &
        -k(20)*n(idx_H2j)  &
        +k(86)*n(idx_H3j)

    !d[C+_dot]/d[H]
    pd(23,5) =  &
        -k(48)*n(idx_Cj)  &
        +k(91)*n(idx_CHj)  &
        -k(200)*n(idx_Cj)

    !d[O+_dot]/d[H]
    pd(24,5) =  &
        -k(44)*n(idx_Oj)

    !d[H3+_dot]/d[H]
    pd(27,5) =  &
        -k(86)*n(idx_H3j)

    !d[CH+_dot]/d[H]
    pd(28,5) =  &
        -k(91)*n(idx_CHj)  &
        +k(94)*n(idx_CH2j)  &
        +k(200)*n(idx_Cj)

    !d[CH2+_dot]/d[H]
    pd(29,5) =  &
        -k(94)*n(idx_CH2j)  &
        +k(97)*n(idx_CH3j)

    !d[CO+_dot]/d[H]
    pd(30,5) =  &
        -k(158)*n(idx_COj)

    !d[CH3+_dot]/d[H]
    pd(31,5) =  &
        -k(97)*n(idx_CH3j)

    !d[Tgas_dot]/d[H]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(5)*1d-3
    if(dnn>0.d0) then
      nn(5) = n(5) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,5) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE]
    pd(1,6) =  &
        -k(4)*n(idx_E)  &
        +2.d0*k(4)*n(idx_E)  &
        +k(247)

    !d[H_dot]/d[HE]
    pd(5,6) =  &
        +k(9)*n(idx_Hj)  &
        +k(10)*n(idx_Hj)  &
        +2.d0*k(11)*n(idx_H2)  &
        -2.d0*k(34)*n(idx_H)*n(idx_H)

    !d[HE_dot]/d[HE]
    pd(6,6) =  &
        -k(4)*n(idx_E)  &
        -k(9)*n(idx_Hj)  &
        -k(10)*n(idx_Hj)  &
        -k(11)*n(idx_H2)  &
        +k(11)*n(idx_H2)  &
        -k(34)*n(idx_H)*n(idx_H)  &
        +k(34)*n(idx_H)*n(idx_H)  &
        -k(247)

    !d[H2_dot]/d[HE]
    pd(7,6) =  &
        -k(11)*n(idx_H2)  &
        +k(34)*n(idx_H)*n(idx_H)

    !d[H+_dot]/d[HE]
    pd(20,6) =  &
        -k(9)*n(idx_Hj)  &
        -k(10)*n(idx_Hj)

    !d[HE+_dot]/d[HE]
    pd(21,6) =  &
        +k(4)*n(idx_E)  &
        +k(9)*n(idx_Hj)  &
        +k(10)*n(idx_Hj)  &
        +k(247)

    !d[Tgas_dot]/d[HE]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(6)*1d-3
    if(dnn>0.d0) then
      nn(6) = n(6) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,6) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2]
    pd(1,7) =  &
        -k(23)*n(idx_E)  &
        +k(23)*n(idx_E)  &
        +k(188)*n(idx_Ck)  &
        +k(191)*n(idx_Ok)  &
        +k(254)  &
        +k(264)

    !d[H-_dot]/d[H2]
    pd(2,7) =  &
        +k(253)

    !d[C-_dot]/d[H2]
    pd(3,7) =  &
        -k(188)*n(idx_Ck)

    !d[O-_dot]/d[H2]
    pd(4,7) =  &
        -k(191)*n(idx_Ok)

    !d[H_dot]/d[H2]
    pd(5,7) =  &
        +2.d0*k(11)*n(idx_HE)  &
        +k(13)*n(idx_HEj)  &
        +2.d0*k(14)*n(idx_HEj)  &
        +k(21)*n(idx_Hj)  &
        +k(22)*n(idx_Hj)  &
        +2.d0*k(23)*n(idx_E)  &
        -k(24)*n(idx_H)  &
        +3.d0*k(24)*n(idx_H)  &
        +4.d0*k(33)*n(idx_H2)  &
        -2.d0*k(36)*n(idx_H)*n(idx_H)  &
        +k(56)*n(idx_C)  &
        +k(58)*n(idx_CH)  &
        +k(70)*n(idx_O)  &
        +k(73)*n(idx_OH)  &
        +k(85)*n(idx_H2j)  &
        +k(90)*n(idx_Cj)  &
        +k(92)*n(idx_CHj)  &
        +k(95)*n(idx_CH2j)  &
        +k(101)*n(idx_Oj)  &
        +k(109)*n(idx_OHj)  &
        +k(110)*n(idx_H2Oj)  &
        +2.d0*k(193)*n(idx_Hj)  &
        +2.d0*k(237)  &
        +2.d0*k(252)  &
        +k(264)

    !d[HE_dot]/d[H2]
    pd(6,7) =  &
        -k(11)*n(idx_HE)  &
        +k(11)*n(idx_HE)  &
        +k(12)*n(idx_HEj)  &
        +k(13)*n(idx_HEj)

    !d[H2_dot]/d[H2]
    pd(7,7) =  &
        -k(11)*n(idx_HE)  &
        -k(12)*n(idx_HEj)  &
        -k(13)*n(idx_HEj)  &
        -k(14)*n(idx_HEj)  &
        -k(21)*n(idx_Hj)  &
        -k(22)*n(idx_Hj)  &
        -k(23)*n(idx_E)  &
        -k(24)*n(idx_H)  &
        -4.d0*k(33)*n(idx_H2)  &
        +2.d0*k(33)*n(idx_H2)  &
        -k(36)*n(idx_H)*n(idx_H)  &
        +2.d0*k(36)*n(idx_H)*n(idx_H)  &
        -k(53)*n(idx_HOCj)  &
        +k(53)*n(idx_HOCj)  &
        -k(56)*n(idx_C)  &
        -k(58)*n(idx_CH)  &
        -k(70)*n(idx_O)  &
        -k(73)*n(idx_OH)  &
        -k(81)*n(idx_O2)  &
        -k(85)*n(idx_H2j)  &
        -k(90)*n(idx_Cj)  &
        -k(92)*n(idx_CHj)  &
        -k(95)*n(idx_CH2j)  &
        -k(101)*n(idx_Oj)  &
        -k(109)*n(idx_OHj)  &
        -k(110)*n(idx_H2Oj)  &
        -k(188)*n(idx_Ck)  &
        -k(191)*n(idx_Ok)  &
        -k(193)*n(idx_Hj)  &
        -k(194)*n(idx_Hj)  &
        -k(197)*n(idx_C)  &
        -k(201)*n(idx_Cj)  &
        -k(237)  &
        -k(252)  &
        -k(253)  &
        -k(254)  &
        -k(264)

    !d[C_dot]/d[H2]
    pd(8,7) =  &
        -k(56)*n(idx_C)  &
        -k(197)*n(idx_C)

    !d[O_dot]/d[H2]
    pd(9,7) =  &
        -k(70)*n(idx_O)

    !d[OH_dot]/d[H2]
    pd(10,7) =  &
        +k(70)*n(idx_O)  &
        -k(73)*n(idx_OH)  &
        +2.d0*k(81)*n(idx_O2)

    !d[CH_dot]/d[H2]
    pd(12,7) =  &
        +k(56)*n(idx_C)  &
        -k(58)*n(idx_CH)

    !d[CH2_dot]/d[H2]
    pd(13,7) =  &
        +k(58)*n(idx_CH)  &
        +k(188)*n(idx_Ck)  &
        +k(197)*n(idx_C)

    !d[H2O_dot]/d[H2]
    pd(16,7) =  &
        +k(73)*n(idx_OH)  &
        +k(191)*n(idx_Ok)

    !d[O2_dot]/d[H2]
    pd(17,7) =  &
        -k(81)*n(idx_O2)

    !d[H+_dot]/d[H2]
    pd(20,7) =  &
        +k(13)*n(idx_HEj)  &
        -k(21)*n(idx_Hj)  &
        -k(22)*n(idx_Hj)  &
        -k(193)*n(idx_Hj)  &
        +k(193)*n(idx_Hj)  &
        -k(194)*n(idx_Hj)  &
        +k(253)  &
        +k(264)

    !d[HE+_dot]/d[H2]
    pd(21,7) =  &
        -k(12)*n(idx_HEj)  &
        -k(13)*n(idx_HEj)  &
        -k(14)*n(idx_HEj)  &
        +k(14)*n(idx_HEj)

    !d[H2+_dot]/d[H2]
    pd(22,7) =  &
        +k(12)*n(idx_HEj)  &
        +k(21)*n(idx_Hj)  &
        +k(22)*n(idx_Hj)  &
        -k(85)*n(idx_H2j)  &
        +k(254)

    !d[C+_dot]/d[H2]
    pd(23,7) =  &
        -k(90)*n(idx_Cj)  &
        -k(201)*n(idx_Cj)

    !d[O+_dot]/d[H2]
    pd(24,7) =  &
        -k(101)*n(idx_Oj)

    !d[HOC+_dot]/d[H2]
    pd(25,7) =  &
        -k(53)*n(idx_HOCj)

    !d[HCO+_dot]/d[H2]
    pd(26,7) =  &
        +k(53)*n(idx_HOCj)

    !d[H3+_dot]/d[H2]
    pd(27,7) =  &
        +k(85)*n(idx_H2j)  &
        +k(194)*n(idx_Hj)

    !d[CH+_dot]/d[H2]
    pd(28,7) =  &
        +k(90)*n(idx_Cj)  &
        -k(92)*n(idx_CHj)

    !d[CH2+_dot]/d[H2]
    pd(29,7) =  &
        +k(92)*n(idx_CHj)  &
        -k(95)*n(idx_CH2j)  &
        +k(201)*n(idx_Cj)

    !d[CH3+_dot]/d[H2]
    pd(31,7) =  &
        +k(95)*n(idx_CH2j)

    !d[OH+_dot]/d[H2]
    pd(32,7) =  &
        +k(101)*n(idx_Oj)  &
        -k(109)*n(idx_OHj)

    !d[H2O+_dot]/d[H2]
    pd(33,7) =  &
        +k(109)*n(idx_OHj)  &
        -k(110)*n(idx_H2Oj)

    !d[H3O+_dot]/d[H2]
    pd(34,7) =  &
        +k(110)*n(idx_H2Oj)

    !d[Tgas_dot]/d[H2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(7)*1d-3
    if(dnn>0.d0) then
      nn(7) = n(7) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,7) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C]
    pd(1,8) =  &
        -k(42)*n(idx_E)  &
        +2.d0*k(42)*n(idx_E)  &
        +k(184)*n(idx_Hk)  &
        +k(192)*n(idx_Ok)  &
        -k(195)*n(idx_E)  &
        +k(217)  &
        +k(255)

    !d[H-_dot]/d[C]
    pd(2,8) =  &
        -k(184)*n(idx_Hk)

    !d[C-_dot]/d[C]
    pd(3,8) =  &
        +k(195)*n(idx_E)

    !d[O-_dot]/d[C]
    pd(4,8) =  &
        -k(192)*n(idx_Ok)

    !d[H_dot]/d[C]
    pd(5,8) =  &
        +k(47)*n(idx_Hj)  &
        +k(56)*n(idx_H2)  &
        +k(59)*n(idx_CH)  &
        +k(74)*n(idx_OH)  &
        +k(75)*n(idx_OH)  &
        +k(87)*n(idx_H2j)  &
        +k(89)*n(idx_H3j)  &
        -k(196)*n(idx_H)

    !d[HE_dot]/d[C]
    pd(6,8) =  &
        +k(49)*n(idx_HEj)  &
        +k(50)*n(idx_HEj)  &
        +k(51)*n(idx_HEj)

    !d[H2_dot]/d[C]
    pd(7,8) =  &
        -k(56)*n(idx_H2)  &
        +k(88)*n(idx_H3j)  &
        +k(117)*n(idx_H3Oj)  &
        -k(197)*n(idx_H2)

    !d[C_dot]/d[C]
    pd(8,8) =  &
        -k(42)*n(idx_E)  &
        -k(47)*n(idx_Hj)  &
        -k(49)*n(idx_HEj)  &
        -k(50)*n(idx_HEj)  &
        -k(51)*n(idx_HEj)  &
        -k(56)*n(idx_H2)  &
        -k(59)*n(idx_CH)  &
        -k(74)*n(idx_OH)  &
        -k(75)*n(idx_OH)  &
        -k(82)*n(idx_O2)  &
        -k(83)*n(idx_O2)  &
        -k(87)*n(idx_H2j)  &
        -k(88)*n(idx_H3j)  &
        -k(89)*n(idx_H3j)  &
        -k(117)*n(idx_H3Oj)  &
        -k(121)*n(idx_O2j)  &
        -k(122)*n(idx_O2j)  &
        -k(127)*n(idx_HCOj)  &
        -k(184)*n(idx_Hk)  &
        -k(192)*n(idx_Ok)  &
        -k(195)*n(idx_E)  &
        -k(196)*n(idx_H)  &
        -k(197)*n(idx_H2)  &
        -4.d0*k(198)*n(idx_C)  &
        -k(199)*n(idx_O)  &
        -k(217)  &
        -k(255)  &
        -4.d0*k(265)*n(idx_C)  &
        -4.d0*k(266)*n(idx_C)  &
        -k(267)*n(idx_O)  &
        -k(268)*n(idx_O)  &
        -k(271)*n(idx_Oj)  &
        -k(272)*n(idx_Oj)

    !d[O_dot]/d[C]
    pd(9,8) =  &
        +k(82)*n(idx_O2)  &
        +k(83)*n(idx_O2)  &
        +k(121)*n(idx_O2j)  &
        -k(199)*n(idx_O)  &
        -k(267)*n(idx_O)  &
        -k(268)*n(idx_O)

    !d[OH_dot]/d[C]
    pd(10,8) =  &
        -k(74)*n(idx_OH)  &
        -k(75)*n(idx_OH)

    !d[CO_dot]/d[C]
    pd(11,8) =  &
        +k(74)*n(idx_OH)  &
        +k(75)*n(idx_OH)  &
        +k(82)*n(idx_O2)  &
        +k(83)*n(idx_O2)  &
        +k(127)*n(idx_HCOj)  &
        +k(192)*n(idx_Ok)  &
        +k(199)*n(idx_O)  &
        +k(267)*n(idx_O)  &
        +k(268)*n(idx_O)

    !d[CH_dot]/d[C]
    pd(12,8) =  &
        +k(56)*n(idx_H2)  &
        -k(59)*n(idx_CH)  &
        +k(184)*n(idx_Hk)  &
        +k(196)*n(idx_H)

    !d[CH2_dot]/d[C]
    pd(13,8) =  &
        +k(197)*n(idx_H2)

    !d[C2_dot]/d[C]
    pd(14,8) =  &
        +k(59)*n(idx_CH)  &
        +2.d0*k(198)*n(idx_C)  &
        +2.d0*k(265)*n(idx_C)  &
        +2.d0*k(266)*n(idx_C)

    !d[O2_dot]/d[C]
    pd(17,8) =  &
        -k(82)*n(idx_O2)  &
        -k(83)*n(idx_O2)  &
        +k(122)*n(idx_O2j)

    !d[H+_dot]/d[C]
    pd(20,8) =  &
        -k(47)*n(idx_Hj)

    !d[HE+_dot]/d[C]
    pd(21,8) =  &
        -k(49)*n(idx_HEj)  &
        -k(50)*n(idx_HEj)  &
        -k(51)*n(idx_HEj)

    !d[H2+_dot]/d[C]
    pd(22,8) =  &
        -k(87)*n(idx_H2j)

    !d[C+_dot]/d[C]
    pd(23,8) =  &
        +k(42)*n(idx_E)  &
        +k(47)*n(idx_Hj)  &
        +k(49)*n(idx_HEj)  &
        +k(50)*n(idx_HEj)  &
        +k(51)*n(idx_HEj)  &
        +k(122)*n(idx_O2j)  &
        +k(217)  &
        +k(255)

    !d[O+_dot]/d[C]
    pd(24,8) =  &
        -k(271)*n(idx_Oj)  &
        -k(272)*n(idx_Oj)

    !d[HCO+_dot]/d[C]
    pd(26,8) =  &
        +k(117)*n(idx_H3Oj)  &
        -k(127)*n(idx_HCOj)

    !d[H3+_dot]/d[C]
    pd(27,8) =  &
        -k(88)*n(idx_H3j)  &
        -k(89)*n(idx_H3j)

    !d[CH+_dot]/d[C]
    pd(28,8) =  &
        +k(87)*n(idx_H2j)  &
        +k(88)*n(idx_H3j)  &
        +k(127)*n(idx_HCOj)

    !d[CH2+_dot]/d[C]
    pd(29,8) =  &
        +k(89)*n(idx_H3j)

    !d[CO+_dot]/d[C]
    pd(30,8) =  &
        +k(121)*n(idx_O2j)  &
        +k(271)*n(idx_Oj)  &
        +k(272)*n(idx_Oj)

    !d[H3O+_dot]/d[C]
    pd(34,8) =  &
        -k(117)*n(idx_H3Oj)

    !d[O2+_dot]/d[C]
    pd(35,8) =  &
        -k(121)*n(idx_O2j)  &
        -k(122)*n(idx_O2j)

    !d[Tgas_dot]/d[C]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(8)*1d-3
    if(dnn>0.d0) then
      nn(8) = n(8) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,8) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O]
    pd(1,9) =  &
        -k(43)*n(idx_E)  &
        +2.d0*k(43)*n(idx_E)  &
        +k(61)*n(idx_CH)  &
        +k(185)*n(idx_Hk)  &
        +k(189)*n(idx_Ck)  &
        -k(204)*n(idx_E)  &
        +k(248)

    !d[H-_dot]/d[O]
    pd(2,9) =  &
        -k(185)*n(idx_Hk)

    !d[C-_dot]/d[O]
    pd(3,9) =  &
        -k(189)*n(idx_Ck)

    !d[O-_dot]/d[O]
    pd(4,9) =  &
        +k(204)*n(idx_E)

    !d[H_dot]/d[O]
    pd(5,9) =  &
        +k(45)*n(idx_Hj)  &
        +k(60)*n(idx_CH)  &
        +2.d0*k(64)*n(idx_CH2)  &
        +k(66)*n(idx_CH2)  &
        +k(70)*n(idx_H2)  &
        +k(76)*n(idx_OH)  &
        +k(77)*n(idx_OH)  &
        +k(93)*n(idx_CHj)  &
        +k(96)*n(idx_CH2j)  &
        +k(102)*n(idx_H2j)  &
        +k(104)*n(idx_H3j)  &
        -k(205)*n(idx_H)  &
        -k(273)*n(idx_H)

    !d[HE_dot]/d[O]
    pd(6,9) =  &
        +k(46)*n(idx_HEj)

    !d[H2_dot]/d[O]
    pd(7,9) =  &
        +k(65)*n(idx_CH2)  &
        -k(70)*n(idx_H2)  &
        +k(98)*n(idx_CH3j)  &
        +k(99)*n(idx_CH3j)  &
        +k(103)*n(idx_H3j)

    !d[C_dot]/d[O]
    pd(8,9) =  &
        +k(62)*n(idx_CH)  &
        +k(68)*n(idx_C2)  &
        +k(69)*n(idx_C2)  &
        -k(199)*n(idx_C)  &
        -k(267)*n(idx_C)  &
        -k(268)*n(idx_C)

    !d[O_dot]/d[O]
    pd(9,9) =  &
        -k(43)*n(idx_E)  &
        -k(45)*n(idx_Hj)  &
        -k(46)*n(idx_HEj)  &
        -k(60)*n(idx_CH)  &
        -k(61)*n(idx_CH)  &
        -k(62)*n(idx_CH)  &
        -k(64)*n(idx_CH2)  &
        -k(65)*n(idx_CH2)  &
        -k(66)*n(idx_CH2)  &
        -k(67)*n(idx_CH2)  &
        -k(68)*n(idx_C2)  &
        -k(69)*n(idx_C2)  &
        -k(70)*n(idx_H2)  &
        -k(76)*n(idx_OH)  &
        -k(77)*n(idx_OH)  &
        -k(93)*n(idx_CHj)  &
        -k(96)*n(idx_CH2j)  &
        -k(98)*n(idx_CH3j)  &
        -k(99)*n(idx_CH3j)  &
        -k(102)*n(idx_H2j)  &
        -k(103)*n(idx_H3j)  &
        -k(104)*n(idx_H3j)  &
        -k(185)*n(idx_Hk)  &
        -k(189)*n(idx_Ck)  &
        -k(199)*n(idx_C)  &
        -k(202)*n(idx_Cj)  &
        -k(203)*n(idx_Cj)  &
        -k(204)*n(idx_E)  &
        -k(205)*n(idx_H)  &
        -4.d0*k(206)*n(idx_O)  &
        -k(248)  &
        -k(267)*n(idx_C)  &
        -k(268)*n(idx_C)  &
        -k(269)*n(idx_Cj)  &
        -k(270)*n(idx_Cj)  &
        -k(273)*n(idx_H)  &
        -4.d0*k(275)*n(idx_O)

    !d[OH_dot]/d[O]
    pd(10,9) =  &
        +k(62)*n(idx_CH)  &
        +k(67)*n(idx_CH2)  &
        +k(70)*n(idx_H2)  &
        -k(76)*n(idx_OH)  &
        -k(77)*n(idx_OH)  &
        +k(185)*n(idx_Hk)  &
        +k(205)*n(idx_H)  &
        +k(273)*n(idx_H)

    !d[CO_dot]/d[O]
    pd(11,9) =  &
        +k(60)*n(idx_CH)  &
        +k(64)*n(idx_CH2)  &
        +k(65)*n(idx_CH2)  &
        +k(68)*n(idx_C2)  &
        +k(69)*n(idx_C2)  &
        +k(189)*n(idx_Ck)  &
        +k(199)*n(idx_C)  &
        +k(267)*n(idx_C)  &
        +k(268)*n(idx_C)

    !d[CH_dot]/d[O]
    pd(12,9) =  &
        -k(60)*n(idx_CH)  &
        -k(61)*n(idx_CH)  &
        -k(62)*n(idx_CH)  &
        +k(67)*n(idx_CH2)

    !d[CH2_dot]/d[O]
    pd(13,9) =  &
        -k(64)*n(idx_CH2)  &
        -k(65)*n(idx_CH2)  &
        -k(66)*n(idx_CH2)  &
        -k(67)*n(idx_CH2)

    !d[C2_dot]/d[O]
    pd(14,9) =  &
        -k(68)*n(idx_C2)  &
        -k(69)*n(idx_C2)

    !d[HCO_dot]/d[O]
    pd(15,9) =  &
        +k(66)*n(idx_CH2)

    !d[O2_dot]/d[O]
    pd(17,9) =  &
        +k(76)*n(idx_OH)  &
        +k(77)*n(idx_OH)  &
        +2.d0*k(206)*n(idx_O)  &
        +2.d0*k(275)*n(idx_O)

    !d[H+_dot]/d[O]
    pd(20,9) =  &
        -k(45)*n(idx_Hj)

    !d[HE+_dot]/d[O]
    pd(21,9) =  &
        -k(46)*n(idx_HEj)

    !d[H2+_dot]/d[O]
    pd(22,9) =  &
        -k(102)*n(idx_H2j)

    !d[C+_dot]/d[O]
    pd(23,9) =  &
        -k(202)*n(idx_Cj)  &
        -k(203)*n(idx_Cj)  &
        -k(269)*n(idx_Cj)  &
        -k(270)*n(idx_Cj)

    !d[O+_dot]/d[O]
    pd(24,9) =  &
        +k(43)*n(idx_E)  &
        +k(45)*n(idx_Hj)  &
        +k(46)*n(idx_HEj)  &
        +k(248)

    !d[HOC+_dot]/d[O]
    pd(25,9) =  &
        +k(98)*n(idx_CH3j)

    !d[HCO+_dot]/d[O]
    pd(26,9) =  &
        +k(61)*n(idx_CH)  &
        +k(96)*n(idx_CH2j)  &
        +k(99)*n(idx_CH3j)

    !d[H3+_dot]/d[O]
    pd(27,9) =  &
        -k(103)*n(idx_H3j)  &
        -k(104)*n(idx_H3j)

    !d[CH+_dot]/d[O]
    pd(28,9) =  &
        -k(93)*n(idx_CHj)

    !d[CH2+_dot]/d[O]
    pd(29,9) =  &
        -k(96)*n(idx_CH2j)

    !d[CO+_dot]/d[O]
    pd(30,9) =  &
        +k(93)*n(idx_CHj)  &
        +k(202)*n(idx_Cj)  &
        +k(203)*n(idx_Cj)  &
        +k(269)*n(idx_Cj)  &
        +k(270)*n(idx_Cj)

    !d[CH3+_dot]/d[O]
    pd(31,9) =  &
        -k(98)*n(idx_CH3j)  &
        -k(99)*n(idx_CH3j)

    !d[OH+_dot]/d[O]
    pd(32,9) =  &
        +k(102)*n(idx_H2j)  &
        +k(103)*n(idx_H3j)

    !d[H2O+_dot]/d[O]
    pd(33,9) =  &
        +k(104)*n(idx_H3j)

    !d[Tgas_dot]/d[O]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(9)*1d-3
    if(dnn>0.d0) then
      nn(9) = n(9) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,9) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[OH]
    pd(1,10) =  &
        +k(186)*n(idx_Hk)  &
        +k(230)

    !d[H-_dot]/d[OH]
    pd(2,10) =  &
        -k(186)*n(idx_Hk)

    !d[H_dot]/d[OH]
    pd(5,10) =  &
        -k(52)*n(idx_H)  &
        +2.d0*k(52)*n(idx_H)  &
        -k(71)*n(idx_H)  &
        -k(72)*n(idx_H)  &
        +k(73)*n(idx_H2)  &
        +k(74)*n(idx_C)  &
        +k(75)*n(idx_C)  &
        +k(76)*n(idx_O)  &
        +k(77)*n(idx_O)  &
        +k(107)*n(idx_Cj)  &
        +k(108)*n(idx_Cj)  &
        +k(141)*n(idx_Hj)  &
        +k(142)*n(idx_Hj)  &
        +k(143)*n(idx_HEj)  &
        +k(144)*n(idx_HEj)  &
        -k(207)*n(idx_H)  &
        +k(229)  &
        +k(259)  &
        -k(274)*n(idx_H)

    !d[HE_dot]/d[OH]
    pd(6,10) =  &
        +k(143)*n(idx_HEj)  &
        +k(144)*n(idx_HEj)

    !d[H2_dot]/d[OH]
    pd(7,10) =  &
        +k(71)*n(idx_H)  &
        +k(72)*n(idx_H)  &
        -k(73)*n(idx_H2)  &
        +k(105)*n(idx_H3j)  &
        +k(106)*n(idx_H3j)

    !d[C_dot]/d[OH]
    pd(8,10) =  &
        -k(74)*n(idx_C)  &
        -k(75)*n(idx_C)

    !d[O_dot]/d[OH]
    pd(9,10) =  &
        +k(52)*n(idx_H)  &
        +k(71)*n(idx_H)  &
        +k(72)*n(idx_H)  &
        -k(76)*n(idx_O)  &
        -k(77)*n(idx_O)  &
        +2.d0*k(78)*n(idx_OH)  &
        +k(229)  &
        +k(259)

    !d[OH_dot]/d[OH]
    pd(10,10) =  &
        -k(52)*n(idx_H)  &
        -k(71)*n(idx_H)  &
        -k(72)*n(idx_H)  &
        -k(73)*n(idx_H2)  &
        -k(74)*n(idx_C)  &
        -k(75)*n(idx_C)  &
        -k(76)*n(idx_O)  &
        -k(77)*n(idx_O)  &
        -4.d0*k(78)*n(idx_OH)  &
        -k(105)*n(idx_H3j)  &
        -k(106)*n(idx_H3j)  &
        -k(107)*n(idx_Cj)  &
        -k(108)*n(idx_Cj)  &
        -k(141)*n(idx_Hj)  &
        -k(142)*n(idx_Hj)  &
        -k(143)*n(idx_HEj)  &
        -k(144)*n(idx_HEj)  &
        -k(186)*n(idx_Hk)  &
        -k(207)*n(idx_H)  &
        -k(229)  &
        -k(230)  &
        -k(259)  &
        -k(274)*n(idx_H)

    !d[CO_dot]/d[OH]
    pd(11,10) =  &
        +k(74)*n(idx_C)  &
        +k(75)*n(idx_C)

    !d[H2O_dot]/d[OH]
    pd(16,10) =  &
        +k(73)*n(idx_H2)  &
        +2.d0*k(78)*n(idx_OH)  &
        +k(186)*n(idx_Hk)  &
        +k(207)*n(idx_H)  &
        +k(274)*n(idx_H)

    !d[O2_dot]/d[OH]
    pd(17,10) =  &
        +k(76)*n(idx_O)  &
        +k(77)*n(idx_O)

    !d[H+_dot]/d[OH]
    pd(20,10) =  &
        -k(141)*n(idx_Hj)  &
        -k(142)*n(idx_Hj)

    !d[HE+_dot]/d[OH]
    pd(21,10) =  &
        -k(143)*n(idx_HEj)  &
        -k(144)*n(idx_HEj)

    !d[C+_dot]/d[OH]
    pd(23,10) =  &
        -k(107)*n(idx_Cj)  &
        -k(108)*n(idx_Cj)

    !d[O+_dot]/d[OH]
    pd(24,10) =  &
        +k(143)*n(idx_HEj)  &
        +k(144)*n(idx_HEj)

    !d[H3+_dot]/d[OH]
    pd(27,10) =  &
        -k(105)*n(idx_H3j)  &
        -k(106)*n(idx_H3j)

    !d[CO+_dot]/d[OH]
    pd(30,10) =  &
        +k(107)*n(idx_Cj)  &
        +k(108)*n(idx_Cj)

    !d[OH+_dot]/d[OH]
    pd(32,10) =  &
        +k(141)*n(idx_Hj)  &
        +k(142)*n(idx_Hj)  &
        +k(230)

    !d[H2O+_dot]/d[OH]
    pd(33,10) =  &
        +k(105)*n(idx_H3j)  &
        +k(106)*n(idx_H3j)

    !d[Tgas_dot]/d[OH]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(10)*1d-3
    if(dnn>0.d0) then
      nn(10) = n(10) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,10) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CO]
    pd(1,11) =  &
        +k(250)

    !d[H_dot]/d[CO]
    pd(5,11) =  &
        -k(84)*n(idx_H)

    !d[HE_dot]/d[CO]
    pd(6,11) =  &
        +k(156)*n(idx_HEj)  &
        +k(157)*n(idx_HEj)

    !d[H2_dot]/d[CO]
    pd(7,11) =  &
        +k(123)*n(idx_H3j)  &
        +k(124)*n(idx_H3j)  &
        +k(125)*n(idx_H3j)  &
        +k(126)*n(idx_H3j)

    !d[C_dot]/d[CO]
    pd(8,11) =  &
        +k(84)*n(idx_H)  &
        +k(157)*n(idx_HEj)  &
        +k(236)  &
        +k(249)

    !d[O_dot]/d[CO]
    pd(9,11) =  &
        +k(156)*n(idx_HEj)  &
        +k(236)  &
        +k(249)

    !d[OH_dot]/d[CO]
    pd(10,11) =  &
        +k(84)*n(idx_H)

    !d[CO_dot]/d[CO]
    pd(11,11) =  &
        -k(54)*n(idx_HOCj)  &
        +k(54)*n(idx_HOCj)  &
        -k(55)*n(idx_HOCj)  &
        +k(55)*n(idx_HOCj)  &
        -k(84)*n(idx_H)  &
        -k(123)*n(idx_H3j)  &
        -k(124)*n(idx_H3j)  &
        -k(125)*n(idx_H3j)  &
        -k(126)*n(idx_H3j)  &
        -k(156)*n(idx_HEj)  &
        -k(157)*n(idx_HEj)  &
        -k(236)  &
        -k(249)  &
        -k(250)

    !d[HE+_dot]/d[CO]
    pd(21,11) =  &
        -k(156)*n(idx_HEj)  &
        -k(157)*n(idx_HEj)

    !d[C+_dot]/d[CO]
    pd(23,11) =  &
        +k(156)*n(idx_HEj)

    !d[O+_dot]/d[CO]
    pd(24,11) =  &
        +k(157)*n(idx_HEj)

    !d[HOC+_dot]/d[CO]
    pd(25,11) =  &
        -k(54)*n(idx_HOCj)  &
        -k(55)*n(idx_HOCj)  &
        +k(125)*n(idx_H3j)  &
        +k(126)*n(idx_H3j)

    !d[HCO+_dot]/d[CO]
    pd(26,11) =  &
        +k(54)*n(idx_HOCj)  &
        +k(55)*n(idx_HOCj)  &
        +k(123)*n(idx_H3j)  &
        +k(124)*n(idx_H3j)

    !d[H3+_dot]/d[CO]
    pd(27,11) =  &
        -k(123)*n(idx_H3j)  &
        -k(124)*n(idx_H3j)  &
        -k(125)*n(idx_H3j)  &
        -k(126)*n(idx_H3j)

    !d[CO+_dot]/d[CO]
    pd(30,11) =  &
        +k(250)

    !d[Tgas_dot]/d[CO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(11)*1d-3
    if(dnn>0.d0) then
      nn(11) = n(11) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,11) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH]
    pd(1,12) =  &
        +k(61)*n(idx_O)  &
        +k(220)

    !d[H_dot]/d[CH]
    pd(5,12) =  &
        -k(57)*n(idx_H)  &
        +k(58)*n(idx_H2)  &
        +k(59)*n(idx_C)  &
        +k(60)*n(idx_O)  &
        +k(130)*n(idx_Hj)  &
        +k(131)*n(idx_Hj)  &
        +k(219)  &
        +k(256)

    !d[H2_dot]/d[CH]
    pd(7,12) =  &
        +k(57)*n(idx_H)  &
        -k(58)*n(idx_H2)

    !d[C_dot]/d[CH]
    pd(8,12) =  &
        +k(57)*n(idx_H)  &
        -k(59)*n(idx_C)  &
        +k(62)*n(idx_O)  &
        +k(219)  &
        +k(256)

    !d[O_dot]/d[CH]
    pd(9,12) =  &
        -k(60)*n(idx_O)  &
        -k(61)*n(idx_O)  &
        -k(62)*n(idx_O)

    !d[OH_dot]/d[CH]
    pd(10,12) =  &
        +k(62)*n(idx_O)

    !d[CO_dot]/d[CH]
    pd(11,12) =  &
        +k(60)*n(idx_O)

    !d[CH_dot]/d[CH]
    pd(12,12) =  &
        -k(57)*n(idx_H)  &
        -k(58)*n(idx_H2)  &
        -k(59)*n(idx_C)  &
        -k(60)*n(idx_O)  &
        -k(61)*n(idx_O)  &
        -k(62)*n(idx_O)  &
        -k(130)*n(idx_Hj)  &
        -k(131)*n(idx_Hj)  &
        -k(219)  &
        -k(220)  &
        -k(256)

    !d[CH2_dot]/d[CH]
    pd(13,12) =  &
        +k(58)*n(idx_H2)

    !d[C2_dot]/d[CH]
    pd(14,12) =  &
        +k(59)*n(idx_C)

    !d[H+_dot]/d[CH]
    pd(20,12) =  &
        -k(130)*n(idx_Hj)  &
        -k(131)*n(idx_Hj)

    !d[HCO+_dot]/d[CH]
    pd(26,12) =  &
        +k(61)*n(idx_O)

    !d[CH+_dot]/d[CH]
    pd(28,12) =  &
        +k(130)*n(idx_Hj)  &
        +k(131)*n(idx_Hj)  &
        +k(220)

    !d[Tgas_dot]/d[CH]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(12)*1d-3
    if(dnn>0.d0) then
      nn(12) = n(12) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,12) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH2]
    pd(1,13) =  &
        +k(223)  &
        +k(260)

    !d[H_dot]/d[CH2]
    pd(5,13) =  &
        -k(63)*n(idx_H)  &
        +2.d0*k(64)*n(idx_O)  &
        +k(66)*n(idx_O)  &
        +k(134)*n(idx_Hj)  &
        +k(135)*n(idx_Hj)  &
        +k(138)*n(idx_HEj)  &
        +k(139)*n(idx_HEj)  &
        +k(222)

    !d[HE_dot]/d[CH2]
    pd(6,13) =  &
        +k(136)*n(idx_HEj)  &
        +k(137)*n(idx_HEj)  &
        +k(138)*n(idx_HEj)  &
        +k(139)*n(idx_HEj)

    !d[H2_dot]/d[CH2]
    pd(7,13) =  &
        +k(63)*n(idx_H)  &
        +k(65)*n(idx_O)  &
        +k(132)*n(idx_Hj)  &
        +k(133)*n(idx_Hj)  &
        +k(136)*n(idx_HEj)  &
        +k(137)*n(idx_HEj)

    !d[O_dot]/d[CH2]
    pd(9,13) =  &
        -k(64)*n(idx_O)  &
        -k(65)*n(idx_O)  &
        -k(66)*n(idx_O)  &
        -k(67)*n(idx_O)

    !d[OH_dot]/d[CH2]
    pd(10,13) =  &
        +k(67)*n(idx_O)

    !d[CO_dot]/d[CH2]
    pd(11,13) =  &
        +k(64)*n(idx_O)  &
        +k(65)*n(idx_O)

    !d[CH_dot]/d[CH2]
    pd(12,13) =  &
        +k(63)*n(idx_H)  &
        +k(67)*n(idx_O)  &
        +k(222)

    !d[CH2_dot]/d[CH2]
    pd(13,13) =  &
        -k(63)*n(idx_H)  &
        -k(64)*n(idx_O)  &
        -k(65)*n(idx_O)  &
        -k(66)*n(idx_O)  &
        -k(67)*n(idx_O)  &
        -k(132)*n(idx_Hj)  &
        -k(133)*n(idx_Hj)  &
        -k(134)*n(idx_Hj)  &
        -k(135)*n(idx_Hj)  &
        -k(136)*n(idx_HEj)  &
        -k(137)*n(idx_HEj)  &
        -k(138)*n(idx_HEj)  &
        -k(139)*n(idx_HEj)  &
        -k(222)  &
        -k(223)  &
        -k(260)

    !d[HCO_dot]/d[CH2]
    pd(15,13) =  &
        +k(66)*n(idx_O)

    !d[H+_dot]/d[CH2]
    pd(20,13) =  &
        -k(132)*n(idx_Hj)  &
        -k(133)*n(idx_Hj)  &
        -k(134)*n(idx_Hj)  &
        -k(135)*n(idx_Hj)

    !d[HE+_dot]/d[CH2]
    pd(21,13) =  &
        -k(136)*n(idx_HEj)  &
        -k(137)*n(idx_HEj)  &
        -k(138)*n(idx_HEj)  &
        -k(139)*n(idx_HEj)

    !d[C+_dot]/d[CH2]
    pd(23,13) =  &
        +k(136)*n(idx_HEj)  &
        +k(137)*n(idx_HEj)

    !d[CH+_dot]/d[CH2]
    pd(28,13) =  &
        +k(132)*n(idx_Hj)  &
        +k(133)*n(idx_Hj)  &
        +k(138)*n(idx_HEj)  &
        +k(139)*n(idx_HEj)

    !d[CH2+_dot]/d[CH2]
    pd(29,13) =  &
        +k(134)*n(idx_Hj)  &
        +k(135)*n(idx_Hj)  &
        +k(223)  &
        +k(260)

    !d[Tgas_dot]/d[CH2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(13)*1d-3
    if(dnn>0.d0) then
      nn(13) = n(13) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,13) = (dn1-dn0)/dnn
    end if

    !d[HE_dot]/d[C2]
    pd(6,14) =  &
        +k(140)*n(idx_HEj)

    !d[C_dot]/d[C2]
    pd(8,14) =  &
        +k(68)*n(idx_O)  &
        +k(69)*n(idx_O)  &
        +k(100)*n(idx_Oj)  &
        +k(140)*n(idx_HEj)  &
        +2.d0*k(227)  &
        +2.d0*k(251)

    !d[O_dot]/d[C2]
    pd(9,14) =  &
        -k(68)*n(idx_O)  &
        -k(69)*n(idx_O)

    !d[CO_dot]/d[C2]
    pd(11,14) =  &
        +k(68)*n(idx_O)  &
        +k(69)*n(idx_O)

    !d[C2_dot]/d[C2]
    pd(14,14) =  &
        -k(68)*n(idx_O)  &
        -k(69)*n(idx_O)  &
        -k(100)*n(idx_Oj)  &
        -k(140)*n(idx_HEj)  &
        -k(227)  &
        -k(251)

    !d[HE+_dot]/d[C2]
    pd(21,14) =  &
        -k(140)*n(idx_HEj)

    !d[C+_dot]/d[C2]
    pd(23,14) =  &
        +k(140)*n(idx_HEj)

    !d[O+_dot]/d[C2]
    pd(24,14) =  &
        -k(100)*n(idx_Oj)

    !d[CO+_dot]/d[C2]
    pd(30,14) =  &
        +k(100)*n(idx_Oj)

    !d[Tgas_dot]/d[C2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(14)*1d-3
    if(dnn>0.d0) then
      nn(14) = n(14) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,14) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HCO]
    pd(1,15) =  &
        +k(263)

    !d[H_dot]/d[HCO]
    pd(5,15) =  &
        +k(262)

    !d[CO_dot]/d[HCO]
    pd(11,15) =  &
        +k(262)

    !d[HCO_dot]/d[HCO]
    pd(15,15) =  &
        -k(262)  &
        -k(263)

    !d[HCO+_dot]/d[HCO]
    pd(26,15) =  &
        +k(263)

    !d[Tgas_dot]/d[HCO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(15)*1d-3
    if(dnn>0.d0) then
      nn(15) = n(15) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,15) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2O]
    pd(1,16) =  &
        +k(233)

    !d[H_dot]/d[H2O]
    pd(5,16) =  &
        -k(79)*n(idx_H)  &
        +k(113)*n(idx_Cj)  &
        +k(114)*n(idx_Cj)  &
        +k(115)*n(idx_Cj)  &
        +k(145)*n(idx_Hj)  &
        +k(146)*n(idx_Hj)  &
        +k(149)*n(idx_HEj)  &
        +k(150)*n(idx_HEj)  &
        +k(232)  &
        +k(261)

    !d[HE_dot]/d[H2O]
    pd(6,16) =  &
        +k(147)*n(idx_HEj)  &
        +k(148)*n(idx_HEj)  &
        +k(149)*n(idx_HEj)  &
        +k(150)*n(idx_HEj)  &
        +k(151)*n(idx_HEj)  &
        +k(152)*n(idx_HEj)

    !d[H2_dot]/d[H2O]
    pd(7,16) =  &
        +k(79)*n(idx_H)  &
        +k(111)*n(idx_H3j)  &
        +k(112)*n(idx_H3j)

    !d[C_dot]/d[H2O]
    pd(8,16) =  &
        +k(116)*n(idx_Cj)

    !d[OH_dot]/d[H2O]
    pd(10,16) =  &
        +k(79)*n(idx_H)  &
        +k(147)*n(idx_HEj)  &
        +k(148)*n(idx_HEj)  &
        +k(232)  &
        +k(261)

    !d[CO_dot]/d[H2O]
    pd(11,16) =  &
        +k(128)*n(idx_HCOj)  &
        +k(129)*n(idx_HCOj)

    !d[H2O_dot]/d[H2O]
    pd(16,16) =  &
        -k(79)*n(idx_H)  &
        -k(111)*n(idx_H3j)  &
        -k(112)*n(idx_H3j)  &
        -k(113)*n(idx_Cj)  &
        -k(114)*n(idx_Cj)  &
        -k(115)*n(idx_Cj)  &
        -k(116)*n(idx_Cj)  &
        -k(128)*n(idx_HCOj)  &
        -k(129)*n(idx_HCOj)  &
        -k(145)*n(idx_Hj)  &
        -k(146)*n(idx_Hj)  &
        -k(147)*n(idx_HEj)  &
        -k(148)*n(idx_HEj)  &
        -k(149)*n(idx_HEj)  &
        -k(150)*n(idx_HEj)  &
        -k(151)*n(idx_HEj)  &
        -k(152)*n(idx_HEj)  &
        -k(232)  &
        -k(233)  &
        -k(261)

    !d[H+_dot]/d[H2O]
    pd(20,16) =  &
        -k(145)*n(idx_Hj)  &
        -k(146)*n(idx_Hj)  &
        +k(147)*n(idx_HEj)  &
        +k(148)*n(idx_HEj)

    !d[HE+_dot]/d[H2O]
    pd(21,16) =  &
        -k(147)*n(idx_HEj)  &
        -k(148)*n(idx_HEj)  &
        -k(149)*n(idx_HEj)  &
        -k(150)*n(idx_HEj)  &
        -k(151)*n(idx_HEj)  &
        -k(152)*n(idx_HEj)

    !d[C+_dot]/d[H2O]
    pd(23,16) =  &
        -k(113)*n(idx_Cj)  &
        -k(114)*n(idx_Cj)  &
        -k(115)*n(idx_Cj)  &
        -k(116)*n(idx_Cj)

    !d[HOC+_dot]/d[H2O]
    pd(25,16) =  &
        +k(113)*n(idx_Cj)

    !d[HCO+_dot]/d[H2O]
    pd(26,16) =  &
        +k(114)*n(idx_Cj)  &
        +k(115)*n(idx_Cj)  &
        -k(128)*n(idx_HCOj)  &
        -k(129)*n(idx_HCOj)

    !d[H3+_dot]/d[H2O]
    pd(27,16) =  &
        -k(111)*n(idx_H3j)  &
        -k(112)*n(idx_H3j)

    !d[OH+_dot]/d[H2O]
    pd(32,16) =  &
        +k(149)*n(idx_HEj)  &
        +k(150)*n(idx_HEj)

    !d[H2O+_dot]/d[H2O]
    pd(33,16) =  &
        +k(116)*n(idx_Cj)  &
        +k(145)*n(idx_Hj)  &
        +k(146)*n(idx_Hj)  &
        +k(151)*n(idx_HEj)  &
        +k(152)*n(idx_HEj)  &
        +k(233)

    !d[H3O+_dot]/d[H2O]
    pd(34,16) =  &
        +k(111)*n(idx_H3j)  &
        +k(112)*n(idx_H3j)  &
        +k(128)*n(idx_HCOj)  &
        +k(129)*n(idx_HCOj)

    !d[Tgas_dot]/d[H2O]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(16)*1d-3
    if(dnn>0.d0) then
      nn(16) = n(16) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,16) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O2]
    pd(1,17) =  &
        +k(234)  &
        +k(258)

    !d[H_dot]/d[O2]
    pd(5,17) =  &
        -k(80)*n(idx_H)  &
        +k(153)*n(idx_Hj)

    !d[HE_dot]/d[O2]
    pd(6,17) =  &
        +k(154)*n(idx_HEj)  &
        +k(155)*n(idx_HEj)

    !d[H2_dot]/d[O2]
    pd(7,17) =  &
        -k(81)*n(idx_H2)

    !d[C_dot]/d[O2]
    pd(8,17) =  &
        -k(82)*n(idx_C)  &
        -k(83)*n(idx_C)

    !d[O_dot]/d[O2]
    pd(9,17) =  &
        +k(80)*n(idx_H)  &
        +k(82)*n(idx_C)  &
        +k(83)*n(idx_C)  &
        +k(118)*n(idx_Cj)  &
        +k(155)*n(idx_HEj)  &
        +2.d0*k(235)  &
        +2.d0*k(257)

    !d[OH_dot]/d[O2]
    pd(10,17) =  &
        +k(80)*n(idx_H)  &
        +2.d0*k(81)*n(idx_H2)  &
        +k(120)*n(idx_CH2j)

    !d[CO_dot]/d[O2]
    pd(11,17) =  &
        +k(82)*n(idx_C)  &
        +k(83)*n(idx_C)  &
        +k(119)*n(idx_Cj)

    !d[O2_dot]/d[O2]
    pd(17,17) =  &
        -k(80)*n(idx_H)  &
        -k(81)*n(idx_H2)  &
        -k(82)*n(idx_C)  &
        -k(83)*n(idx_C)  &
        -k(118)*n(idx_Cj)  &
        -k(119)*n(idx_Cj)  &
        -k(120)*n(idx_CH2j)  &
        -k(153)*n(idx_Hj)  &
        -k(154)*n(idx_HEj)  &
        -k(155)*n(idx_HEj)  &
        -k(234)  &
        -k(235)  &
        -k(257)  &
        -k(258)

    !d[H+_dot]/d[O2]
    pd(20,17) =  &
        -k(153)*n(idx_Hj)

    !d[HE+_dot]/d[O2]
    pd(21,17) =  &
        -k(154)*n(idx_HEj)  &
        -k(155)*n(idx_HEj)

    !d[C+_dot]/d[O2]
    pd(23,17) =  &
        -k(118)*n(idx_Cj)  &
        -k(119)*n(idx_Cj)

    !d[O+_dot]/d[O2]
    pd(24,17) =  &
        +k(119)*n(idx_Cj)  &
        +k(155)*n(idx_HEj)

    !d[HCO+_dot]/d[O2]
    pd(26,17) =  &
        +k(120)*n(idx_CH2j)

    !d[CH2+_dot]/d[O2]
    pd(29,17) =  &
        -k(120)*n(idx_CH2j)

    !d[CO+_dot]/d[O2]
    pd(30,17) =  &
        +k(118)*n(idx_Cj)

    !d[O2+_dot]/d[O2]
    pd(35,17) =  &
        +k(153)*n(idx_Hj)  &
        +k(154)*n(idx_HEj)  &
        +k(234)  &
        +k(258)

    !d[Tgas_dot]/d[O2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(17)*1d-3
    if(dnn>0.d0) then
      nn(17) = n(17) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,17) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[CO_total]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(18)*1d-3
    if(dnn>0.d0) then
      nn(18) = n(18) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,18) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[H2O_total]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(19)*1d-3
    if(dnn>0.d0) then
      nn(19) = n(19) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,19) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H+]
    pd(1,20) =  &
        -k(2)*n(idx_E)  &
        -k(3)*n(idx_E)  &
        +k(29)*n(idx_Hk)

    !d[H-_dot]/d[H+]
    pd(2,20) =  &
        -k(28)*n(idx_Hk)  &
        -k(29)*n(idx_Hk)

    !d[C-_dot]/d[H+]
    pd(3,20) =  &
        -k(159)*n(idx_Ck)

    !d[O-_dot]/d[H+]
    pd(4,20) =  &
        -k(160)*n(idx_Ok)

    !d[H_dot]/d[H+]
    pd(5,20) =  &
        +k(2)*n(idx_E)  &
        +k(3)*n(idx_E)  &
        +k(9)*n(idx_HE)  &
        +k(10)*n(idx_HE)  &
        -k(18)*n(idx_H)  &
        -k(19)*n(idx_H)  &
        +k(21)*n(idx_H2)  &
        +k(22)*n(idx_H2)  &
        +2.d0*k(28)*n(idx_Hk)  &
        +k(45)*n(idx_O)  &
        +k(47)*n(idx_C)  &
        +k(130)*n(idx_CH)  &
        +k(131)*n(idx_CH)  &
        +k(134)*n(idx_CH2)  &
        +k(135)*n(idx_CH2)  &
        +k(141)*n(idx_OH)  &
        +k(142)*n(idx_OH)  &
        +k(145)*n(idx_H2O)  &
        +k(146)*n(idx_H2O)  &
        +k(153)*n(idx_O2)  &
        +k(159)*n(idx_Ck)  &
        +k(160)*n(idx_Ok)  &
        +2.d0*k(193)*n(idx_H2)

    !d[HE_dot]/d[H+]
    pd(6,20) =  &
        -k(9)*n(idx_HE)  &
        -k(10)*n(idx_HE)

    !d[H2_dot]/d[H+]
    pd(7,20) =  &
        -k(21)*n(idx_H2)  &
        -k(22)*n(idx_H2)  &
        +k(132)*n(idx_CH2)  &
        +k(133)*n(idx_CH2)  &
        -k(193)*n(idx_H2)  &
        -k(194)*n(idx_H2)

    !d[C_dot]/d[H+]
    pd(8,20) =  &
        -k(47)*n(idx_C)  &
        +k(159)*n(idx_Ck)

    !d[O_dot]/d[H+]
    pd(9,20) =  &
        -k(45)*n(idx_O)  &
        +k(160)*n(idx_Ok)

    !d[OH_dot]/d[H+]
    pd(10,20) =  &
        -k(141)*n(idx_OH)  &
        -k(142)*n(idx_OH)

    !d[CH_dot]/d[H+]
    pd(12,20) =  &
        -k(130)*n(idx_CH)  &
        -k(131)*n(idx_CH)

    !d[CH2_dot]/d[H+]
    pd(13,20) =  &
        -k(132)*n(idx_CH2)  &
        -k(133)*n(idx_CH2)  &
        -k(134)*n(idx_CH2)  &
        -k(135)*n(idx_CH2)

    !d[H2O_dot]/d[H+]
    pd(16,20) =  &
        -k(145)*n(idx_H2O)  &
        -k(146)*n(idx_H2O)

    !d[O2_dot]/d[H+]
    pd(17,20) =  &
        -k(153)*n(idx_O2)

    !d[H+_dot]/d[H+]
    pd(20,20) =  &
        -k(2)*n(idx_E)  &
        -k(3)*n(idx_E)  &
        -k(9)*n(idx_HE)  &
        -k(10)*n(idx_HE)  &
        -k(18)*n(idx_H)  &
        -k(19)*n(idx_H)  &
        -k(21)*n(idx_H2)  &
        -k(22)*n(idx_H2)  &
        -k(28)*n(idx_Hk)  &
        -k(29)*n(idx_Hk)  &
        -k(45)*n(idx_O)  &
        -k(47)*n(idx_C)  &
        -k(130)*n(idx_CH)  &
        -k(131)*n(idx_CH)  &
        -k(132)*n(idx_CH2)  &
        -k(133)*n(idx_CH2)  &
        -k(134)*n(idx_CH2)  &
        -k(135)*n(idx_CH2)  &
        -k(141)*n(idx_OH)  &
        -k(142)*n(idx_OH)  &
        -k(145)*n(idx_H2O)  &
        -k(146)*n(idx_H2O)  &
        -k(153)*n(idx_O2)  &
        -k(159)*n(idx_Ck)  &
        -k(160)*n(idx_Ok)  &
        -k(193)*n(idx_H2)  &
        +k(193)*n(idx_H2)  &
        -k(194)*n(idx_H2)

    !d[HE+_dot]/d[H+]
    pd(21,20) =  &
        +k(9)*n(idx_HE)  &
        +k(10)*n(idx_HE)

    !d[H2+_dot]/d[H+]
    pd(22,20) =  &
        +k(18)*n(idx_H)  &
        +k(19)*n(idx_H)  &
        +k(21)*n(idx_H2)  &
        +k(22)*n(idx_H2)  &
        +k(29)*n(idx_Hk)

    !d[C+_dot]/d[H+]
    pd(23,20) =  &
        +k(47)*n(idx_C)

    !d[O+_dot]/d[H+]
    pd(24,20) =  &
        +k(45)*n(idx_O)

    !d[H3+_dot]/d[H+]
    pd(27,20) =  &
        +k(194)*n(idx_H2)

    !d[CH+_dot]/d[H+]
    pd(28,20) =  &
        +k(130)*n(idx_CH)  &
        +k(131)*n(idx_CH)  &
        +k(132)*n(idx_CH2)  &
        +k(133)*n(idx_CH2)

    !d[CH2+_dot]/d[H+]
    pd(29,20) =  &
        +k(134)*n(idx_CH2)  &
        +k(135)*n(idx_CH2)

    !d[OH+_dot]/d[H+]
    pd(32,20) =  &
        +k(141)*n(idx_OH)  &
        +k(142)*n(idx_OH)

    !d[H2O+_dot]/d[H+]
    pd(33,20) =  &
        +k(145)*n(idx_H2O)  &
        +k(146)*n(idx_H2O)

    !d[O2+_dot]/d[H+]
    pd(35,20) =  &
        +k(153)*n(idx_O2)

    !d[Tgas_dot]/d[H+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(20)*1d-3
    if(dnn>0.d0) then
      nn(20) = n(20) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,20) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE+]
    pd(1,21) =  &
        -k(5)*n(idx_E)  &
        -k(6)*n(idx_E)  &
        -k(7)*n(idx_E)  &
        +2.d0*k(7)*n(idx_E)

    !d[H-_dot]/d[HE+]
    pd(2,21) =  &
        -k(161)*n(idx_Hk)

    !d[H_dot]/d[HE+]
    pd(5,21) =  &
        -k(8)*n(idx_H)  &
        +k(13)*n(idx_H2)  &
        +2.d0*k(14)*n(idx_H2)  &
        +k(138)*n(idx_CH2)  &
        +k(139)*n(idx_CH2)  &
        +k(143)*n(idx_OH)  &
        +k(144)*n(idx_OH)  &
        +k(149)*n(idx_H2O)  &
        +k(150)*n(idx_H2O)  &
        +k(161)*n(idx_Hk)

    !d[HE_dot]/d[HE+]
    pd(6,21) =  &
        +k(5)*n(idx_E)  &
        +k(6)*n(idx_E)  &
        +k(8)*n(idx_H)  &
        +k(12)*n(idx_H2)  &
        +k(13)*n(idx_H2)  &
        +k(46)*n(idx_O)  &
        +k(49)*n(idx_C)  &
        +k(50)*n(idx_C)  &
        +k(51)*n(idx_C)  &
        +k(136)*n(idx_CH2)  &
        +k(137)*n(idx_CH2)  &
        +k(138)*n(idx_CH2)  &
        +k(139)*n(idx_CH2)  &
        +k(140)*n(idx_C2)  &
        +k(143)*n(idx_OH)  &
        +k(144)*n(idx_OH)  &
        +k(147)*n(idx_H2O)  &
        +k(148)*n(idx_H2O)  &
        +k(149)*n(idx_H2O)  &
        +k(150)*n(idx_H2O)  &
        +k(151)*n(idx_H2O)  &
        +k(152)*n(idx_H2O)  &
        +k(154)*n(idx_O2)  &
        +k(155)*n(idx_O2)  &
        +k(156)*n(idx_CO)  &
        +k(157)*n(idx_CO)  &
        +k(161)*n(idx_Hk)

    !d[H2_dot]/d[HE+]
    pd(7,21) =  &
        -k(12)*n(idx_H2)  &
        -k(13)*n(idx_H2)  &
        -k(14)*n(idx_H2)  &
        +k(136)*n(idx_CH2)  &
        +k(137)*n(idx_CH2)

    !d[C_dot]/d[HE+]
    pd(8,21) =  &
        -k(49)*n(idx_C)  &
        -k(50)*n(idx_C)  &
        -k(51)*n(idx_C)  &
        +k(140)*n(idx_C2)  &
        +k(157)*n(idx_CO)

    !d[O_dot]/d[HE+]
    pd(9,21) =  &
        -k(46)*n(idx_O)  &
        +k(155)*n(idx_O2)  &
        +k(156)*n(idx_CO)

    !d[OH_dot]/d[HE+]
    pd(10,21) =  &
        -k(143)*n(idx_OH)  &
        -k(144)*n(idx_OH)  &
        +k(147)*n(idx_H2O)  &
        +k(148)*n(idx_H2O)

    !d[CO_dot]/d[HE+]
    pd(11,21) =  &
        -k(156)*n(idx_CO)  &
        -k(157)*n(idx_CO)

    !d[CH2_dot]/d[HE+]
    pd(13,21) =  &
        -k(136)*n(idx_CH2)  &
        -k(137)*n(idx_CH2)  &
        -k(138)*n(idx_CH2)  &
        -k(139)*n(idx_CH2)

    !d[C2_dot]/d[HE+]
    pd(14,21) =  &
        -k(140)*n(idx_C2)

    !d[H2O_dot]/d[HE+]
    pd(16,21) =  &
        -k(147)*n(idx_H2O)  &
        -k(148)*n(idx_H2O)  &
        -k(149)*n(idx_H2O)  &
        -k(150)*n(idx_H2O)  &
        -k(151)*n(idx_H2O)  &
        -k(152)*n(idx_H2O)

    !d[O2_dot]/d[HE+]
    pd(17,21) =  &
        -k(154)*n(idx_O2)  &
        -k(155)*n(idx_O2)

    !d[H+_dot]/d[HE+]
    pd(20,21) =  &
        +k(8)*n(idx_H)  &
        +k(13)*n(idx_H2)  &
        +k(147)*n(idx_H2O)  &
        +k(148)*n(idx_H2O)

    !d[HE+_dot]/d[HE+]
    pd(21,21) =  &
        -k(5)*n(idx_E)  &
        -k(6)*n(idx_E)  &
        -k(7)*n(idx_E)  &
        -k(8)*n(idx_H)  &
        -k(12)*n(idx_H2)  &
        -k(13)*n(idx_H2)  &
        -k(14)*n(idx_H2)  &
        +k(14)*n(idx_H2)  &
        -k(46)*n(idx_O)  &
        -k(49)*n(idx_C)  &
        -k(50)*n(idx_C)  &
        -k(51)*n(idx_C)  &
        -k(136)*n(idx_CH2)  &
        -k(137)*n(idx_CH2)  &
        -k(138)*n(idx_CH2)  &
        -k(139)*n(idx_CH2)  &
        -k(140)*n(idx_C2)  &
        -k(143)*n(idx_OH)  &
        -k(144)*n(idx_OH)  &
        -k(147)*n(idx_H2O)  &
        -k(148)*n(idx_H2O)  &
        -k(149)*n(idx_H2O)  &
        -k(150)*n(idx_H2O)  &
        -k(151)*n(idx_H2O)  &
        -k(152)*n(idx_H2O)  &
        -k(154)*n(idx_O2)  &
        -k(155)*n(idx_O2)  &
        -k(156)*n(idx_CO)  &
        -k(157)*n(idx_CO)  &
        -k(161)*n(idx_Hk)

    !d[H2+_dot]/d[HE+]
    pd(22,21) =  &
        +k(12)*n(idx_H2)

    !d[C+_dot]/d[HE+]
    pd(23,21) =  &
        +k(49)*n(idx_C)  &
        +k(50)*n(idx_C)  &
        +k(51)*n(idx_C)  &
        +k(136)*n(idx_CH2)  &
        +k(137)*n(idx_CH2)  &
        +k(140)*n(idx_C2)  &
        +k(156)*n(idx_CO)

    !d[O+_dot]/d[HE+]
    pd(24,21) =  &
        +k(46)*n(idx_O)  &
        +k(143)*n(idx_OH)  &
        +k(144)*n(idx_OH)  &
        +k(155)*n(idx_O2)  &
        +k(157)*n(idx_CO)

    !d[CH+_dot]/d[HE+]
    pd(28,21) =  &
        +k(138)*n(idx_CH2)  &
        +k(139)*n(idx_CH2)

    !d[OH+_dot]/d[HE+]
    pd(32,21) =  &
        +k(149)*n(idx_H2O)  &
        +k(150)*n(idx_H2O)

    !d[H2O+_dot]/d[HE+]
    pd(33,21) =  &
        +k(151)*n(idx_H2O)  &
        +k(152)*n(idx_H2O)

    !d[O2+_dot]/d[HE+]
    pd(35,21) =  &
        +k(154)*n(idx_O2)

    !d[HE++_dot]/d[HE+]
    pd(36,21) =  &
        +k(7)*n(idx_E)

    !d[Tgas_dot]/d[HE+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(21)*1d-3
    if(dnn>0.d0) then
      nn(21) = n(21) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,21) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2+]
    pd(1,22) =  &
        -k(30)*n(idx_E)  &
        -k(31)*n(idx_E)

    !d[H-_dot]/d[H2+]
    pd(2,22) =  &
        -k(32)*n(idx_Hk)

    !d[H_dot]/d[H2+]
    pd(5,22) =  &
        -k(20)*n(idx_H)  &
        +2.d0*k(30)*n(idx_E)  &
        +2.d0*k(31)*n(idx_E)  &
        +k(32)*n(idx_Hk)  &
        +k(85)*n(idx_H2)  &
        +k(87)*n(idx_C)  &
        +k(102)*n(idx_O)  &
        +k(214)

    !d[H2_dot]/d[H2+]
    pd(7,22) =  &
        +k(20)*n(idx_H)  &
        +k(32)*n(idx_Hk)  &
        -k(85)*n(idx_H2)

    !d[C_dot]/d[H2+]
    pd(8,22) =  &
        -k(87)*n(idx_C)

    !d[O_dot]/d[H2+]
    pd(9,22) =  &
        -k(102)*n(idx_O)

    !d[H+_dot]/d[H2+]
    pd(20,22) =  &
        +k(20)*n(idx_H)  &
        +k(214)

    !d[H2+_dot]/d[H2+]
    pd(22,22) =  &
        -k(20)*n(idx_H)  &
        -k(30)*n(idx_E)  &
        -k(31)*n(idx_E)  &
        -k(32)*n(idx_Hk)  &
        -k(85)*n(idx_H2)  &
        -k(87)*n(idx_C)  &
        -k(102)*n(idx_O)  &
        -k(214)

    !d[H3+_dot]/d[H2+]
    pd(27,22) =  &
        +k(85)*n(idx_H2)

    !d[CH+_dot]/d[H2+]
    pd(28,22) =  &
        +k(87)*n(idx_C)

    !d[OH+_dot]/d[H2+]
    pd(32,22) =  &
        +k(102)*n(idx_O)

    !d[Tgas_dot]/d[H2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(22)*1d-3
    if(dnn>0.d0) then
      nn(22) = n(22) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,22) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C+]
    pd(1,23) =  &
        -k(37)*n(idx_E)  &
        -k(38)*n(idx_E)  &
        -k(39)*n(idx_E)

    !d[H_dot]/d[C+]
    pd(5,23) =  &
        -k(48)*n(idx_H)  &
        +k(90)*n(idx_H2)  &
        +k(107)*n(idx_OH)  &
        +k(108)*n(idx_OH)  &
        +k(113)*n(idx_H2O)  &
        +k(114)*n(idx_H2O)  &
        +k(115)*n(idx_H2O)  &
        -k(200)*n(idx_H)

    !d[H2_dot]/d[C+]
    pd(7,23) =  &
        -k(90)*n(idx_H2)  &
        -k(201)*n(idx_H2)

    !d[C_dot]/d[C+]
    pd(8,23) =  &
        +k(37)*n(idx_E)  &
        +k(38)*n(idx_E)  &
        +k(39)*n(idx_E)  &
        +k(48)*n(idx_H)  &
        +k(116)*n(idx_H2O)

    !d[O_dot]/d[C+]
    pd(9,23) =  &
        +k(118)*n(idx_O2)  &
        -k(202)*n(idx_O)  &
        -k(203)*n(idx_O)  &
        -k(269)*n(idx_O)  &
        -k(270)*n(idx_O)

    !d[OH_dot]/d[C+]
    pd(10,23) =  &
        -k(107)*n(idx_OH)  &
        -k(108)*n(idx_OH)

    !d[CO_dot]/d[C+]
    pd(11,23) =  &
        +k(119)*n(idx_O2)

    !d[H2O_dot]/d[C+]
    pd(16,23) =  &
        -k(113)*n(idx_H2O)  &
        -k(114)*n(idx_H2O)  &
        -k(115)*n(idx_H2O)  &
        -k(116)*n(idx_H2O)

    !d[O2_dot]/d[C+]
    pd(17,23) =  &
        -k(118)*n(idx_O2)  &
        -k(119)*n(idx_O2)

    !d[H+_dot]/d[C+]
    pd(20,23) =  &
        +k(48)*n(idx_H)

    !d[C+_dot]/d[C+]
    pd(23,23) =  &
        -k(37)*n(idx_E)  &
        -k(38)*n(idx_E)  &
        -k(39)*n(idx_E)  &
        -k(48)*n(idx_H)  &
        -k(90)*n(idx_H2)  &
        -k(107)*n(idx_OH)  &
        -k(108)*n(idx_OH)  &
        -k(113)*n(idx_H2O)  &
        -k(114)*n(idx_H2O)  &
        -k(115)*n(idx_H2O)  &
        -k(116)*n(idx_H2O)  &
        -k(118)*n(idx_O2)  &
        -k(119)*n(idx_O2)  &
        -k(200)*n(idx_H)  &
        -k(201)*n(idx_H2)  &
        -k(202)*n(idx_O)  &
        -k(203)*n(idx_O)  &
        -k(269)*n(idx_O)  &
        -k(270)*n(idx_O)

    !d[O+_dot]/d[C+]
    pd(24,23) =  &
        +k(119)*n(idx_O2)

    !d[HOC+_dot]/d[C+]
    pd(25,23) =  &
        +k(113)*n(idx_H2O)

    !d[HCO+_dot]/d[C+]
    pd(26,23) =  &
        +k(114)*n(idx_H2O)  &
        +k(115)*n(idx_H2O)

    !d[CH+_dot]/d[C+]
    pd(28,23) =  &
        +k(90)*n(idx_H2)  &
        +k(200)*n(idx_H)

    !d[CH2+_dot]/d[C+]
    pd(29,23) =  &
        +k(201)*n(idx_H2)

    !d[CO+_dot]/d[C+]
    pd(30,23) =  &
        +k(107)*n(idx_OH)  &
        +k(108)*n(idx_OH)  &
        +k(118)*n(idx_O2)  &
        +k(202)*n(idx_O)  &
        +k(203)*n(idx_O)  &
        +k(269)*n(idx_O)  &
        +k(270)*n(idx_O)

    !d[H2O+_dot]/d[C+]
    pd(33,23) =  &
        +k(116)*n(idx_H2O)

    !d[Tgas_dot]/d[C+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(23)*1d-3
    if(dnn>0.d0) then
      nn(23) = n(23) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,23) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O+]
    pd(1,24) =  &
        -k(40)*n(idx_E)  &
        -k(41)*n(idx_E)

    !d[H_dot]/d[O+]
    pd(5,24) =  &
        -k(44)*n(idx_H)  &
        +k(101)*n(idx_H2)

    !d[H2_dot]/d[O+]
    pd(7,24) =  &
        -k(101)*n(idx_H2)

    !d[C_dot]/d[O+]
    pd(8,24) =  &
        +k(100)*n(idx_C2)  &
        -k(271)*n(idx_C)  &
        -k(272)*n(idx_C)

    !d[O_dot]/d[O+]
    pd(9,24) =  &
        +k(40)*n(idx_E)  &
        +k(41)*n(idx_E)  &
        +k(44)*n(idx_H)

    !d[C2_dot]/d[O+]
    pd(14,24) =  &
        -k(100)*n(idx_C2)

    !d[H+_dot]/d[O+]
    pd(20,24) =  &
        +k(44)*n(idx_H)

    !d[O+_dot]/d[O+]
    pd(24,24) =  &
        -k(40)*n(idx_E)  &
        -k(41)*n(idx_E)  &
        -k(44)*n(idx_H)  &
        -k(100)*n(idx_C2)  &
        -k(101)*n(idx_H2)  &
        -k(271)*n(idx_C)  &
        -k(272)*n(idx_C)

    !d[CO+_dot]/d[O+]
    pd(30,24) =  &
        +k(100)*n(idx_C2)  &
        +k(271)*n(idx_C)  &
        +k(272)*n(idx_C)

    !d[OH+_dot]/d[O+]
    pd(32,24) =  &
        +k(101)*n(idx_H2)

    !d[Tgas_dot]/d[O+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(24)*1d-3
    if(dnn>0.d0) then
      nn(24) = n(24) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,24) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HOC+]
    pd(1,25) =  &
        -k(183)*n(idx_E)

    !d[H_dot]/d[HOC+]
    pd(5,25) =  &
        +k(183)*n(idx_E)

    !d[H2_dot]/d[HOC+]
    pd(7,25) =  &
        -k(53)*n(idx_H2)  &
        +k(53)*n(idx_H2)

    !d[CO_dot]/d[HOC+]
    pd(11,25) =  &
        -k(54)*n(idx_CO)  &
        +k(54)*n(idx_CO)  &
        -k(55)*n(idx_CO)  &
        +k(55)*n(idx_CO)  &
        +k(183)*n(idx_E)

    !d[HOC+_dot]/d[HOC+]
    pd(25,25) =  &
        -k(53)*n(idx_H2)  &
        -k(54)*n(idx_CO)  &
        -k(55)*n(idx_CO)  &
        -k(183)*n(idx_E)

    !d[HCO+_dot]/d[HOC+]
    pd(26,25) =  &
        +k(53)*n(idx_H2)  &
        +k(54)*n(idx_CO)  &
        +k(55)*n(idx_CO)

    !d[Tgas_dot]/d[HOC+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(25)*1d-3
    if(dnn>0.d0) then
      nn(25) = n(25) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,25) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HCO+]
    pd(1,26) =  &
        -k(181)*n(idx_E)  &
        -k(182)*n(idx_E)

    !d[H_dot]/d[HCO+]
    pd(5,26) =  &
        +k(181)*n(idx_E)

    !d[C_dot]/d[HCO+]
    pd(8,26) =  &
        -k(127)*n(idx_C)  &
        +k(182)*n(idx_E)

    !d[OH_dot]/d[HCO+]
    pd(10,26) =  &
        +k(182)*n(idx_E)

    !d[CO_dot]/d[HCO+]
    pd(11,26) =  &
        +k(127)*n(idx_C)  &
        +k(128)*n(idx_H2O)  &
        +k(129)*n(idx_H2O)  &
        +k(181)*n(idx_E)

    !d[H2O_dot]/d[HCO+]
    pd(16,26) =  &
        -k(128)*n(idx_H2O)  &
        -k(129)*n(idx_H2O)

    !d[HCO+_dot]/d[HCO+]
    pd(26,26) =  &
        -k(127)*n(idx_C)  &
        -k(128)*n(idx_H2O)  &
        -k(129)*n(idx_H2O)  &
        -k(181)*n(idx_E)  &
        -k(182)*n(idx_E)

    !d[CH+_dot]/d[HCO+]
    pd(28,26) =  &
        +k(127)*n(idx_C)

    !d[H3O+_dot]/d[HCO+]
    pd(34,26) =  &
        +k(128)*n(idx_H2O)  &
        +k(129)*n(idx_H2O)

    !d[Tgas_dot]/d[HCO+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(26)*1d-3
    if(dnn>0.d0) then
      nn(26) = n(26) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,26) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H3+]
    pd(1,27) =  &
        -k(162)*n(idx_E)  &
        -k(163)*n(idx_E)

    !d[H_dot]/d[H3+]
    pd(5,27) =  &
        -k(86)*n(idx_H)  &
        +k(89)*n(idx_C)  &
        +k(104)*n(idx_O)  &
        +k(162)*n(idx_E)  &
        +3.d0*k(163)*n(idx_E)  &
        +k(216)

    !d[H2_dot]/d[H3+]
    pd(7,27) =  &
        +k(86)*n(idx_H)  &
        +k(88)*n(idx_C)  &
        +k(103)*n(idx_O)  &
        +k(105)*n(idx_OH)  &
        +k(106)*n(idx_OH)  &
        +k(111)*n(idx_H2O)  &
        +k(112)*n(idx_H2O)  &
        +k(123)*n(idx_CO)  &
        +k(124)*n(idx_CO)  &
        +k(125)*n(idx_CO)  &
        +k(126)*n(idx_CO)  &
        +k(162)*n(idx_E)  &
        +k(215)

    !d[C_dot]/d[H3+]
    pd(8,27) =  &
        -k(88)*n(idx_C)  &
        -k(89)*n(idx_C)

    !d[O_dot]/d[H3+]
    pd(9,27) =  &
        -k(103)*n(idx_O)  &
        -k(104)*n(idx_O)

    !d[OH_dot]/d[H3+]
    pd(10,27) =  &
        -k(105)*n(idx_OH)  &
        -k(106)*n(idx_OH)

    !d[CO_dot]/d[H3+]
    pd(11,27) =  &
        -k(123)*n(idx_CO)  &
        -k(124)*n(idx_CO)  &
        -k(125)*n(idx_CO)  &
        -k(126)*n(idx_CO)

    !d[H2O_dot]/d[H3+]
    pd(16,27) =  &
        -k(111)*n(idx_H2O)  &
        -k(112)*n(idx_H2O)

    !d[H+_dot]/d[H3+]
    pd(20,27) =  &
        +k(215)

    !d[H2+_dot]/d[H3+]
    pd(22,27) =  &
        +k(86)*n(idx_H)  &
        +k(216)

    !d[HOC+_dot]/d[H3+]
    pd(25,27) =  &
        +k(125)*n(idx_CO)  &
        +k(126)*n(idx_CO)

    !d[HCO+_dot]/d[H3+]
    pd(26,27) =  &
        +k(123)*n(idx_CO)  &
        +k(124)*n(idx_CO)

    !d[H3+_dot]/d[H3+]
    pd(27,27) =  &
        -k(86)*n(idx_H)  &
        -k(88)*n(idx_C)  &
        -k(89)*n(idx_C)  &
        -k(103)*n(idx_O)  &
        -k(104)*n(idx_O)  &
        -k(105)*n(idx_OH)  &
        -k(106)*n(idx_OH)  &
        -k(111)*n(idx_H2O)  &
        -k(112)*n(idx_H2O)  &
        -k(123)*n(idx_CO)  &
        -k(124)*n(idx_CO)  &
        -k(125)*n(idx_CO)  &
        -k(126)*n(idx_CO)  &
        -k(162)*n(idx_E)  &
        -k(163)*n(idx_E)  &
        -k(215)  &
        -k(216)

    !d[CH+_dot]/d[H3+]
    pd(28,27) =  &
        +k(88)*n(idx_C)

    !d[CH2+_dot]/d[H3+]
    pd(29,27) =  &
        +k(89)*n(idx_C)

    !d[OH+_dot]/d[H3+]
    pd(32,27) =  &
        +k(103)*n(idx_O)

    !d[H2O+_dot]/d[H3+]
    pd(33,27) =  &
        +k(104)*n(idx_O)  &
        +k(105)*n(idx_OH)  &
        +k(106)*n(idx_OH)

    !d[H3O+_dot]/d[H3+]
    pd(34,27) =  &
        +k(111)*n(idx_H2O)  &
        +k(112)*n(idx_H2O)

    !d[Tgas_dot]/d[H3+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(27)*1d-3
    if(dnn>0.d0) then
      nn(27) = n(27) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,27) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH+]
    pd(1,28) =  &
        -k(164)*n(idx_E)

    !d[H_dot]/d[CH+]
    pd(5,28) =  &
        -k(91)*n(idx_H)  &
        +k(92)*n(idx_H2)  &
        +k(93)*n(idx_O)  &
        +k(164)*n(idx_E)

    !d[H2_dot]/d[CH+]
    pd(7,28) =  &
        +k(91)*n(idx_H)  &
        -k(92)*n(idx_H2)

    !d[C_dot]/d[CH+]
    pd(8,28) =  &
        +k(164)*n(idx_E)  &
        +k(221)

    !d[O_dot]/d[CH+]
    pd(9,28) =  &
        -k(93)*n(idx_O)

    !d[H+_dot]/d[CH+]
    pd(20,28) =  &
        +k(221)

    !d[C+_dot]/d[CH+]
    pd(23,28) =  &
        +k(91)*n(idx_H)

    !d[CH+_dot]/d[CH+]
    pd(28,28) =  &
        -k(91)*n(idx_H)  &
        -k(92)*n(idx_H2)  &
        -k(93)*n(idx_O)  &
        -k(164)*n(idx_E)  &
        -k(221)

    !d[CH2+_dot]/d[CH+]
    pd(29,28) =  &
        +k(92)*n(idx_H2)

    !d[CO+_dot]/d[CH+]
    pd(30,28) =  &
        +k(93)*n(idx_O)

    !d[Tgas_dot]/d[CH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(28)*1d-3
    if(dnn>0.d0) then
      nn(28) = n(28) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,28) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH2+]
    pd(1,29) =  &
        -k(165)*n(idx_E)  &
        -k(166)*n(idx_E)  &
        -k(167)*n(idx_E)

    !d[H_dot]/d[CH2+]
    pd(5,29) =  &
        -k(94)*n(idx_H)  &
        +k(95)*n(idx_H2)  &
        +k(96)*n(idx_O)  &
        +k(165)*n(idx_E)  &
        +2.d0*k(167)*n(idx_E)  &
        +k(224)

    !d[H2_dot]/d[CH2+]
    pd(7,29) =  &
        +k(94)*n(idx_H)  &
        -k(95)*n(idx_H2)  &
        +k(166)*n(idx_E)

    !d[C_dot]/d[CH2+]
    pd(8,29) =  &
        +k(166)*n(idx_E)  &
        +k(167)*n(idx_E)

    !d[O_dot]/d[CH2+]
    pd(9,29) =  &
        -k(96)*n(idx_O)

    !d[OH_dot]/d[CH2+]
    pd(10,29) =  &
        +k(120)*n(idx_O2)

    !d[CH_dot]/d[CH2+]
    pd(12,29) =  &
        +k(165)*n(idx_E)

    !d[O2_dot]/d[CH2+]
    pd(17,29) =  &
        -k(120)*n(idx_O2)

    !d[HCO+_dot]/d[CH2+]
    pd(26,29) =  &
        +k(96)*n(idx_O)  &
        +k(120)*n(idx_O2)

    !d[CH+_dot]/d[CH2+]
    pd(28,29) =  &
        +k(94)*n(idx_H)  &
        +k(224)

    !d[CH2+_dot]/d[CH2+]
    pd(29,29) =  &
        -k(94)*n(idx_H)  &
        -k(95)*n(idx_H2)  &
        -k(96)*n(idx_O)  &
        -k(120)*n(idx_O2)  &
        -k(165)*n(idx_E)  &
        -k(166)*n(idx_E)  &
        -k(167)*n(idx_E)  &
        -k(224)

    !d[CH3+_dot]/d[CH2+]
    pd(31,29) =  &
        +k(95)*n(idx_H2)

    !d[Tgas_dot]/d[CH2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(29)*1d-3
    if(dnn>0.d0) then
      nn(29) = n(29) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,29) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CO+]
    pd(1,30) =  &
        -k(180)*n(idx_E)

    !d[H_dot]/d[CO+]
    pd(5,30) =  &
        -k(158)*n(idx_H)

    !d[C_dot]/d[CO+]
    pd(8,30) =  &
        +k(180)*n(idx_E)

    !d[O_dot]/d[CO+]
    pd(9,30) =  &
        +k(180)*n(idx_E)

    !d[CO_dot]/d[CO+]
    pd(11,30) =  &
        +k(158)*n(idx_H)

    !d[H+_dot]/d[CO+]
    pd(20,30) =  &
        +k(158)*n(idx_H)

    !d[CO+_dot]/d[CO+]
    pd(30,30) =  &
        -k(158)*n(idx_H)  &
        -k(180)*n(idx_E)

    !d[Tgas_dot]/d[CO+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(30)*1d-3
    if(dnn>0.d0) then
      nn(30) = n(30) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,30) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH3+]
    pd(1,31) =  &
        -k(168)*n(idx_E)  &
        -k(169)*n(idx_E)  &
        -k(170)*n(idx_E)

    !d[H_dot]/d[CH3+]
    pd(5,31) =  &
        -k(97)*n(idx_H)  &
        +k(168)*n(idx_E)  &
        +2.d0*k(170)*n(idx_E)  &
        +k(225)

    !d[H2_dot]/d[CH3+]
    pd(7,31) =  &
        +k(97)*n(idx_H)  &
        +k(98)*n(idx_O)  &
        +k(99)*n(idx_O)  &
        +k(169)*n(idx_E)  &
        +k(226)

    !d[O_dot]/d[CH3+]
    pd(9,31) =  &
        -k(98)*n(idx_O)  &
        -k(99)*n(idx_O)

    !d[CH_dot]/d[CH3+]
    pd(12,31) =  &
        +k(169)*n(idx_E)  &
        +k(170)*n(idx_E)

    !d[CH2_dot]/d[CH3+]
    pd(13,31) =  &
        +k(168)*n(idx_E)

    !d[HOC+_dot]/d[CH3+]
    pd(25,31) =  &
        +k(98)*n(idx_O)

    !d[HCO+_dot]/d[CH3+]
    pd(26,31) =  &
        +k(99)*n(idx_O)

    !d[CH+_dot]/d[CH3+]
    pd(28,31) =  &
        +k(226)

    !d[CH2+_dot]/d[CH3+]
    pd(29,31) =  &
        +k(97)*n(idx_H)  &
        +k(225)

    !d[CH3+_dot]/d[CH3+]
    pd(31,31) =  &
        -k(97)*n(idx_H)  &
        -k(98)*n(idx_O)  &
        -k(99)*n(idx_O)  &
        -k(168)*n(idx_E)  &
        -k(169)*n(idx_E)  &
        -k(170)*n(idx_E)  &
        -k(225)  &
        -k(226)

    !d[Tgas_dot]/d[CH3+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(31)*1d-3
    if(dnn>0.d0) then
      nn(31) = n(31) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,31) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[OH+]
    pd(1,32) =  &
        -k(171)*n(idx_E)

    !d[H_dot]/d[OH+]
    pd(5,32) =  &
        +k(109)*n(idx_H2)  &
        +k(171)*n(idx_E)

    !d[H2_dot]/d[OH+]
    pd(7,32) =  &
        -k(109)*n(idx_H2)

    !d[O_dot]/d[OH+]
    pd(9,32) =  &
        +k(171)*n(idx_E)  &
        +k(231)

    !d[H+_dot]/d[OH+]
    pd(20,32) =  &
        +k(231)

    !d[OH+_dot]/d[OH+]
    pd(32,32) =  &
        -k(109)*n(idx_H2)  &
        -k(171)*n(idx_E)  &
        -k(231)

    !d[H2O+_dot]/d[OH+]
    pd(33,32) =  &
        +k(109)*n(idx_H2)

    !d[Tgas_dot]/d[OH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(32)*1d-3
    if(dnn>0.d0) then
      nn(32) = n(32) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,32) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2O+]
    pd(1,33) =  &
        -k(172)*n(idx_E)  &
        -k(173)*n(idx_E)  &
        -k(174)*n(idx_E)

    !d[H_dot]/d[H2O+]
    pd(5,33) =  &
        +k(110)*n(idx_H2)  &
        +k(173)*n(idx_E)  &
        +2.d0*k(174)*n(idx_E)  &
        +k(241)

    !d[H2_dot]/d[H2O+]
    pd(7,33) =  &
        -k(110)*n(idx_H2)  &
        +k(172)*n(idx_E)  &
        +k(240)

    !d[O_dot]/d[H2O+]
    pd(9,33) =  &
        +k(172)*n(idx_E)  &
        +k(174)*n(idx_E)  &
        +k(238)

    !d[OH_dot]/d[H2O+]
    pd(10,33) =  &
        +k(173)*n(idx_E)  &
        +k(239)

    !d[H+_dot]/d[H2O+]
    pd(20,33) =  &
        +k(239)

    !d[H2+_dot]/d[H2O+]
    pd(22,33) =  &
        +k(238)

    !d[O+_dot]/d[H2O+]
    pd(24,33) =  &
        +k(240)

    !d[OH+_dot]/d[H2O+]
    pd(32,33) =  &
        +k(241)

    !d[H2O+_dot]/d[H2O+]
    pd(33,33) =  &
        -k(110)*n(idx_H2)  &
        -k(172)*n(idx_E)  &
        -k(173)*n(idx_E)  &
        -k(174)*n(idx_E)  &
        -k(238)  &
        -k(239)  &
        -k(240)  &
        -k(241)

    !d[H3O+_dot]/d[H2O+]
    pd(34,33) =  &
        +k(110)*n(idx_H2)

    !d[Tgas_dot]/d[H2O+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(33)*1d-3
    if(dnn>0.d0) then
      nn(33) = n(33) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,33) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H3O+]
    pd(1,34) =  &
        -k(175)*n(idx_E)  &
        -k(176)*n(idx_E)  &
        -k(177)*n(idx_E)  &
        -k(178)*n(idx_E)

    !d[H_dot]/d[H3O+]
    pd(5,34) =  &
        +2.d0*k(175)*n(idx_E)  &
        +k(176)*n(idx_E)  &
        +k(177)*n(idx_E)  &
        +k(244)

    !d[H2_dot]/d[H3O+]
    pd(7,34) =  &
        +k(117)*n(idx_C)  &
        +k(176)*n(idx_E)  &
        +k(178)*n(idx_E)  &
        +k(245)

    !d[C_dot]/d[H3O+]
    pd(8,34) =  &
        -k(117)*n(idx_C)

    !d[O_dot]/d[H3O+]
    pd(9,34) =  &
        +k(176)*n(idx_E)

    !d[OH_dot]/d[H3O+]
    pd(10,34) =  &
        +k(175)*n(idx_E)  &
        +k(178)*n(idx_E)  &
        +k(243)

    !d[H2O_dot]/d[H3O+]
    pd(16,34) =  &
        +k(177)*n(idx_E)  &
        +k(242)

    !d[H+_dot]/d[H3O+]
    pd(20,34) =  &
        +k(242)

    !d[H2+_dot]/d[H3O+]
    pd(22,34) =  &
        +k(243)

    !d[HCO+_dot]/d[H3O+]
    pd(26,34) =  &
        +k(117)*n(idx_C)

    !d[OH+_dot]/d[H3O+]
    pd(32,34) =  &
        +k(245)

    !d[H2O+_dot]/d[H3O+]
    pd(33,34) =  &
        +k(244)

    !d[H3O+_dot]/d[H3O+]
    pd(34,34) =  &
        -k(117)*n(idx_C)  &
        -k(175)*n(idx_E)  &
        -k(176)*n(idx_E)  &
        -k(177)*n(idx_E)  &
        -k(178)*n(idx_E)  &
        -k(242)  &
        -k(243)  &
        -k(244)  &
        -k(245)

    !d[Tgas_dot]/d[H3O+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(34)*1d-3
    if(dnn>0.d0) then
      nn(34) = n(34) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,34) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O2+]
    pd(1,35) =  &
        -k(179)*n(idx_E)

    !d[C_dot]/d[O2+]
    pd(8,35) =  &
        -k(121)*n(idx_C)  &
        -k(122)*n(idx_C)

    !d[O_dot]/d[O2+]
    pd(9,35) =  &
        +k(121)*n(idx_C)  &
        +2.d0*k(179)*n(idx_E)

    !d[O2_dot]/d[O2+]
    pd(17,35) =  &
        +k(122)*n(idx_C)

    !d[C+_dot]/d[O2+]
    pd(23,35) =  &
        +k(122)*n(idx_C)

    !d[CO+_dot]/d[O2+]
    pd(30,35) =  &
        +k(121)*n(idx_C)

    !d[O2+_dot]/d[O2+]
    pd(35,35) =  &
        -k(121)*n(idx_C)  &
        -k(122)*n(idx_C)  &
        -k(179)*n(idx_E)

    !d[Tgas_dot]/d[O2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(35)*1d-3
    if(dnn>0.d0) then
      nn(35) = n(35) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,35) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE++]
    pd(1,36) =  &
        -k(15)*n(idx_E)

    !d[HE+_dot]/d[HE++]
    pd(21,36) =  &
        +k(15)*n(idx_E)

    !d[HE++_dot]/d[HE++]
    pd(36,36) =  &
        -k(15)*n(idx_E)

    !d[Tgas_dot]/d[HE++]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(36)*1d-3
    if(dnn>0.d0) then
      nn(36) = n(36) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,36) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[CR]
    pd(39,37) = 0.d0

    !d[Tgas_dot]/d[g]
    pd(39,38) = 0.d0

    !d[E_dot]/d[Tgas]
    pd(1,39) = 0.d0

    !d[H-_dot]/d[Tgas]
    pd(2,39) = 0.d0

    !d[C-_dot]/d[Tgas]
    pd(3,39) = 0.d0

    !d[O-_dot]/d[Tgas]
    pd(4,39) = 0.d0

    !d[H_dot]/d[Tgas]
    pd(5,39) = 0.d0

    !d[HE_dot]/d[Tgas]
    pd(6,39) = 0.d0

    !d[H2_dot]/d[Tgas]
    pd(7,39) = 0.d0

    !d[C_dot]/d[Tgas]
    pd(8,39) = 0.d0

    !d[O_dot]/d[Tgas]
    pd(9,39) = 0.d0

    !d[OH_dot]/d[Tgas]
    pd(10,39) = 0.d0

    !d[CO_dot]/d[Tgas]
    pd(11,39) = 0.d0

    !d[CH_dot]/d[Tgas]
    pd(12,39) = 0.d0

    !d[CH2_dot]/d[Tgas]
    pd(13,39) = 0.d0

    !d[C2_dot]/d[Tgas]
    pd(14,39) = 0.d0

    !d[HCO_dot]/d[Tgas]
    pd(15,39) = 0.d0

    !d[H2O_dot]/d[Tgas]
    pd(16,39) = 0.d0

    !d[O2_dot]/d[Tgas]
    pd(17,39) = 0.d0

    !d[CO_total_dot]/d[Tgas]
    pd(18,39) = 0.d0

    !d[H2O_total_dot]/d[Tgas]
    pd(19,39) = 0.d0

    !d[H+_dot]/d[Tgas]
    pd(20,39) = 0.d0

    !d[HE+_dot]/d[Tgas]
    pd(21,39) = 0.d0

    !d[H2+_dot]/d[Tgas]
    pd(22,39) = 0.d0

    !d[C+_dot]/d[Tgas]
    pd(23,39) = 0.d0

    !d[O+_dot]/d[Tgas]
    pd(24,39) = 0.d0

    !d[HOC+_dot]/d[Tgas]
    pd(25,39) = 0.d0

    !d[HCO+_dot]/d[Tgas]
    pd(26,39) = 0.d0

    !d[H3+_dot]/d[Tgas]
    pd(27,39) = 0.d0

    !d[CH+_dot]/d[Tgas]
    pd(28,39) = 0.d0

    !d[CH2+_dot]/d[Tgas]
    pd(29,39) = 0.d0

    !d[CO+_dot]/d[Tgas]
    pd(30,39) = 0.d0

    !d[CH3+_dot]/d[Tgas]
    pd(31,39) = 0.d0

    !d[OH+_dot]/d[Tgas]
    pd(32,39) = 0.d0

    !d[H2O+_dot]/d[Tgas]
    pd(33,39) = 0.d0

    !d[H3O+_dot]/d[Tgas]
    pd(34,39) = 0.d0

    !d[O2+_dot]/d[Tgas]
    pd(35,39) = 0.d0

    !d[HE++_dot]/d[Tgas]
    pd(36,39) = 0.d0

    !d[CR_dot]/d[Tgas]
    pd(37,39) = 0.d0

    !d[g_dot]/d[Tgas]
    pd(38,39) = 0.d0

    !d[Tgas_dot]/d[Tgas]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(39)*1d-3
    if(dnn>0.d0) then
      nn(39) = n(39) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,39) = (dn1-dn0)/dnn
    end if

    !d[dummy_dot]/d[Tgas]
    pd(40,39) = 0.d0

    !d[Tgas_dot]/d[dummy]
    pd(39,40) = 0.d0

  end subroutine jex

end module krome_ode

!############### MODULE ##############
module krome_user
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  integer,parameter::KROME_idx_E = 1	!E
  integer,parameter::KROME_idx_Hk = 2	!H-
  integer,parameter::KROME_idx_Ck = 3	!C-
  integer,parameter::KROME_idx_Ok = 4	!O-
  integer,parameter::KROME_idx_H = 5	!H
  integer,parameter::KROME_idx_HE = 6	!HE
  integer,parameter::KROME_idx_H2 = 7	!H2
  integer,parameter::KROME_idx_C = 8	!C
  integer,parameter::KROME_idx_O = 9	!O
  integer,parameter::KROME_idx_OH = 10	!OH
  integer,parameter::KROME_idx_CO = 11	!CO
  integer,parameter::KROME_idx_CH = 12	!CH
  integer,parameter::KROME_idx_CH2 = 13	!CH2
  integer,parameter::KROME_idx_C2 = 14	!C2
  integer,parameter::KROME_idx_HCO = 15	!HCO
  integer,parameter::KROME_idx_H2O = 16	!H2O
  integer,parameter::KROME_idx_O2 = 17	!O2
  integer,parameter::KROME_idx_CO_total = 18	!CO_total
  integer,parameter::KROME_idx_H2O_total = 19	!H2O_total
  integer,parameter::KROME_idx_Hj = 20	!H+
  integer,parameter::KROME_idx_HEj = 21	!HE+
  integer,parameter::KROME_idx_H2j = 22	!H2+
  integer,parameter::KROME_idx_Cj = 23	!C+
  integer,parameter::KROME_idx_Oj = 24	!O+
  integer,parameter::KROME_idx_HOCj = 25	!HOC+
  integer,parameter::KROME_idx_HCOj = 26	!HCO+
  integer,parameter::KROME_idx_H3j = 27	!H3+
  integer,parameter::KROME_idx_CHj = 28	!CH+
  integer,parameter::KROME_idx_CH2j = 29	!CH2+
  integer,parameter::KROME_idx_COj = 30	!CO+
  integer,parameter::KROME_idx_CH3j = 31	!CH3+
  integer,parameter::KROME_idx_OHj = 32	!OH+
  integer,parameter::KROME_idx_H2Oj = 33	!H2O+
  integer,parameter::KROME_idx_H3Oj = 34	!H3O+
  integer,parameter::KROME_idx_O2j = 35	!O2+
  integer,parameter::KROME_idx_HEjj = 36	!HE++
  integer,parameter::KROME_idx_CR = 37	!CR
  integer,parameter::KROME_idx_g = 38	!g
  integer,parameter::KROME_idx_Tgas = 39	!Tgas
  integer,parameter::KROME_idx_dummy = 40	!dummy

  integer,parameter::krome_idx_cool_h2 = 1
  integer,parameter::krome_idx_cool_h2gp = 2
  integer,parameter::krome_idx_cool_atomic = 3
  integer,parameter::krome_idx_cool_cen = 3
  integer,parameter::krome_idx_cool_hd = 4
  integer,parameter::krome_idx_cool_z = 5
  integer,parameter::krome_idx_cool_metal = 5
  integer,parameter::krome_idx_cool_dh = 6
  integer,parameter::krome_idx_cool_enthalpic = 6
  integer,parameter::krome_idx_cool_dust = 7
  integer,parameter::krome_idx_cool_compton = 8
  integer,parameter::krome_idx_cool_cie = 9
  integer,parameter::krome_idx_cool_continuum = 10
  integer,parameter::krome_idx_cool_cont = 10
  integer,parameter::krome_idx_cool_exp = 11
  integer,parameter::krome_idx_cool_expansion = 11
  integer,parameter::krome_idx_cool_ff = 12
  integer,parameter::krome_idx_cool_bss = 12
  integer,parameter::krome_idx_cool_custom = 13
  integer,parameter::krome_idx_cool_co = 14
  integer,parameter::krome_idx_cool_zcie = 15
  integer,parameter::krome_idx_cool_zcienouv = 16
  integer,parameter::krome_idx_cool_zextend = 17
  integer,parameter::krome_idx_cool_gh = 18
  integer,parameter::krome_idx_cool_oh = 19
  integer,parameter::krome_idx_cool_h2o = 20
  integer,parameter::krome_idx_cool_hcn = 21
  integer,parameter::krome_ncools = 21

  integer,parameter::krome_idx_heat_chem = 1
  integer,parameter::krome_idx_heat_compress = 2
  integer,parameter::krome_idx_heat_compr = 2
  integer,parameter::krome_idx_heat_photo = 3
  integer,parameter::krome_idx_heat_dh = 4
  integer,parameter::krome_idx_heat_enthalpic = 4
  integer,parameter::krome_idx_heat_photoav = 5
  integer,parameter::krome_idx_heat_av = 5
  integer,parameter::krome_idx_heat_cr = 6
  integer,parameter::krome_idx_heat_dust = 7
  integer,parameter::krome_idx_heat_xray = 8
  integer,parameter::krome_idx_heat_visc = 9
  integer,parameter::krome_idx_heat_viscous = 9
  integer,parameter::krome_idx_heat_custom = 10
  integer,parameter::krome_idx_heat_zcie = 11
  integer,parameter::krome_nheats = 11

  integer,parameter::krome_nrea=275
  integer,parameter::krome_nmols=36
  integer,parameter::krome_nspec=40
  integer,parameter::krome_natoms=5
  integer,parameter::krome_ndust=0
  integer,parameter::krome_ndustTypes=0
  integer,parameter::krome_nPhotoBins=0
  integer,parameter::krome_nPhotoRates=0

  real*8,parameter::krome_boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::krome_boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::krome_boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::krome_iboltzmann_eV = 1d0/krome_boltzmann_eV !K / eV
  real*8,parameter::krome_iboltzmann_erg = 1d0/krome_boltzmann_erg !K / erg
  real*8,parameter::krome_planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::krome_planck_J = 6.62606957d-34 !J s
  real*8,parameter::krome_planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::krome_iplanck_eV = 1d0/krome_planck_eV !1 / eV / s
  real*8,parameter::krome_iplanck_J = 1d0/krome_planck_J !1 / J / s
  real*8,parameter::krome_iplanck_erg = 1d0/krome_planck_erg !1 / erg / s
  real*8,parameter::krome_gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::krome_e_mass = 9.10938188d-28 !g
  real*8,parameter::krome_p_mass = 1.67262158d-24 !g
  real*8,parameter::krome_n_mass = 1.674920d-24 !g
  real*8,parameter::krome_ip_mass = 1d0/krome_p_mass !1/g
  real*8,parameter::krome_clight = 2.99792458e10 !cm/s
  real*8,parameter::krome_pi = 3.14159265359d0 !#
  real*8,parameter::krome_eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::krome_ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::krome_ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::krome_seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::krome_km_to_cm = 1d5 !km -> cm
  real*8,parameter::krome_cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::krome_kvgas_erg = 8.d0*krome_boltzmann_erg/krome_pi/krome_p_mass !
  real*8,parameter::krome_pre_kvgas_sqrt = sqrt(8.d0*krome_boltzmann_erg/krome_pi) !
  real*8,parameter::krome_pre_planck = 2.d0*krome_planck_erg/krome_clight**2 !erg/cm2*s3
  real*8,parameter::krome_exp_planck = krome_planck_erg / krome_boltzmann_erg !s*K
  real*8,parameter::krome_stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::krome_N_avogadro = 6.0221d23 !#
  real*8,parameter::krome_Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::krome_Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::krome_hubble = 0.704d0 !dimensionless
  real*8,parameter::krome_Omega0 = 1.0d0 !dimensionless
  real*8,parameter::krome_Omegab = 0.0456d0 !dimensionless
  real*8,parameter::krome_Hubble0 = 1.d2*krome_hubble*krome_km_to_cm*krome_cm_to_Mpc !1/s

contains

  !*******************
  subroutine krome_set_user_crate(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_crate = argset
  end subroutine krome_set_user_crate

  !*******************
  function krome_get_user_crate()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_crate
    krome_get_user_crate = user_crate
  end function krome_get_user_crate

  !*******************
  subroutine krome_set_user_Av(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_Av = argset
  end subroutine krome_set_user_Av

  !*******************
  function krome_get_user_Av()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_Av
    krome_get_user_Av = user_Av
  end function krome_get_user_Av

  !************************
  !returns the Tdust averaged over the number density
  ! as computed in the tables
  function krome_get_table_Tdust(x,Tgas)
    use krome_commons
    use krome_grfuncs
    implicit none
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_table_Tdust
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_table_Tdust = get_table_Tdust(n(:))

  end function krome_get_table_Tdust

  !**********************
  !convert from MOCASSIN abundances to KROME
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  ! imap: matrix position index map, integer
  ! returns KROME abundances (cm-3, real*8)
  function krome_convert_xmoc(xmoc,imap) result(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*4,intent(in):: xmoc(:,:)
    real*8::x(nmols),n(nspec)
    integer,intent(in)::imap(:)

    x(:) = 0d0

    x(idx_H) = xmoc(imap(1), 1)
    x(idx_HE) = xmoc(imap(2), 1)
    x(idx_C) = xmoc(imap(6), 1)
    x(idx_O) = xmoc(imap(8), 1)
    x(idx_Hj) = xmoc(imap(1), 2)
    x(idx_HEj) = xmoc(imap(2), 2)
    x(idx_Cj) = xmoc(imap(6), 2)
    x(idx_Oj) = xmoc(imap(8), 2)
    x(idx_HEjj) = xmoc(imap(2), 3)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    x(idx_e) = get_electrons(n(:))

  end function krome_convert_xmoc

  !*************************
  !convert from KROME abundances to MOCASSIN
  ! x: KROME abuances (cm-3, real*8)
  ! imap: matrix position index map, integer
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  subroutine krome_return_xmoc(x,imap,xmoc)
    use krome_commons
    implicit none
    real*8,intent(in)::x(nmols)
    real*4,intent(out)::xmoc(:,:)
    integer,intent(in)::imap(:)

    xmoc(:,:) = 0d0

    xmoc(imap(1), 1) = x(idx_H)
    xmoc(imap(2), 1) = x(idx_HE)
    xmoc(imap(6), 1) = x(idx_C)
    xmoc(imap(8), 1) = x(idx_O)
    xmoc(imap(1), 2) = x(idx_Hj)
    xmoc(imap(2), 2) = x(idx_HEj)
    xmoc(imap(6), 2) = x(idx_Cj)
    xmoc(imap(8), 2) = x(idx_Oj)
    xmoc(imap(2), 3) = x(idx_HEjj)

  end subroutine krome_return_xmoc

  !**********************
  !convert number density (cm-3) into column
  ! density (cm-2) using the specific density
  ! column method (see help for option
  ! -columnDensityMethod)
  ! num is the number density, x(:) is the species
  ! array, Tgas is the gas temperature
  ! If the method is not JEANS, x(:) and Tgas
  ! are dummy variables
  function krome_num2col(num,x,Tgas)
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: x(nmols),krome_num2col
    real*8 :: Tgas,num
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_num2col = num2col(num,n(:))

  end function krome_num2col

  !***********************
  !print on screen the current values of all phys variables
  subroutine krome_print_phys_variables()
    use krome_commons
    implicit none

    print *, "Tcmb:", phys_Tcmb
    print *, "zredshift:", phys_zredshift
    print *, "orthoParaRatio:", phys_orthoParaRatio
    print *, "metallicity:", phys_metallicity
    print *, "Tfloor:", phys_Tfloor

  end subroutine krome_print_phys_variables

  !*******************
  subroutine krome_set_Tcmb(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tcmb = arg
  end subroutine krome_set_Tcmb

  !*******************
  function krome_get_Tcmb()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tcmb
    krome_get_Tcmb = phys_Tcmb
  end function krome_get_Tcmb

  !*******************
  subroutine krome_set_zredshift(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_zredshift = arg
  end subroutine krome_set_zredshift

  !*******************
  function krome_get_zredshift()
    use krome_commons
    implicit none
    real*8 :: krome_get_zredshift
    krome_get_zredshift = phys_zredshift
  end function krome_get_zredshift

  !*******************
  subroutine krome_set_orthoParaRatio(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_orthoParaRatio = arg
  end subroutine krome_set_orthoParaRatio

  !*******************
  function krome_get_orthoParaRatio()
    use krome_commons
    implicit none
    real*8 :: krome_get_orthoParaRatio
    krome_get_orthoParaRatio = phys_orthoParaRatio
  end function krome_get_orthoParaRatio

  !*******************
  subroutine krome_set_metallicity(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_metallicity = arg
  end subroutine krome_set_metallicity

  !*******************
  function krome_get_metallicity()
    use krome_commons
    implicit none
    real*8 :: krome_get_metallicity
    krome_get_metallicity = phys_metallicity
  end function krome_get_metallicity

  !*******************
  subroutine krome_set_Tfloor(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tfloor = arg
  end subroutine krome_set_Tfloor

  !*******************
  function krome_get_Tfloor()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tfloor
    krome_get_Tfloor = phys_Tfloor
  end function krome_get_Tfloor

  !*****************************
  !dump the data for restart (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_store(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="replace")
    !dump temperature
    write(nfile,*) Tgas
    write(nfile,*) dt
    !dump species
    do i=1,nmols
      write(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_store

  !*****************************
  !restore the data from a dump (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_restore(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="old")
    !restore temperature
    read(nfile,*) Tgas
    read(nfile,*) dt
    !restore species
    do i=1,nmols
      read(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_restore

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on()
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off()
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

  !***************************
  !alias for coe in krome_subs
  ! returns the coefficient array of size krome_nrea
  ! for a given Tgas
  function krome_get_coef(Tgas,x)
    use krome_commons
    use krome_subs
    use krome_tabs
    real*8 :: krome_get_coef(nrea),x(nmols)
    real*8,value:: Tgas
    real*8::n(nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    call makeStoreOnceRates(n(:))

    krome_get_coef(:) = coe(n(:))

  end function krome_get_coef

  !****************************
  !get the mean molecular weight from
  ! mass fractions
  function krome_get_mu_x(xin)
    use krome_commons
    implicit none
    real*8 :: xin(nmols), krome_get_mu_x
    real*8::n(nmols)
    n(:) = krome_x2n(xin(:),1d0)
    krome_get_mu_x = krome_get_mu(n(:))
  end function krome_get_mu_x

  !****************************
  !return the adiabatic index from mass fractions
  ! and temperature in K
  function krome_get_gamma_x(xin,inTgas)
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: xin(nmols), krome_get_gamma_x
    real*8::x(nmols),Tgas,rhogas

    Tgas = inTgas
    x(:) = krome_x2n(xin(:),1d0)
    krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

  end function krome_get_gamma_x

  !***************************
  !normalize mass fractions and
  ! set charge to zero
  subroutine krome_consistent_x(x)
    use krome_commons
    use krome_constants
    implicit none
    real*8 :: x(nmols)
    real*8::isumx,sumx,xerr,imass(nmols),ee

    !1. charge consistency
    imass(:) = krome_get_imass()

    x(idx_e) = 0.d0

    ee = sum(krome_get_charges()*x(:)*imass(:))
    ee = max(ee*e_mass,0d0)
    x(idx_e) = ee

    !2. mass fraction consistency
    sumx = sum(x)

    !NOTE: uncomment here if you want some additional control
    !conservation error threshold: rise an error if above xerr
    !xerr = 1d-2
    !if(abs(sum-1d0)>xerr) then
    !   print *,"ERROR: some problem with conservation!"
    !   print *,"|sum(x)-1|=",abs(sum-1d0)
    !   stop
    !end if

    isumx = 1d0/sumx
    x(:) = x(:) * isumx

  end subroutine krome_consistent_x

  !*********************
  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    use krome_commons
    implicit none
    real*8 :: n(nmols),krome_n2x(nmols)
    real*8,value :: rhogas

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    use krome_commons
    implicit none
    real*8 :: x(nmols),krome_x2n(nmols)
    real*8,value :: rhogas

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !******************
  !returns free-fall time using the number density
  ! abundances of array x(:)
  function krome_get_free_fall_time(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time
    real*8::x(:),n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_free_fall_time = get_free_fall_time(n(:))

  end function krome_get_free_fall_time

  !******************
  !returns free-fall time using the total mass density
  !  of gas, rhogas (g/cm3)
  function krome_get_free_fall_time_rho(rhogas)
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time_rho
    real*8::rhogas

    krome_get_free_fall_time_rho = get_free_fall_time_rho(rhogas)

  end function krome_get_free_fall_time_rho

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt)
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    use krome_gadiab
    implicit none
    real*8 :: x(nmols), Tgas
    real*8 :: dt
    real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas

  end subroutine krome_thermo

  !*************************
  !get heating (erg/cm3/s) for a given species
  ! array x(:) and Tgas
  function krome_get_heating(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_heating
    real*8::Tgas,k(nrea),nH2dust,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    k(:) = coe(n(:))
    nH2dust = 0d0
    krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
  end function krome_get_heating

  !*****************************
  ! get an array containing individual heatings (erg/cm3/s)
  ! the array has size krome_nheats. see heatcool.gps
  ! for index list
  function krome_get_heating_array(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k(nrea),nH2dust
    real*8 :: x(nmols),krome_get_heating_array(nheats)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    k(:) = coe(n(:))
    Tgas = inTgas
    nH2dust = 0d0
    krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

  end function krome_get_heating_array

  !*************************
  !get cooling (erg/cm3/s) for x(:) species array
  ! and Tgas
  function krome_get_cooling(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_cooling
    real*8::Tgas,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    krome_get_cooling = cooling(n,Tgas)
  end function krome_get_cooling

  !*****************************
  ! get an array containing individual coolings (erg/cm3/s)
  ! the array has size krome_ncools. see heatcool.gps
  ! for index list
  function krome_get_cooling_array(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas
    real*8 :: x(nmols),krome_get_cooling_array(ncools)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    Tgas = inTgas
    krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

  end function krome_get_cooling_array

  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n)
    use krome_cooling
    implicit none
    real*8 :: n(krome_nmols)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::x(nspec)
    integer, optional :: nfile_in
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
    if(present(nfile_in)) nfile = nfile_in
    call dump_cool(x(:),Tgas,nfile)

  end subroutine krome_dump_cooling

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  subroutine krome_conserveLin_x(x,ref)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),ref(natoms)

    call conserveLin_x(x(:),ref(:))

  end subroutine krome_conserveLin_x

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  function krome_conserveLinGetRef_x(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

    krome_conserveLinGetRef_x(:) = &
        conserveLinGetRef_x(x(:))

  end function krome_conserveLinGetRef_x

  !*************************
  !force conservation to array x(:)
  !using xi(:) as initial abundances.
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
    real*8::n(krome_nspec),ni(krome_nspec)

    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !get the adiabatic index for x(:) species abundances
  ! and Tgas.
  ! alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas)
    use krome_subs
    use krome_commons
    use krome_gadiab
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_gamma
    real*8::n(nspec)
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies.
  ! alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    integer :: krome_get_zatoms(nmols)
    integer::zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)

  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from
  ! number density and mass density.
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols), krome_get_mu
    real*8::n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    krome_get_mu = get_mu(n(:))
  end function krome_get_mu

  !***************************
  !get the names of the reactions as a
  ! character*50 array of krome_nrea
  ! elements
  !! !! cannot yet be called from C
  function krome_get_rnames()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    character*50 :: krome_get_rnames(nrea)

    krome_get_rnames(:) = get_rnames()

  end function krome_get_rnames

  !*****************
  !get an array of double containing the masses in g
  ! of the species.
  ! alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(nspec)
    real*8 :: krome_get_Hnuclei, x(nmols)
    n(:) = 0d0
    n(1:nmols) = x(:)

    krome_get_Hnuclei = get_Hnuclei(n(:))

  end function krome_get_Hnuclei

  !*****************
  !get an array of size krome_nmols containing the
  ! charges of the species.
  ! alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 and size krome_nmols
  ! containing the names of all the species.
  ! alias for get_names
  !!  !! cannot yet be called from C
  function krome_get_names()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    character*16 :: krome_get_names(nmols)
    character*16::tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !********************
  !get space-separated header of chemical species
  function krome_get_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*141::krome_get_names_header
    character*16::tmp(nspec)
    integer::i

    tmp(:) = get_names()

    krome_get_names_header = ""
    do i=1,nmols
      krome_get_names_header = trim(krome_get_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_names_header

  !********************
  !get space-separated header of coolings
  function krome_get_cooling_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*141::krome_get_cooling_names_header
    character*16::tmp(ncools)
    integer::i

    tmp(:) = get_cooling_names()

    krome_get_cooling_names_header = ""
    do i=1,ncools
      if(trim(tmp(i))=="") cycle
      krome_get_cooling_names_header = trim(krome_get_cooling_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_cooling_names_header

  !********************
  !get space-separated header of heatings
  function krome_get_heating_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*87::krome_get_heating_names_header
    character*16::tmp(nheats)
    integer::i

    tmp(:) = get_heating_names()

    krome_get_heating_names_header = ""
    do i=1,nheats
      if(trim(tmp(i))=="") cycle
      krome_get_heating_names_header = trim(krome_get_heating_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_heating_names_header

  !*****************
  !get the index of the species with name name.
  ! alias for get_index
  !! !! cannot yet be called from C
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer :: krome_get_index
    character*(*) :: name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n)
    use krome_commons
    real*8 :: krome_get_rho, n(nmols)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho

  !*************************
  !scale the abundances of the metals contained in n(:)
  ! to Z according to Asplund+2009.
  ! note that this applies only to neutral atoms.
  subroutine krome_scale_Z(x,Z)
    use krome_commons
    use krome_getphys
    real*8 :: x(nmols)
    real*8 :: Z
    real*8::Htot,n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0

    Htot = get_Hnuclei(n(:))
    x(idx_C) = max(Htot * 1d1**(Z+(-3.5700000000000003)), 1d-40)
    x(idx_O) = max(Htot * 1d1**(Z+(-3.3100000000000005)), 1d-40)

  end subroutine krome_scale_Z

  !*************************
  !set the total metallicity
  ! in terms of Z/Z_solar
  subroutine krome_set_Z(xarg)
    use krome_commons
    real*8 :: xarg

    total_Z = xarg

  end subroutine krome_set_Z

  !*************************
  !set D is in terms of D_solar (D/D_sol).
  subroutine krome_set_dust_to_gas(xarg)
    use krome_commons
    real*8 :: xarg

    dust2gas_ratio = xarg

  end subroutine

  !*************************
  !set the clumping factor
  subroutine krome_set_clump(xarg)
    use krome_commons
    real*8 :: xarg

    clump_factor = xarg

  end subroutine krome_set_clump

  !***********************
  !get the number of electrons assuming
  ! total neutral charge (cations-anions)
  function krome_get_electrons(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8 :: x(nmols), krome_get_electrons
    real*8::n(nspec)
    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_electrons = get_electrons(n(:))
  end function krome_get_electrons

  !**********************
  !print on screen the first nbest highest reaction fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::x(nmols),n(nspec)
    integer :: nbest
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)

  end subroutine krome_print_best_flux

  !*********************
  !print only the highest fluxes greater than a fraction frac
  ! of the maximum flux
  subroutine krome_print_best_flux_frac(xin,Tgas,frac)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas,frac
    real*8::n(nspec)
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_frac(n,Tgas,frac)

  end subroutine krome_print_best_flux_frac

  !**********************
  !print the highest nbest fluxes for reactions involving
  !a given species using the index idx_find (e.g. krome_idx_H2)
  subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::n(nspec)
    integer :: nbest,idx_find
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_spec(n,Tgas,nbest,idx_find)
  end subroutine krome_print_best_flux_spec

  !*******************************
  !get an array of size krome_nrea with
  ! the fluxes of all the reactions in cm-3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8 :: krome_get_flux(nrea),n(nmols)
    real*8,value :: Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*****************************
  !store the fluxes to the file unit ifile
  ! using the chemical composition x(:), and the
  ! gas temperature Tgas. xvar is th value of an
  ! user-defined independent variable that
  ! can be employed for plots.
  ! the file columns are as follow
  ! rate number, xvar, absolute flux,
  !  flux/maxflux, flux fraction wrt total,
  !  reaction name (*50 string)
  subroutine krome_explore_flux(x,Tgas,ifile,xvar)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols)
    real*8 :: Tgas,xvar
    real*8::flux(nrea),fluxmax,sumflux,n(nspec)
    integer :: ifile
    integer::i
    character*50::rname(nrea)

    !get reaction names
    rname(:) = get_rnames()
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    !get fluxes
    flux(:) = get_flux(n(:), Tgas)
    fluxmax = maxval(flux) !maximum flux
    sumflux = sum(flux) !sum of all the fluxes
    !loop on reactions
    do i=1,nrea
      write(ifile,'(I8,5E17.8e3,a3,a50)') i,xvar,Tgas,flux(i),&
          flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
    end do
    write(ifile,*)

  end subroutine krome_explore_flux

  !*********************
  !get nulcear qeff for the reactions
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: krome_get_qeff(nrea)

    krome_get_qeff(:) = get_qeff()

  end function krome_get_qeff

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile)
    use krome_commons
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::flux(nrea)
    integer :: nfile
    integer::i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
      write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  !dump all the evaluation of the coefficient rates in
  ! the file funit, in the range inTmin, inTmax, using
  ! imax points
  subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer :: funit,imax
    real*8 :: inTmin,inTmax
    real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

    Tmin = log10(inTmin)
    Tmax = log10(inTmax)

    n(:) = 1d-40
    do i=1,imax
      Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
      n(idx_Tgas) = Tgas
      k(:) = coe(n(:))
      do j=1,nrea
        write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
      end do
      write(funit,*)
    end do

  end subroutine krome_dump_rates

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer::i,charges(nspec)
    real*8 :: x(nmols)
    real*8 :: Tgas
    real*8::masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
      print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info

  !*****************************
  subroutine krome_set_mpi_rank(xarg)
    use krome_commons
    implicit none
    integer :: xarg
    krome_mpi_rank=xarg
  end subroutine krome_set_mpi_rank

  !**************************
  function krome_get_jacobian(j,x,Tgas)
    use krome_ode
    use krome_commons
    implicit none
    integer, value :: j
    real*8,value :: Tgas
    real*8 :: x(nmols),krome_get_jacobian(nspec)
    integer::ian, jan, i
    real*8::tt, n(nspec)
    real*8::pdj(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = tgas

    tt = 0d0
    ian = 0
    jan = 0

    call jes(nspec, tt, n, j, ian, jan, pdj)
    krome_get_jacobian(:) = pdj(:)

  end function krome_get_jacobian

end module krome_user

!############### MODULE ##############
module krome_reduction
contains

  !**************************
  function fex_check(n,Tgas)
    use krome_commons
    use krome_tabs
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe_tab(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
      rrmax = max(rrmax, arr_flux(i))
    end do
    fex_check = rrmax

  end function fex_check

end module krome_reduction

!############### MODULE ##############
module krome_main

  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2026-01-22 10:17:20
  !  Changeset cd85309
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !*******************************
  !KROME main (interface to the solver library)

  subroutine krome(x,Tgas,dt  )
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    use krome_getphys
    use krome_tabs
    implicit none
    real*8 :: Tgas,dt
    real*8 :: x(nmols)
    real*8 :: rhogas

    real*8::mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
    integer::icount,i,icount_max
    integer:: ierr

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
    integer::iwork(641)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(3972)
    logical::got_error,equil

    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1.000000d-04 !relative tolerance
    atol(:) = 1.000000d-10 !absolute tolerance
    icount_max = 100 !maximum number of iterations

    itask = 1
    iopt = 0

    !MF=
    !  = 222 internal-generated JAC and sparsity
    !  = 121 user-provided JAC and internal generated sparsity
    !  =  22 internal-generated JAC but sparsity user-provided
    !  =  21 user-provided JAC and sparsity
    MF = 222
    !end init DLSODES
    !****************************

    ierr = 0 !error flag, zero==OK!
    n(:) = 0d0 !initialize densities

    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas !put temperature in the input array

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    !store initial values
    ni(:) = n(:)
    n_global(:) = n(:)

    call makeStoreOnceRates(n(:))

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
      icount = icount + 1
      !solve ODE
      CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, &
          ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
          LIW, JES, MF)

      krome_call_to_fex = krome_call_to_fex + IWORK(12)
      !check DLSODES exit status
      if(istate==2) then
        exit !sucsessful integration
      elseif(istate==-1) then
        istate = 1 !exceeded internal max iterations
      elseif(istate==-5 .or. istate==-4) then
        istate = 3 !wrong sparsity recompute
      elseif(istate==-3) then
        n(:) = ni(:)
        istate = 1
      else
        got_error = .true.
      end if

      if(got_error.or.icount>icount_max) then
        if (krome_mpi_rank>0) then
          print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
          print *,krome_mpi_rank,"istate:",istate
          print *,krome_mpi_rank,"iter count:",icount
          print *,krome_mpi_rank,"max iter count:",icount_max
          print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
        else
          print *,"ERROR: wrong solver exit status!"
          print *,"istate:",istate
          print *,"iter count:",icount
          print *,"max iter count:",icount_max
          print *,"SEE KROME_ERROR_REPORT file"
        end if
        call krome_dump(n(:), rwork(:), iwork(:), ni(:))
        stop
      end if

    end do

    !avoid negative species
    do i=1,nspec
      n(i) = max(n(i),0d0)
    end do

    !returns to user array
    x(:) = n(1:nmols)

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !*********************************
  !integrates to equilibrium using constant temperature
  subroutine krome_equilibrium(x,Tgas,verbosity)
    use krome_ode
    use krome_subs
    use krome_commons
    use krome_constants
    use krome_getphys
    use krome_tabs
    implicit none
    integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
    integer::i,imax
    integer,optional::verbosity
    integer::verbose
    real*8 :: Tgas
    real*8 :: x(nmols)
    real*8 :: rhogas
    real*8::tloc,n(nspec),mass(nspec),ni(nspec)
    real*8::dt,xin
    integer::iwork(641)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(3972)
    real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
    logical::converged

    integer, save :: ncall=0
    integer, parameter :: ncall_print_frequency=20000
    integer :: ncallp
    integer::charges(nspec)
    real*8::masses(nspec)
    character*16::names(nspec)

    !set verbosity from argument
    verbose = 1 !default is verbose
    if(present(verbosity)) verbose = verbosity

    call XSETF(0)!toggle solver verbosity
    meth = 2
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1d-6 !relative tolerance
    atol(:) = 1d-20 !absolute tolerance

    ! Switches to decide when equilibrium has been reached
    ertol = 1d-5  ! relative min change in a species
    eatol = 1d-12 ! absolute min change in a species
    max_time=seconds_per_year*5d8 ! max time we will be integrating for

    !for DLSODES options see its manual
    iopt = 0
    itask = 1
    istate = 1

    mf = 222 !internally evaluated sparsity and jacobian
    tloc = 0d0 !initial time

    n(:) = 0d0 !initialize densities
    !copy into array
    n(nmols+1:) = 0d0
    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas

    !store previous values
    ni(:) = n(:)
    n_global(:) = ni(:)

    call makeStoreOnceRates(n(:))

    imax = 1000

    dt = seconds_per_year * 1d2
    t_tot = dt
    converged = .false.
    do while (.not. converged)
      do i=1,imax
        !solve ODE
        CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
        if(istate==2) then
          exit
        else
          istate=1
        end if
      end do
      !check errors
      if(istate.ne.2) then
        print *,"ERROR: no equilibrium found!"
        stop
      end if

      !avoid negative species
      do i=1,nspec
        n(i) = max(n(i),0d0)
      end do

      ! check if we have converged by comparing the error in any species with an relative abundance above eatol
      converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
          .or. t_tot .gt. max_time

      ! Increase integration time by a reasonable factor
      if(.not. converged) then
        dt = dt * 3.
        t_tot = t_tot + dt
        ni = n
        n_global = n
      endif
    enddo
    !returns to user array
    x(:) = n(1:nmols)

    if(t_tot > max_time .and. &
        maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
    print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
    print *, 'Tgas :', Tgas
    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
    do i=1,nmols
      print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
    end do
    print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
    print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
    print *, 'Absolute and relative floors:', eatol, ertol
  end if

  ! Print info ever so often
  !$omp critical
  ncall=ncall+1
  ncallp = ncall
  !$omp end critical

  if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
    print *, 'Found equilibrium for ', ncallp, ' cells.'
  end if

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
  implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
  use krome_commons
  use krome_ode
  implicit none
  integer::n,ierr
  real*8::x(n),f(n),tt
  call fex(n,tt,x(:),f(:))
  f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
  use krome_commons
  use krome_subs
  use krome_tabs
  use krome_reduction
  use krome_ode
  use krome_getphys
  integer::fnum,i,iwork(:),idx(nrea),j
  real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
  character*16::names(nspec),FMTi,FMTr
  character*50::rnames(nrea),fname,prex
  integer,save::mx_dump=1000 ! max nr of reports before terminating
  fnum = 99
  if (krome_mpi_rank>0) then
    write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
  else
    fname = "KROME_ERROR_REPORT"
  endif
  open(fnum,FILE=trim(fname),status="replace")
  tt = 0d0
  names(:) = get_names()
  rnames(:) = get_rnames()
  call fex(nspec,tt,n(:),dn(:))

  write(fnum,*) "KROME ERROR REPORT"
  write(fnum,*)
  !SPECIES
  write(fnum,*) "Species abundances"
  write(fnum,*) "**********************"
  write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
  end do
  write(fnum,*) "**********************"

  !F90 FRIENDLY RESTART
  write(fnum,*)
  write(fnum,*) "**********************"
  write(fnum,*) "F90-friendly species"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(prex,'(a,i3,a)') "x(",i,") = "
    write(fnum,*) trim(prex),ni(i),"!"//names(i)
  end do

  write(fnum,*) "**********************"

  !RATE COEFFIECIENTS
  k(:) = coe_tab(n(:))
  idx(:) = idx_sort(k(:))
  kmax = maxval(k)
  write(fnum,*)
  write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    kperc = 0.d0
    if(kmax>0.d0) kperc = k(i)*1d2/kmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !FLUXES
  call load_arrays
  rrmax = fex_check(n(:), n(idx_Tgas))
  idx(:) = idx_sort(arr_flux(:))
  write(fnum,*)
  write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    rperc = 0.d0
    if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !SOLVER
  FMTr = "(a30,E16.7e3)"
  FMTi = "(a30,I10)"
  write(fnum,*) "Solver-related information:"
  write(fnum,FMTr) "step size last",rwork(11)
  write(fnum,FMTr) "step size attempt",rwork(12)
  write(fnum,FMTr) "time current",rwork(13)
  write(fnum,FMTr) "tol scale factor",rwork(14)
  write(fnum,FMTi) "numeber of steps",iwork(11)
  write(fnum,FMTi) "call to fex",iwork(12)
  write(fnum,FMTi) "call to jex",iwork(13)
  write(fnum,FMTi) "last order used",iwork(14)
  write(fnum,FMTi) "order attempt",iwork(15)
  write(fnum,FMTi) "idx largest error",iwork(16)
  write(fnum,FMTi) "RWORK size required",iwork(17)
  write(fnum,FMTi) "IWORK size required",iwork(18)
  write(fnum,FMTi) "NNZ in Jac",iwork(19)
  write(fnum,FMTi) "extra fex to compute jac",iwork(20)
  write(fnum,FMTi) "number of LU decomp",iwork(21)
  write(fnum,FMTi) "base address in RWORK",iwork(22)
  write(fnum,FMTi) "base address of IAN",iwork(23)
  write(fnum,FMTi) "base address of JAN",iwork(24)
  write(fnum,FMTi) "NNZ in lower LU",iwork(25)
  write(fnum,FMTi) "NNZ in upper LU",iwork(21)
  write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
  write(fnum,*)
  write(fnum,*) "END KROME ERROR REPORT"
  write(fnum,*)
  close(fnum)

  mx_dump = mx_dump - 1
  if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
  use krome_commons
  use krome_tabs
  use krome_subs
  use krome_reduction
  use krome_dust
  use krome_cooling
  use krome_photo
  use krome_fit
  use krome_getphys

  !init phys common variables
  !$omp parallel
  phys_Tcmb = 2.73d0
  phys_zredshift = 0d0
  phys_orthoParaRatio = 3d0
  phys_metallicity = 0d0
  phys_Tfloor = 2.73d0
  !$omp end parallel

  !init metallicity default
  !assuming solar
  total_Z = 1d0

  !default D/D_sol = Z/Z_sol
  !assuming linear scaling
  dust2gas_ratio = total_Z

  !default broadening turubulence velocity
  broadeningVturb2 = 0d0

  !default clumping factor for
  ! H2 formation on dust by Jura/Gnedin
  clump_factor = 1d0

  !default for thermo and chem toggle is ON
  !$omp parallel
  krome_thermo_toggle = 1
  krome_chemo_toggle = 1
  !$omp end parallel

  !load arrays with ractants/products indexes
  call load_arrays()

  !initialize the table for exp(-a/T) function
  call init_exp_table()

  call load_parts()

  !init photo reactants indexes

  !get machine precision
  krome_epsilon = epsilon(0d0)

  !load verbatim reactions
  call loadReactionsVerbatim()

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
  !krome_get_coe: public interface to obtain rate coefficients
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coe(nrea), x(nmols), Tgas
  real*8::n(nspec)

  n(:) = 0d0
  n(1:nmols) = x(:)
  n(idx_Tgas) = Tgas
  krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
  !krome_get_coeT: public interface to obtain rate coefficients
  ! with argument Tgas only
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coeT(nrea),Tgas
  real*8::n(nspec)
  n(idx_Tgas) = Tgas
  krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
