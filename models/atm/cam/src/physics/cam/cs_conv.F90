module cs_conv
!---------------------------------------------------------------------------------
! Purpose:
!
! Interface for Chikira-Sugiyama convection scheme 
!
! Author: Minoru Chikira
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use abortutils,    only: endrun
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog
  
!  to be delt with later
!  use physconst,     only: GRAV  => gravit, CP   => cpair, EL   => latvap, &
!                           EMELT => latice, rair         , RVAP => rh20

  implicit none

  include 'mpif.h'

  private                ! Make default type private to the module
!
! Private data
!
  integer , parameter, public :: nctp  = 14            ! number of cloud types
  real(r8), parameter         :: unset_r8 = -999._r8   ! missing value
!
! Physical constants (temporarily following MIROC)
!
  real(r8), parameter :: &
    GRAV  = 9.8_r8   ,   & ! gravity
    CP    = 1004.6_r8,   & ! specific heat of air
    EL    = 2.5e6_r8,    & ! latent heat of condensation
    EMELT = 3.4e5_r8,    & ! latent heat of fusion
    RAIR  = 287.04_r8,   & ! gas constant of air
    RVAP  = 461._r8,     & ! gas constant of vapor
    TMELT = 273.15_r8,   & ! melting point of water
    ES0   = 611._r8,     & ! saturation e at 0 deg C (Pa)
    TQICE = 273.15_r8      ! T threshold for ice QSAT
!
  real(r8), save :: EPSV, EPSVT
!
! Shared variables
!
  integer, save :: ITL          ! index of liquid water
  integer, save :: ITI          ! index of ice water
  integer, save :: ICHNK        ! chunk identifier
!
  integer, save :: irank, ierror   ! to obtain RANK
!
! Tuning parameters set from namelist
!
  real(r8), save :: &
    CLMD,   & ! entrainment efficiency
    PA,     & ! factor for buoyancy to affect updraft velocity
    CPRES,  & ! pressure factor for momentum transport
    ALP0      ! alpha parameter in prognostic closure
!
! PUBLIC: interfaces
!
  public csconv_readnl   ! read csconv_nl namelist
  public cs_convi        ! CS scheme initialization
  public cs_convr        ! CS scheme main driver
  
contains
!---------------------------------------------------------------------------------
! temporarily following MIROC
function FQSAT( T, P )   ! calculate saturation water vapor 

  implicit none
  
  real(r8) :: FQSAT           ! saturation water vapor
  real(r8), intent(in) :: T   ! temperature [K]
  real(r8), intent(in) :: P   ! pressure [Pa]
  
  FQSAT = EPSV * ES0 / P &
        * EXP( (EL+EMELT/2._r8*(1._r8-SIGN(1._r8,T-TQICE))) &
               /RVAP *( 1._r8/TMELT - 1._r8/T )           )

end function FQSAT
!---------------------------------------------------------------------------------
! temporarily following MIROC
function FDQSAT( T, QS )   ! calculate d(qs)/dT

  implicit none
  
  real(r8) :: FDQSAT           ! d(QSAT)/d(T)
  real(r8), intent(in) :: T    ! temperature [K]
  real(r8), intent(in) :: QS   ! saturation water vapor [kg/kg]
  
  FDQSAT = (EL+EMELT/2._r8*(1._r8-SIGN(1._r8,T-TMELT))) &
         * QS / ( RVAP * T*T )

end function FDQSAT
!---------------------------------------------------------------------------------
subroutine csconv_readnl(nlfile)

   use namelist_utils, only: find_group_name
   use units,          only: getunit, freeunit
   use mpishorthand

   implicit none

   character(len=*), intent(in) :: nlfile   ! filepath for file containing namelist input

!
! Local variables
!
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'csconv_readnl'
!
   real(r8) :: &
     csconv_clmd  = unset_r8,   &
     csconv_pa    = unset_r8,   &
     csconv_cpres = unset_r8,   &
     csconv_alp0  = unset_r8

   namelist /csconv_nl/ csconv_clmd, csconv_pa, csconv_cpres, csconv_alp0

   if ( masterproc ) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name( unitn, 'csconv_nl', status=ierr )
      if ( ierr == 0 ) then
         read( unitn, csconv_nl, iostat=ierr)
         if ( ierr /= 0 ) then
            call endrun( subname // ':: ERROR reading namelist' )
         end if
      end if
      close( unitn )
      call freeunit( unitn )

      CLMD  = csconv_clmd
      PA    = csconv_pa
      CPRES = csconv_cpres
      ALP0  = csconv_alp0
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast( CLMD , 1, mpir8, 0, mpicom)
   call mpibcast( PA   , 1, mpir8, 0, mpicom)
   call mpibcast( CPRES, 1, mpir8, 0, mpicom) 
   call mpibcast( ALP0 , 1, mpir8, 0, mpicom) 
#endif

end subroutine csconv_readnl
!---------------------------------------------------------------------------------
subroutine cs_convi(limcnv_in, no_deep_pbl_in)

   use constituents, only: cnst_get_ind

   implicit none

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection (not supported yet)
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates CS convection entirely within PBL (not supported yet)
   
   EPSV  = RAIR / RVAP
   EPSVT = 1._r8 / EPSV - 1._r8
   
   call cnst_get_ind('CLDLIQ', ITL)
   call cnst_get_ind('CLDICE', ITI)

   if ( CLMD  == unset_r8 ) &
      call endrun( 'cs_convi: csconv_clmd must be set in the namelist.' )
   if ( PA    == unset_r8 ) &
      call endrun( 'cs_convi: csconv_pa must be set in the namelist.' )
   if ( CPRES == unset_r8 ) &
      call endrun( 'cs_convi: csconv_cpres must be set in the namelist.' )
   if ( ALP0  == unset_r8 ) &
      call endrun( 'cs_convi: csconv_alp0 must be set in the namelist.' )

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters cs_convi: CLMD' , CLMD
      write(iulog,*) 'tuning parameters cs_convi: PA'   , PA
      write(iulog,*) 'tuning parameters cs_convi: CPRES', CPRES
      write(iulog,*) 'tuning parameters cs_convi: ALP0' , ALP0
   end if

end subroutine cs_convi
!---------------------------------------------------------------------------------
subroutine cs_convr(lchnk   ,ncol    , &
                    t       ,q       ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    , &
                    delta   ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     , &
                    rliq    ,landfrac, &
                    u       ,v       ,utnd    ,vtnd    ,snow    , &
                    flxprec ,flxsnow ,cbmfx, &
                    sigmai, sigma)                            !DDsigma
!---------------------------------------------------------------------------------
! Purpose:
!
! Main driver for Chikira-Sugiyama convective scheme 
!
! Author: Minoru Chikira
!
!---------------------------------------------------------------------------------
   use ppgrid,        only: pcols, pver, pverp
   use cam_history,   only: outfld
   use constituents,  only: pcnst

   implicit none
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! temperature at mid-layer (K)
   real(r8), intent(in) :: q(pcols,pver,pcnst)    ! tracer array including moisture (kg/kg)
   real(r8), intent(in) :: pap(pcols,pver)        ! pressure at mid-layer (Pa)
   real(r8), intent(in) :: paph(pcols,pver+1)     ! pressure at boundaries (Pa)
   real(r8), intent(in) :: zm(pcols,pver)         ! height from surface at mid-layer (m)
   real(r8), intent(in) :: geos(pcols)            ! geopotential height at surface (m^2/s^2)
   real(r8), intent(in) :: zi(pcols,pver+1)       ! height from surface at boundaries (m)
   real(r8), intent(in) :: pblh(pcols)            ! PBL height (not used)
   real(r8), intent(in) :: tpert(pcols)           ! Theta perturbation in PBL (not used)
   real(r8), intent(in) :: landfrac(pcols)        ! RBN Landfrac
! added for cs_convr
   real(r8), intent(in) :: u(pcols,pver)          ! zonal wind at mid-layer (m/s)
   real(r8), intent(in) :: v(pcols,pver)          ! meridional wind at mid-layer (m/s)
   
   real(r8), intent(in) :: DELTA                  ! 2 delta t (model time increment in seconds)
!
! modified arguments
!
   real(r8), intent(inout) :: CBMFX(pcols,nctp)   ! cloud base mass flux (kg/m2/s)
!
! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver,pcnst)   ! tracer tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)         ! dry static energy tendency (W/kg)
   real(r8), intent(out) :: utnd(pcols,pver)         ! zonal wind tendency (m/s^2)
   real(r8), intent(out) :: vtnd(pcols,pver)         ! meridional wind tendency (m/s^2)
   
   real(r8), intent(out) :: mcon(pcols,pverp)      ! convective mass flux (kg/m2/s)
   real(r8), intent(out) :: dlf(pcols,pver)        ! scattered version of the detraining cld h2o tend (kg/kg/s)
   real(r8), intent(out) :: pflx(pcols,pverp)      ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)        ! condensation - evaporation
   real(r8), intent(out) :: cape(pcols)            ! convective available potential energy (J/kg)
   real(r8), intent(out) :: zdu(pcols,pver)        ! detraining mass flux from deep convection
   real(r8), intent(out) :: jctop(pcols)           ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot(pcols)           ! o row of base of cloud indices passed out.
   
   real(r8), intent(out) :: prec(pcols)            ! precipitation at surface (including snowfall) (m/s)
   real(r8), intent(out) :: snow(pcols)            ! snowfall at surface (m/s)
   real(r8), intent(out) :: rliq(pcols)            ! reserved liquid (not yet in cldliq) for energy integrals (m/s)
   real(r8), intent(out) :: flxprec(pcols,pverp)   ! precipitation flux (including snowfall) at interfaces (kg/m2/s)
   real(r8), intent(out) :: flxsnow(pcols,pverp)   ! snowfall flux at interfaces (kg/m2/s)

   !DDsigma - output added for sigma diagnostics
   real(r8), intent(out)   :: sigmai(pcols,pverp,nctp)  !DDsigma  sigma by cloud type - on interfaces (1=sfc)
   real(r8), intent(out)   :: sigma(pcols,pverp)        !DDsigma  sigma totaled over cloud type - on interfaces (1=sfc)
   real(r8)                :: sigmacs(pcols,pverp)        !DDsigma 
   !DDsigma

!
! output arguments of CS_CUMLUS
!
   real(r8) GTT(pcols,pver)                       ! temperature tendency [K/s]
   real(r8) GTQ(pcols,pver,pcnst)                 ! tracer tendency [kg/kg/s]
   real(r8) GTU(pcols,pver)                       ! zonal velocity tendency [m/s2]
   real(r8) GTV(pcols,pver)                       ! meridional velocity tendency [m/s2]
   real(r8) CMDET(pcols,pver)                     ! detrainment mass flux [kg/m2/s]
   real(r8) GTLDET(pcols,pver)                    ! cloud liquid tendency by detrainment [1/s]
   real(r8) GTIDET(pcols,pver)                    ! cloud ice tendency by detrainment [1/s]
   real(r8) GTPRP(pcols,pver+1)                   ! precipitation (including snowfall) flux at interfaces [kg/m2/s]
   real(r8) GSNWP(pcols,pver+1)                   ! snowfall flux at interfaces [kg/m2/s]
   real(r8) GMFX0(pcols,pver+1)                   ! updraft mass flux [kg/m2/s]
   integer  KT(pcols,nctp)                        ! cloud top index for each cloud type
!
! input arguments of CS_CUMLUS
!
   real(r8) GDT(pcols,pver)                       ! temperature [K]
   real(r8) GDQ(pcols,pver,pcnst)                 ! tracers including moisture [kg/kg]
   real(r8) GDU(pcols,pver)                       ! zonal wind [m/s]
   real(r8) GDV(pcols,pver)                       ! meridional wind [m/s]
   real(r8) GDTM(pcols,pver+1)                    ! temperature at boundaries of layers [K]
   real(r8) GDP(pcols,pver)                       ! pressure [Pa]
   real(r8) GDPM(pcols,pver+1)                    ! pressure at boundaries of layers [Pa]
   real(r8) GDZ(pcols,pver)                       ! altitude [m]
   real(r8) GDZM(pcols,pver+1)                    ! altitude at boundaries of layers [m]
   integer ISTS, IENS
!
! local variables
!
   real(r8) :: zs(pcols)                          ! surface height [m]
   real(r8) :: DELTI                              ! 1 delta t
   real(r8) :: ftem(pcols,pver)
   integer KTMAX(pcols)                           ! max of KT
   integer i, k, n

   call mpi_comm_rank( mpi_comm_world, irank, ierror )
!
! convert CAM input variables to MIROC counterparts
!
   ICHNK = lchnk
   DELTI = 0.5_r8*DELTA
   ISTS = 1
   IENS = ncol

   do i =1, ncol
      zs(i) = geos(i)/GRAV
   end do
   do k = 1, pver
      do i = 1, ncol
         GDT (i,pver-k+1) = t(i,k)
         GDU (i,pver-k+1) = u(i,k)
         GDV (i,pver-k+1) = v(i,k)
         GDZ (i,pver-k+1) = zm(i,k) + zs(i)
         GDP (i,pver-k+1) = pap(i,k)
      end do
   end do
   do k = 1, pver+1
      do i = 1, ncol
         GDZM(i,pver-k+2) = zi(i,k) + zs(i)
         GDPM(i,pver-k+2) = paph(i,k)
      end do
   end do
   do n = 1, pcnst
      do k = 1, pver
         do i = 1, ncol
            GDQ(i,pver-k+1,n) = q(i,k,n)
         end do
      end do
   end do
!
! calculate temperature at interfaces
!
   call TINTP( GDTM,               & ! output
               GDT, GDP, GDPM,     & ! input
               ISTS, IENS      )     ! active array size

   !DDsigma - initialize the sigma diagnostics
   sigmai = 0.0 !DDsigma
   sigma  = 0.0 !DDsigma
   sigmacs  = 0.0 !DDsigma
!
! call main routine
!
   call CS_CUMLUS &
      ( GTT   , GTQ   , GTU   , GTV   ,   & ! output
        CMDET , GTLDET, GTIDET,           & ! output
        GTPRP , GSNWP , GMFX0 ,           & ! output
        cape  , KT    ,                   & ! output
        CBMFX ,                           & ! modified
        GDT   , GDQ   , GDU   , GDV   ,   & ! input
        GDTM  ,                           & ! input
        GDP   , GDPM  , GDZ   , GDZM  ,   & ! input
        DELTA , DELTI , ISTS  , IENS  ,   & ! input
        sigmai, sigmacs)                      ! output !DDsigma
!
! convert MIROC output variables to CAM counterparts
!
  do n = 1, pcnst
     do k = 1, pver
        do i = 1, ncol
           qtnd(i,pver-k+1,n) = GTQ(i,k,n)
        end do
     end do
  end do
!
  do k = 1, pver
     do i = 1, ncol
        heat(i,pver-k+1) = CP*GTT(i,k) - EMELT*GTIDET(i,k)
        utnd(i,pver-k+1) = GTU(i,k)
        vtnd(i,pver-k+1) = GTV(i,k)
!        zdu (i,pver-k+1) = CMDET(i,k)
        dlf (i,pver-k+1) = GTLDET(i,k) + GTIDET(i,k)
        rliq(i) = ( GTLDET(i,k)+GTIDET(i,k) )*( GDPM(i,k+1)-GDPM(i,k) )/GRAV

        sigma(i,pver-k+1) = sigmacs(i,k)   !DDsigma

     end do
  end do
!
  do k = 1, pver+1
     do i = 1, ncol
        flxprec(i,pver-k+2) = GTPRP(i,k)
        flxsnow(i,pver-k+2) = GSNWP(i,k)
! Only updraft mass flux (excluding downdraft) is given to mcon
! since mcon is used to estimate cloud fraction.
        mcon   (i,pver-k+2) = GMFX0(i,k)
     end do
  end do
!
  KTMAX = 1
  do n = 1, nctp
     do i = 1, ncol
        KTMAX(i) = max( KTMAX(i), KT(i,n) )
     end do
  end do
!
  do i = 1, ncol
     jctop(i) = pver-KTMAX(i)+1
     prec(i) = GTPRP(i,1)/1000._r8   ! kg/m2/s => m/s
     snow(i) = GSNWP(i,1)/1000._r8   ! kg/m2/s => m/s
     rliq(i) = rliq(i)/1000._r8      ! kg/m2/s => m/s
  end do
  
  cme   = 0._r8   ! temporarily set to be zero
  pflx  = 0._r8   ! temporarily set to be zero
  zdu   = 0._r8   ! temporarily set to be zero
  jcbot = pver    ! set to be the lowest layer

  call outfld('CAPE'    , cape         , pcols, lchnk)
  ftem(:ncol,:pver) = heat(:ncol,:pver)/CP
  call outfld('CSDT'    , ftem         , pcols, lchnk)
  call outfld('CSDQ'    , qtnd(1,1,1)  , pcols, lchnk)
  call outfld('CSDLIQ'  , qtnd(1,1,ITL), pcols, lchnk)
  call outfld('CSDICE'  , qtnd(1,1,ITI), pcols, lchnk)
  call outfld('CSMTU'   , utnd         , pcols, lchnk)
  call outfld('CSMTV'   , vtnd         , pcols, lchnk)
  call outfld('CSFLXPRC', flxprec      , pcols, lchnk)
  call outfld('CSFLXSNW', flxsnow      , pcols, lchnk)
  call outfld('CSMU'    , mcon         , pcols, lchnk)
  call outfld('PRECCDCS', prec         , pcols, lchnk)

   !DDsigma 
  call outfld('SIGMA'   , sigma        , pcols, lchnk)
   !DDsigma 

end subroutine cs_convr
!************************************************************************
!* Original source code in MIROC5
!*
!* PACKAGE PCUMC  !!  physics: cumulus parameterization with
!*                             state-dependent entrainment rate
!*                             developed by Minoru Chikira
!* [Note]
!* -This routine works as the prognostic Arakawa-Schubert scheme
!*  if OPT_ASMODE is specified.
!* -Specify OPT_NS02 to use entrainment rate of Neggers et al. (2002)
!* -Specify OPT_CUMBGT to check water and energy budget.
!* -Specify OPT_CUMCHK to check range of output values.
!*
!*   [HIS] 08/09/19(chikira)   MIROC4.1
!*         08/10/30(hiro)      CMT modified
!*         08/11/11(chikira)   Neggers et al. (2002)
!*         08/12/3 (chikira)   downdraft detrainment modified
!*         08/12/3 (chikira)   COSP output
!*         09/02/24(chikira)   fix convective inhibition
!*         09/04/16(hiro)      CMIP5 output (cbasep,ctopp)
!*         09/09/03(yokohata)  COSP
!*         10/11/19(toshi)     small bug fix
!*         14/02/07(chikira)   CUMDWN bug fix, CMT modified
!************************************************************************
       SUBROUTINE CS_CUMLUS   & !! cumulus main routine
                ( GTT   , GTQ   , GTU   , GTV   ,   & ! output
                  CMDET , GTLDET, GTIDET,           & ! output
                  GTPRP , GSNWP , GMFX0 ,           & ! output
                  CAPE  , KT    ,                   & ! output
!                  CUMCLW, CUMFRC,
!*COSP
!                  QLIQC , QICEC , GPRCPF, GSNWPF,
!
!                  GTCFRC, FLIQC ,
!#ifdef OPT_CHASER
!     O            RFXC  , SFXC  , LEVCUM, LNFRC , REVC  , ! <<CHEM>>
!#endif
                  CBMFX ,                           & ! modified
                  GDT   , GDQ   , GDU   , GDV   ,   & ! input
                  GDTM  ,                           & ! input
                  GDP   , GDPM  , GDZ   , GDZM  ,   & ! input
!                  GDCFRC,
                  DELTA , DELTI , ISTS  , IENS  ,   & ! input
                  sigmai, sigma )                     ! output !DDsigma
! 
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [OUTPUT]
      REAL(r8), INTENT(OUT) :: GTT   ( IJSDIM, KMAX      ) !! heating rate
      REAL(r8), INTENT(OUT) :: GTQ   ( IJSDIM, KMAX, NTR ) !! change in q
      REAL(r8), INTENT(OUT) :: GTU   ( IJSDIM, KMAX      ) !! tendency of u
      REAL(r8), INTENT(OUT) :: GTV   ( IJSDIM, KMAX      ) !! tendency of v
      REAL(r8), INTENT(OUT) :: CMDET ( IJSDIM, KMAX      ) !! detrainment mass flux
      REAL(r8), INTENT(OUT) :: GTLDET( IJSDIM, KMAX      ) !! cloud liquid tendency by detrainment
      REAL(r8), INTENT(OUT) :: GTIDET( IJSDIM, KMAX      ) !! cloud ice tendency by detrainment
      REAL(r8), INTENT(OUT) :: GTPRP ( IJSDIM, KMAX+1    ) !! rain+snow flux
      REAL(r8), INTENT(OUT) :: GSNWP ( IJSDIM, KMAX+1    ) !! snowfall flux
      REAL(r8), INTENT(OUT) :: GMFX0 ( IJSDIM, KMAX+1    ) !! updraft mass flux
      REAL(r8), INTENT(OUT) :: CAPE  ( IJSDIM            )
      INTEGER , INTENT(OUT) :: KT    ( IJSDIM, NCTP      ) !! cloud top
!
!   [MODIFIED]
      REAL(r8), INTENT(INOUT) :: CBMFX ( IJSDIM, NCTP      ) !! cloud base mass flux
      
   !DDsigma - output added for sigma diagnostics
   real(r8), intent(out)   :: sigmai(IJSDIM,KMAX+1,nctp)  !DDsigma  sigma by cloud type - on interfaces (1=sfc)
   real(r8), intent(out)   :: sigma(IJSDIM,KMAX+1)        !DDsigma  sigma totaled over cloud type - on interfaces (1=sfc)
   !DDsigma

!
!   [INPUT]
      REAL(r8), INTENT(IN) :: GDT   ( IJSDIM, KMAX      ) !! temperature T
      REAL(r8), INTENT(IN) :: GDQ   ( IJSDIM, KMAX, NTR ) !! humidity, tracer
      REAL(r8), INTENT(IN) :: GDU   ( IJSDIM, KMAX      ) !! westerly u
      REAL(r8), INTENT(IN) :: GDV   ( IJSDIM, KMAX      ) !! southern wind v
      REAL(r8), INTENT(IN) :: GDTM  ( IJSDIM, KMAX+1    ) !! temperature T
      REAL(r8), INTENT(IN) :: GDP   ( IJSDIM, KMAX      ) !! pressure P
      REAL(r8), INTENT(IN) :: GDPM  ( IJSDIM, KMAX+1    ) !! pressure (half lev)
      REAL(r8), INTENT(IN) :: GDZ   ( IJSDIM, KMAX      ) !! altitude
      REAL(r8), INTENT(IN) :: GDZM  ( IJSDIM, KMAX+1    ) !! altitude
      REAL(r8), INTENT(IN) :: DELTA                       !! delta(t) (dynamics)
      REAL(r8), INTENT(IN) :: DELTI                       !! delta(t) (internal variable)
      INTEGER, INTENT(IN) :: ISTS, IENS   !! array range
!
!   [INTERNAL WORK]
      REAL(r8)     GPRCC ( IJSDIM, NTR       ) !! rainfall
      REAL(r8)     GSNWC ( IJSDIM            ) !! snowfall
      REAL(r8)     CUMCLW( IJSDIM, KMAX      ) !! cloud water in cumulus
      REAL(r8)     CUMFRC( IJSDIM            ) !! cumulus cloud fraction
!COSP
      REAL(r8)     QLIQC ( IJSDIM, KMAX   )    !! cumulus cloud liquid water [kg/kg]
      REAL(r8)     QICEC ( IJSDIM, KMAX   )    !! cumulus cloud ice [kg/kg]
      REAL(r8)     GPRCPF( IJSDIM, KMAX   )    !! rainfall flux at full level
      REAL(r8)     GSNWPF( IJSDIM, KMAX   )    !! snowfall flux at full level
!
      REAL(r8)     GTCFRC( IJSDIM, KMAX      ) !! change in cloud fraction
      REAL(r8)     FLIQC ( IJSDIM, KMAX      ) !! liquid ratio in cumulus
!
!#ifdef OPT_CHASER
!      REAL(r8)     RFXC  ( IJSDIM, KMAX+1    ) !! precipi. flx [kg/m2/s]
!      REAL(r8)     SFXC  ( IJSDIM, KMAX+1    ) !! ice/snow flx [kg/m2/s]
!      INTEGER      LEVCUM( IJSDIM, KMAX      ) !! flag for cum. cloud top
!      REAL(r8)     LNFRC ( IJSDIM, KMAX      ) !! areal rates of clouds
!      REAL(r8)     REVC  ( IJSDIM, KMAX      ) !! evaporation rates
!#endif
!
      REAL(r8)     GDCFRC( IJSDIM, KMAX      ) !! cloud fraction
!
      REAL(r8)     GDQI  ( IJSDIM, KMAX )      !! cloud ice
      REAL(r8)     GTQI  ( IJSDIM, KMAX )      !! tendency of cloud ice
      REAL(r8)     GTQL  ( IJSDIM, KMAX )      !! tendency of cloud liquid
!
      REAL(r8)     GDW   ( IJSDIM, KMAX )      !! total water
      REAL(r8)     DELP  ( IJSDIM, KMAX )
      REAL(r8)     GDQS  ( IJSDIM, KMAX )      !! saturate moisture
      REAL(r8)     FDQS  ( IJSDIM, KMAX )
      REAL(r8)     GAM   ( IJSDIM, KMAX )
      REAL(r8)     GDS   ( IJSDIM, KMAX )      !! dry static energy
      REAL(r8)     GDH   ( IJSDIM, KMAX )      !! moist static energy
      REAL(r8)     GDHS  ( IJSDIM, KMAX )      !! saturate MSE
!
      REAL(r8)     GCYM  ( IJSDIM, KMAX )      !! norm. mass flux (half lev)
      REAL(r8)     GCHB  ( IJSDIM )            !! cloud base MSE-Li*Qi
      REAL(r8)     GCWB  ( IJSDIM )            !! cloud base total water
      REAL(r8)     GCUB  ( IJSDIM )            !! cloud base U
      REAL(r8)     GCVB  ( IJSDIM )            !! cloud base V
      REAL(r8)     GCIB  ( IJSDIM )            !! cloud base ice
      REAL(r8)     ELAM  ( IJSDIM, KMAX, NCTP )   !! entrainment (rate*massflux)
      REAL(r8)     GCYT  ( IJSDIM, NCTP )      !! norm. mass flux @top
      REAL(r8)     GCHT  ( IJSDIM, NCTP )      !! cloud top MSE
      REAL(r8)     GCQT  ( IJSDIM, NCTP )      !! cloud top q
      REAL(r8)     GCUT  ( IJSDIM, NCTP )      !! cloud top U
      REAL(r8)     GCVT  ( IJSDIM, NCTP )      !! cloud top V
      REAL(r8)     GCLT  ( IJSDIM, NCTP )      !! cloud top cloud water
      REAL(r8)     GCIT  ( IJSDIM, NCTP )      !! cloud top cloud ice
      REAL(r8)     GTPRT ( IJSDIM, NCTP )      !! precipitation/M
      REAL(r8)     GCLZ  ( IJSDIM, KMAX )      !! cloud liquid for each CTP
      REAL(r8)     GCIZ  ( IJSDIM, KMAX )      !! cloud ice for each CTP
!
      REAL(r8)     ACWF  ( IJSDIM, NCTP )      !! cloud work function
      REAL(r8)     GPRCIZ( IJSDIM, KMAX+1 )    !! precipitation
      REAL(r8)     GSNWIZ( IJSDIM, KMAX+1 )    !! snowfall
      REAL(r8)     GTPRC0( IJSDIM       )      !! precip. before evap.
!
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 )    !! mass flux (updraft+downdraft)
      REAL(r8)     QLIQ  ( IJSDIM, KMAX   )    !! total cloud liquid
      REAL(r8)     QICE  ( IJSDIM, KMAX   )    !! total cloud ice
      REAL(r8)     GPRCI ( IJSDIM, KMAX   )    !! rainfall generation
      REAL(r8)     GSNWI ( IJSDIM, KMAX   )    !! snowfall generation
!
      REAL(r8)     GPRCP ( IJSDIM, KMAX+1 )    !! rainfall flux
!
      REAL(r8)     GTEVP ( IJSDIM, KMAX   )    !! evaporation+sublimation
      REAL(r8)     GMDD  ( IJSDIM, KMAX+1 )    !! downdraft mass flux
!
      REAL(r8)     CUMHGT( IJSDIM, NCTP   )    !! cloud top height
      REAL(r8)     CTOPP ( IJSDIM         )    !! cloud top pressure

      REAL(r8)     GDZTR ( IJSDIM         )   !! tropopause height
      REAL(r8)     FLIQOU( IJSDIM, KMAX   )   !! liquid ratio in cumulus
!#ifdef OPT_CHASER
!      REAL(r8)     TOPFLX( IJSDIM, NCTP   )    !! flux at each cloud top
!#endif
      INTEGER    KB    ( IJSDIM )
      INTEGER    KSTRT ( IJSDIM ) !! tropopause level
      REAL(r8)   GAMX
      REAL(r8)   CIN   ( IJSDIM )
      INTEGER    JBUOY ( IJSDIM )
      REAL(r8)   DELZ, BUOY, DELWC, DELER
      REAL(r8)   WCB   ( NCTP )                !! updraft velocity**2 @base
      SAVE       WCB
      REAL(r8)   WCBX
      REAL(r8)   ERMR  ( NCTP )                !! entrainment rate (ASMODE)
      SAVE       ERMR
      INTEGER    KTMX  ( NCTP )                !! max of cloud top
      INTEGER    KTMXT                         !! max of cloud top
      REAL(r8)   TIMED
      REAL(r8)   GDCLDX, GDMU2X, GDMU3X
      INTEGER    IERR
!
      LOGICAL    OOUT1, OOUT2
      INTEGER    KBMX, i
      INTEGER    IJ, K, CTP

      REAL(r8)     HBGT ( IJSDIM )               !! imbalance in column heat
      REAL(r8)     WBGT ( IJSDIM )               !! imbalance in column water
      
!DDsigma begin local work variables - all on model interfaces (sfc=1)
   REAL(r8)     lamdai         !! lamda for cloud type ctp
   REAL(r8)     lamdaprod( IJSDIM, KMAX+1   )   !! product of (1+lamda) through cloud type ctp
   REAL(r8)     gdrhom         !!  density
   REAL(r8)     gdtvm          !!  virtual temperature
   REAL(r8)     gdqm           !!  water vaper

   ! the following are new arguments to cumup to get them out 
   REAL(r8)     wcv( IJSDIM, KMAX+1   )        !! in-cloud vertical velocity
!DDsigma end local work variables

!
!   [INTERNAL PARM]
      REAL(r8) :: WCBMIN = 0._r8       !! min. of updraft velocity at cloud base
      REAL(r8) :: WCBMAX = 1.4_r8      !! max. of updraft velocity at cloud base
      REAL(r8) :: WCBAS  = 2._r8       !! updraft velocity**2 at cloud base (ASMODE)
      REAL(r8) :: ERAMIN = 1.e-5_r8    !! min. of entrainment rate
                                       !! used only in OPT_ASMODE
      REAL(r8) :: ERAMAX = 2.e-3_r8    !! max. of entrainment rate
                                       !! used only in OPT_ASMODE

      LOGICAL  :: OINICB = .false.     !! set 0.d0 to CBMFX

      REAL(r8) :: VARMIN = 1.e-13_r8   !! minimum of PDF variance
      REAL(r8) :: VARMAX = 5.e-7_r8    !! maximum of PDF variance
      REAL(r8) :: SKWMAX = 0.566_r8    !! maximum of PDF skewness

      REAL(r8) :: PSTRMX = 400.e2_r8   !! max P of tropopause
      REAL(r8) :: PSTRMN = 50.e2_r8    !! min P of tropopause
      REAL(r8) :: GCRSTR = 1.e-4_r8    !! crit. dT/dz tropopause

      LOGICAL, SAVE :: OTSPT1( NTR )   ! tracer transport by updraft, downdraft on/off
                                       ! should not include subgrid PDF and turbulence
      LOGICAL, SAVE :: OTSPT2( NTR )   ! tracer transport by subsidence on/off
                                       ! should include subgrid PDF and turbulence
      INTEGER, SAVE :: IMFXR ( NTR )
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
!
      LOGICAL, SAVE :: OFIRST = .TRUE.   !! called first time?
!
!   [ONCE]
      IF ( OFIRST ) THEN
         IF ( masterproc ) &
            WRITE ( iulog,* ) ' @@@ CHIKIRA-SUGIYAMA CUMULUS SCHEME'

         OFIRST = .FALSE.
!
         IF ( OINICB ) THEN
            IF ( masterproc ) &
               WRITE ( iulog,* ) &
                  ' ### PCUMC: OINICB=T - DEFAULT USED: CBMFX', 0.D0
            CBMFX = 0.D0
         END IF
!
         IF ( NCTP < 1 ) &
            CALL endrun( 'CS_CUMLUS: NCTP must be positive.' )
         DELWC  = ( WCBMAX-WCBMIN ) / DBLE( NCTP )
         DO CTP = 1, NCTP
            WCB ( CTP ) = ( CTP*DELWC )**2
         END DO
#ifdef OPT_ASMODE
         IF ( NCTP >= 2 ) THEN
            IF ( ERAMIN <= 0._r8 ) &
               CALL endrun( 'CS_CUMLUS: ERAMIN must be positive.' )
            DELER  = LOG10( ERAMAX / ERAMIN ) / DBLE( NCTP-1 )
            DO CTP = 1, NCTP
               ERMR( CTP ) = ERAMAX*( 10.D0**( -DELER*( CTP-1 ) ) )
            END DO
         ELSE
            ERMR( 1 ) = ERAMIN
         END IF
         IF ( masterproc ) &
            WRITE( iulog,* ) &
              ' ### PCUMC: ERMR =', ( ERMR( CTP ), CTP=1, NCTP )
#else
         ERMR = 0._r8
         IF ( masterproc ) &
            WRITE( iulog,* ) &
              ' ### PCUMC: WCB =', ( WCB( CTP ), CTP=1, NCTP )
#endif
         OTSPT1 = .false.
         OTSPT2 = .true.
         OTSPT2(1) = .false.
         OTSPT2(ITL) = .false.
         OTSPT2(ITI) = .false.

         IMFXR( :   ) = 0
         IMFXR( 1   ) = 1
         IMFXR( ITL ) = 1
         IMFXR( ITI ) = 1
      END IF
!
      GPRCC = 0.D0
      GTQ  = 0.D0
      GTQI = 0.D0
      GTQI = 0.D0
      GTT  = 0.D0
      GTU  = 0.D0
      GTV  = 0.D0
      GMFLX = 0.D0
      GMFX0 = 0.D0
      GPRCI = 0.D0
      GSNWI = 0.D0
      QLIQ  = 0.D0
      QICE  = 0.D0
      GTPRC0 = 0.D0
      GTCFRC = 0.D0
      CUMCLW = 0.D0
      FLIQC  = 0.D0
      FLIQOU = -999.D0
      HBGT  = 0.D0
      WBGT  = 0.D0
      GPRCPF = 0.D0
      GSNWPF = 0.D0
      GDZTR  = 0.D0
      KSTRT = KMAX
!#ifdef OPT_CHASER
!      TOPFLX = 0.D0
!      LNFRC  = 0.D0
!      REVC   = 0.D0
!      LEVCUM = 0
!#endif
      !GDQI = GDQ( :,:,ITI )
      !GDW  = GDQ( :,:,1 ) + GDQ( :,:,ITL ) + GDQI( :,: )
      GDQI(:IENS,:) = GDQ( :IENS,:,ITI )
      GDW(:IENS,:)  = GDQ( :IENS,:,1 ) + GDQ( :IENS,:,ITL ) + GDQI( :IENS,: )
!
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            DELP ( IJ,K ) = GDPM( IJ,K ) - GDPM( IJ,K+1 )
            GDQS ( IJ,K ) = FQSAT( GDT( IJ,K ), GDP( IJ,K ) )
            FDQS ( IJ,K ) = FDQSAT( GDT( IJ,K ), GDQS( IJ,K ) )
            GAM  ( IJ,K ) = EL/CP*FDQS( IJ,K )
            GDS  ( IJ,K ) = CP*GDT( IJ,K ) + GRAV*GDZ( IJ,K )
            GDH  ( IJ,K ) = GDS( IJ,K ) + EL*GDQ( IJ,K,1 )
            GDHS ( IJ,K ) = GDS( IJ,K ) + EL*GDQS( IJ,K )
         END DO
      END DO
!
!        < tropopause >
!
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            GAMX = ( GDTM( IJ,K+1 )-GDTM( IJ,K ) ) &
                  /( GDZM( IJ,K+1 )-GDZM( IJ,K ) )
            IF (   ( GDP(IJ,K).LT.PSTRMX .AND. GAMX.GT.GCRSTR ) &
                .OR. GDP(IJ,K) .LT. PSTRMN                     ) THEN
               KSTRT( IJ ) = MIN( K, KSTRT(IJ) )
            END IF
         END DO
      END DO
      DO IJ = ISTS, IENS
         K = KSTRT( IJ )
         GDZTR( IJ ) = GDZM( IJ,K )
      END DO
!
      CALL CUMBAS   & !! Cloud Base properties
         ( KB    , GCYM  , KBMX  ,           & ! output
           GCHB  , GCWB  , GCUB  , GCVB  ,   & ! output
           GCIB  ,                           & ! output
           GDH   , GDW   , GDHS  , GDQS  ,   & ! input
           GDQI  , GDU   , GDV   , GDZM  ,   & ! input
           GDPM  , FDQS  , GAM   ,           & ! input
           ISTS  , IENS                    )   ! input
           
!DDsigma some initialization
      lamdaprod = 1.0
!DDsigma end local work variables

!
      DO CTP = 1, NCTP
!
#ifdef OPT_ASMODE
         WCBX = WCBAS
#else
         WCBX = WCB( CTP )
#endif
         CALL CUMUP   & !! In-cloud Properties
            ( ACWF(1,CTP) , ELAM(1,1,CTP),                        & ! output
              GCLZ        , GCIZ        , GPRCIZ      , GSNWIZ,   & ! output
              GCYT(1,CTP) , GCHT(1,CTP) , GCQT (1,CTP),           & ! output
              GCLT(1,CTP) , GCIT(1,CTP) , GTPRT(1,CTP),           & ! output
              GCUT(1,CTP) , GCVT(1,CTP) ,                         & ! output
              KT  (1,CTP) , KTMX(CTP)   ,                         & ! output
              GCYM  ,                                             & ! modified
              wcv   ,                                             & ! !DD-sigma new output
              GCHB  , GCWB  , GCUB  , GCVB  ,                     & ! input
              GCIB  ,                                             & ! input
              GDU   , GDV   , GDH   , GDW   ,                     & ! input
              GDHS  , GDQS  , GDT   , GDTM  ,                     & ! input
              GDQ   , GDQI  , GDZ   , GDZM  ,                     & ! input
              GDPM  , FDQS  , GAM   , GDZTR ,                     & ! input
              CPRES , WCBX  , ERMR(CTP),                          & ! input
              KB    , CTP   , ISTS  , IENS   )                      ! input
!
         CALL CUMBMX   & !! Cloud Base Mass Flux
            ( CBMFX(1,CTP),                           & ! modified
              ACWF (1,CTP), GCYT(1,CTP), GDZM     ,   & ! input
              GDW         , GDQS       , DELP     ,   & ! input
              KT   (1,CTP), KTMX(CTP)  , KB       ,   & ! input
              DELTI       ,ISTS        , IENS       )
!
         CALL CUMFLX   & !! Cloud Mass Flux & Precip.
            ( GMFX0 , GPRCI , GSNWI ,                                 & ! output
              QLIQ  , QICE  , GTPRC0,                                 & ! output
!#ifdef OPT_CHASER
!     M        TOPFLX(1,CTP),                   ! <<CHEM>>
!#endif
              CBMFX(1,CTP), GCYM        , GPRCIZ     , GSNWIZ     ,   & ! input
              GTPRT(1,CTP), GCLZ        , GCIZ       ,                & ! input
              KB          , KT   (1,CTP), KTMX (CTP) ,                & ! input
              ISTS        , IENS                                   )    ! input

!DDsigma -  begin sigma computation
      do i = ISTS, IENS
        do k=kb(i),kt(i,ctp)                           !loop from cloud base to cloud top
          if(cbmfx(i,ctp) > 0.0) then                  ! this should avoid zero wcv in the denominator
            GDQM  = 0.5D0*( GDQ( I,K,1 )     + GDQ( I,K-1,1 ) )  ! as computed in cumup
            gdtvm = gdtm(i,k) * (1 + epsvt * gdqm)
            gdrhom = gdpm(i,k) / (rair * gdtvm)        ! gas law
            lamdai = gcym(i,k) * cbmfx(i,ctp) / (gdrhom*wcv(i,k))
            lamdaprod(i,k) = lamdaprod(i,k) * (1.+lamdai)
            sigmai(i,k,ctp) = lamdai / lamdaprod(i,k)
            sigma(i,k) = sigma(i,k) + sigmai(i,k,ctp)
          endif
        enddo
      enddo
!DDsigma -  end sigma computation

      END DO
!
      GMFLX( ISTS:IENS,: ) = GMFX0( ISTS:IENS,: )
      KTMXT = 3
      DO CTP = 1, NCTP
         IF ( KTMX( CTP ) .GT. KTMXT ) KTMXT = KTMX( CTP )
      END DO
      DO K = 1, KTMXT
         DO IJ = ISTS, IENS
            CUMCLW( IJ,K ) = QLIQ( IJ,K )+QICE( IJ,K )
            IF ( CUMCLW( IJ,K ) .GT. 0.D0 ) THEN
               FLIQC( IJ,K ) = QLIQ( IJ,K )/CUMCLW( IJ,K )
               FLIQOU( IJ,K ) = FLIQC( IJ,K )
            END IF
         END DO
      END DO
!
      CALL CUMCLD   & !! Cumulus Cloudiness
         ( CUMCLW, QLIQ  , QICE  , FLIQC  ,   & ! modified
           CUMFRC,                            & ! output
!#ifdef OPT_CHASER
!     M     LEVCUM, LNFRC ,           ! <<CHEM>>
!     I     TOPFLX,                   ! <<CHEM>>
!#endif
           GMFLX , KTMXT , ISTS  , IENS    )    ! input
!
      CALL CUMDET   & !! Cloud Detrainment Heating
         ( CMDET , GTLDET, GTIDET,                   & ! output
           GTT   , GTQ   , GTCFRC, GTU   , GTV   ,   & ! modified
           GTQI  ,                                   & ! modified
           GDH   , GDQ   , GDCFRC, GDU   , GDV   ,   & ! input
           CBMFX , GCYT  , DELP  , GCHT  , GCQT  ,   & ! input
           GCLT  , GCIT  , GCUT  , GCVT  , GDQI  ,   & ! input
           KT    , ISTS  , IENS                    )   ! input
!
      CALL CUMDWN   & !! Melt & Freeze & Evaporation
         ( GTT   , GTQ   , GTU   , GTV   ,   & ! modified
           GTQI  , GMFLX ,                   & ! modified
           GPRCP , GSNWP , GTEVP , GMDD  ,   & ! output
!#ifdef OPT_CHASER
!     O     REVC  ,                   ! <<CHEM>>
!#endif
           GPRCI , GSNWI ,                   & ! input
           GDH   , GDW   , GDQ   , GDQI  ,   & ! input
           GDQS  , GDS   , GDHS  , GDT   ,   & ! input
           GDU   , GDV   , GDZ   ,           & ! input
           GDZM  , GCYM  , FDQS  , DELP  ,   & ! input
           KB    , KTMXT , ISTS  , IENS    )   ! input
!
      GPRCC( ISTS:IENS,1 ) = GPRCP( ISTS:IENS,1 )
      GSNWC( ISTS:IENS   ) = GSNWP( ISTS:IENS,1 )
      GTPRP( ISTS:IENS,: ) = GPRCP( ISTS:IENS,: ) &
                           + GSNWP( ISTS:IENS,: )
!
      CALL CUMSBH   & !! Cloud Subsidence Heating
         ( GTT   , GTQ   , GTQI  ,        & ! modified
           GTU   , GTV   ,                & ! modified
           GDH   , GDQ   , GDQI  ,        & ! input
           GDU   , GDV   ,                & ! input
           DELP  , GMFLX , GMFX0 ,        & ! input
           KTMXT , CPRES , ISTS  , IENS )   ! input
!
      CALL CUMUPR   & !! Tracer Updraft
         ( GTQ   , GPRCC ,                           & ! modified
           GDQ   , CBMFX , ELAM  , GDZ   , GDZM  ,   & ! input
           GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
           GTPRT , GTEVP , GTPRC0,                   & ! input
           KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
           DELP  , OTSPT1, ISTS  , IENS            )   ! input
!
      CALL CUMDNR   & !! Tracer Downdraft
         ( GTQ   ,                        & ! modified
           GDQ   , GMDD  , DELP  ,        & ! input
           KTMXT , OTSPT1, ISTS  , IENS )   ! input
!
      CALL CUMSBR   & !! Tracer Subsidence
         ( GTQ   ,                   & ! modified
           GDQ   , DELP  ,           & ! input
           GMFLX , KTMXT , OTSPT2,   & ! input
           ISTS  , IENS            )   ! input
!
      GTQ( ISTS:IENS,:,ITI ) = GTQI( ISTS:IENS,: )
!
      CALL CUMFXR   & !! Tracer mass fixer without detrainment
         ( GTQ   ,                                   & ! modified
           GDQ   , DELP  , DELTA , KTMXT , IMFXR,    & ! input
           ISTS  , IENS                            )   ! input
!
      GTQL( ISTS:IENS,: ) = GTQ( ISTS:IENS,:,ITL ) &
                          + GTLDET( ISTS:IENS,: ) &
                          + GTIDET( ISTS:IENS,: )
!
      CALL CUMFXR1   & !! Tracer mass fixer with detrainment
         ( GTQL        ,                                     & ! modified
           GDQ(1,1,ITL), DELP, DELTA, KTMXT, IMFXR( ITL ),   & ! input
           ISTS        , IENS                              )   ! input
!
      GTLDET( ISTS:IENS,: ) = GTQL( ISTS:IENS,: ) &
                            - GTQ( ISTS:IENS,:,ITL ) &
                            - GTIDET( ISTS:IENS,: )
!
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
! tendencies of subgrid PDF (turned off)
!            GDCLDX = GDCFRC( IJ,K ) + GTCFRC( IJ,K )*DELTA
!            GDCLDX = MIN( MAX( GDCLDX, 0.D0 ), 1.D0 )
!            GTCFRC( IJ,K ) = ( GDCLDX - GDCFRC( IJ,K ) )/DELTA
!
!            GDMU2X = GDQ( IJ,K,IMU2 ) + GTQ( IJ,K,IMU2 )*DELTA
!            GDMU2X = MIN( MAX( GDMU2X,VARMIN ),VARMAX )
!            GDMU3X = GDQ( IJ,K,IMU3 ) + GTQ( IJ,K,IMU3 )*DELTA
!            GDMU3X = MIN( MAX( GDMU3X,-SKWMAX ),SKWMAX )
!            GTQ( IJ,K,IMU2 ) = ( GDMU2X - GDQ( IJ,K,IMU2 ))/DELTA
!            GTQ( IJ,K,IMU3 ) = ( GDMU3X - GDQ( IJ,K,IMU3 ))/DELTA
!
            HBGT( IJ ) = HBGT( IJ ) &
                       + ( CP*GTT(IJ,K) + EL*GTQ(IJ,K,1) &
                         - EMELT*( GTQ(IJ,K,ITI)+GTIDET(IJ,K) ) &
                         )*DELP(IJ,K)/GRAV
            WBGT( IJ ) = WBGT( IJ ) &
                       + ( GTQ(IJ,K,1) + GTQ(IJ,K,ITL) &
                         + GTQ(IJ,K,ITI) + GTLDET(IJ,K) &
                         + GTIDET(IJ,K) )*DELP(IJ,K)/GRAV
         END DO
      END DO
!
      DO IJ = ISTS, IENS
         HBGT( IJ ) = HBGT( IJ ) - EMELT*GSNWC( IJ )
         WBGT( IJ ) = WBGT( IJ ) + GPRCC( IJ,1 ) + GSNWC( IJ )
      END DO
!
      CTOPP( ISTS:IENS ) = 1.D6
      DO CTP = 1, NCTP
         DO IJ = ISTS, IENS
            IF ( KT( IJ,CTP ) .GT. KB( IJ ) ) THEN
               CUMHGT ( IJ,CTP ) = GDZ( IJ,KT( IJ,CTP ) )
               CTOPP( IJ ) = MIN( CTOPP( IJ ),GDP( IJ,KT( IJ,CTP ) ))
            ELSE
               CUMHGT ( IJ,CTP ) = -999.D0
            END IF
         END DO
      END DO
      DO IJ = ISTS, IENS
         IF( CTOPP( IJ ) .GE. 1.D6 ) THEN
            CTOPP( IJ ) = -999.D0
         END IF
      END DO
!
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            GPRCPF( IJ,K ) = 0.5D0*( GPRCP( IJ,K )+GPRCP( IJ,K+1 ) )
            GSNWPF( IJ,K ) = 0.5D0*( GSNWP( IJ,K )+GSNWP( IJ,K+1 ) )
         END DO
      END DO
!COSP
!necessary?
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            QLIQC( IJ,K ) = QLIQ( IJ,K )
            QICEC( IJ,K ) = QICE( IJ,K )
         END DO
      END DO
!
      CALL OUTFLD_CS('CSDLDET' , GTLDET       , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSDIDET' , GTIDET       , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSES'    , GDS          , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSEH'    , GDH          , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSEHS'   , GDHS         , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSEQS'   , GDQS         , IJSDIM, KMAX  , ICHNK ) 
      CALL OUTFLD_CS('CSCBMF01', CBMFX(1,NCTP), IJSDIM, 1     , ICHNK ) 
      CALL OUTFLD_CS('CSCWF01' , ACWF(1,NCTP) , IJSDIM, 1     , ICHNK ) 
      CALL OUTFLD_CS('CSHBGT'  , HBGT         , IJSDIM, 1     , ICHNK )
      CALL OUTFLD_CS('CSWBGT'  , WBGT         , IJSDIM, 1     , ICHNK )
!
!      CALL HISTNN( GCYT, 'GCYT', 'norm mass-flux at top',
!     &             'kg/m**2/s', 'A', HCLAS, NCTP )
!      CALL HISTNN( ACWF,  'CWF',   'cloud work function',
!     &             'J/kg'     , 'A', HCLAS, NCTP )
!      CALL HISTNN( CBMFX, 'CBMFX', 'cloud-base mass flux',
!     &             'kg/m**2/s', 'A', HCLAS, NCTP )
!      CALL HISTNN( CUMHGT, 'CUMHGT', 'cloud top height',
!     &             'm'        , 'A', HCLAS, NCTP )
!      CALL HISTIN( CTOPP, 'CTOPP', 'cloud top pressure',
!     &                                      'Pa'       , 'ASFC', HCLAS )
!      CALL HISTIN( GMFLX, 'CMFLX', 'cloud mass flux',
!     &                                      'kg/m**2/s', 'AMLV', HCLAS )
!      CALL HISTIN( GMDD,  'CDFLX', 'downdraft mass flux',
!     &                                      'kg/m**2/s', 'AMLV', HCLAS )
!      CALL HISTIN( GPRCP, 'FRANC', 'cumulus rain flux',
!     &                                      'kg/m**2/s', 'AMLV', HCLAS )
!      CALL HISTIN( GPRCPF, 'FRANCF', 'cumulus rain flux',
!     &                                      'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( GSNWP, 'FSNWC', 'cumulus snow flux',
!     &                                      'kg/m**2/s', 'AMLV', HCLAS )
!      CALL HISTIN( GSNWPF, 'FSNWCF', 'cumulus snow flux',
!     &                                      'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( GPRCP, 'FPRCC', 'cumulus rain+snow flux',
!     &                                      'kg/m**2/s', 'AMLV', HCLAS )
!      CALL HISTAD( GSNWP, 'FPRCC', 1.D0 )
!      CALL HISTIN( GTEVP, 'EPRCC', 'cumulus rain+snow evap',
!     &                                      'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( GDS,   'GDS',   'dry static energy',
!     &                                      'm**2/s**2', 'ALEV', HCLAS )
!      CALL HISTIN( GDH,   'GDH',   'moist static energy',
!     &                                      'm**2/s**2', 'ALEV', HCLAS )
!      CALL HISTIN( GDHS,  'GDHS',  'saturate moist static energy',
!     &                                      'm**2/s**2', 'ALEV', HCLAS )
!      CALL HISTRG( OOUT1, 'CAPE', 'conv. available potential energy',
!     &                                      'm**2/s**2', 'ASFC', HCLAS )
!      CALL HISTRG( OOUT2, 'CIN',  'conv. inhibition',
!     &                                      'm**2/s**2', 'ASFC', HCLAS )
!      CALL HISTIN( QLIQ,  'QLIQC', 'cumulus cloud liquid water',
!     &                                      'kg/kg',     'ALEV', HCLAS )
!      CALL HISTIN( QICE,  'QICEC', 'cumulus cloud ice',
!     &                                      'kg/kg',     'ALEV', HCLAS )
!      CALL HISTIN( FLIQOU, 'FLIQC', 'cumulus cloud liquid fraction',
!     &                                      '     ',     'ALEV', HCLAS )
!      CALL HISTIN( HBGT, 'DHBGTC', 'inbalance in column heat',
!     &                                    'J/m**2/s','ASFC', HCLAS )
!      CALL HISTIN( WBGT, 'DWBGTC', 'inbalance in column water',
!     &                                    'kg/m**2/s','ASFC', HCLAS )
!      CALL HISTIN( GDZTR, 'ZTROP', 'tropopause height',
!     &                                    'm', 'ASFC', HCLAS )
!
!      IF ( OOUT1 .OR. OOUT2 ) THEN
         CAPE  = 0.D0
         CIN   = 0.D0
         JBUOY = 0
         DO K = 2, KMAX
            DO IJ = ISTS, IENS
               IF ( K .GE. KB( IJ ) ) THEN
                  BUOY = ( GDH( IJ,1 )-GDHS( IJ,K ) ) &
                       / ( 1.D0+EL/CP*FDQS( IJ,K ) )  &
                       / ( CP*GDT( IJ,K ) )
               ELSE
                  BUOY = ( GDS( IJ,1 )-GDS( IJ,K ) ) &
                       / ( CP*GDT( IJ,K ) )
               END IF
               DELZ = GDZM( IJ,K+1 )-GDZM( IJ,K )
               IF ( BUOY .GT. 0.D0 .AND. &
                    JBUOY( IJ ) .NE. 0   ) THEN
                  CAPE( IJ ) = CAPE( IJ ) + BUOY*GRAV*DELZ
                  JBUOY( IJ ) = 2
               ELSE IF ( BUOY .LT. 0.D0 .AND. &
                         JBUOY( IJ ) .NE. 2 ) THEN
                  CIN( IJ ) = CIN( IJ ) - BUOY*GRAV*DELZ
                  JBUOY( IJ ) = 1
               END IF
            END DO
         END DO
         DO IJ = ISTS, IENS
            IF ( JBUOY(IJ) .NE. 2 ) CIN( IJ ) = -999.D0
         END DO
!         CALL HISTAX( CAPE, 'CAPE', 1.D0, .FALSE. )
!         CALL HISTAX( CIN,  'CIN',  1.D0, .FALSE. )
!      END IF
!
#ifdef OPT_CUMCHK
      CALL CUMCHK   & !! check range of output values
         ( GTT   , GTQ   , GTU   , GTV   ,   & ! input
           GTCFRC, GPRCC , GSNWC , CUMCLW,   & ! input
           CUMFRC, FLIQC , GTPRP ,           & ! input
           ISTS  , IENS                    )   ! input
#endif
!
      END SUBROUTINE CS_CUMLUS
!***********************************************************************
      SUBROUTINE CUMBAS   & !! cloud base
               ( KB    , GCYM  , KBMX  ,           & ! output
                 GCHB  , GCWB  , GCUB  , GCVB  ,   & ! output
                 GCIB  ,                           & ! output
                 GDH   , GDW   , GDHS  , GDQS  ,   & ! input
                 GDQI  , GDU   , GDV   , GDZM  ,   & ! input
                 GDPM  , FDQS  , GAM   ,           & ! input
                 ISTS  , IENS                    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!
!   [OUTPUT]
      INTEGER    KB    ( IJSDIM )         !! cloud base
      REAL(r8)   GCYM  ( IJSDIM, KMAX )   !! norm. mass flux (half lev)
      INTEGER    KBMX
      REAL(r8)   GCHB  ( IJSDIM )         !! cloud base MSE
      REAL(r8)   GCWB  ( IJSDIM )         !! cloud base total water
      REAL(r8)   GCUB  ( IJSDIM )         !! cloud base U
      REAL(r8)   GCVB  ( IJSDIM )         !! cloud base V
      REAL(r8)   GCIB  ( IJSDIM )         !! cloud base ice
!
!   [INPUT]
      REAL(r8)   GDH   ( IJSDIM, KMAX )        !! moist static energy
      REAL(r8)   GDW   ( IJSDIM, KMAX )        !! total water
      REAL(r8)   GDHS  ( IJSDIM, KMAX )        !! saturate MSE
      REAL(r8)   GDQS  ( IJSDIM, KMAX )        !! saturate humidity
      REAL(r8)   GDQI  ( IJSDIM, KMAX )        !! cloud ice
      REAL(r8)   GDU   ( IJSDIM, KMAX )        !! u-velocity
      REAL(r8)   GDV   ( IJSDIM, KMAX )        !! v-velocity
      REAL(r8)   GDZM  ( IJSDIM, KMAX+1 )      !! Altitude (half lev)
      REAL(r8)   GDPM  ( IJSDIM, KMAX+1 )      !! pressure (half lev)
      REAL(r8)   FDQS  ( IJSDIM, KMAX )
      REAL(r8)   GAM   ( IJSDIM, KMAX )
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)   CBASE ( IJSDIM )            !! cloud base height
      REAL(r8)   CBASEP( IJSDIM )            !! cloud base pressure
      REAL(r8)   DELZ, QSL, GAMX
      INTEGER    IJ, K
!
!   [INTERNAL PARM]
      INTEGER, PARAMETER :: KMAXM1 = KMAX-1
      INTEGER :: KLCLB = 1                   !! LCL base level
      INTEGER :: KCB   = 0                   !! fix cloud bottom
      INTEGER :: KBMAX = KMAXM1              !! cloud base max
      INTEGER :: KBOFS = 0                   !! cloud base offset
!
      GCYM( ISTS:IENS,: ) = 0.D0
!
      IF ( KCB .GT. 0 ) THEN
         DO IJ = ISTS, IENS
            KB( IJ ) = KCB
         END DO
      ELSE
         DO IJ = ISTS, IENS
            KB( IJ ) = KBMAX
         END DO
         DO K = KBMAX-1, KLCLB+1, -1
            DO IJ = ISTS, IENS
               GAMX = FDQS( IJ,K )/( 1.D0+GAM( IJ,K ) )/CP
               QSL = GDQS( IJ,K ) &
                   + GAMX*( GDH( IJ,KLCLB )-GDHS( IJ,K ) )
               IF ( GDW( IJ,KLCLB ) .GE. QSL ) THEN
                  KB( IJ ) = K + KBOFS
               END IF
            END DO
         END DO
      END IF
!
      KBMX = 1
      DO IJ = ISTS, IENS
         KBMX = MAX( KBMX, KB( IJ ) )
         CBASE ( IJ ) = GDZM( IJ,KB( IJ ) )-GDZM( IJ,1 )
         CBASEP( IJ ) = GDPM( IJ,KB( IJ ) )
      END DO
!
      DO K = 1, KBMX
         DO IJ = ISTS, IENS
            IF ( K .LE. KB( IJ ) ) THEN
               GCYM( IJ,K ) = ( GDZM( IJ,K      ) - GDZM( IJ,1 ) ) &
                             /( GDZM( IJ,KB(IJ) ) - GDZM( IJ,1 ) )
               GCYM( IJ,K ) = SQRT( GCYM( IJ,K ) )
            END IF
         END DO
      END DO
!
      GCHB( ISTS:IENS ) = 0.D0
      GCWB( ISTS:IENS ) = 0.D0
      GCUB( ISTS:IENS ) = 0.D0
      GCVB( ISTS:IENS ) = 0.D0
      GCIB( ISTS:IENS ) = 0.D0
!
      DO K = 1, KBMX
         DO IJ = ISTS, IENS
            IF ( K .LT. KB( IJ ) ) THEN
               DELZ       = GCYM( IJ,K+1 ) - GCYM( IJ,K )
               GCHB( IJ ) = GCHB( IJ ) + DELZ * GDH ( IJ,K )
               GCWB( IJ ) = GCWB( IJ ) + DELZ * GDW ( IJ,K )
               GCUB( IJ ) = GCUB( IJ ) + DELZ * GDU ( IJ,K )
               GCVB( IJ ) = GCVB( IJ ) + DELZ * GDV ( IJ,K )
               GCIB( IJ ) = GCIB( IJ ) + DELZ * GDQI( IJ,K )
            END IF
         END DO
      END DO
!
!      CALL HISTIN( CBASE , 'CBASE' , 'cloud base height',
!     &             'm'  , 'ASFC' , HCLAS )
!      CALL HISTIN( CBASEP, 'CBASEP', 'cloud base pressure',
!     &             'Pa' , 'ASFC' , HCLAS )
!
      END SUBROUTINE CUMBAS
!***********************************************************************
      SUBROUTINE CUMUP   & !! in-cloud properties
               ( ACWF  , ELAM  ,                   & ! output
                 GCLZ  , GCIZ  , GPRCIZ, GSNWIZ,   & ! output
                 GCYT  , GCHT  , GCQT  ,           & ! output
                 GCLT  , GCIT  , GTPRT ,           & ! output
                 GCUT  , GCVT  ,                   & ! output
                 KT    , KTMX  ,                   & ! output
                 GCYM  ,                           & ! modified
                 wcv   ,                           & ! !DDsigma new output
                 GCHB  , GCWB  , GCUB  , GCVB  ,   & ! input
                 GCIB  ,                           & ! input
                 GDU   , GDV   , GDH   , GDW   ,   & ! input
                 GDHS  , GDQS  , GDT   , GDTM  ,   & ! input
                 GDQ   , GDQI  , GDZ   , GDZM  ,   & ! input
                 GDPM  , FDQS  , GAM   , GDZTR ,   & ! input
                 CPRES , WCB   , ERMR  ,           & ! input
                 KB    , CTP   , ISTS  , IENS    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [OUTPUT]
      REAL(r8)   ACWF  ( IJSDIM         )   !! cloud work function
      REAL(r8)   ELAM  ( IJSDIM, KMAX   )   !! entrainment (rate*massflux)
      REAL(r8)   GCLZ  ( IJSDIM, KMAX   )   !! cloud liquid water*eta
      REAL(r8)   GCIZ  ( IJSDIM, KMAX   )   !! cloud ice*eta
      REAL(r8)   GPRCIZ( IJSDIM, KMAX   )   !! rain generation*eta
      REAL(r8)   GSNWIZ( IJSDIM, KMAX   )   !! snow generation*eta
      REAL(r8)   GCYT  ( IJSDIM         )   !! norm. mass flux @top
      REAL(r8)   GCHT  ( IJSDIM         )   !! cloud top MSE*eta
      REAL(r8)   GCQT  ( IJSDIM         )   !! cloud top moisture*eta
      REAL(r8)   GCLT  ( IJSDIM         )   !! cloud top liquid water*eta
      REAL(r8)   GCIT  ( IJSDIM         )   !! cloud top ice*eta
      REAL(r8)   GTPRT ( IJSDIM         )   !! cloud top (rain+snow)*eta
      REAL(r8)   GCUT  ( IJSDIM         )   !! cloud top u*eta
      REAL(r8)   GCVT  ( IJSDIM         )   !! cloud top v*eta
      INTEGER    KT    ( IJSDIM         )   !! cloud top
      INTEGER    KTMX                       !! max of cloud top
!
!   [MODIFIED]
      REAL(r8)   GCYM  ( IJSDIM, KMAX   )   !! norm. mass flux
!
!   [INPUT]
      REAL(r8)   GCHB  ( IJSDIM         )   !! MSE at cloud base
      REAL(r8)   GCWB  ( IJSDIM         )   !! total water @cloud base
      REAL(r8)   GCUB  ( IJSDIM         )   !! U at cloud base
      REAL(r8)   GCVB  ( IJSDIM         )   !! V at cloud base
      REAL(r8)   GCIB  ( IJSDIM         )   !! cloud ice at cloud base
      REAL(r8)   GDU   ( IJSDIM, KMAX   )   !! U
      REAL(r8)   GDV   ( IJSDIM, KMAX   )   !! V
      REAL(r8)   GDH   ( IJSDIM, KMAX   )   !! moist static energy
      REAL(r8)   GDW   ( IJSDIM, KMAX   )   !! total water
      REAL(r8)   GDHS  ( IJSDIM, KMAX   )   !! saturation MSE
      REAL(r8)   GDQS  ( IJSDIM, KMAX   )   !! saturation q
      REAL(r8)   GDT   ( IJSDIM, KMAX   )   !! T
      REAL(r8)   GDTM  ( IJSDIM, KMAX+1 )   !! T (half lev)
      REAL(r8)   GDQ   ( IJSDIM, KMAX, NTR )   !! q
      REAL(r8)   GDQI  ( IJSDIM, KMAX   )   !! cloud ice
      REAL(r8)   GDZ   ( IJSDIM, KMAX   )   !! z
      REAL(r8)   GDZM  ( IJSDIM, KMAX+1 )   !! z (half lev)
      REAL(r8)   GDPM  ( IJSDIM, KMAX+1 )   !! p (half lev)
      REAL(r8)   FDQS  ( IJSDIM, KMAX   )
      REAL(r8)   GAM   ( IJSDIM, KMAX   )
      REAL(r8)   GDZTR ( IJSDIM         )   !! tropopause height
      REAL(r8)   CPRES                      !! pres. fac. for cum. fric.
      REAL(r8)   WCB                        !! updraft velocity**2 @base
      REAL(r8)   ERMR                       !! entrainment rate (ASMODE)
      INTEGER    KB    ( IJSDIM         )
      INTEGER    CTP
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GCHMZ ( IJSDIM, KMAX+1 )   !! cloud h *eta (half lev)
      REAL(r8)     GCWMZ ( IJSDIM, KMAX+1 )   !! cloud Qt*eta (half lev)
      REAL(r8)     GCUMZ ( IJSDIM, KMAX+1 )   !! cloud U *eta (half lev)
      REAL(r8)     GCVMZ ( IJSDIM, KMAX+1 )   !! cloud V *eta (half lev)
      REAL(r8)     GCIMZ ( IJSDIM, KMAX+1 )   !! cloud Qi*eta (half lev)
      REAL(r8)     GTPRMZ( IJSDIM, KMAX+1 )   !! rain+snow *eta (half lev)
!
      REAL(r8)     BUOY  ( IJSDIM, KMAX   )   !! buoyancy
      REAL(r8)     BUOYM ( IJSDIM, KMAX+1 )   !! buoyancy (half lev)
      REAL(r8)     WCM   ( IJSDIM, KMAX+1 )   !! updraft velocity**2 (half lev)
      REAL(r8)     WCV   ( IJSDIM, KMAX+1 )   !! updraft velocity (half lev)
      REAL(r8)     GCY   ( IJSDIM, KMAX   )   !! norm. mass flux
      REAL(r8)     ELAR  ( IJSDIM, KMAX   )   !! entrainment rate
!
      REAL(r8)     GCHM  ( IJSDIM, KMAX+1 )   !! cloud MSE (half lev)
      REAL(r8)     GCWM  ( IJSDIM, KMAX+1 )   !! cloud Qt  (half lev)
      REAL(r8)     GCTM  ( IJSDIM, KMAX+1 )   !! cloud T (half lev)
      REAL(r8)     GCQM  ( IJSDIM, KMAX+1 )   !! cloud q (half lev)
      REAL(r8)     GCLM  ( IJSDIM, KMAX+1 )   !! cloud liquid ( half lev)
      REAL(r8)     GCIM  ( IJSDIM, KMAX+1 )   !! cloud ice (half lev)
      REAL(r8)     GCUM  ( IJSDIM, KMAX+1 )   !! cloud U (half lev)
      REAL(r8)     GCVM  ( IJSDIM, KMAX+1 )   !! cloud V (half lev)
!
      REAL(r8)     WCM_  ( IJSDIM         )
      REAL(r8)     ELARM1( IJSDIM         )
      REAL(r8)     GDZMKB( IJSDIM         )
      REAL(r8)     GDQSM, GDHSM, GDQM, GDCM, FDQSM, GCCM
      REAL(r8)     DELZ, ELADZ, DCTM , CPGM, DELC, FICE, ELARM2
      REAL(r8)     GCQMZ, GCCMZ, PRECR, GTPRIZ, DELZL
      REAL(r8)     GCWT, GCCT, DCT, WCVX
      REAL(r8)     PRCZH
      INTEGER      K, IJ
      CHARACTER    CTNUM*2
!
#ifdef OPT_CUMBGT
      REAL(r8)     HBGT  ( IJSDIM )           !! heat budget
      REAL(r8)     WBGT  ( IJSDIM )           !! water budget
      REAL(r8)     PBGT  ( IJSDIM )           !! precipitation budget
      REAL(r8)     MBGT  ( IJSDIM )           !! mass budget
      REAL(r8)     GTPRX ( IJSDIM )           !! (rain+snow)*eta at top
      REAL(r8)     GSNWT ( IJSDIM )           !! cloud top snow*eta
      REAL(r8)     HBMX, WBMX, PBMX, MBMX
      SAVE       HBMX, WBMX, PBMX, MBMX
#endif
!
!   [INTERNAL PARAM]

      REAL(r8), SAVE :: CLMP
      REAL(r8) ::  PRECZ0 = 1.5e3_r8
      REAL(r8) ::  PRECZH = 4.e3_r8
      REAL(r8) ::  ZTREF  = 1._r8
      REAL(r8) ::  PB     = 1._r8
      REAL(r8) ::  TAUZ   = 1.e4_r8
      REAL(r8) ::  ELMD   = 2.4e-3     !! for Neggers and Siebesma (2002)
      REAL(r8) ::  ELAMIN = 0._r8      !! min. of entrainment rate
      REAL(r8) ::  ELAMAX = 4.e-3      !! max. of entrainment rate
      REAL(r8) ::  WCCRT  = 0._r8
      REAL(r8) ::  TSICE  = 268.15_r8  !! compatible with macrop_driver
      REAL(r8) ::  TWICE  = 238.15_r8  !! compatible with macrop_driver
      REAL(r8) ::  EPS    = 1.e-10

      LOGICAL, SAVE :: OFIRST = .TRUE.
!
!   [INTERNAL FUNC]
      REAL(r8)     FPREC   !! precipitation ratio in condensate
      REAL(r8)     FRICE   !! ice ratio in cloud water
      REAL(r8)     Z       !! altitude
      REAL(r8)     ZH      !! scale height
      REAL(r8)     T       !! temperature
!
      FPREC( Z,ZH ) = MIN( MAX(1.D0-EXP(-(Z-PRECZ0)/ZH), 0.D0), 1.D0)
      FRICE( T ) = MIN( MAX( (TSICE-T)/(TSICE-TWICE), 0.D0 ), 1.D0 )
!
! Note: iteration is not made to diagnose cloud ice for simplicity
!
      IF ( OFIRST ) THEN
         CLMP = 2.D0*(1.D0-CLMD)*PA
         OFIRST = .FALSE.
      END IF

      ACWF  ( ISTS:IENS   ) = 0.D0
      ELAM  ( ISTS:IENS,: ) = unset_r8
      GCLZ  ( ISTS:IENS,: ) = 0.D0
      GCIZ  ( ISTS:IENS,: ) = 0.D0
      GPRCIZ( ISTS:IENS,: ) = 0.D0
      GSNWIZ( ISTS:IENS,: ) = 0.D0
      GCYT  ( ISTS:IENS   ) = 0.D0
      GCHT  ( ISTS:IENS   ) = 0.D0
      GCQT  ( ISTS:IENS   ) = 0.D0
      GCLT  ( ISTS:IENS   ) = 0.D0
      GCIT  ( ISTS:IENS   ) = 0.D0
      GTPRT ( ISTS:IENS   ) = 0.D0
      GCUT  ( ISTS:IENS   ) = 0.D0
      GCVT  ( ISTS:IENS   ) = 0.D0
!
      GCHMZ ( ISTS:IENS,: ) = 0.D0
      GCWMZ ( ISTS:IENS,: ) = 0.D0
      GCIMZ ( ISTS:IENS,: ) = 0.D0
      GCUMZ ( ISTS:IENS,: ) = 0.D0
      GCVMZ ( ISTS:IENS,: ) = 0.D0
      GTPRMZ( ISTS:IENS,: ) = 0.D0
!
      BUOY  ( ISTS:IENS,: ) = unset_r8
      BUOYM ( ISTS:IENS,: ) = unset_r8
      WCM   ( ISTS:IENS,: ) = unset_r8
      WCV   ( ISTS:IENS,: ) = unset_r8
      GCY   ( ISTS:IENS,: ) = unset_r8
      ELAR  ( ISTS:IENS,: ) = unset_r8
!
      GCHM  ( ISTS:IENS,: ) = unset_r8
      GCWM  ( ISTS:IENS,: ) = unset_r8
      GCTM  ( ISTS:IENS,: ) = unset_r8
      GCQM  ( ISTS:IENS,: ) = unset_r8
      GCLM  ( ISTS:IENS,: ) = unset_r8
      GCIM  ( ISTS:IENS,: ) = unset_r8
      GCUM  ( ISTS:IENS,: ) = unset_r8
      GCVM  ( ISTS:IENS,: ) = unset_r8
!#ifdef SYS_SX
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            IF ( K .GT. KB( IJ ) ) THEN
               GCYM( IJ,K ) = 0.D0
            END IF
         END DO
      END DO
!#else
!      DO IJ = ISTS, IENS
!         GCYM( IJ,KB(IJ)+1:KMAX ) = 0.D0
!      END DO
!#endif
      DO IJ = ISTS, IENS
         GDZMKB( IJ ) = GDZM( IJ,KB(IJ) )
      END DO
!
!     < cloud base properties >
!
      DO IJ = ISTS, IENS
         K = KB( IJ )
         GCHM( IJ,K ) = GCHB( IJ )
         GCWM( IJ,K ) = GCWB( IJ )
         WCM ( IJ,K ) = WCB
         GCUM( IJ,K ) = GCUB( IJ )
         GCVM( IJ,K ) = GCVB( IJ )
!
         GDQSM = FQSAT( GDTM( IJ,K ), GDPM( IJ,K ) )
         GDHSM = CP*GDTM( IJ,K ) + GRAV*GDZMKB( IJ ) + EL*GDQSM
         FDQSM = FDQSAT( GDTM( IJ,K ), GDQSM )
!
         DCTM  = ( GCHM( IJ,K ) - GDHSM )/( CP+EL*FDQSM )
         GCTM( IJ,K )  = GDT( IJ,K ) + DCTM
         GCQM( IJ,K )  = GDQSM + FDQSM*DCTM
         GCQM( IJ,K )  = MIN( GCQM( IJ,K ), GCWM( IJ,K ) )
         GCCM          = MAX( GCWM( IJ,K )-GCQM( IJ,K ), 0.D0 )
!
         GCIM( IJ,K ) = FRICE( GCTM( IJ,K ) )*GCCM
         GCLM( IJ,K ) = MAX( GCCM-GCIM( IJ,K ), 0.D0 )
         GCHM( IJ,K ) = GCHM( IJ,K )+EMELT*( GCIM( IJ,K )-GCIB( IJ ) )
         DCTM  = ( GCHM( IJ,K ) - GDHSM )/( CP+EL*FDQSM )
         GCTM( IJ,K ) = GDT( IJ,K ) + DCTM
!
         GDQM  = 0.5D0*( GDQ( IJ,K,1 ) + GDQ( IJ,K-1,1 ) )
         GDCM  = 0.5D0*( GDQ( IJ,K,ITL ) + GDQI( IJ,K ) &
                       + GDQ( IJ,K-1,ITL ) + GDQI( IJ,K-1 ) )
!
         BUOYM( IJ,K ) = ( DCTM/GDTM( IJ,K ) &
           + EPSVT*( GCQM(IJ,K)-GDQM )-GCCM+GDCM )*GRAV
!
#ifdef OPT_ASMODE
         ELARM1( IJ ) = ERMR
#elif defined OPT_NS02
         ELARM1( IJ ) = ELMD/SQRT( WCM( IJ,K ) )
#else
         ELARM1( IJ ) = CLMD*PA*BUOYM( IJ,K )/WCM( IJ,K )
#endif
         ELARM1( IJ ) = MIN( MAX( ELARM1( IJ ), ELAMIN ), ELAMAX )
!
         GCHMZ ( IJ,K ) = GCHM( IJ,K )
         GCWMZ ( IJ,K ) = GCWM( IJ,K )
         GCUMZ ( IJ,K ) = GCUM( IJ,K )
         GCVMZ ( IJ,K ) = GCVM( IJ,K )
         GCIMZ ( IJ,K ) = GCIM( IJ,K )
         WCM_( IJ )  = WCM( IJ,K )
      END DO
!
!     < in-cloud properties >
!
      DO K = 3, KMAX
         DO IJ = ISTS, IENS
            IF ( K .GT. KB( IJ ) .AND. WCM_( IJ ) .GT. WCCRT ) THEN
               WCV( IJ,K-1 ) = SQRT( MAX( WCM_( IJ ), 0.D0 ) )
               DELZ  = GDZM( IJ,K ) - GDZM( IJ,K-1 )
               GCYM( IJ,K ) = EXP( ELARM1( IJ )*DELZ )*GCYM( IJ,K-1 )
               ELADZ = GCYM( IJ,K ) - GCYM( IJ,K-1 )
!
               GCHMZ( IJ,K ) = GCHMZ( IJ,K-1 ) + GDH( IJ,K-1 )*ELADZ
               GCWMZ( IJ,K ) = GCWMZ( IJ,K-1 ) + GDW( IJ,K-1 )*ELADZ
!
               GDQSM = FQSAT( GDTM( IJ,K ), GDPM( IJ,K ) )
               GDHSM = CP*GDTM( IJ,K )+GRAV*GDZM( IJ,K )+EL*GDQSM
               FDQSM = FDQSAT( GDTM( IJ,K ), GDQSM )
               CPGM  = CP+EL*FDQSM
               PRCZH = PRECZH * MIN( GDZTR( IJ ) / ZTREF, 1.D0 )
               PRECR = FPREC( GDZM( IJ,K )-GDZMKB( IJ ),PRCZH )
!
               DCTM  = ( GCHMZ( IJ,K )/GCYM( IJ,K )-GDHSM )/CPGM
               GCQMZ = ( GDQSM+FDQSM*DCTM )*GCYM( IJ,K )
               GCQMZ = MIN( GCQMZ, GCWMZ( IJ,K ) )
               GTPRMZ( IJ,K ) = PRECR*( GCWMZ( IJ,K )-GCQMZ )
               GTPRMZ( IJ,K ) = MAX( GTPRMZ(IJ,K), GTPRMZ(IJ,K-1) )
               GCCMZ = GCWMZ( IJ,K )-GCQMZ-GTPRMZ( IJ,K )
               DELC  = MIN( GCCMZ, 0.D0 )
               GCCMZ = GCCMZ - DELC
               GCQMZ = GCQMZ + DELC
!
               FICE  = FRICE( GDTM( IJ,K )+DCTM )
               GCIMZ( IJ,K ) = FICE*GCCMZ
               GSNWIZ( IJ,K-1 ) = FICE*( GTPRMZ(IJ,K)-GTPRMZ(IJ,K-1) )
               GCHMZ( IJ,K ) = GCHMZ( IJ,K ) &
                 + EMELT*( GCIMZ( IJ,K   ) + GSNWIZ( IJ,K-1 ) &
                         - GCIMZ( IJ,K-1 ) - GDQI( IJ,K-1 )*ELADZ )
               DCTM  = ( GCHMZ( IJ,K )/GCYM( IJ,K )-GDHSM )/CPGM
!
               GDQM  = 0.5D0*( GDQ( IJ,K,1 ) + GDQ( IJ,K-1,1 ) )
               GDCM  = 0.5D0*( GDQ( IJ,K,ITL )+GDQI( IJ,K ) &
                             + GDQ( IJ,K-1,ITL )+GDQI( IJ,K-1 ) )
               GCQM( IJ,K ) = GCQMZ/GCYM( IJ,K )
               GCCM         = GCCMZ/GCYM( IJ,K )
!
               BUOYM( IJ,K ) = ( DCTM/GDTM( IJ,K ) &
                 + EPSVT*( GCQM(IJ,K)-GDQM )-GCCM+GDCM )*GRAV
               BUOY( IJ,K-1 ) = 0.5D0*( BUOYM( IJ,K )+BUOYM( IJ,K-1 ) )
!
#ifdef OPT_ASMODE
               WCM( IJ,K ) &
                 = ( WCM_( IJ ) + 2.D0*PA*DELZ*BUOY( IJ,K-1 ) ) &
                 / ( 1.D0 + 2.D0*PB*DELZ*ERMR )
#elif OPT_NS02
               WCM( IJ,K ) = WCM_( IJ ) &
                 + 2.D0*DELZ*( PA*BUOYM( IJ,K-1 )-ELMD*WCV( IJ,K-1 ) )
               WCM( IJ,K ) = MAX( WCM( IJ,K ), 0.D0 )
               WCVX = SQRT( 0.5D0*( WCM( IJ,K )+WCM_( IJ ) ) )
               WCM( IJ,K ) = WCM_( IJ ) &
                 + 2.D0*DELZ*( PA*BUOY( IJ,K-1 )-ELMD*WCVX )
#else
               IF ( BUOY( IJ,K-1 ) .GT. 0.D0 ) THEN
                  WCM ( IJ,K ) &
                    = ( WCM_( IJ ) + CLMP*DELZ*BUOY( IJ,K-1 ) ) &
                    / ( 1.D0 + DELZ/TAUZ )
               ELSE
                  WCM ( IJ,K ) &
                    = ( WCM_( IJ ) + 2.D0*PA*DELZ*BUOY( IJ,K-1 ) ) &
                    / ( 1.D0 + DELZ/TAUZ + 2.D0*DELZ*ELAMIN )
               END IF
#endif
!
#ifdef OPT_ASMODE
               ELARM2 = ERMR
#elif OPT_NS02
               ELARM2 = ELMD/SQRT( MAX( WCM( IJ,K ), EPS ) )
#else
               ELARM2 = CLMD*PA*BUOYM( IJ,K )/MAX( WCM( IJ,K ), EPS )
#endif
               ELARM2 = MIN( MAX( ELARM2, ELAMIN ), ELAMAX )
               ELAR( IJ,K-1 ) = 0.5D0*( ELARM1( IJ ) + ELARM2 )
               GCYM( IJ,K ) = EXP( ELAR(IJ,K-1)*DELZ )*GCYM( IJ,K-1 )
               ELADZ  = GCYM( IJ,K ) - GCYM( IJ,K-1 )
               ELAM( IJ,K-1 ) = ELADZ/DELZ
!
               GCHMZ( IJ,K ) = GCHMZ( IJ,K-1 ) + GDH( IJ,K-1 )*ELADZ
               GCWMZ( IJ,K ) = GCWMZ( IJ,K-1 ) + GDW( IJ,K-1 )*ELADZ
               GCUMZ( IJ,K ) = GCUMZ( IJ,K-1 ) + GDU( IJ,K-1 )*ELADZ
               GCVMZ( IJ,K ) = GCVMZ( IJ,K-1 ) + GDV( IJ,K-1 )*ELADZ
!
               DCTM  = ( GCHMZ( IJ,K )/GCYM( IJ,K )-GDHSM )/CPGM
               GCQMZ = ( GDQSM+FDQSM*DCTM )*GCYM( IJ,K )
               GCQMZ = MIN( GCQMZ, GCWMZ( IJ,K ) )
               GTPRMZ( IJ,K ) = PRECR*( GCWMZ( IJ,K )-GCQMZ )
               GTPRMZ( IJ,K ) = MAX( GTPRMZ(IJ,K), GTPRMZ(IJ,K-1) )
               GCCMZ = GCWMZ( IJ,K )-GCQMZ-GTPRMZ( IJ,K )
               DELC  = MIN( GCCMZ, 0.D0 )
               GCCMZ = GCCMZ - DELC
               GCQMZ = GCQMZ + DELC
               GCCM         = GCCMZ/GCYM( IJ,K )
               GCQM( IJ,K ) = GCQMZ/GCYM( IJ,K )
!
               FICE  = FRICE( GDTM( IJ,K )+DCTM )
               GCIMZ( IJ,K ) = FICE*GCCMZ
               GCIM( IJ,K ) = GCIMZ( IJ,K )/GCYM( IJ,K )
               GCLM( IJ,K ) = MAX( GCCM-GCIM( IJ,K ), 0.D0 )
               GTPRIZ = GTPRMZ( IJ,K ) - GTPRMZ( IJ,K-1 )
               GSNWIZ( IJ,K-1 ) = FICE*GTPRIZ
               GPRCIZ( IJ,K-1 ) = ( 1.D0-FICE )*GTPRIZ
               GCHMZ( IJ,K ) = GCHMZ( IJ,K ) &
                 + EMELT*( GCIMZ( IJ,K   ) + GSNWIZ( IJ,K-1 ) &
                         - GCIMZ( IJ,K-1 ) - GDQI( IJ,K-1 )*ELADZ )
               GCHM( IJ,K ) = GCHMZ( IJ,K )/GCYM( IJ,K )
               DCTM  = ( GCHM( IJ,K )-GDHSM )/CPGM
               GCTM( IJ,K ) = GDTM( IJ,K ) + DCTM
!
               GCWM( IJ,K ) = GCWMZ( IJ,K )/GCYM( IJ,K )
               GCUM( IJ,K ) = GCUMZ( IJ,K )/GCYM( IJ,K )
               GCVM( IJ,K ) = GCVMZ( IJ,K )/GCYM( IJ,K )
               DELZL = GDZ( IJ,K-1 )-GDZM( IJ,K-1 )
               GCY ( IJ,K-1 ) = GCYM(IJ,K-1)*EXP( ELAR(IJ,K-1)*DELZL )
               GCLZ( IJ,K-1 ) = 0.5D0*( GCLM( IJ,K ) + GCLM( IJ,K-1 ) ) &
                              * GCY( IJ,K-1 )
               GCIZ( IJ,K-1 ) = 0.5D0*( GCIM( IJ,K ) + GCIM( IJ,K-1 ) ) &
                              * GCY( IJ,K-1 )
               IF ( BUOY( IJ,K-1 ) .GT. 0.D0 ) THEN
                  ACWF( IJ ) = ACWF( IJ ) &
                             + BUOY( IJ,K-1 )*GCY( IJ,K-1 )*DELZ
               END IF
!
               ELARM1( IJ ) = ELARM2
               WCM_( IJ ) = WCM( IJ,K )
            END IF
         END DO
      END DO
!
!     < find cloud top >
!
      KT( ISTS:IENS ) = -1
      DO K = KMAX, 2, -1
         DO IJ = ISTS, IENS
            IF ( K .GT. KB( IJ ) &
                 .AND. KT( IJ ) .EQ. -1 &
                 .AND. BUOYM( IJ,K ) .GT. 0.D0 &
                 .AND. WCM  ( IJ,K ) .GT. WCCRT ) THEN
               KT( IJ ) = K
            END IF
         END DO
      END DO
!
      KTMX = 2
      DO IJ = ISTS, IENS
         IF ( KT( IJ ) .GT. KTMX ) KTMX = KT( IJ )
      END DO
!
      DO IJ = ISTS, IENS
         IF ( KT( IJ ) .GT. 0 ) THEN
            GCYM  ( IJ,KT( IJ )+1:KMAX ) = 0.D0
            GCLZ  ( IJ,KT( IJ ):KMAX ) = 0.D0
            GCIZ  ( IJ,KT( IJ ):KMAX ) = 0.D0
            GPRCIZ( IJ,KT( IJ ):KMAX ) = 0.D0
            GSNWIZ( IJ,KT( IJ ):KMAX ) = 0.D0
         ELSE
            GCYM  ( IJ,KB( IJ )+1:KMAX ) = 0.D0
            GCLZ  ( IJ,: ) = 0.D0
            GCIZ  ( IJ,: ) = 0.D0
            GPRCIZ( IJ,: ) = 0.D0
            GSNWIZ( IJ,: ) = 0.D0
         END IF
      END DO
!
!     < cloud top properties >
!
      DO IJ = ISTS, IENS
         IF ( KT( IJ ) .GT. 0 ) THEN
            K = KT( IJ )
            GCYT( IJ ) = GCY( IJ,K )
            ELADZ = GCYT( IJ ) - GCYM( IJ,K )
            ELAM( IJ,K ) = ELADZ/( GDZ( IJ,K )-GDZM( IJ,K ) )
!
            GCHT( IJ ) = GCHMZ( IJ,K ) + GDH( IJ,K )*ELADZ
            GCWT       = GCWMZ( IJ,K ) + GDW( IJ,K )*ELADZ
            GCUT( IJ ) = GCUMZ( IJ,K ) + GDU( IJ,K )*ELADZ
            GCVT( IJ ) = GCVMZ( IJ,K ) + GDV( IJ,K )*ELADZ
!
            DCT  = ( GCHT( IJ )/GCYT( IJ ) - GDHS( IJ,K ) ) &
                 / ( CP*( 1.D0 + GAM( IJ,K ) ) )
            GCQT( IJ ) = ( GDQS( IJ,K ) + FDQS( IJ,K )*DCT ) &
                       * GCYT( IJ )
            GCQT( IJ ) = MIN( GCQT( IJ ), GCWT )
            PRCZH      = PRECZH * MIN( GDZTR( IJ ) / ZTREF, 1.D0 )
            GTPRT( IJ ) = FPREC( GDZ( IJ,K )-GDZMKB( IJ ), &
                                 PRCZH ) &
                        *( GCWT-GCQT( IJ ) )
            GTPRT( IJ ) = MAX( GTPRT( IJ ), GTPRMZ( IJ,K ) )
            GCCT  = GCWT-GCQT( IJ )-GTPRT( IJ )
            DELC = MIN( GCCT, 0.D0 )
            GCCT = GCCT - DELC
            GCQT( IJ ) = GCQT( IJ ) + DELC
!
            FICE = FRICE( GDT( IJ,K )+DCT )
            GCIT( IJ ) = FICE*GCCT
            GCLT( IJ ) = ( 1.D0-FICE )*GCCT
            GTPRIZ = GTPRT( IJ ) - GTPRMZ( IJ,K )
            GPRCIZ( IJ,K ) = ( 1.D0-FICE )*GTPRIZ
            GSNWIZ( IJ,K ) = FICE*GTPRIZ
            GCHT( IJ ) = GCHT( IJ ) &
              + EMELT*( GCIT( IJ ) + GSNWIZ( IJ,K ) &
                      - GCIMZ( IJ,K ) - GDQI( IJ,K )*ELADZ )
!
            GCUT( IJ ) = GCUT( IJ )*( 1.D0-CPRES ) &
                       + GCY( IJ,K )*GDU( IJ,K )*CPRES
            GCVT( IJ ) = GCVT( IJ )*( 1.D0-CPRES ) &
                       + GCY( IJ,K )*GDV( IJ,K )*CPRES
            GCLZ( IJ,K ) = GCLT( IJ )
            GCIZ( IJ,K ) = GCIT( IJ )
         END IF
      END DO
!
#ifdef OPT_CUMBGT   /* budget check */
      HBGT ( ISTS:IENS ) = 0.D0
      WBGT ( ISTS:IENS ) = 0.D0
      PBGT ( ISTS:IENS ) = 0.D0
      MBGT ( ISTS:IENS ) = 0.D0
      GTPRX( ISTS:IENS ) = 0.D0
      GSNWT( ISTS:IENS ) = 0.D0
!
      IF ( CTP .EQ. 1 ) THEN
         HBMX = 0.D0
         WBMX = 0.D0
         PBMX = 0.D0
         MBMX = 0.D0
      END IF
!
      DO K = 2, KMAX
         DO IJ = ISTS, IENS
            IF ( K .GE. KB( IJ ) .AND. K .LT. KT( IJ ) ) THEN
               ELADZ = GCYM( IJ,K+1 ) - GCYM( IJ,K )
               DELZ  = GDZM( IJ,K+1 ) - GDZM( IJ,K )
               HBGT( IJ ) = HBGT( IJ ) &
                          + ( GDH( IJ,K )-EMELT*GDQI( IJ,K ) )*ELADZ
               WBGT( IJ ) = WBGT( IJ ) + GDW( IJ,K )*ELADZ
               MBGT( IJ ) = MBGT( IJ ) + ELAM( IJ,K )*DELZ
               GTPRX( IJ ) = GTPRX( IJ ) &
                           + GPRCIZ( IJ,K ) + GSNWIZ( IJ,K )
               GSNWT( IJ ) = GSNWT( IJ ) + GSNWIZ( IJ,K )
            END IF
         END DO
      END DO
!
      DO IJ = ISTS, IENS
         IF ( KT( IJ ) .GT. KB( IJ ) ) THEN
            ELADZ = GCYT( IJ ) - GCYM( IJ,KT(IJ) )
            DELZ  = GDZ( IJ,KT(IJ) )-GDZM( IJ,KT(IJ) )
            GTPRX( IJ ) = GTPRX( IJ ) &
                        + GPRCIZ( IJ,KT(IJ) ) + GSNWIZ( IJ,KT(IJ) )
            GSNWT( IJ ) = GSNWT( IJ ) + GSNWIZ( IJ,KT(IJ) )
            HBGT( IJ ) = HBGT( IJ ) &
                       + GCHB( IJ ) - EMELT*GCIB( IJ ) &
                       + ( GDH( IJ,KT(IJ) )-EMELT*GDQI( IJ,KT(IJ) ) ) &
                         *ELADZ &
                       - ( GCHT(IJ)-EMELT*( GCIT(IJ)+GSNWT(IJ) ) )
            WBGT( IJ ) = WBGT( IJ ) &
                       + GCWB( IJ ) + GDW( IJ,KT(IJ) )*ELADZ &
                       - GCQT( IJ ) - GCLT( IJ ) - GCIT( IJ ) &
                       - GTPRT( IJ )
            MBGT( IJ ) = MBGT( IJ ) + 1.D0 + ELAM( IJ,KT(IJ) )*DELZ &
                       - GCYT( IJ )
            PBGT( IJ ) = GTPRT( IJ ) - GTPRX( IJ )
!
            IF ( ABS( HBGT(IJ) ) .GT. ABS( HBMX ) ) HBMX = HBGT(IJ)
            IF ( ABS( WBGT(IJ) ) .GT. ABS( WBMX ) ) WBMX = WBGT(IJ)
            IF ( ABS( PBGT(IJ) ) .GT. ABS( PBMX ) ) PBMX = PBGT(IJ)
            IF ( ABS( MBGT(IJ) ) .GT. ABS( MBMX ) ) MBMX = MBGT(IJ)
         END IF
      END DO
!
      IF ( CTP .EQ. NCTP ) THEN
         WRITE( iulog,* ) &
            '### CUMUP(rank=',irank,'): energy imbalance =', HBMX
         WRITE( iulog,* ) &
            '### CUMUP(rank=',irank,'): water imbalance =', WBMX
         WRITE( iulog,* ) &
            '### CUMUP(rank=',irank,'): precipitation imbalance =', PBMX
         WRITE( iulog,* ) &
            '### CUMUP(rank=',irank,'): mass imbalance =', MBMX
      END IF
#endif
!
      CALL OUTFLD_CS('CSCH01', GCHM, IJSDIM, KMAX, ICHNK )
!
!      WRITE( CTNUM, '(I2.2)' ) CTP
!
!      CALL HISTIN( BUOY, 'BUOY'//CTNUM, 'cumulus buoyancy '//CTNUM,
!     &             'm/s**2', 'ALEV', HCLAS )
!      CALL HISTIN( ELAR, 'ELAR'//CTNUM, 'entrainment rate '//CTNUM,
!     &             '1/m', 'ALEV', HCLAS )
!      CALL HISTIN( GCHM, 'CUMH'//CTNUM,
!     &             'cum. moist static energy '//CTNUM,
!     &             'm**2/s**2', 'AMLV', HCLAS )
!      CALL HISTIN( GCWM, 'CUMW'//CTNUM, 'cumulus total water '//CTNUM,
!     &             'kg/kg', 'AMLV', HCLAS )
!      CALL HISTIN( WCV, 'WC'//CTNUM, 'updraft velocity '//CTNUM,
!     &             'm/s', 'AMLV', HCLAS )
!      CALL HISTIN( GCTM, 'CUMT'//CTNUM, 'cumulus temperature '//CTNUM,
!     &             'K', 'AMLV', HCLAS )
!      CALL HISTIN( GCQM, 'CUMQ'//CTNUM, 'cumulus water vapor '//CTNUM,
!     &             'kg/kg', 'AMLV', HCLAS )
!      CALL HISTIN( GCLM, 'CUML'//CTNUM, 'cumulus liquid '//CTNUM,
!     &             'kg/kg', 'AMLV', HCLAS )
!      CALL HISTIN( GCIM, 'CUMI'//CTNUM, 'cumulus ice '//CTNUM,
!     &             'kg/kg', 'AMLV', HCLAS )
!      CALL HISTIN( GCLM, 'CUMC'//CTNUM, 'cumulus cloud water '//CTNUM,
!     &             'kg/kg', 'AMLV', HCLAS )
!      CALL HISTAD( GCIM, 'CUMC'//CTNUM, 1.D0 )
!      CALL HISTIN( GCUM, 'CUMU'//CTNUM, 'cumulus U '//CTNUM,
!     &             'm/s', 'AMLV', HCLAS )
!      CALL HISTIN( GCVM, 'CUMV'//CTNUM, 'cumulus V '//CTNUM,
!     &             'm/s', 'AMLV', HCLAS )
!      CALL HISTIN( GPRCIZ, 'CUMRI'//CTNUM,
!     &             'cum. rain generation '//CTNUM,
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( GSNWIZ, 'CUMSI'//CTNUM,
!     &             'cum. snow generation '//CTNUM,
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( GPRCIZ, 'CUMPI'//CTNUM,
!     &             'cum. rain+snow generation '//CTNUM,
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTAD( GSNWIZ, 'CUMPI'//CTNUM, 1.D0 )
!      CALL HISTIN( GCYM, 'CUMY'//CTNUM, 'norm. mass flux '//CTNUM,
!     &             '', 'AMLV', HCLAS )
!
      END SUBROUTINE CUMUP
!***********************************************************************
      SUBROUTINE CUMBMX   & !! cloud base mass flux
               ( CBMFX ,                   & ! modified
                 ACWF  , GCYT  , GDZM  ,   & ! input
                 GDW   , GDQS  , DELP  ,   & ! input
                 KT    , KTMX  , KB    ,   & ! input
                 DELT  , ISTS  , IENS    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     CBMFX ( IJSDIM )            !! cloud base mass flux
!
!   [INPUT]
      REAL(r8)     ACWF  ( IJSDIM )            !! cloud work function
      REAL(r8)     GCYT  ( IJSDIM )            !! norm mass flux @top
      REAL(r8)     GDZM  ( IJSDIM, KMAX+1 )    !! height
      REAL(r8)     GDW   ( IJSDIM, KMAX   )    !! total water
      REAL(r8)     GDQS  ( IJSDIM, KMAX   )    !! saturate humidity
      REAL(r8)     DELP  ( IJSDIM, KMAX   )    !! delt pressure
      INTEGER      KT    ( IJSDIM )            !! cloud top
      INTEGER      KTMX                        !! max. of cloud top
      INTEGER      KB    ( IJSDIM )            !! cloud base
      REAL(r8)     DELT                        !! time step
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     QX    ( IJSDIM )
      REAL(r8)     QSX   ( IJSDIM )
      REAL(r8)     RHM   ( IJSDIM )
      INTEGER      IJ, K
      REAL(r8)     ALP
      REAL(r8)     FMAX1
!
!   [INTERNAL PARAM]
      REAL(r8) :: FMAX   = 1.5e-2_r8          !! maximum flux
      REAL(r8) :: RHMCRT = 0._r8              !! critical val. of RH@ all could
      REAL(r8) :: ALP1   = 0._r8
      REAL(r8) :: TAUD   = 1.e3_r8
      REAL(r8) :: ZFMAX  = 3.5e3_r8
      REAL(r8) :: ZDFMAX = 5.e2_r8
      REAL(r8) :: FMAXP  = 2._r8
      REAL(r8) :: EPS    = 1.e-10_r8
!
      QX ( ISTS:IENS )  = 0.D0
      QSX( ISTS:IENS )  = 0.D0
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            IF ( K .GE. KB( IJ ) .AND. K .LE. KT( IJ ) ) THEN
               QX ( IJ ) = QX ( IJ ) + GDW ( IJ,K ) * DELP( IJ,K )
               QSX( IJ ) = QSX( IJ ) + GDQS( IJ,K ) * DELP( IJ,K )
            END IF
         END DO
      END DO
      DO IJ = ISTS, IENS
         RHM( IJ ) = QX( IJ ) / MAX( QSX( IJ ), EPS )
      END DO
!
      DO IJ = ISTS, IENS
         IF ( KT( IJ ) .GT. KB( IJ ) &
              .AND. RHM( IJ ) .GE. RHMCRT ) THEN
            ALP = ALP0 + ALP1*( GDZM( IJ,KT(IJ) )-GDZM( IJ,KB(IJ) ) )
            FMAX1 = &
                   ( 1.D0 - TANH( (GDZM(IJ,1)-ZFMAX)/ZDFMAX  ) ) / 2.D0
            FMAX1 = FMAX * FMAX1**FMAXP
            CBMFX( IJ ) = CBMFX( IJ ) &
                        + MAX( ACWF( IJ ), 0.D0 )/( ALP*2.D0 )*DELT
            CBMFX( IJ ) = CBMFX( IJ )/( 1.D0 + DELT/(TAUD*2.D0) )
            CBMFX( IJ ) = MAX( CBMFX( IJ ), 0.D0 )
            CBMFX( IJ ) = MIN( CBMFX( IJ ), FMAX1/GCYT( IJ ) )
         ELSE
            CBMFX( IJ ) = 0.D0
         END IF
      END DO
!
      END SUBROUTINE CUMBMX
!***********************************************************************
      SUBROUTINE CUMFLX   & !! cloud mass flux
               ( GMFLX , GPRCI , GSNWI ,           & ! output
                 QLIQ  , QICE  , GTPRC0,           & ! output
!#ifdef OPT_CHASER
!     M           TOPFLX,                         ! <<CHEM>>
!#endif
                 CBMFX , GCYM  , GPRCIZ, GSNWIZ,   & ! input
                 GTPRT , GCLZ  , GCIZ  ,           & ! input
                 KB    , KT    , KTMX  ,           & ! input
                 ISTS  , IENS                    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!
!   [OUTPUT]
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 )     !! mass flux
      REAL(r8)     GPRCI ( IJSDIM, KMAX   )     !! rainfall generation
      REAL(r8)     GSNWI ( IJSDIM, KMAX   )     !! snowfall generation
      REAL(r8)     QLIQ  ( IJSDIM, KMAX   )     !! cloud liquid
      REAL(r8)     QICE  ( IJSDIM, KMAX   )     !! cloud ice
      REAL(r8)     GTPRC0( IJSDIM         )     !! precip. before evap.
!
!   [MODIFY]
!#ifdef OPT_CHASER
!      REAL(r8)     TOPFLX( IJSDIM         )     !! mass flux at cloud top
!#endif
!
!   [INPUT]
      REAL(r8)     CBMFX ( IJSDIM         )     !! cloud base mass flux
      REAL(r8)     GCYM  ( IJSDIM, KMAX   )     !! normalized mass flux
      REAL(r8)     GPRCIZ( IJSDIM, KMAX   )     !! precipitation/M
      REAL(r8)     GSNWIZ( IJSDIM, KMAX   )     !! snowfall/M
      REAL(r8)     GTPRT ( IJSDIM         )     !! rain+snow @top
      REAL(r8)     GCLZ  ( IJSDIM, KMAX   )     !! cloud liquid/M
      REAL(r8)     GCIZ  ( IJSDIM, KMAX   )     !! cloud ice/M
      INTEGER      KB    ( IJSDIM         )     !! cloud base
      INTEGER      KT    ( IJSDIM         )     !! cloud top
      INTEGER      KTMX                         !! max of cloud top
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER    IJ, K
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            GMFLX( IJ,K ) = GMFLX( IJ,K ) &
                          + GCYM( IJ,K )*CBMFX( IJ )
         END DO
      END DO
!
!#ifdef OPT_CHASER
!      DO IJ = ISTS, IENS
!         IF ( KT( IJ ) .GT. KB( IJ ) ) THEN
!           TOPFLX( IJ ) = GCYM( IJ,KT(IJ) )*CBMFX( IJ )
!         ELSE
!           TOPFLX( IJ ) = 0.D0
!         END IF
!      END DO
!#endif
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            GPRCI( IJ,K ) = GPRCI( IJ,K ) + GPRCIZ( IJ,K )*CBMFX( IJ )
            GSNWI( IJ,K ) = GSNWI( IJ,K ) + GSNWIZ( IJ,K )*CBMFX( IJ )
         END DO
      END DO
!
      DO IJ = ISTS, IENS
         GTPRC0( IJ ) = GTPRC0( IJ ) + GTPRT( IJ )*CBMFX( IJ )
      END DO
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            QLIQ( IJ,K ) = QLIQ( IJ,K ) + GCLZ( IJ,K )*CBMFX( IJ )
            QICE( IJ,K ) = QICE( IJ,K ) + GCIZ( IJ,K )*CBMFX( IJ )
         END DO
      END DO
!
      END SUBROUTINE CUMFLX
!***********************************************************************
      SUBROUTINE CUMDET   & !! detrainment
               ( CMDET , GTLDET, GTIDET,                   & ! output
                 GTT   , GTQ   , GTCFRC, GTU   , GTV   ,   & ! modified
                 GTQI  ,                                   & ! modified
                 GDH   , GDQ   , GDCFRC, GDU   , GDV   ,   & ! input
                 CBMFX , GCYT  , DELP  , GCHT  , GCQT  ,   & ! input
                 GCLT  , GCIT  , GCUT  , GCVT  , GDQI  ,   & ! input
                 KT    , ISTS  , IENS                    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [OUTPUT]
      REAL(r8)     CMDET ( IJSDIM, KMAX )   !! detrainment mass flux
      REAL(r8)     GTLDET( IJSDIM, KMAX )   !! cloud liquid tendency by detrainment
      REAL(r8)     GTIDET( IJSDIM, KMAX )   !! cloud ice tendency by detrainment
!
!   [MODIFY]
      REAL(r8)     GTT   ( IJSDIM, KMAX )   !! temperature tendency
      REAL(r8)     GTQ   ( IJSDIM, KMAX, NTR )   !! moisture tendency
      REAL(r8)     GTCFRC( IJSDIM, KMAX )   !! cloud fraction tendency
      REAL(r8)     GTU   ( IJSDIM, KMAX )   !! u tendency
      REAL(r8)     GTV   ( IJSDIM, KMAX )   !! v tendency
      REAL(r8)     GTQI  ( IJSDIM, KMAX )   !! cloud ice tendency
!
!   [INPUT]
      REAL(r8)     GDH   ( IJSDIM, KMAX )   !! moist static energy
      REAL(r8)     GDQ   ( IJSDIM, KMAX, NTR ) !! humidity qv
      REAL(r8)     GDCFRC( IJSDIM, KMAX )   !! cloud fraction
      REAL(r8)     GDU   ( IJSDIM, KMAX )
      REAL(r8)     GDV   ( IJSDIM, KMAX )
      REAL(r8)     DELP  ( IJSDIM, KMAX )
      REAL(r8)     CBMFX ( IJSDIM, NCTP )   !! cloud base mass flux
      REAL(r8)     GCYT  ( IJSDIM, NCTP )   !! detraining mass flux
      REAL(r8)     GCHT  ( IJSDIM, NCTP )   !! detraining MSE
      REAL(r8)     GCQT  ( IJSDIM, NCTP )   !! detraining qv
      REAL(r8)     GCLT  ( IJSDIM, NCTP )   !! detraining ql
      REAL(r8)     GCIT  ( IJSDIM, NCTP )   !! detraining qi
      REAL(r8)     GCUT  ( IJSDIM, NCTP )   !! detraining u
      REAL(r8)     GCVT  ( IJSDIM, NCTP )   !! detraining v
      REAL(r8)     GDQI  ( IJSDIM, KMAX )   !! cloud ice
      INTEGER      KT    ( IJSDIM, NCTP )   !! cloud top
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GTHCI, GTQVCI, GTQLCI, GTQICI
      REAL(r8)     GTCCI
      REAL(r8)     GTUCI, GTVCI
      REAL(r8)     GTXCI
      INTEGER      IJ, K, CTP
!
!
!PARALLEL_FORBID
      CMDET ( ISTS:IENS,: ) = 0.D0
      GTLDET( ISTS:IENS,: ) = 0.D0
      GTIDET( ISTS:IENS,: ) = 0.D0
!PARALLEL_FORBID
      DO CTP = 1, NCTP
         DO IJ = ISTS, IENS
            IF ( KT( IJ,CTP ) .GT. 0 ) THEN
               K = KT( IJ,CTP )
               GTXCI = GRAV/DELP( IJ,K )*CBMFX( IJ,CTP )
               GTHCI &
                = GTXCI &
                * ( GCHT( IJ,CTP ) - GCYT( IJ,CTP )*GDH( IJ,K ) )
               GTQVCI &
                = GTXCI &
                * ( GCQT( IJ,CTP ) - GCYT( IJ,CTP )*GDQ( IJ,K,1 ) )
               GTQLCI &
                = GTXCI * ( GCLT( IJ,CTP ) &
                           -GCYT( IJ,CTP )*GDQ( IJ,K,ITL ) )
               GTQICI &
                = GTXCI * ( GCIT( IJ,CTP ) &
                           -GCYT( IJ,CTP )*GDQI( IJ,K ) )
               GTCCI &
                = GTXCI &
                * ( GCYT( IJ,CTP ) - GCYT( IJ,CTP )*GDCFRC( IJ,K ) )
               GTUCI &
                = GTXCI &
                * ( GCUT( IJ,CTP ) - GCYT( IJ,CTP )*GDU( IJ,K ) )
               GTVCI &
                = GTXCI &
                * ( GCVT( IJ,CTP ) - GCYT( IJ,CTP )*GDV( IJ,K ) )
!
               GTQ( IJ,K,1   ) = GTQ( IJ,K,1 ) + GTQVCI
               GTT( IJ,K     ) = GTT( IJ,K   ) &
                               + ( GTHCI - EL*GTQVCI )/CP
! ql tendency by detrainment is treated by stratiform scheme
!               GTQ( IJ,K,ITL ) = GTQ( IJ,K,ITL ) + GTQLCI
               GTLDET( IJ,K )  = GTLDET( IJ,K ) + GTQLCI
! qi tendency by detrainment is treated by stratiform scheme
!               GTQI  ( IJ,K )  = GTQI( IJ,K ) + GTQICI
               GTIDET( IJ,K )  = GTIDET( IJ,K ) + GTQICI
               GTCFRC( IJ,K )  = GTCFRC( IJ,K )  + GTCCI
               GTU( IJ,K )     = GTU( IJ,K ) + GTUCI
               GTV( IJ,K )     = GTV( IJ,K ) + GTVCI
!
               CMDET( IJ,K ) = CMDET( IJ,K ) &
                             + GCYT( IJ,CTP )*CBMFX( IJ,CTP )
            END IF
         END DO
      END DO
!
!      CALL HISTIN( CMDET, 'CMDET', 'detrainment',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!
      END SUBROUTINE CUMDET
!***********************************************************************
      SUBROUTINE CUMSBH   & !! adiabat. descent
               ( GTT   , GTQ   , GTQI  ,        & ! modified
                 GTU   , GTV   ,                & ! modified
                 GDH   , GDQ   , GDQI  ,        & ! input
                 GDU   , GDV   ,                & ! input
                 DELP  , GMFLX , GMFX0 ,        & ! input
                 KTMX  , CPRES , ISTS  , IENS )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTT   ( IJSDIM, KMAX )   !! Temperature tendency
      REAL(r8)     GTQ   ( IJSDIM, KMAX, NTR )   !! Moisture etc tendency
      REAL(r8)     GTQI  ( IJSDIM, KMAX )
      REAL(r8)     GTU   ( IJSDIM, KMAX )   !! u tendency
      REAL(r8)     GTV   ( IJSDIM, KMAX )   !! v tendency
!
!   [INPUT]
      REAL(r8)     GDH   ( IJSDIM, KMAX )
      REAL(r8)     GDQ   ( IJSDIM, KMAX, NTR )   !! humidity etc
      REAL(r8)     GDQI  ( IJSDIM, KMAX )
      REAL(r8)     GDU   ( IJSDIM, KMAX )
      REAL(r8)     GDV   ( IJSDIM, KMAX )
      REAL(r8)     DELP  ( IJSDIM, KMAX )
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 )   !! mass flux (updraft+downdraft)
      REAL(r8)     GMFX0 ( IJSDIM, KMAX+1 )   !! mass flux (updraft only)
      INTEGER      KTMX
      REAL(r8)     CPRES                    !! pressure factor for cumulus friction
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      IJ, K, KM, KP
      REAL(r8)     SBH0, SBQ0, SBL0, SBI0, SBC0, SBS0
      REAL(r8)     SBH1, SBQ1, SBL1, SBI1, SBC1, SBS1, FX1
      REAL(r8)     SBU0, SBV0, SBU1, SBV1
      REAL(r8)     GTHCI, GTQVCI, GTQLCI, GTQICI
      REAL(r8)     GTM2CI, GTM3CI
      REAL(r8)     GTUCI, GTVCI
      REAL(r8)     FX    ( ISTS:IENS )

      REAL(r8) :: GTLSBH( IJSDIM, KMAX )
      REAL(r8) :: GTISBH( IJSDIM, KMAX )
!
!
      FX = 0.D0
      GTLSBH = 0._r8
      GTISBH = 0._r8
!
      DO 1100 K = KTMX, 1, -1
         KM = MAX( K-1, 1    )
         KP = MIN( K+1, KMAX )
         DO 1100 IJ = ISTS, IENS
            SBH0 = GMFLX( IJ,K+1 )*( GDH( IJ,KP )-GDH( IJ,K  ) )
            SBQ0 = GMFLX( IJ,K+1 )*( GDQ( IJ,KP,1 )-GDQ( IJ,K,1  ) )
            SBL0 = GMFLX( IJ,K+1 )*( GDQ( IJ,KP,ITL )-GDQ( IJ,K,ITL ) )
            SBI0 = GMFLX( IJ,K+1 )*( GDQI( IJ,KP )-GDQI( IJ,K ) )
            SBU0 = GMFLX( IJ,K+1 )*( GDU( IJ,KP )-GDU( IJ,K  ) ) &
                 - GMFX0( IJ,K+1 )*( GDU( IJ,KP )-GDU( IJ,K  ) )*CPRES
            SBV0 = GMFLX( IJ,K+1 )*( GDV( IJ,KP )-GDV( IJ,K  ) ) &
                 - GMFX0( IJ,K+1 )*( GDV( IJ,KP )-GDV( IJ,K  ) )*CPRES
!
            SBH1 = GMFLX( IJ,K   )*( GDH( IJ,K  )-GDH( IJ,KM ) )
            SBQ1 = GMFLX( IJ,K   )*( GDQ( IJ,K,1   )-GDQ( IJ,KM,1 ) )
            SBL1 = GMFLX( IJ,K   )*( GDQ( IJ,K,ITL )-GDQ( IJ,KM,ITL ) )
            SBI1 = GMFLX( IJ,K   )*( GDQI( IJ,K )-GDQI( IJ,KM ) )
            SBU1 = GMFLX( IJ,K   )*( GDU( IJ,K  )-GDU( IJ,KM ) ) &
                 - GMFX0( IJ,K   )*( GDU( IJ,K  )-GDU( IJ,KM ) )*CPRES
            SBV1 = GMFLX( IJ,K   )*( GDV( IJ,K  )-GDV( IJ,KM ) ) &
                 - GMFX0( IJ,K   )*( GDV( IJ,K  )-GDV( IJ,KM ) )*CPRES
!
!#ifndef SYS_SX   /* original */
            IF ( GMFLX( IJ,K ) .GT. GMFLX( IJ,K+1 ) ) THEN
               FX1 = 0.5D0
            ELSE
               FX1 = 0.D0
            END IF
!#else            /* optimized for NEC SX series */
!            FX1 = 0.25D0 - SIGN(0.25D0,GMFLX(IJ,K+1)-GMFLX(IJ,K)) !! 0.5 or 0.
!#endif
!
            GTHCI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBH0 + FX1 *SBH1 )
            GTQVCI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBQ0 + FX1 *SBQ1 )
            GTQLCI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBL0 + FX1 *SBL1 )
            GTQICI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBI0 + FX1 *SBI1 )
            GTUCI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBU0 + FX1 *SBU1 )
            GTVCI = GRAV/DELP( IJ,K ) &
                   *( ( 1.D0-FX( IJ ) )*SBV0 + FX1 *SBV1 )
!
            GTT ( IJ,K     ) = GTT( IJ,K )+( GTHCI-EL*GTQVCI )/CP
            GTQ ( IJ,K,1   ) = GTQ( IJ,K,1 ) + GTQVCI
            GTQ ( IJ,K,ITL ) = GTQ( IJ,K,ITL  ) + GTQLCI
            GTQI( IJ,K ) = GTQI( IJ,K ) + GTQICI
            GTU ( IJ,K ) = GTU( IJ,K ) + GTUCI
            GTV ( IJ,K ) = GTV( IJ,K ) + GTVCI

            GTLSBH( IJ,K ) = GTQLCI
            GTISBH( IJ,K ) = GTQICI
!
!            SBC0 = GMFLX( IJ,K+1 )*( GDQ( IJ,KP,IMU2 )-GDQ( IJ,K,IMU2 ))
!            SBS0 = GMFLX( IJ,K+1 )*( GDQ( IJ,KP,IMU3 )-GDQ( IJ,K,IMU3 ))
!            SBC1 = GMFLX( IJ,K   )*( GDQ( IJ,K,IMU2 )-GDQ( IJ,KM,IMU2 ))
!            SBS1 = GMFLX( IJ,K   )*( GDQ( IJ,K,IMU3 )-GDQ( IJ,KM,IMU3 ))
!            GTM2CI = GRAV/DELP( IJ,K )
!     &             *( ( 1.D0-FX( IJ ) )*SBC0 + FX1 *SBC1 )
!            GTM3CI = GRAV/DELP( IJ,K )
!     &             *( ( 1.D0-FX( IJ ) )*SBS0 + FX1 *SBS1 )
!            GTQ( IJ,K,IMU2 ) = GTQ( IJ,K,IMU2 ) + GTM2CI
!            GTQ( IJ,K,IMU3 ) = GTQ( IJ,K,IMU3 ) + GTM3CI
!
            FX ( IJ )   = FX1
 1100 CONTINUE
!
      CALL OUTFLD_CS('CSDLSB', GTLSBH, IJSDIM, KMAX, ICHNK )
      CALL OUTFLD_CS('CSDISB', GTISBH, IJSDIM, KMAX, ICHNK )
!
      END SUBROUTINE CUMSBH
!***********************************************************************
      SUBROUTINE CUMDWN   & !! Freeze & Melt & Evaporation
               ( GTT   , GTQ   , GTU   , GTV   ,   & ! modified
                 GTQI  , GMFLX ,                   & ! modified
                 GPRCP , GSNWP , GTEVP , GMDD  ,   & ! output
!#ifdef OPT_CHASER
!     O           REVC  ,                   ! <<CHEM>>
!#endif
                 GPRCI , GSNWI ,                   & ! input
                 GDH   , GDW   , GDQ   , GDQI  ,   & ! input
                 GDQS  , GDS   , GDHS  , GDT   ,   & ! input
                 GDU   , GDV   , GDZ   ,           & ! input
                 GDZM  , GCYM  , FDQS  , DELP  ,   & ! input
                 KB    , KTMX  , ISTS  , IENS    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTT   ( IJSDIM, KMAX )   !! Temperature tendency
      REAL(r8)     GTQ   ( IJSDIM, KMAX, NTR )  !! Moisture etc tendency
      REAL(r8)     GTU   ( IJSDIM, KMAX )   !! u tendency
      REAL(r8)     GTV   ( IJSDIM, KMAX )   !! v tendency
      REAL(r8)     GTQI  ( IJSDIM, KMAX )   !! cloud ice tendency
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 ) !! mass flux
!
!   [OUTPUT]
      REAL(r8)     GPRCP ( IJSDIM, KMAX+1 ) !! rainfall flux
      REAL(r8)     GSNWP ( IJSDIM, KMAX+1 ) !! snowfall flux
      REAL(r8)     GTEVP ( IJSDIM, KMAX   ) !! evaporation+sublimation
      REAL(r8)     GMDD  ( IJSDIM, KMAX+1 ) !! downdraft mass flux
!#ifdef OPT_CHASER
!      REAL(r8)     REVC  ( IJSDIM, KMAX )   !! evapo. rate <<CHEM>>
!#endif
!
!   [INPUT]
      REAL(r8)     GPRCI ( IJSDIM, KMAX )   !! rainfall generation
      REAL(r8)     GSNWI ( IJSDIM, KMAX )   !! snowfall generation
      REAL(r8)     GDH   ( IJSDIM, KMAX )   !! moist static energy
      REAL(r8)     GDW   ( IJSDIM, KMAX )   !! total water
      REAL(r8)     GDQ   ( IJSDIM, KMAX, NTR )   !! humidity etc
      REAL(r8)     GDQI  ( IJSDIM, KMAX )   !! cloud ice
      REAL(r8)     GDQS  ( IJSDIM, KMAX )   !! saturate humidity
      REAL(r8)     GDS   ( IJSDIM, KMAX )   !! dry static energy
      REAL(r8)     GDHS  ( IJSDIM, KMAX ) !! saturate moist static energy
      REAL(r8)     GDT   ( IJSDIM, KMAX )   !! air temperature T
      REAL(r8)     GDU   ( IJSDIM, KMAX )   !! u-velocity
      REAL(r8)     GDV   ( IJSDIM, KMAX )   !! v-velocity
      REAL(r8)     GDZ   ( IJSDIM, KMAX )   !! altitude
      REAL(r8)     GDZM  ( IJSDIM, KMAX+1 ) !! altitude (half lev)
      REAL(r8)     GCYM  ( IJSDIM, KMAX )   !! norm. mass flux
      REAL(r8)     FDQS  ( IJSDIM, KMAX )
      REAL(r8)     DELP  ( IJSDIM, KMAX )
      INTEGER      KB    ( IJSDIM )
      INTEGER      KTMX
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
! Note: Some variables have 3-dimensions for the purpose of budget check.
      REAL(r8)     EVAPD ( IJSDIM, KMAX )      !! evap. in downdraft
      REAL(r8)     SUBLD ( IJSDIM, KMAX )      !! subl. in downdraft
      REAL(r8)     EVAPE ( IJSDIM, KMAX )      !! evap. in environment
      REAL(r8)     SUBLE ( IJSDIM, KMAX )      !! subl. in environment
      REAL(r8)     EVAPX ( IJSDIM, KMAX )      !! evap. env. to DD
      REAL(r8)     SUBLX ( IJSDIM, KMAX )      !! subl. env. to DD
      REAL(r8)     GMDDE ( IJSDIM, KMAX )      !! downdraft entrainment
      REAL(r8)     SNMLT ( IJSDIM, KMAX )      !! melt - freeze
      REAL(r8)     GCHDD ( IJSDIM, KMAX )      !! MSE detrainment
      REAL(r8)     GCWDD ( IJSDIM, KMAX )      !! water detrainment
      REAL(r8)     GTTEV ( IJSDIM, KMAX )      !! T tendency by evaporation
      REAL(r8)     GTQEV ( IJSDIM, KMAX )      !! q tendency by evaporation
      REAL(r8)     GCHD  ( ISTS:IENS )         !! downdraft MSE
      REAL(r8)     GCWD  ( ISTS:IENS )         !! downdraft q
      REAL(r8)     GCUD  ( ISTS:IENS )         !! downdraft u
      REAL(r8)     GCVD  ( ISTS:IENS )         !! downdraft v
      REAL(r8)     FSNOW ( ISTS:IENS )
      REAL(r8)     GMDDD ( ISTS:IENS )
      INTEGER      IJ, K
      REAL(r8)     FMELT( IJSDIM,KMAX )
      REAL(r8)     FEVP ( IJSDIM,KMAX )
      REAL(r8)     GDTW
      REAL(r8)     GCHX, GCTX, GCQSX, GTPRP, EVSU, GTEVE, LVIC
      REAL(r8)     DQW, DTW, GDQW, DZ, GCSD, FDET, GDHI
      REAL(r8)     GMDDX, GMDDMX
      REAL(r8)     GCHDX, GCWDX
      REAL(r8)     GCUDD, GCVDD
      REAL(r8)     GTHCI, GTQVCI, GTQLCI, GTQICI, GTUCI, GTVCI
#ifdef OPT_CUMBGT
      REAL(r8)     WBGT  ( ISTS:IENS )         !! water budget
      REAL(r8)     HBGT  ( ISTS:IENS )         !! energy budget
      REAL(r8)     DDWBGT( ISTS:IENS )         !! downdraft water budget
      REAL(r8)     DDHBGT( ISTS:IENS )         !! downdraft energy budget
      REAL(r8)     WMX, HMX, DDWMX, DDHMX
#endif
!
!   [INTERNAL PARM]
      REAL(r8) :: TWSNOW = 273.15_r8   !! wet-bulb temp. rain/snow
      REAL(r8) :: FTMLT  = 4._r8       !! temp. factor for melt
      REAL(r8) :: GMFLXC = 5.e-2_r8    !! critical mass flux
      REAL(r8) :: VTERMS = 2._r8       !! terminal velocity of snowflake
      REAL(r8) :: MELTAU = 10._r8      !! melting timescale
!
      REAL(r8) :: EVAPR  = 0.3_r8      !! evaporation factor
      REAL(r8) :: REVPDD = 1._r8       !! max rate of DD to evapolation
      REAL(r8) :: RDDR   = 5.e-4_r8    !! DD rate (T0 R0 W0)^-1
      REAL(r8) :: RDDMX  = 0.5_r8      !! norm. flux of downdraft
      REAL(r8) :: VTERM  = 5._r8       !! term. vel. of precip.
      REAL(r8) :: EVATAU = 2._r8    !! evaporation/sublimation timescale
      REAL(r8) :: ZDMIN  = 5.e2_r8     !! min altitude of downdraft detrainment
!
! Note: It is assumed that condensate evaporates in downdraft air.
!
      GPRCP( ISTS:IENS,: ) = 0.D0
      GSNWP( ISTS:IENS,: ) = 0.D0
      GMDD ( ISTS:IENS,: ) = 0.D0
      GTEVP( ISTS:IENS,: ) = 0.D0
      EVAPD( ISTS:IENS,: ) = 0.D0
      SUBLD( ISTS:IENS,: ) = 0.D0
      EVAPE( ISTS:IENS,: ) = 0.D0
      SUBLE( ISTS:IENS,: ) = 0.D0
      EVAPX( ISTS:IENS,: ) = 0.D0
      SUBLX( ISTS:IENS,: ) = 0.D0
      GMDDE( ISTS:IENS,: ) = 0.D0
      SNMLT( ISTS:IENS,: ) = 0.D0
      GCHDD( ISTS:IENS,: ) = 0.D0
      GCWDD( ISTS:IENS,: ) = 0.D0
      GTTEV( ISTS:IENS,: ) = 0.D0
      GTQEV( ISTS:IENS,: ) = 0.D0
      GCHD ( ISTS:IENS ) = 0.D0
      GCWD ( ISTS:IENS ) = 0.D0
      GCUD ( ISTS:IENS ) = 0.D0
      GCVD ( ISTS:IENS ) = 0.D0
      FMELT( ISTS:IENS,: ) = 0.D0
      FEVP ( ISTS:IENS,: ) = 0.D0
!#ifdef OPT_CHASER
!      REVC ( ISTS:IENS,: )  = 0.D0
!#endif
!
      DO K = KTMX, 1, -1   ! loop A
         DO IJ = ISTS, IENS
            DZ   = GDZM( IJ,K+1 ) - GDZM( IJ,K )
            FEVP( IJ,K ) = &
                  ( 1.D0 - TANH( EVATAU*VTERM/DZ ) )
         END DO
!
!     < precipitation melt & freeze >
!
         DO IJ = ISTS, IENS
            GTPRP = GPRCP( IJ,K+1 )+GSNWP( IJ,K+1 )
            IF ( GTPRP .GT. 0.D0 ) THEN
               FSNOW( IJ ) = GSNWP( IJ,K+1 )/GTPRP
            ELSE
               FSNOW( IJ ) = 0.D0
            END IF
            LVIC  = ( EL+EMELT*FSNOW( IJ ) )/CP
            GDTW  = GDT( IJ,K ) &
                  - LVIC*( GDQS( IJ,K ) - GDQ( IJ,K,1 ) ) &
                        /( 1.D0 + LVIC*FDQS( IJ,K ) )
            IF ( GDTW .LT. TWSNOW ) THEN
               GSNWP( IJ,K ) = GSNWP( IJ,K+1 ) &
                             + GPRCI( IJ,K ) + GSNWI( IJ,K )
               GTTEV( IJ,K ) = EMELT/CP*GPRCI( IJ,K ) &
                                       *GRAV/DELP( IJ,K )
               SNMLT( IJ,K ) = -GPRCI( IJ,K )
            ELSE
               DZ   = GDZM( IJ,K+1 ) - GDZM( IJ,K )
               FMELT( IJ,K ) = ( 1.D0 + FTMLT*( GDTW - TWSNOW ) ) &
                     * ( 1.D0 - TANH( GMFLX( IJ,K+1 )/GMFLXC ) ) &
                     * ( 1.D0 - TANH( VTERMS*MELTAU/DZ ) )
               FMELT( IJ,K ) = MAX( MIN( FMELT(IJ,K), 1.D0 ), 0.D0 )
               SNMLT( IJ,K ) = GSNWP( IJ,K+1 )*FMELT( IJ,K )
               GSNWP( IJ,K ) = GSNWP( IJ,K+1 )+GSNWI( IJ,K ) &
                             - SNMLT( IJ,K )
               GPRCP( IJ,K ) = GPRCP( IJ,K+1 )+GPRCI( IJ,K ) &
                             + SNMLT( IJ,K )
               GTTEV( IJ,K ) = -EMELT/CP*SNMLT( IJ,K ) &
                                        *GRAV/DELP( IJ,K )
            END IF
         END DO
!
!     < downdraft >
!
         DO IJ = ISTS, IENS   ! loop B
            IF ( GMDD( IJ,K+1 ) .GT. 0.D0 ) THEN
               GCHX  = GCHD( IJ )/GMDD( IJ,K+1 )
               GCTX  = GDT( IJ,K ) &
                     + ( GCHX-GDHS( IJ,K ) )/( CP+EL*FDQS( IJ,K ) )
               GCQSX = GDQS( IJ,K ) &
                     + FDQS( IJ,K )*( GCTX-GDT( IJ,K ) )
               GCQSX = GCQSX*GMDD( IJ,K+1 )
               EVSU  = MAX( GCQSX - GCWD( IJ ),0.D0 ) * FEVP( IJ,K )
               GTPRP = GPRCP( IJ,K ) + GSNWP( IJ,K )
               IF ( GTPRP .GT. 0.D0 ) THEN
                  FSNOW( IJ ) = GSNWP( IJ,K )/GTPRP
               ELSE
                  FSNOW( IJ ) = 0.D0
               END IF
               EVAPD( IJ,K ) = EVSU*( 1.D0-FSNOW( IJ ) )
               EVAPD( IJ,K ) = MIN( EVAPD( IJ,K ), GPRCP( IJ,K ) )
               SUBLD( IJ,K ) = EVSU*FSNOW( IJ )
               SUBLD( IJ,K ) = MIN( SUBLD( IJ,K ), GSNWP( IJ,K ) )
               GPRCP( IJ,K ) = GPRCP( IJ,K ) - EVAPD( IJ,K )
               GSNWP( IJ,K ) = GSNWP( IJ,K ) - SUBLD( IJ,K )
               GCWD( IJ ) = GCWD( IJ ) + EVAPD( IJ,K ) + SUBLD( IJ,K )
               GCHD( IJ ) = GCHD( IJ ) - EMELT*SUBLD( IJ,K )
            END IF

            GMDD( IJ,K ) = GMDD( IJ,K+1 )
!
            LVIC = ( EL + EMELT*FSNOW( IJ ) )/CP
            DQW  = ( GDQS( IJ,K ) - GDW( IJ,K ) ) &
                 / ( 1.D0 + LVIC*FDQS( IJ,K ) )
            DQW  = MAX( DQW, 0.D0 )
            DTW  = LVIC*DQW
            GDQW = GDW( IJ,K ) + DQW*FEVP( IJ,K )
            DZ   = GDZM( IJ,K+1 ) - GDZM( IJ,K )
!
            EVSU = EVAPR/VTERM*DQW*DZ * FEVP( IJ,K )
            EVAPE( IJ,K ) = EVSU*GPRCP( IJ,K )
            EVAPE( IJ,K ) = MIN( EVAPE( IJ,K ), GPRCP( IJ,K ) )
            SUBLE( IJ,K ) = EVSU*GSNWP( IJ,K )
            SUBLE( IJ,K ) = MIN( SUBLE( IJ,K ), GSNWP( IJ,K ) )
            GTEVP( IJ,K ) = EVAPD( IJ,K ) + SUBLD( IJ,K ) &
                          + EVAPE( IJ,K ) + SUBLE( IJ,K )
!
!#ifdef OPT_CHASER
!            GTPRP = GPRCP( IJ,K+1 )+GSNWP( IJ,K+1 )
!            IF ( GTPRP .GT. 0.D0 ) THEN
!               REVC ( IJ,K )  = GTEVP( IJ,K )/GTPRP
!     $                          *GRAV/DELP( IJ,K )        ! <<CHEM>>
!            END IF
!#endif
!
            GTPRP = GPRCP( IJ,K ) + GSNWP( IJ,K )
            GPRCP( IJ,K ) = GPRCP( IJ,K ) - EVAPE( IJ,K )
            GSNWP( IJ,K ) = GSNWP( IJ,K ) - SUBLE( IJ,K )
!
            GMDDD( IJ ) = 0.D0
            IF ( GDZ( IJ,K )-GDZM( IJ,1 ) .GT. ZDMIN ) THEN
               GTEVE  = EVAPE( IJ,K )+SUBLE( IJ,K )
               GMDDMX = REVPDD*GTEVE/MAX( DQW, 1.D-10 )
               GMDDE( IJ,K ) = RDDR*( DTW*GTPRP*DELP( IJ,K ) )
               GMDDE( IJ,K ) = MAX( MIN( GMDDE(IJ,K), GMDDMX ), 0.D0 )
               GMDDX  = GMDD( IJ,K+1 ) + GMDDE( IJ,K )
               EVSU   = GMDDE( IJ,K )*DQW*FEVP( IJ,K )
               IF ( GTEVE .GT. 0.D0 ) THEN
                  FSNOW( IJ ) = SUBLE( IJ,K )/GTEVE
               ELSE
                  FSNOW( IJ ) = 0.D0
               END IF
               EVAPX( IJ,K ) = ( 1.D0-FSNOW( IJ ) )*EVSU
               SUBLX( IJ,K ) = FSNOW( IJ )*EVSU
!
               IF ( GMDDX .GT. 0.D0 ) THEN
                  GDHI  = GDH( IJ,K ) - EMELT*GDQI( IJ,K )
                  GCHDX = GCHD( IJ ) + GDHI*GMDDE(IJ,K) &
                        - EMELT*SUBLX( IJ,K )
                  GCWDX = GCWD( IJ ) + GDQW*GMDDE(IJ,K)
                  GCSD  = ( GCHDX - EL*GCWDX )/GMDDX
                  IF ( GCSD .LT. GDS( IJ,K ) ) THEN
                     GCHD( IJ   ) = GCHDX
                     GCWD( IJ   ) = GCWDX
                     GCUD( IJ   ) = GCUD( IJ ) &
                                  + GDU( IJ,K )*GMDDE( IJ,K )
                     GCVD( IJ   ) = GCVD( IJ ) &
                                  + GDV( IJ,K )*GMDDE( IJ,K )
                     GMDD( IJ,K ) = GMDDX
                     EVAPE( IJ,K ) = EVAPE( IJ,K ) - EVAPX( IJ,K )
                     SUBLE( IJ,K ) = SUBLE( IJ,K ) - SUBLX( IJ,K )
                     EVAPD( IJ,K ) = EVAPD( IJ,K ) + EVAPX( IJ,K )
                     SUBLD( IJ,K ) = SUBLD( IJ,K ) + SUBLX( IJ,K )
                     GMDDD( IJ )  = 0.D0
                  ELSE
                     GMDDE( IJ,K ) = 0.D0
                     GMDDD( IJ   )  = GMDD( IJ,K+1 )
                  END IF
               END IF
            ELSE
               GMDDD( IJ ) = ( GDZM( IJ,K+1 )-GDZM( IJ,K ) ) &
                            /( GDZM( IJ,K+1 )-GDZM( IJ,1 ) ) &
                            *GMDD( IJ,K+1 )
            END IF
!
            GMDDD( IJ ) = MAX( GMDDD(IJ), GMDD(IJ,K)-RDDMX*GMFLX(IJ,K) )
!
            IF ( GMDDD( IJ ) .GT. 0.D0 ) THEN
               FDET = GMDDD( IJ )/GMDD( IJ,K )
               GCHDD( IJ,K ) = FDET*GCHD( IJ )
               GCWDD( IJ,K ) = FDET*GCWD( IJ )
               GCUDD = FDET*GCUD( IJ )
               GCVDD = FDET*GCVD( IJ )
!
               GTHCI  = GRAV/DELP( IJ,K ) &
                      *( GCHDD( IJ,K ) - GMDDD( IJ )*GDH( IJ,K ) )
               GTQVCI = GRAV/DELP( IJ,K ) &
                      *( GCWDD( IJ,K ) - GMDDD( IJ )*GDQ(IJ,K,1) )
               GTQLCI = -GRAV/DELP( IJ,K ) &
                      * GMDDD( IJ )*GDQ( IJ,K,ITL )
               GTQICI = -GRAV/DELP( IJ,K ) &
                      * GMDDD( IJ )*GDQI( IJ,K )
               GTUCI  = GRAV/DELP( IJ,K ) &
                      *( GCUDD - GMDDD( IJ )*GDU( IJ,K ) )
               GTVCI  = GRAV/DELP( IJ,K ) &
                      *( GCVDD - GMDDD( IJ )*GDV( IJ,K ) )
!
               GTT ( IJ,K ) = GTT( IJ,K ) &
                            + ( GTHCI - EL*GTQVCI )/CP
               GTQ ( IJ,K,1   ) = GTQ( IJ,K,1   ) + GTQVCI
               GTQ ( IJ,K,ITL ) = GTQ( IJ,K,ITL ) + GTQLCI
               GTQI( IJ,K ) = GTQI( IJ,K ) + GTQICI
               GTU ( IJ,K ) = GTU( IJ,K ) + GTUCI
               GTV ( IJ,K ) = GTV( IJ,K ) + GTVCI
!
               GCHD( IJ   ) = GCHD( IJ   ) - GCHDD( IJ,K )
               GCWD( IJ   ) = GCWD( IJ   ) - GCWDD( IJ,K )
               GCUD( IJ   ) = GCUD( IJ   ) - GCUDD
               GCVD( IJ   ) = GCVD( IJ   ) - GCVDD
               GMDD( IJ,K ) = GMDD( IJ,K ) - GMDDD( IJ )
            END IF
         END DO   ! loop B
!
      END DO   ! loop A
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            GTTEV( IJ,K ) = GTTEV( IJ,K ) &
              - ( EL*EVAPE( IJ,K )+( EL+EMELT )*SUBLE( IJ,K ) ) &
                *GRAV/DELP( IJ,K )/CP
            GTT( IJ,K ) = GTT( IJ,K ) + GTTEV( IJ,K )
!
            GTQEV( IJ,K ) = GTQEV( IJ,K ) &
              + ( EVAPE( IJ,K )+SUBLE( IJ,K ) )*GRAV/DELP( IJ,K )
            GTQ( IJ,K,1 ) = GTQ( IJ,K,1 ) + GTQEV( IJ,K )
!
            GMFLX( IJ,K ) = GMFLX( IJ,K ) - GMDD( IJ,K )
         END DO
      END DO
!
!      CALL HISTIN( FMELT, 'FMELT', 'melting rate in cumdown',
!     &             '   ', 'ALEV', HCLAS )
!      CALL HISTIN( FEVP, 'FEVP', 'evap/subl factor',
!     &             '   ', 'ALEV', HCLAS )
!      CALL HISTIN( EVAPD, 'EVAPD', 'cum. rain evap. into DD',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( SUBLD, 'SUBLD', 'cum. snow subl. into DD',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( EVAPE, 'EVAPE', 'cum. rain evap. into env.',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( SUBLE, 'SUBLE', 'cum. snow subl. into env.',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( EVAPX, 'EVAPX', 'cum. rain evap. from env. into DD.',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!      CALL HISTIN( SUBLX, 'SUBLX', 'cum. snow subl. from env. into DD.',
!     &             'kg/m**2/s', 'ALEV', HCLAS )
!
#ifdef OPT_CUMBGT   /* budget check */
      WBGT( ISTS:IENS ) = 0.D0
      HBGT( ISTS:IENS ) = 0.D0
      DDWBGT( ISTS:IENS ) = 0.D0
      DDHBGT( ISTS:IENS ) = 0.D0
      WMX = 0.D0
      HMX = 0.D0
      DDWMX = 0.D0
      DDHMX = 0.D0
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            WBGT( IJ ) = WBGT( IJ ) + GPRCI( IJ,K ) + GSNWI( IJ,K ) &
                       - EVAPD( IJ,K ) - SUBLD( IJ,K ) &
                       - GTQEV( IJ,K )*DELP( IJ,K )/GRAV
            HBGT( IJ ) = HBGT( IJ ) &
                       + EL*EVAPE( IJ,K ) + EMELT*SNMLT( IJ,K ) &
                       + ( EL+EMELT )*SUBLE( IJ,K ) &
                       + CP*GTTEV( IJ,K )*DELP( IJ,K )/GRAV
            DDWBGT( IJ ) = DDWBGT( IJ ) &
                         + EVAPD( IJ,K ) + SUBLD( IJ,K ) &
                         + GDW( IJ,K )*GMDDE( IJ,K ) &
                         - GCWDD( IJ,K )
            DDHBGT( IJ ) = DDHBGT( IJ ) &
                         + ( GDH(IJ,K)-EMELT*GDQI(IJ,K) )*GMDDE(IJ,K) &
                         - EMELT*SUBLD( IJ,K ) &
                         - GCHDD( IJ,K )
         END DO
      END DO
      DO IJ = ISTS, IENS
         WBGT( IJ ) = WBGT( IJ ) - GPRCP( IJ,1 ) - GSNWP( IJ,1 )
         IF ( ABS( WBGT(IJ) ) .GT. ABS( WMX ) ) WMX   = WBGT( IJ )
         IF ( ABS( HBGT(IJ) ) .GT. ABS( HMX ) ) HMX   = HBGT( IJ )
         IF ( ABS(DDWBGT(IJ)) .GT. ABS(DDWMX) ) DDWMX = DDWBGT(IJ)
         IF ( ABS(DDHBGT(IJ)) .GT. ABS(DDHMX) ) DDHMX = DDHBGT(IJ)
      END DO
!
      WRITE( iulog,* ) &
         '### CUMDWN(rank=',irank,'): water imbalance =', WMX
      WRITE( iulog,* ) &
         '### CUMDWN(rank=',irank,'): energy imbalance =', HMX
      WRITE( iulog,* ) &
         '### CUMDWN(rank=',irank,'): downdraft water imbalance =', DDWMX
      WRITE( iulog,* ) &
         '### CUMDWN(rank=',irank,'): downdraft energy imbalance =', DDHMX
#endif
!
      END SUBROUTINE CUMDWN
!***********************************************************************
      SUBROUTINE CUMCLD   & !! cloudiness
               ( CUMCLW, QLIQ  , QICE  , FLIQC  ,   & ! modified
                 CUMFRC,                            & ! output
!#ifdef OPT_CHASER
!     M           LEVCUM, LNFRC ,           ! <<CHEM>>
!     I           TOPFLX,                   ! <<CHEM>>
!#endif
                 GMFLX , KTMX  , ISTS, IENS   )       ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!
!   [OUTPUT]
      REAL(r8)     CUMFRC( IJSDIM )            !! cumulus cloud fraction
!
!   [MODIFY]
      REAL(r8)     CUMCLW( IJSDIM, KMAX   )    !! cloud water in cumulus
      REAL(r8)     QLIQ  ( IJSDIM, KMAX   )    !! cloud liquid
      REAL(r8)     QICE  ( IJSDIM, KMAX   )    !! cloud ice
      REAL(r8)     FLIQC ( IJSDIM, KMAX   )    !! liquid ratio in cumulus
!#ifdef OPT_CHASER
!      INTEGER      LEVCUM( IJSDIM, KMAX )      !!
!      REAL(r8)     LNFRC ( IJSDIM, KMAX )      !! fraction of each cloud
!#endif
!
!   [INPUT]
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 ) !! cumulus mass flux
      INTEGER      KTMX
!#ifdef OPT_CHASER
!      REAL(r8)     TOPFLX( IJSDIM, NCTP   ) !! mass flux at each cloud top
!#endif
      INTEGER      ISTS, IENS
!
!   [WORK]
      INTEGER      IJ, K
      REAL(r8)     CUMF
      REAL(r8)     QC
!#ifdef OPT_CHASER
!      REAL(r8)     SNFRC( ISTS:IENS   )
!      REAL(r8)     LNF
!#endif
      LOGICAL, SAVE :: OFIRST = .TRUE.
!
!   [INTERNAL PARAM]
      REAL(r8) :: FACLW  = 0.1_r8     !! Mc->CLW
      REAL(r8) :: CMFMIN = 2.e-3_r8   !! Mc->cloudiness
      REAL(r8) :: CMFMAX = 3.e-1_r8   !! Mc->cloudiness
      REAL(r8) :: CLMIN  = 1.e-3_r8   !! cloudiness Min.
      REAL(r8) :: CLMAX  = 0.1_r8     !! cloudiness Max.
      REAL(r8), SAVE :: FACLF
!
      IF ( OFIRST ) THEN
         FACLF = (CLMAX-CLMIN)/LOG(CMFMAX/CMFMIN)
         OFIRST = .FALSE.
      END IF
!
      CUMFRC( ISTS:IENS ) = 0.D0
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            CUMFRC( IJ ) = MAX( CUMFRC( IJ ), GMFLX( IJ,K ) )
         END DO
      END DO
      DO IJ = ISTS, IENS
         IF ( CUMFRC( IJ ) .GT. 0.D0 ) THEN
            CUMF         = LOG( MAX( CUMFRC( IJ ),CMFMIN )/CMFMIN )
            CUMFRC( IJ ) = MIN( FACLF*CUMF+CLMIN, CLMAX )
         END IF
      END DO
!
      DO K = 1, KTMX
         DO IJ = ISTS, IENS
            IF ( GMFLX( IJ,K ) .GT. 0.D0 ) THEN
               QLIQ  ( IJ,K ) = FACLW*QLIQ  ( IJ,K )/GMFLX( IJ,K ) &
                               *CUMFRC( IJ )
               QICE  ( IJ,K ) = FACLW*QICE  ( IJ,K )/GMFLX( IJ,K ) &
                               *CUMFRC( IJ )
               CUMCLW( IJ,K ) = FACLW*CUMCLW( IJ,K )/GMFLX( IJ,K ) &
                               *CUMFRC( IJ )
               QC = QLIQ( IJ,K ) + QICE( IJ,K )
               IF ( QC .GT. 0.D0 ) THEN
                  FLIQC( IJ,K ) = QLIQ( IJ,K ) / QC
               END IF
            END IF
         END DO
      END DO
!
!#ifdef OPT_CHASER
! =  <<CHEM>> ========================================================
!
!      DO 3100 K  = 1    , KTMX
!      DO 3100 IJ = ISTS , IENS
!         IF ( TOPFLX( IJ,K ) .GT. CMFMIN*1.D-2 ) THEN
!            LEVCUM( IJ,K ) = 1
!         END IF
!
!         IF ( TOPFLX( IJ,K ) .GT. 0. ) THEN
!            LNF           = LOG( MAX(TOPFLX( IJ,K ),CMFMIN)/CMFMIN )
!            LNFRC( IJ,K ) = MIN( FACLF*LNF+CLMIN, CLMAX )
!         END IF
! 3100 CONTINUE
!
!      DO IJ = ISTS, IENS
!         SNFRC( IJ ) = 0.D0
!      END DO
!
!      DO K  = 1    , KTMX
!      DO IJ = ISTS , IENS
!         SNFRC( IJ ) = SNFRC( IJ ) + LNFRC( IJ,K )
!      END DO
!      END DO
!
!      DO K  = 1    , KTMX
!      DO IJ = ISTS , IENS
!         IF ( SNFRC( IJ ) .GT. 0.D0 ) THEN
!            LNFRC( IJ,K ) = LNFRC( IJ,K )*CUMFRC( IJ )/ SNFRC( IJ )
!         END IF
!      END DO
!      END DO
!#endif /* OPT_CHASER */
!
      END SUBROUTINE CUMCLD
!***********************************************************************
      SUBROUTINE CUMUPR   & !! Tracer Updraft
               ( GTR   , GPRCC ,                           & ! modified
                 GDR   , CBMFX , ELAM  , GDZ   , GDZM  ,   & ! input
                 GCYM  , GCYT  , GCQT  , GCLT  , GCIT  ,   & ! input
                 GTPRT , GTEVP , GTPRC0,                   & ! input
                 KB    , KBMX  , KT    , KTMX  , KTMXT ,   & ! input
                 DELP  , OTSPT , ISTS  , IENS            )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTR   ( IJSDIM, KMAX, NTR )
      REAL(r8)     GPRCC ( IJSDIM, NTR  )
!
!   [INPUT]
      REAL(r8)     GDR   ( IJSDIM, KMAX, NTR )
      REAL(r8)     CBMFX ( IJSDIM, NCTP   )
      REAL(r8)     ELAM  ( IJSDIM, KMAX, NCTP )
      REAL(r8)     GDZ   ( IJSDIM, KMAX   )
      REAL(r8)     GDZM  ( IJSDIM, KMAX+1 )
      REAL(r8)     GCYM  ( IJSDIM, KMAX   )
      REAL(r8)     GCYT  ( IJSDIM, NCTP   )
      REAL(r8)     GCQT  ( IJSDIM, NCTP   )
      REAL(r8)     GCLT  ( IJSDIM, NCTP   )
      REAL(r8)     GCIT  ( IJSDIM, NCTP   )
      REAL(r8)     GTPRT ( IJSDIM, NCTP   )
      REAL(r8)     GTEVP ( IJSDIM, KMAX   )
      REAL(r8)     GTPRC0( IJSDIM         )   !! precip. before evap.
      INTEGER      KB    ( IJSDIM )
      INTEGER      KBMX
      INTEGER      KT    ( IJSDIM, NCTP   )
      INTEGER      KTMX  ( NCTP           )
      INTEGER      KTMXT
      REAL(r8)     DELP  ( IJSDIM, KMAX   )
      LOGICAL      OTSPT ( NTR )              !! transport with this routine?
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      IJ, K, LT, TP, CTP
      REAL(r8)     GCRTD
      REAL(r8)     SCAV
      REAL(r8)     GCWT, GPRCR
      REAL(r8)     GCRB  ( ISTS:IENS )
      REAL(r8)     GCRT  ( ISTS:IENS )
      REAL(r8)     DR    ( ISTS:IENS )
      REAL(r8)     DGCB  ( ISTS:IENS, KMAX )
      REAL(r8)     DZ    ( ISTS:IENS, KMAX )
      REAL(r8)     DZT   ( ISTS:IENS, NCTP )
      REAL(r8)     RGCWT ( ISTS:IENS, NCTP )
      REAL(r8)     RDZM  ( ISTS:IENS, KMAX )
      REAL(r8)     EVPF  ( ISTS:IENS, KMAX )
      REAL(r8)     MASK1 ( ISTS:IENS, NCTP )
      REAL(r8)     MASK2 ( ISTS:IENS, NCTP )
!
!   [INTERNAL PARM]
      REAL(r8) :: FSCAV ( NTR ) = 0._r8
      REAL(r8) :: FSWTR ( NTR ) = 0._r8
!
      DO K = 1, KBMX
         DO IJ = ISTS, IENS
            DGCB( IJ,K ) = GCYM( IJ,K+1 ) - GCYM( IJ,K )
         END DO
      END DO
      DO K = 1, KTMXT
         DO IJ = ISTS, IENS
            DZ  ( IJ,K ) = GDZM( IJ,K+1 ) - GDZM( IJ,K )
            RDZM( IJ,K ) = GRAV / DELP( IJ,K )
            EVPF( IJ,K ) = 0.D0
            IF ( GTPRC0( IJ ) .GT. 0.D0 ) THEN
               EVPF( IJ,K ) = GTEVP( IJ,K ) / GTPRC0( IJ )
            END IF
         END DO
      END DO
      DO CTP = 1, NCTP
         DO IJ = ISTS, IENS
            K = KT( IJ, CTP )
!
            GCWT = GCQT( IJ,CTP ) + GCLT( IJ,CTP ) + GCIT( IJ,CTP )
            RGCWT( IJ,CTP ) = 0.D0
            IF ( GCWT .GT. 0.D0 ) THEN
               RGCWT( IJ,CTP ) = 1.D0 / GCWT
            END IF
!
            MASK1( IJ,CTP ) = 0.D0
            DZT  ( IJ,CTP ) = 0.D0
            IF ( K .GT. KB( IJ ) ) THEN
               MASK1( IJ,CTP ) = 1.D0
               DZT  ( IJ,CTP ) = GDZ( IJ,K ) - GDZM( IJ,K )
            END IF
            MASK2( IJ,CTP ) = 0.D0
            IF ( CBMFX( IJ,CTP ) .GT. 0.D0 ) then
               MASK2( IJ,CTP ) = 1.D0
            END IF
         END DO
      END DO
!
      DO LT = 1, NTR   ! outermost loop
!
         IF ( OTSPT( LT ) ) THEN
            GCRB = 0.D0
            DO K = 1, KBMX
               DO IJ = ISTS, IENS
                  IF ( K .LT. KB( IJ ) ) THEN
                     GCRB( IJ ) = GCRB( IJ ) &
                                + DGCB( IJ,K ) * GDR( IJ,K,LT )
                  END IF
               END DO
            END DO
!
            DO CTP = 1, NCTP
               DR = 0.D0
               DO K = 2, KTMX( CTP )
                  DO IJ = ISTS, IENS
                     IF ( K .GE. KB( IJ     ) .AND. &
                          K .LT. KT( IJ,CTP ) ) THEN
                        DR( IJ ) = DR( IJ ) &
                                 + DZ( IJ,K ) * ELAM( IJ,K,CTP ) &
                                              * GDR ( IJ,K,LT  )
                     END IF
                  END DO
               END DO
!
               DO IJ = ISTS, IENS
                  K = MAX( KT( IJ,CTP ), 1 )
                  DR( IJ ) = DR( IJ ) &
                           + DZT( IJ,CTP ) * ELAM( IJ,K,CTP ) &
                                           * GDR ( IJ,K,LT  ) &
                                           * MASK1( IJ,CTP )
                  GCRT( IJ ) = ( GCRB( IJ ) + DR( IJ ) ) * MASK1( IJ,CTP )
!
                  SCAV = FSCAV( LT ) * GTPRT( IJ,CTP ) &
                       + FSWTR( LT ) &
                       * GTPRT( IJ,CTP ) * RGCWT( IJ,CTP )
                  SCAV  = MIN( SCAV, 1.D0 )
                  GCRTD = GCRT( IJ ) * ( 1.D0 - SCAV )
                  GPRCR = SCAV * GCRT( IJ ) * CBMFX( IJ,CTP )

                  GTR( IJ,K,LT ) = GTR( IJ,K,LT ) &
                       + RDZM( IJ,K ) * CBMFX( IJ,CTP ) &
                       * ( GCRTD - GCYT( IJ,CTP ) * GDR( IJ,K,LT ) ) &
                       * MASK2( IJ,CTP )
                  GPRCC( IJ,LT ) = GPRCC( IJ,LT ) + GPRCR &
                                 * MASK2( IJ,CTP )
               END DO
            END DO
!
            DO K = KTMXT, 1, -1
               DO IJ = ISTS, IENS
                  GTR( IJ,K,LT ) = GTR( IJ,K,LT ) &
                                 + RDZM( IJ,K ) &
                                 * GPRCC( IJ,LT ) * EVPF( IJ,K )
                  GPRCC( IJ,LT ) = GPRCC( IJ,LT ) &
                                 * ( 1.D0 - EVPF( IJ,K ) )
               END DO
            END DO
!
         END IF
!
      END DO   ! outermost loop
!
      END SUBROUTINE CUMUPR
!***********************************************************************
      SUBROUTINE CUMDNR   & !! Tracer Downdraft
               ( GTR   ,                         & ! modified
                 GDR   , GMDD  , DELP  ,         & ! input
                 KTMX  , OTSPT , ISTS  , IENS )    ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTR   ( IJSDIM, KMAX, NTR )   !! Temperature tendency
!
!   [INPUT]
      REAL(r8)     GDR   ( IJSDIM, KMAX, NTR )
      REAL(r8)     GMDD  ( IJSDIM, KMAX+1 )      !! downdraft mass flux
      REAL(r8)     DELP  ( IJSDIM, KMAX   )
      INTEGER      KTMX
      LOGICAL      OTSPT ( NTR )
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GCRD  ( ISTS:IENS )           !! downdraft q
      REAL(r8)     GMDDE, GMDDD, GCRDD
      INTEGER      IJ, K, LT
!
!
      DO LT = 1, NTR
!
         IF ( OTSPT( LT ) ) THEN
            GCRD = 0.D0
            DO K = KTMX, 1, -1
               DO IJ = ISTS, IENS
                  GMDDE = GMDD( IJ,K ) - GMDD( IJ,K+1 )
                  IF ( GMDDE .GE. 0.D0 ) THEN
                     GCRD( IJ ) = GCRD( IJ ) + GDR( IJ,K,LT )*GMDDE
                  ELSE IF ( GMDD( IJ,K+1 ) .GT. 0.D0 ) THEN
                     GMDDD = - GMDDE
                     GCRDD = GMDDD/GMDD( IJ,K+1 ) * GCRD( IJ )
                     GTR( IJ,K,LT ) = GTR( IJ,K,LT ) &
                                    + GRAV/DELP( IJ,K ) &
                                     *( GCRDD - GMDDD*GDR( IJ,K,LT ) )
                     GCRD( IJ ) = GCRD( IJ ) - GCRDD
                  END IF
               END DO
            END DO
!
         END IF
!
      END DO
!
      END SUBROUTINE CUMDNR
!***********************************************************************
      SUBROUTINE CUMSBR   & !! Tracer Subsidence
               ( GTR   ,                              & ! modified
                 GDR   , DELP  ,                      & ! input
                 GMFLX , KTMX  , OTSPT , ISTS, IENS )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTR   ( IJSDIM, KMAX, NTR )   !! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   ( IJSDIM, KMAX, NTR )   !! tracer
      REAL(r8)     DELP  ( IJSDIM, KMAX   )
      REAL(r8)     GMFLX ( IJSDIM, KMAX+1 )      !! mass flux
      INTEGER      KTMX
      LOGICAL      OTSPT ( NTR )                 !! tracer transport on/off
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER      IJ, K, KM, KP, LT
      REAL(r8)     SBR0, SBR1, FX1
      REAL(r8)     FX    ( ISTS:IENS )
!
      DO LT = 1, NTR
         IF ( OTSPT( LT ) ) THEN
            FX( ISTS:IENS ) = 0.D0
            DO K = KTMX, 1, -1
               KM = MAX( K-1, 1    )
               KP = MIN( K+1, KMAX )
               DO IJ = ISTS, IENS
                  SBR0 = GMFLX( IJ,K+1 ) &
                         *( GDR( IJ,KP,LT )-GDR( IJ,K ,LT ) )
                  SBR1 = GMFLX( IJ,K   ) &
                         *( GDR( IJ,K ,LT )-GDR( IJ,KM,LT ) )
!
                  IF ( GMFLX( IJ,K ) .GT. GMFLX( IJ,K+1 ) ) THEN
                     FX1 = 0.5D0
                  ELSE
                     FX1 = 0.D0
                  END IF
!
                  GTR( IJ,K,LT ) &
                        = GTR( IJ,K,LT ) &
                        + GRAV/DELP( IJ,K ) &
                         *( ( 1.D0-FX( IJ ) )*SBR0 + FX1 *SBR1 )
!
                  FX( IJ ) = FX1
               END DO
            END DO
!
         END IF
      END DO
!
      END SUBROUTINE CUMSBR
!*********************************************************************
      SUBROUTINE CUMFXR   & !! Tracer mass fixer
               ( GTR   ,                                   & ! modified
                 GDR   , DELP  , DELTA , KTMX  , IMFXR ,   & ! input
                 ISTS  , IENS                            )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTR   ( IJSDIM, KMAX, NTR )   !! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   ( IJSDIM, KMAX, NTR )   !! tracer
      REAL(r8)     DELP  ( IJSDIM, KMAX      )
      REAL(r8)     DELTA                         !! time step
      INTEGER      KTMX
      INTEGER      IMFXR ( NTR )
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GDR1
      REAL(r8)     GDR2  ( ISTS:IENS, KMAX )
      REAL(r8)     TOT0  ( ISTS:IENS )
      REAL(r8)     TOT1  ( ISTS:IENS )
      REAL(r8)     TRAT  ( ISTS:IENS )
      REAL(r8)     FWAT
      INTEGER      IJ, K, LT
#ifdef OPT_CUMBGT
      REAL(r8)     GTWB  ( ISTS:IENS )           !! Qt tendency (before)
      REAL(r8)     GTWA  ( ISTS:IENS )           !! Qt tendency (after)
      REAL(r8)     WBGT, WBMX
#endif
!
! Attention: tracers are forced to be positive unless IMFXR=0.
!
#ifdef OPT_CUMBGT
      GTWB( ISTS:IENS ) = 0.D0
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            GTWB( IJ ) = GTWB( IJ ) &
                       + ( GTR( IJ,K,1 )+GTR( IJ,K,ITL ) &
                         + GTR( IJ,K,ITI ) )*DELP( IJ,K )/GRAV
         END DO
      END DO
#endif
!
      DO LT = 1, NTR
         SELECT CASE ( IMFXR( LT ) )
            CASE (0)
               CYCLE
            CASE (1)
               FWAT = 1.D0
            CASE (2)
               FWAT = 0.D0
            CASE DEFAULT
               EXIT
         END SELECT
!
         TOT0( ISTS:IENS ) = 0.D0
         TOT1( ISTS:IENS ) = 0.D0
!
         DO K = KTMX, 1, -1
            DO IJ = ISTS, IENS
               IF ( GTR( IJ,K,LT ) .NE. 0.D0 ) THEN
                  GDR1         = GDR( IJ,K,LT ) &
                               + DELTA*GTR( IJ,K,LT )
                  GDR2( IJ,K ) = MAX( GDR1, 0.D0 )
                  GDR1         = GDR1          *         FWAT &
                               + GDR( IJ,K,LT )*( 1.D0 - FWAT )
                  TOT0( IJ )   = TOT0( IJ ) &
                               + GDR1        *(DELP( IJ,K )/GRAV)
                  TOT1( IJ )   = TOT1( IJ ) &
                               + GDR2( IJ,K )*(DELP( IJ,K )/GRAV)
               END IF
            END DO
         END DO
!
         DO IJ = ISTS, IENS
            IF ( TOT1( IJ ) .GT. 0.D0 ) THEN
               TRAT( IJ ) = MAX( TOT0( IJ ), 0.D0 )/TOT1( IJ )
            ELSE
               TRAT( IJ ) = 1.D0
            END IF
         END DO
!
         DO K = KTMX, 1, -1
            DO IJ = ISTS, IENS
               IF ( GTR( IJ,K,LT ) .NE. 0.D0 ) THEN
                  GDR2( IJ,K    ) = GDR2( IJ,K )*TRAT( IJ )
                  GTR ( IJ,K,LT ) = ( GDR2( IJ,K )-GDR( IJ,K,LT ) ) &
                                  / DELTA
               END IF
            END DO
         END DO
!
      END DO   ! LT-loop
!
#ifdef OPT_CUMBGT   /* budget check */
      GTWA( ISTS:IENS ) = 0.D0
      WBMX = 0.D0
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            GTWA( IJ ) = GTWA( IJ ) &
                       + ( GTR( IJ,K,1 )+GTR( IJ,K,ITL ) &
                         + GTR( IJ,K,ITI ) )*DELP( IJ,K )/GRAV
         END DO
      END DO
      DO IJ = ISTS, IENS
         WBGT = GTWA( IJ )-GTWB( IJ )
         IF ( ABS( WBGT ) .GT. ABS( WBMX ) ) WBMX = WBGT
      END DO
      WRITE( iulog,* ) &
         '### CUMFXR(rank=',irank,'): water imbalance =', WBMX
#endif
!
      END SUBROUTINE CUMFXR
!*********************************************************************
      SUBROUTINE CUMFXR1   & !! Tracer mass fixer
               ( GTR   ,                                   & ! modified
                 GDR   , DELP  , DELTA , KTMX  , IMFXR ,   & ! input
                 ISTS  , IENS                            )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!
!   [MODIFY]
      REAL(r8)     GTR   ( IJSDIM, KMAX )   !! tracer tendency
!
!   [INPUT]
      REAL(r8)     GDR   ( IJSDIM, KMAX )   !! tracer
      REAL(r8)     DELP  ( IJSDIM, KMAX )
      REAL(r8)     DELTA                         !! time step
      INTEGER      KTMX
      INTEGER      IMFXR
        ! 0: mass fixer is not applied
        !    tracers which may become negative values
        !    e.g. subgrid-PDFs
        ! 1: mass fixer is applied, total mass may change through cumulus scheme
        !    e.g. moisture, liquid cloud, ice cloud, aerosols
        ! 2: mass fixer is applied, total mass never change through cumulus scheme
        !    e.g. CO2
      INTEGER      ISTS, IENS
!
!   [INTERNAL WORK]
      REAL(r8)     GDR1
      REAL(r8)     GDR2  ( ISTS:IENS, KMAX )
      REAL(r8)     TOT0  ( ISTS:IENS )
      REAL(r8)     TOT1  ( ISTS:IENS )
      REAL(r8)     TRAT  ( ISTS:IENS )
      REAL(r8)     FWAT
      INTEGER      IJ, K
!
! Attention: tracers are forced to be positive unless IMFXR=0.
!
      SELECT CASE ( IMFXR )
         CASE (0)
            RETURN
         CASE (1)
            FWAT = 1.D0
         CASE (2)
            FWAT = 0.D0
         CASE DEFAULT
            RETURN
      END SELECT
!
      TOT0( ISTS:IENS ) = 0.D0
      TOT1( ISTS:IENS ) = 0.D0
!
      DO K = KTMX, 1, -1
         DO IJ = ISTS, IENS
            IF ( GTR( IJ,K ) .NE. 0.D0 ) THEN
               GDR1         = GDR( IJ,K ) &
                            + DELTA*GTR( IJ,K )
               GDR2( IJ,K ) = MAX( GDR1, 0.D0 )
               GDR1         = GDR1*FWAT &
                            + GDR( IJ,K )*( 1.D0 - FWAT )
               TOT0( IJ )   = TOT0( IJ ) &
                            + GDR1        *(DELP( IJ,K )/GRAV)
               TOT1( IJ )   = TOT1( IJ ) &
                            + GDR2( IJ,K )*(DELP( IJ,K )/GRAV)
            END IF
         END DO
      END DO
!
      DO IJ = ISTS, IENS
         IF ( TOT1( IJ ) .GT. 0.D0 ) THEN
            TRAT( IJ ) = MAX( TOT0( IJ ), 0.D0 )/TOT1( IJ )
         ELSE
            TRAT( IJ ) = 1.D0
         END IF
      END DO
!
      DO K = KTMX, 1, -1
         DO IJ = ISTS, IENS
            IF ( GTR( IJ,K ) .NE. 0.D0 ) THEN
               GDR2( IJ,K ) = GDR2( IJ,K )*TRAT( IJ )
               GTR ( IJ,K ) = ( GDR2( IJ,K )-GDR( IJ,K ) ) &
                            / DELTA
            END IF
         END DO
      END DO
!
      END SUBROUTINE CUMFXR1
!*********************************************************************
      SUBROUTINE CUMCHK   & !! check range of output values
               ( GTT   , GTQ   , GTU   , GTV   ,   & ! input
                 GTCFRC, GPRCC , GSNWC , CUMCLW,   & ! input
                 CUMFRC, FLIQC , GTPRP ,           & ! input
                 ISTS  , IENS                    )   ! input
!
      use ppgrid      , only: IJSDIM => pcols, KMAX => pver
      use constituents, only: NTR => pcnst
!
      IMPLICIT NONE
!
!   [INPUT]
      REAL(r8)     GTT   ( IJSDIM, KMAX      ) !! heating rate
      REAL(r8)     GTQ   ( IJSDIM, KMAX, NTR ) !! change in q
      REAL(r8)     GTU   ( IJSDIM, KMAX      ) !! tendency of u
      REAL(r8)     GTV   ( IJSDIM, KMAX      ) !! tendency of v
      REAL(r8)     GPRCC ( IJSDIM, NTR       ) !! rainfall
      REAL(r8)     GSNWC ( IJSDIM            ) !! snowfall
      REAL(r8)     CUMCLW( IJSDIM, KMAX      ) !! cloud water in cumulus
      REAL(r8)     CUMFRC( IJSDIM            ) !! cumulus cloud fraction
      REAL(r8)     GTCFRC( IJSDIM, KMAX      ) !! change in cloud fraction
      REAL(r8)     FLIQC ( IJSDIM, KMAX      ) !! liquid ratio in cumulus
      REAL(r8)     GTPRP ( IJSDIM, KMAX+1    ) !! rain+snow flux
!
      INTEGER    ISTS, IENS
!
!   [INTERNAL WORK]
      INTEGER    IJ, K
!
!   [INTERNAL PARM]
      REAL(r8) :: GTTMAX  = 1.e-2_r8
      REAL(r8) :: GTQVMAX = 1.e-4_r8
      REAL(r8) :: GTQLMAX = 1.e-5_r8
      REAL(r8) :: GTUMAX  = 1.e-2_r8
      REAL(r8) :: GTVMAX  = 1.e-2_r8
      REAL(r8) :: GTCFMAX = 1.e-3_r8
      REAL(r8) :: PRCCMAX = 1.e-2_r8
      REAL(r8) :: SNWCMAX = 1.e-2_r8
      REAL(r8) :: CLWMAX  = 1.e-3_r8
      REAL(r8) :: TPRPMAX = 1.e-2_r8
      REAL(r8) :: GTQIMAX = 1.e-5_r8
      REAL(r8) :: GTM2MAX = 1._r8
      REAL(r8) :: GTM3MAX = 1._r8
!
      DO K = 1, KMAX
         DO IJ = ISTS, IENS
            IF ( ABS( GTT( IJ,K ) ) .GT. GTTMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTT(',IJ,',',K,')=',GTT( IJ,K )
            END IF
            IF ( ABS( GTQ( IJ,K,1 ) ) .GT. GTQVMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTQ(',IJ,',',K,',1 )=', &
                          GTQ( IJ,K,1 )
            END IF
            IF ( ABS( GTQ( IJ,K,ITL ) ) .GT. GTQLMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTQ(',IJ,',',K,',ITL )=', &
                          GTQ( IJ,K,ITL )
            END IF
            IF ( ABS( GTU( IJ,K ) ) .GT. GTUMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTU(',IJ,',',K,')=',GTU( IJ,K )
            END IF
            IF ( ABS( GTV( IJ,K ) ) .GT. GTVMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTV(',IJ,',',K,')=',GTV( IJ,K )
            END IF
            IF ( ABS( GTCFRC( IJ,K ) ) .GT. GTCFMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTCFRC(',IJ,',',K,')=', &
                          GTCFRC( IJ,K )
            END IF
            IF ( CUMCLW( IJ,K ) .GT. CLWMAX .OR. &
                 CUMCLW( IJ,K ) .LT. 0.D0       ) THEN
               WRITE(iulog,*) '### CUMCHK: CUMCLW(',IJ,',',K,')=', &
                          CUMCLW( IJ,K )
            END IF
            IF ( FLIQC( IJ,K ) .GT. 1.D0 .OR. &
                 FLIQC( IJ,K ) .LT. 0.D0      ) THEN
               WRITE(iulog,*) '### CUMCHK: FLIQC(',IJ,',',K,')=', &
                          FLIQC( IJ,K )
            END IF
            IF ( GTPRP( IJ,K ) .GT. TPRPMAX .OR. &
                 GTPRP( IJ,K ) .LT. 0.D0       ) THEN
               WRITE(iulog,*) '### CUMCHK: GTPRP(',IJ,',',K,')=', &
                          GTPRP( IJ,K )
            END IF
            IF ( ABS( GTQ( IJ,K,ITI ) ) .GT. GTQIMAX ) THEN
               WRITE(iulog,*) '### CUMCHK: GTQ(',IJ,',',K,',ITI )=', &
                          GTQ( IJ,K,ITI )
            END IF
!            IF ( ABS( GTQ( IJ,K,IMU2 ) ) .GT. GTM2MAX ) THEN
!               WRITE(iulog,*) '### CUMCHK: GTQ(',IJ,',',K,',IMU2 )=', &
!                          GTQ( IJ,K,IMU2 )
!            END IF
!            IF ( ABS( GTQ( IJ,K,IMU3 ) ) .GT. GTM3MAX ) THEN
!               WRITE(iulog,*) '### CUMCHK: GTQ(',IJ,',',K,',IMU3 )=', &
!                          GTQ( IJ,K,IMU3 )
!            END IF
         END DO
      END DO
!
      DO IJ = ISTS, IENS
         IF ( GPRCC( IJ,1 ) .GT. PRCCMAX .OR. &
              GPRCC( IJ,1 ) .LT. 0.D0       ) THEN
            WRITE(iulog,*) '### CUMCHK: GPRCC(',IJ,')=',GPRCC( IJ,1 )
         END IF
         IF ( GSNWC( IJ ) .GT. SNWCMAX .OR. &
              GSNWC( IJ ) .LT. 0.D0       ) THEN
            WRITE(iulog,*) '### CUMCHK: GSNWC(',IJ,')=',GSNWC( IJ )
         END IF
         IF ( CUMFRC( IJ ) .GT. 1.D0 .OR. &
              CUMFRC( IJ ) .LT. 0.D0      ) THEN
            WRITE(iulog,*) '### CUMCHK: CUMFRC(',IJ,')=',CUMFRC( IJ )
         END IF
      END DO
!
      END SUBROUTINE CUMCHK
!***********************************************************************
      SUBROUTINE TINTP & !! vertical interpolation of temperature
               ( GDTM  ,                  & ! output
                 GDT   , GDP   , GDPM ,   & ! input
                 ISTS  , IENS           )   ! input

      use ppgrid, only: IJSDIM => pcols, KMAX => pver
!
      IMPLICIT NONE
!*
!*   [OUTPUT]
      REAL(r8)     GDTM  ( IJSDIM, KMAX+1 )     !! temperature (half lev)
!*
!*   [INPUT]
      REAL(r8)     GDT   ( IJSDIM, KMAX )       !! temperature (full lev)
      REAL(r8)     GDP   ( IJSDIM, KMAX )       !! pressure (full lev)
      REAL(r8)     GDPM  ( IJSDIM, KMAX+1 )     !! pressure (half lev)
      INTEGER      ISTS, IENS                   !! range of active grids
!*
!*   [INTERNAL WORK]
      REAL(r8)     FTINT ( KMAX )               !! intrp. coef.
      REAL(r8)     FTINTM( KMAX )               !! intrp. coef.

      INTEGER    IJ, K
!*
!*          < interp. temp. >
!*
      DO K = 2, KMAX
         DO IJ = ISTS, IENS
            FTINTM( K ) = LOG( GDPM(IJ, K ) / GDP(IJ, K)  ) &
                        / LOG( GDP(IJ, K-1) / GDP(IJ, K)  )
            FTINT ( K ) = LOG( GDP(IJ, K-1) / GDPM(IJ, K) ) &
                        / LOG( GDP(IJ, K-1) / GDP(IJ, K)  )
            GDTM( IJ,K ) = FTINTM( K ) * GDT( IJ,K-1 ) &
                         + FTINT ( K ) * GDT( IJ,K   )
         END DO
      END DO

      DO IJ = ISTS, IENS
            GDTM( IJ,KMAX+1 ) = GDT( IJ,KMAX )
            GDTM( IJ,1      ) = GDT( IJ,1 )
      END DO

      RETURN
      END SUBROUTINE TINTP
!***********************************************************************
      subroutine outfld_cs & !! outfld for MIROC variables
               ( name, var_cs, idim, kdim, lchnk )

      use cam_history,   only: outfld

      implicit none

      character(len=*), intent(in) :: name
      integer, intent(in) :: idim, kdim, lchnk
      real(r8), intent(in) :: var_cs(idim,kdim)

      real(r8) :: var_cam(idim,kdim)
      integer :: i, k

      do k = 1, kdim
         do i = 1, idim
            var_cam(i,kdim-k+1) = var_cs(i,k)
         end do
      end do

      call outfld( name, var_cam, idim, lchnk )

      end subroutine outfld_cs
!***********************************************************************

end module cs_conv
