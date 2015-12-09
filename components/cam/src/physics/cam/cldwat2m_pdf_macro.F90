!-------------------------------------------------------------------------- !
! Purpose: Compute the stratiform cloud fraction and liquid macrophysical   !
!          tendency (ie condensation minus evaporation in the non-convective!
!          portion of each grid cell using a bivariate Gaussian probability !
!          distribution for subgrid-scale variations in theta_l and qt      !
!                                                                           ! 
! Authors: Peter Caldwell (caldwell19@llnl.gov), Sungsu Park, Steve Klein   !
!                                                                           !
! Date:    2/3/10                                                           !
!-------------------------------------------------------------------------- !

module cldwat2m_pdf_macro

!STUFF USED BY ALL SUBROUTINES BELOW:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use shr_spfn_mod,     only: erf => shr_spfn_erf

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use cam_abortutils,       only: endrun
use cam_logfile,      only: iulog
use physconst,        only: pi
use ref_pres,         only: top_lev => trop_cloud_top_lev


use physconst,    only: latvap,&      !latent heat of vaporization
			cpair         !specific heat of dry air

use wv_saturation, only: qsat_water !subroutine calculates LIQ sat mix ratio

implicit none
private
save

public :: ini_pdf_macro, pdf_mmacro_pcond

!STUFF FOR COMPUTING CRITICAL RH IN get_pdf_width FUNCTION:
!----------------------------------------------------------
real(r8), private            :: premib                       ! Bottom height for mid-level liquid stratus fraction
real(r8), private            :: premit                       ! Top    height for mid-level liquid stratus fraction
real(r8), private            :: rhminl                       ! Critical RH for low-level  liquid stratus clouds
real(r8), private            :: rhminh                       ! Critical RH for high-level liquid stratus clouds
real(r8), private            :: rhminl_adj_land              ! rhminl adjustment for snowfree land
real(r8), private            :: trunc_macro                  ! if >0, use truncated Gaussian w/ this sigma cutoff. Else untruncated.

real(r8), parameter          :: qsmall       = 1.e-18_r8     ! Smallest mixing ratio considered in the macrophysics

real(r8)                     :: erfN2

!END OF STUFF USED BY ALL SUBROUTINES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


contains

!=======================================================================
  subroutine ini_pdf_macro()

    !--------------------------------------------------------------------- !
    !                                                                      ! 
    ! Purpose: Initialize constants for the liquid stratiform macrophysics !
    !          called from stratiform_init in stratiform.F90.              !  
    ! Author:  Peter Caldwell                                              !
    ! Reference: Peter's latex document                                     !
    ! History: created 2/3/10                                              !
    !                                                                      !
    !--------------------------------------------------------------------- !
    use cloud_fraction, only: cldfrc_getparams
    use phys_control, only: phys_getopts

    implicit none

    ! THESE THINGS ARE NEEDED FOR CALCULATING CRITICAL RH IN get_pdf_width
    !----------------------------------------------------------
    call cldfrc_getparams(rhminl_out = rhminl, rhminl_adj_land_out = rhminl_adj_land,    &
         rhminh_out = rhminh, premit_out = premit, premib_out = premib)

    call phys_getopts(trunc_macro_out = trunc_macro)


    if( masterproc ) then
       write(iulog,*) 'ini_macro: tuning parameters : rhminl = ', rhminl, &
            'rhminl_adj_land = ', rhminl_adj_land, & 
            'rhminh = ', rhminh, & 
            'premit = ', premit, & 
            'premib = ',  premib, &
            'trunc_macro = ',trunc_macro
    end if

    if (trunc_macro.gt.0._r8) then
       erfN2=erf(trunc_macro/sqrt(2._r8))
    else
       erfN2=-9999._r8 !should never be used
    end if

    return
  end subroutine ini_pdf_macro


!=======================================================================
subroutine pdf_mmacro_pcond(ncol  , dtime, nl     , qv,ql,     &
     T, pres,a_cu, landfrac, snowh,     & !inputs
     s_tend, qv_tend, ql_tend, nl_tend, & !output tendencies
     cme   , alst   , qlst   )            !other outputs

  !---------------------------------------------------------------------------!
  !                                                                           ! 
  ! Purpose: Compute the stratiform cloud fraction and liquid macrophysical   !
  !          tendency (ie condensation minus evaporation in the non-convective!
  !          portion of each grid cell using a bivariate Gaussian probability !
  !          distribution for subgrid-scale variations in theta_l and qt      !  
  !                                                                           !
  ! Author:  Peter Caldwell                                                   !
  ! Reference: Peter's latex document                                         !
  ! History: created 2/3/10                                                   !
  !                                                                           !
  !---------------------------------------------------------------------------!

  implicit none

  ! INPUTS:
  !-------------------------
  integer,  intent(in)    :: ncol               !number of columns to act on.
  real(r8), intent(in)    :: dtime              ! Timestep					      [ sec ]
  real(r8), intent(in)    :: nl(pcols,pver)     !input Cell-mean liq droplet concentration              [ #/kg ]
  real(r8), intent(in)    :: qv(pcols,pver)     !input cell-mean water vapor content                    [ kg/kg ]
  real(r8), intent(in)    :: ql(pcols,pver)     !input cell-mean liq water content                      [ kg/kg ]
  real(r8), intent(in)    :: T(pcols,pver)      !input cell-mean temperature                            [ K ]
  real(r8), intent(in)    :: pres(pcols,pver)   !input cell-mean pressure                               [ Pa ]
  real(r8), intent(in)    :: a_cu(pcols,pver)   !Convective Cloud Fraction                              [ frac ]
  real(r8), intent(in)    :: landfrac(pcols)    !land fraction                                        [ frac ]
  real(r8), intent(in)    :: snowh(pcols)       !snow depth (water equivalent)                         [ m ]

  ! OUTPUT:
  !-------------------------
  real(r8), intent(out)   :: s_tend(pcols,pver)   !Cell-mean dry static energy tendency due to macrophys[ W/kg ]
  real(r8), intent(out)   :: qv_tend(pcols,pver)  !Cell-mean water vapor content tendency due to macro  [ kg/kg/s ]
  real(r8), intent(out)   :: ql_tend(pcols,pver)  !Cell-mean tendency of liq water due to macrophys     [ kg/kg/s ]
  real(r8), intent(out)   :: nl_tend(pcols,pver)  !Cell-mean liq water droplet tendency due to macro    [ #/kg/s ]
  real(r8), intent(out)   :: cme(pcols,pver)      !Net condensation rate (for this ver, = ql_tend)      [ kg/kg/s ]
  real(r8), intent(out)   :: alst(pcols,pver)     !Cell-mean liquid stratus fraction                    [ frac ]
  real(r8), intent(out)   :: qlst(pcols,pver)     !in-stratiform-cloud liquid water content             [ kg/kg ]
!removed these from function call - just implementing pdf MACRO, not MICRO right now
!  real(r8), intent(out)   :: ql_pk(pcols,pver)    !peak of the Gaussian used for the ql PDF             [ kg/kg ]
!  real(r8), intent(out)   :: ql_std(pcols,pver)   !standard dev of the Gaussian used for the ql PDF     [ kg/kg ]

  ! INTERNAL VARIABLES
  !-------------------------
  real(r8)                :: Tl(pcols,pver)       ! Cell-mean liquid water temperature                  [ K ]
  real(r8)                :: qw(pcols,pver)       ! qv+ql (aka qt - qi)                                 [ kg/kg ]
  real(r8)                :: qs(pcols,pver)       ! saturation mixing ratio at Tl                       [ kg/kg ]
  real(r8)                :: ql_new(pcols,pver)   ! ql at end of macrophysics                           [ kg/kg ]
  real(r8)                :: es(pcols,pver)       ! saturation vapor pressure at mT (not used)          [ Pa ]
  real(r8)                :: Q(pcols,pver)        ! Cell-mean saturation excess (mean(qt-qi)-qs(meanT)) [ kg/kg ]
  real(r8)                :: stdS(pcols,pver)     ! S distribution standard deviation                   [ kg/kg ]
  real(r8)                :: dqs_dT(pcols,pver)   ! dqs_dT evaluated at Tl                              [ kg/(kg K) ]
  real(r8)                :: A(pcols,pver)        ! condensational heating correction                   [ kg/(kg K) ]
  real(r8)                :: erf_part(pcols,pver) ! The error-function part of the cldfrac calculation  [ ??? ]
  real(r8)                :: lcldfrac(pcols,pver) ! fraction of stratiform area of cell having liq cld  [ frac ]
  integer i,k

  ! The PDF is defined over the non-convective region only, so input/output variables should be non-convective region
  ! averages (as opposed to cell-averages).  Sungsu assures me that the input thermodynamic variables are already 
  ! non-convective region averages so I shouldn't have to normalize by Cu cldfrac...

  !COMPUTE LIQUID WATER TEMPERATURE:
  !------------------------
  !PMC not using top_lev here to avoid get_pdf_width operating on
  !uninitialized elements.
  Tl(:ncol,:) = T(:ncol,:) - latvap/cpair*ql(:ncol,:)

  !COMPUTE ICE-FREE TOTAL WATER
  !------------------------
  qw(:ncol,:)=qv(:ncol,:) + ql(:ncol,:) !this is total water mixing ratio minus ice component. 

  !DEAL WITH PDF WIDTH
  !----------------------------------------------------
  !PMC note: arrays in routines always declared w/ pcol size -> must pass pcol sized variables into them...
  stdS(:,:) = get_pdf_width( ncol, qw(:,:), pres(:,:), landfrac(:), snowh(:)  )
!  stdS(:,:) = get_pdf_width_T(ncol, T(:,:), pres(:,:), landfrac(:), snowh(:)  )

  !GET SATURATION MIXING RATIO FOR TL
  !------------------------
  ! Note I'm computing qs(Tl) here, not qs(T). This is because I'm approximating qs(T)=qs(Tl)+dqs/dT*dT
  ! below to avoid iteration.
  !call qsat_water(Tl,pres,es,qs,dqsdt=dqs_dT) !BSINGH - Commented out as it was causing problems with Intel compiler
  do k=top_lev,pver !BSINGH - Added this loop to replace the above call
     call qsat_water(Tl(:ncol,k),pres(:ncol,k),es(:ncol,k),qs(:ncol,k),dqsdt=dqs_dT(:ncol,k))
  enddo
  !this is like total cell-mean water content minus cell-mean saturation, but with ice neglected.
  Q(:ncol,top_lev:) = qw(:ncol,top_lev:) - qs(:ncol,top_lev:)

  !GET A = correction to account for latent heat release:
  !------------------------
  A(:ncol,top_lev:) =  (1._r8+latvap*dqs_dT(:ncol,top_lev:)/cpair)**(-1._r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !NON-TRUNCATED PDF VERSION:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trunc_macro.lt.0._r8) then

     !COMPUTE cldfrac
     !------------------------
     ! Need to loop because erf (the error function) is a only available as a scalar.

     do k=top_lev,pver
        do i=1,ncol
           erf_part(i,k) = erf(Q(i,k) / (sqrt(2._r8)*stdS(i,k) ) )
        end do !for i
     end do !for k

     lcldfrac(:ncol,top_lev:) = 0.5_r8*(1._r8+erf_part(:ncol,top_lev:) ) 

     ! The PDF is only applied over the non-convective region.  Thus, the cell-average cld fraction
     ! is the non-convective region stratiform cld frac (above) times the fraction of the cell that 
     ! is non-convective plus (stratiform cld fraction in convective region... which is zero)*the 
     ! fraction of the cell that is convective...

     alst(:ncol,top_lev:) = (1._r8 - a_cu(:ncol,top_lev:))*lcldfrac(:ncol,top_lev:)

     !COMPUTE CELL-AVE QL
     !------------------------
     ! I've found that rounding can make mql<0. The max() below fixes this. 

     do k=top_lev,pver
        do i=1,ncol

           !PMC 8/1/11 - bugfix: made this *cell-ave* ql
           ql_new(i,k) = (1._r8-a_cu(i,k))*max(0._r8, &
                A(i,k)*( Q(i,k)*lcldfrac(i,k)+stdS(i,k)/sqrt(2._r8*pi)*exp(-Q(i,k)**2._r8/(2._r8*stdS(i,k)**2._r8))))

           !qlst = ql within stratus cloud (IE OUTPUT IS IN-CLOUD QL). 
           !PMC 8/1/11 checked in-cld ql = (non-convect ql)/(non-convect cldfrac)=(all-cell ql)/(all-cell cldfrac)
           if (alst(i,k).gt.qsmall) then 
              qlst(i,k) = ql_new(i,k)/alst(i,k)
           else
              qlst(i,k) = 0._r8
           endif

        end do !for i
     end do !for k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TRUNCATED PDF VERSION:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else !if truncated pdf
     do k=1,pver
        do i=1,ncol    
           if (Q(i,k).le.-trunc_macro*stdS(i,k)) then
              lcldfrac(i,k) = 0._r8
              alst(i,k) = 0._r8
              ql_new(i,k) = 0._r8
              qlst(i,k) = 0._r8
           else if (Q(i,k).ge.trunc_macro*stdS(i,k)) then 
              lcldfrac(i,k) = 1._r8
              alst(i,k) = 1._r8 - a_cu(i,k)
              !A*expected value over whole distn is = A*mean (which is A*Q)
              ql_new(i,k) = (1._r8-a_cu(i,k)) * A(i,k)*Q(i,k)
              qlst(i,k) = ql_new(i,k)/alst(i,k)
           else
              erf_part = erf(Q(i,k) / (sqrt(2._r8)*stdS(i,k) ) )
              lcldfrac(i,k) = 0.5_r8*(1._r8+erf_part(i,k)/erfN2 ) 
              alst(i,k) = (1._r8 - a_cu(i,k))*lcldfrac(i,k)
              ql_new(i,k) = (1._r8-a_cu(i,k)) * A(i,k) *(  &
                   Q(i,k)*lcldfrac(i,k) +stdS(i,k)/(erfN2*sqrt(2._r8*pi))* &
                   (  exp(-Q(i,k)**2._r8/(2._r8*stdS(i,k)**2._r8)) - exp(-trunc_macro**2._r8 / 2._r8) ) )
              if (alst(i,k).gt.qsmall) then !prevent blow up if alst very small?
                 qlst(i,k) = ql_new(i,k)/alst(i,k) 
              else
                 qlst(i,k) = 0._r8
              end if
           end if

        end do !for i
     end do !for k

  end if

  !COMPUTE MACROPHYSICAL TENDENCIES:
  !------------------------
  ql_tend(:ncol,top_lev:) = (ql_new(:ncol,top_lev:) - ql(:ncol,top_lev:))/dtime

  ! in this routine, qv + ql = constant (qi doesn't change), so dqv/dt = -dql/dt.
  qv_tend(:ncol,top_lev:) = -ql_tend(:ncol,top_lev:) 

  ! This routine conserves sl=s-latvap*ql, so ds/dt = +latvap*dql/dt
  s_tend(:ncol,top_lev:)  = latvap*ql_tend(:ncol,top_lev:)

  !according to Sungsu, cme=ql_tend + fixes to keep within bounds. Since I don't apply
  !fixes, cme=ql_tend for this routine. 

  cme(:ncol,top_lev:)=ql_tend(:ncol,top_lev:)

  !HANDLE DROPLET CONCENTRATION CHANGE:
  !-------------------------
  ! If condensation occurs, leave nl alone. Else reduce it proportional to ql decrease.
  ! Since the PDF implementation prevents ql from going negative, the max is probably
  ! not needed.

  do k=top_lev,pver
     do i=1,ncol

        if (ql_tend(i,k) .lt. 0._r8) then
           !Note ql_tend<0 only possible if ql_new<ql => ql must be >0. Thus no check on ql is needed.
           nl_tend(i,k) = max( -nl(i,k)/dtime, nl(i,k)*ql_tend(i,k)/ql(i,k) )
        else
           nl_tend(i,k) = 0._r8
        end if

     end do !for i
  end do !for k

  !COMPUTE PEAK AND STANDARD DEV FOR GAUSSIAN QL DISTN (NEEDED BY MICROPHYS)
  !-------------------------
 ! ql_pk(:ncol,top_lev:)  = A(:ncol,top_lev:)*Q(:ncol,top_lev:)
 ! ql_std(:ncol,top_lev:) = A(:ncol,top_lev:)*stdS(:ncol,top_lev:)

  return
end subroutine pdf_mmacro_pcond


!=======================================================================
function get_pdf_width(ncol, qw, p, landfrac, snowh )
  !---------------------------------------------------------------------------!
  !                                                                           ! 
  ! Purpose: Compute the standard deviation of s=qt-qi-qs(T,pres) for use in  !
  !          pdf_mmacro_pcond.  As an initial cut, this is done by setting T  !
  !          variations to zero and assuming qt variations are proportional   !
  !          to the cell-mean qt.  The constant of proportionality is chosen  !
  !          to give a variance similar to that implied by CAM4's Slingo-style!
  !          cloud fraction parameterization.                                 !
  !                                                                           !
  ! Author:  Peter Caldwell, 2/3/10                                           !
  ! Reference: Peter's gauss_approx.tex document                              !
  ! History: created 2/3/10                                                   !
  !                                                                           !
  ! Future Work: This routine needs major improvement.  The current           !
  !              implementation is just a placeholder for future development. !
  !              The call to this routine may eventually move into            !
  !              macrop_driver or tphysbc in order to have better access to   !
  !              needed quantities.                                           !
  !---------------------------------------------------------------------------!
  implicit none

  ! INPUTS:
  !-----------------------------------------------------------------
  integer , intent(in)    :: ncol                 !number of columns to act on.
  real(r8), intent(in)    :: qw(pcols,pver)       !cell-mean qt - qi                 [ kg/kg ]
  real(r8), intent(in)    :: p(pcols,pver)        !cell-center pressure              [ Pa ]
  real(r8), intent(in)    :: landfrac(pcols)      !land fraction                     [ frac ]
  real(r8), intent(in)    :: snowh(pcols)         !snowdepth (water equivalent)      [ m ]

  !INTERNAL VARIABLES:
  !-----------------------------------------------------------------
  real(r8)                 :: rhcrit(pcols,pver)        ! critical RH                                       [ frac ]
  real(r8)                 :: rhwght                    ! RH parameter inherited from cldwat2m_macro        [ ???  ]
  real(r8)                 :: ds(pcols,pver)            ! half-width of triangular distribution s           [ kg/kg ]
  real(r8), parameter      :: cldrh=1._r8		      ! Relative Humidity in cloud                        [ frac ]
  integer i,k

  !OUTPUT:
  !-----------------------------------------------------------------
  real(r8)                 :: get_pdf_width(pcols,pver) ! Standard deviation of s=qt-qi-qs(T,pres)          [ kg/kg ]

  !GET CRITICAL RH FOLLOWING FORMULA FROM CAM5:
  !--------------------------
  ! This loop is taken more or less directly from cldwat2m_macro.F90

  do k=1,pver
     do i=1,ncol
        if( p(i,k) .ge. premib ) then
           if( (nint(landfrac(i)).eq.1) .and. (snowh(i).le.0.000001_r8) ) then
              rhcrit(i,k) = rhminl - rhminl_adj_land
           else
              rhcrit(i,k) = rhminl
           endif
        elseif( p(i,k) .lt. premit ) then
           rhcrit(i,k)   = rhminh
        else
           rhwght = (premib-(max(p(i,k),premit)))/(premib-premit)
           rhcrit(i,k) = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
        endif
     end do
  end do

  !FIND THE WIDTH OF A TRIANGULAR S DISTRIBUTION CONSISTENT WITH THIS CRITICAL RH:
  !--------------------------
  ! EXPLANATION: Let all quantities refer to cell averages in this note. Define
  ! s=qt-qs(Tl) and ds = the half-width of the truncated Gaussian. If cloud 
  ! is just beginning to form in the cell, RH=qt/qs(T)=rhcrit and s+ds = 0. Use the
  ! definition of s to rewrite the 2nd eq as qs = qt+ds [noting that qs(Tl)=qs(T) 
  ! since cld just beginning to form]. Substitute into the 1st eq to get 
  ! rhcrit*(qt+ds)=qt => ds = (1-rhcrit)/rhcrit*qt

  ! Now assume 2*sigma cutoff => sigma=(1-rhcrit)/(2.*rhcrit)*qt.

  get_pdf_width(:ncol,:) = (1._r8-rhcrit(:ncol,:))/(2._r8*rhcrit(:ncol,:))*qw(:ncol,:)

  ! Note this derivation is really only appropriate for qt, but I'm applying it
  ! to qw because qi is generally a very small portion of qt. 

  return
end function get_pdf_width

!=======================================================================
function get_pdf_width_T(ncol, T, p, landfrac, snowh )
  !---------------------------------------------------------------------------!
  !                                                                           ! 
  ! Purpose: Compute the standard deviation of s=qt-qi-qs(T,pres) for use in  !
  !          pdf_mmacro_pcond following the approach of default CAM5 where    !
  !          width is proportional to qs(T).                                  !
  !                                                                           !
  ! Author:  Peter Caldwell, 6/22/12                                          !
  ! Reference: Peter's gauss_approx.tex document                              !
  ! History: spawned from get_pdf_width() on 6/22/12 to test whether high-cld !
  !          decrease is due to my different formulation for pdf_width        !
  !                                                                           !
  !---------------------------------------------------------------------------!
  implicit none

  ! INPUTS:
  !-----------------------------------------------------------------
  integer , intent(in)    :: ncol                 !number of columns to act on.
  real(r8), intent(in)    :: T(pcols,pver)        !cell-mean T                       [ K ]
  real(r8), intent(in)    :: p(pcols,pver)        !cell-center pressure              [ Pa ]
  real(r8), intent(in)    :: landfrac(pcols)      !land fraction                     [ frac ]
  real(r8), intent(in)    :: snowh(pcols)         !snowdepth (water equivalent)      [ m ]

  !INTERNAL VARIABLES:
  !-----------------------------------------------------------------
  real(r8)                 :: rhcrit(pcols,pver)        ! critical RH                                       [ frac ]
  real(r8)                 :: rhwght                    ! RH parameter inherited from cldwat2m_macro        [ ???  ]
  real(r8)                 :: ds(pcols,pver)            ! half-width of triangular distribution s           [ kg/kg ]
  real(r8), parameter      :: cldrh=1._r8		! Relative Humidity in cloud                        [ frac ]
  real(r8)                 :: qs(pcols,pver)            ! saturation specific humidity                      [ kg/kg ]
  real(r8)                 :: es(pcols,pver)            ! saturation vapor pressure (unneeded, but part of qsat output)
  integer i,k

  !OUTPUT:
  !-----------------------------------------------------------------
  real(r8)                 :: get_pdf_width_T(pcols,pver) ! Standard deviation of s=qt-qi-qs(T,pres)          [ kg/kg ]

  !GET CRITICAL RH FOLLOWING FORMULA FROM CAM5:
  !--------------------------
  ! This loop is taken more or less directly from cldwat2m_macro.F90

  do k=1,pver
     do i=1,ncol
        if( p(i,k) .ge. premib ) then
           if( (nint(landfrac(i)).eq.1) .and. (snowh(i).le.0.000001_r8) ) then
              rhcrit(i,k) = rhminl - rhminl_adj_land
           else
              rhcrit(i,k) = rhminl
           endif
        elseif( p(i,k) .lt. premit ) then
           rhcrit(i,k)   = rhminh
        else
           rhwght = (premib-(max(p(i,k),premit)))/(premib-premit)
           rhcrit(i,k) = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
        endif
     end do
  end do

  !GET QS ARRAY
  !------------------------
  !call qsat_water(T,p,es,qs) !BSINGH - Commented out as it was causing problems with Intel compiler
  do k=top_lev,pver !BSINGH - Added this loop to replace the above call
     call qsat_water(T(:ncol,k),p(:ncol,k),es(:ncol,k),qs(:ncol,k))
  enddo

  !FIND THE WIDTH OF A TRIANGULAR S DISTRIBUTION CONSISTENT WITH THIS CRITICAL RH:
  !--------------------------
  ! EXPLANATION: Let all quantities refer to cell averages in this note. Define
  ! s=qt-qs(T)=(RH-1)*qs and ds = the half-width of a trunc Gauss distribution for s. 
  ! If cloud is just beginning to form in the cell, RH=qt/qs(T)=rhcrit and s+ds = 0
  ! so ds=(1-rhcrit)*qs. Assume 2*sigma cutoff-> ds=2*sigma -> sigma=(1-rhcrit)/2*qs.

  get_pdf_width_T(:ncol,:) = (1._r8-rhcrit(:ncol,:))/2._r8*qs(:ncol,:)


  return
end function get_pdf_width_T

end module cldwat2m_pdf_macro


