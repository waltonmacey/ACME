subroutine massborrow(subnam,lchnk,ncol,ndim,mbeg,mend,q,pdel) 

!!....................................................................... 
!! Borrow mass from layer above or below to conserve mass 
!! 
!! Author : Kai Zhang (kai.zhang@pnnl.gov) 
!!....................................................................... 

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pver, begchunk, endchunk
  use spmd_utils,      only: masterproc
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice
  use constituents,    only: cnst_get_ind, pcnst
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog

  implicit none

!! interface 
!!....................................................................... 

  character*(*), intent(in) :: subnam                 ! name of calling routine
  integer, intent(in) :: lchnk                        ! chunk identifier
  integer, intent(in) :: ncol                         ! number of atmospheric columns
  integer, intent(in) :: ndim                         ! number of dim members
  integer, intent(in) :: mbeg
  integer, intent(in) :: mend
  real(r8), intent(inout) :: q(ndim,pver,mbeg:mend)  ! moisture/tracer field
  real(r8), intent(in) :: pdel(ndim,pver)            ! 

!! local 
!!....................................................................... 

  integer :: i, k, m, j 
  real(r8):: nmass, zeps
  real(r8):: bmass(ndim)
  integer :: ixcldliq
  integer :: ixcldice
  integer :: ixrain
  integer :: ixsnow
  integer :: tind(5)

  !! init
  !!....................................................................... 

  zeps = epsilon(1.0_r8)

  !! loop over tracers
  !!....................................................................... 

  do m = mbeg, mend

     !!if(masterproc) write(iulog,*) '### tracer index : ', m

     bmass(1:ncol) = 0.0_r8
     
     !! top to bottom
     !!....................................................................... 

     do k = 1, pver
        do i = 1, ncol

           !! new mass in the current layer
           !!....................................................................... 

           nmass = q(i,k,m) + bmass(i)/pdel(i,k)

           if ( nmass > 0.0_r8 ) then

              !! if new mass in the current layer is positive, don't borrow mass any more 
              !!....................................................................... 

              q(i,k,m) = nmass
              bmass(i) = 0.0_r8

           else

              !! set mass to zero in the current layer, and save bmass
              !!....................................................................... 

              bmass(i) = nmass * pdel(i,k)

              q(i,k,m) = 0._r8 

           end if !! nmass > 0.0_r8 

        end do !! i 
     end do !! k 


     !!....................................................................... 
     !! bottom to top
     !!....................................................................... 
     
     do k = pver, 1, -1 

        do i = 1, ncol

           !! if the surface layer still needs to borrow mass 
           !!....................................................................... 

           if (bmass(i) < 0._r8 ) then

              !! new mass in the current layer
              !!....................................................................... 

              nmass = q(i,k,m) + bmass(i)/pdel(i,k)

              if ( nmass > 0.0_r8 ) then

                 !! if new mass in the current layer is positive, don't borrow mass any more 
                 !!....................................................................... 

                 q(i,k,m) = nmass 
                 bmass(i) = 0.0_r8

              else

                 !! if new mass in the current layer is negative, continue to borrow mass
                 !!....................................................................... 

                 bmass(i) = nmass*pdel(i,k)
                 q(i,k,m) = 0._r8 

              end if !! nmass > 0.0_r8 

           end if !! bmass(i) < -zeps 

        end do !! i 

     end do !! k 

  end do !! m

  return 
  end subroutine massborrow

