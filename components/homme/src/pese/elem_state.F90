!
! element-state for PESE primitive-equation vertical spectral-elements
!_______________________________________________________________________
#if 0
type, public :: elem_state_t

  real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)                ! horizontal velocity
  real (kind=real_kind) :: T   (np,np,nlev,timelevels)                  ! temperature
  real (kind=real_kind) :: lnps(np,np,timelevels)                       ! log surface pressure
  real (kind=real_kind) :: ps_v(np,np,timelevels)                       ! surface pressure
  real (kind=real_kind) :: phis(np,np)                                  ! surface geopotential (prescribed)
  real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                     ! tracer concentration
  real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)                   ! tracer mass
  real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)                  ! delta p on levels

  real (kind=real_kind) :: O1  (np,np,nlev,timelevels)                  ! diagnostic output fields
  real (kind=real_kind) :: O2  (np,np,nlev,timelevels)                  ! diagnostic output fields
  real (kind=real_kind) :: O3  (np,np,nlev,timelevels)                  ! diagnostic output fields
  real (kind=real_kind) :: O4  (np,np,nlev,timelevels)                  ! diagnostic output fields

  real (kind=real_kind) :: omega(np,np,nlev)                            ! dp/dt

end type elem_state_t

! num vars in state struct used for standalone HOMME restart file
integer(kind=int_kind),public,parameter :: StateComponents=8

#endif
