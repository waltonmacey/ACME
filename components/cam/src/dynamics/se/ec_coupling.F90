!
!-------------------------------------------------------------------------------
! dynamics - physics coupling for elevation classes
!-------------------------------------------------------------------------------
module ec_coupling

  use shr_kind_mod,   only: r8 => shr_kind_r8, SHR_KIND_CL
  use constituents,   only: pcnst, cnst_name
  use physics_types,  only: physics_state

  implicit none
  private
  save

  ! phys_column_t represents all elevation class info for one physics column
  !       A physics column corresponds to the grid cell scale
  ! elevation classes are packed in physics-column order
  ! ec_loc specifies the chunk and column for all the elevation classes
  !       which make up this physics column (elevation class set).
  !    An elevation class set is all the elevation classes making up one
  !          physics grid cell (typically the same as a GLL cell).
  !    The first dimension of ec_loc matches area and elevation
  type, public :: phys_column_t
    integer                          :: num_elevation_classes
    real(r8),            allocatable :: area(:)                !  fraction?
    real(r8),            allocatable :: elevation(:)           ! m
    ! The second dimension of ec_loc are for the chunk index and column
    integer,             allocatable :: ec_loc(:,:)            ! (#ecs, 2)
    type(physics_state), allocatable :: grid_phys_state
  end type phys_column_t

  ! ec_state_t is used to hold information on a physics state or tendency or
  !            a dynamics state or tendency
  !            The first dimension is usually the number of points per
  !            block (dynamics) or chunk (physics) while the last dimension
  !            is the block or chunk number.
  !            Middle dimensions are number of levels and tracer number.
  type, public :: ec_state_t
    real(r8), allocatable :: u(:,:,:)
    real(r8), allocatable :: v(:,:,:)
    real(r8), allocatable :: t(:,:,:)
    real(r8), allocatable :: omega(:,:,:)
    real(r8), allocatable :: ps(:,:)
    real(r8), allocatable :: phis(:,:)
    real(r8), allocatable :: q(:,:,:,:)
  contains
    procedure          :: allocate      => ec_state_allocate
    procedure          :: copy          => ec_state_copy
  end type ec_state_t


  ! Public module variables
  logical, public              :: ec_active = .false.
  integer, public              :: max_elevation_classes = 1
  ! ec_sets holds information about the elevation classes on a PE (task)
  !    A single task can hold one or more sets of elevation classes.
  !    A set of elevation classes is all the elevation classes making up one
  !          physics grid cell (typically the same as a GLL cell).
  !    
  type(phys_column_t), public, allocatable :: ec_sets(:,:) ! # EC sets on this PE (task)

  ! Private module variables
  character(len=SHR_KIND_CL)   :: elevation_classes
!This should move to dyn_grid
!  type(phys_column_t), pointer :: phys_columns(:,:) ! (nphys_pts, nelemd)

  ! Public interface functions
  public avg_elevation_classes_to_phys_state
  public dyn_state_to_elevation_classes
  public elevation_classes_to_dyn_tend
  public elevation_classes_readnl
  public elevation_classes_init

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !---------------------------------------------------------------------------
  !
  !  ec_state_allocate
  !
  !  Allocate the variables inside an ec_state object
  !
  !---------------------------------------------------------------------------
  subroutine ec_state_allocate(this, nphys, nlev, nelem, pcnst)
    class(ec_state_t), intent(inout)    :: this
    integer,           intent(in)       :: nphys
    integer,           intent(in)       :: nlev
    integer,           intent(in)       :: nelem
    integer,           intent(in)       :: pcnst

    ! Check u to see if the shape is correct
    if (allocated(this%u)) then
      if ( (size(this%u, 1) /= nphys) .or. (size(this%u, 2) /= nlev) .or.     &
           (size(this%u, 3) /= nelem)) then
        deallocate(this%u)
      end if
    end if

    ! Check v to see if the shape is correct
    if (allocated(this%v)) then
      if ( (size(this%v, 1) /= nphys) .or. (size(this%v, 2) /= nlev) .or.     &
           (size(this%v, 3) /= nelem)) then
        deallocate(this%v)
      end if
    end if

    ! Check T to see if the shape is correct
    if (allocated(this%T)) then
      if ( (size(this%T, 1) /= nphys) .or. (size(this%T, 2) /= nlev) .or.     &
           (size(this%T, 3) /= nelem)) then
        deallocate(this%T)
      end if
    end if

    ! Check omega to see if the shape is correct
    if (allocated(this%omega)) then
      if ( (size(this%omega, 1) /= nphys) .or.                                &
           (size(this%omega, 2) /= nlev) .or.                                 &
           (size(this%omega, 3) /= nelem)) then
        deallocate(this%omega)
      end if
    end if

    ! Check ps to see if the shape is correct
    if (allocated(this%ps)) then
      if ( (size(this%ps, 1) /= nphys) .or. (size(this%ps, 2) /= nelem)) then
        deallocate(this%ps)
      end if
    end if

    ! Check phis to see if the shape is correct
    if (allocated(this%phis)) then
      if ( (size(this%phis, 1) /= nphys) .or. (size(this%phis, 2) /= nelem)) then
        deallocate(this%phis)
      end if
    end if

    ! Check q to see if the shape is correct
    if (allocated(this%Q)) then
      if ( (size(this%Q, 1) /= nphys) .or. (size(this%Q, 2) /= nlev) .or.     &
           (size(this%Q, 3) /= pcnst) .or. (size(this%Q, 4) /= nelem)) then
        deallocate(this%Q)
      end if
    end if

    if (.not. allocated(this%u)) then
      allocate(this%u(nphys, nlev, nelem))
    end if

    if (.not. allocated(this%v)) then
      allocate(this%v(nphys, nlev, nelem))
    end if

    if (.not. allocated(this%T)) then
      allocate(this%T(nphys, nlev, nelem))
    end if

    if (.not. allocated(this%omega)) then
      allocate(this%omega(nphys, nlev, nelem))
    end if

    if (.not. allocated(this%ps)) then
      allocate(this%ps(nphys, nelem))
    end if

    if (.not. allocated(this%phis)) then
      allocate(this%phis(nphys, nelem))
    end if

    if (.not. allocated(this%Q)) then
      allocate(this%Q(nphys, nlev, pcnst, nelem))
    end if

  end subroutine ec_state_allocate
  
  !---------------------------------------------------------------------------
  !
  !  ec_state_copy
  !
  !  Copy the contents of one elevation class state into another
  !        Allocate variables, if necessary
  !
  !---------------------------------------------------------------------------
  subroutine ec_state_copy(this, ec_state)
    ! Dummy arguments
    class(ec_state_t), intent(inout) :: this
    class(ec_state_t), intent(in)    :: ec_state

    ! Local variables
    integer                          :: nphys
    integer                          :: nlev
    integer                          :: nelem
    integer                          :: pcnst

    if ( (.not. allocated(this%u))  .or. (.not. allocated(this%v))     .or.   &
         (.not. allocated(this%t))  .or. (.not. allocated(this%omega)) .or.   &
         (.not. allocated(this%ps)) .or. (.not. allocated(this%phis))  .or.   &
         (.not. allocated(this%Q))) then
      ! Assume we can use any ec_state variable for size information
      nphys = size(ec_state%u, 1)
      nlev  = size (ec_state%u, 2)
      nelem = size(ec_state%u, 3)
      pcnst = size(ec_state%Q, 3)
      call this%allocate(nphys, nlev, nelem, pcnst)
    end if

    this%u     = ec_state%u
    this%v     = ec_state%v
    this%t     = ec_state%t
    this%omega = ec_state%omega
    this%ps    = ec_state%ps
    this%phis  = ec_state%phis
    this%Q     = ec_state%Q

  end subroutine ec_state_copy

  !---------------------------------------------------------------------------
  !
  !  avg_elevation_classes_to_phys_state
  !
  !  Average a set of elevation class state or tendencies to the grid cell mean
  !
  !---------------------------------------------------------------------------
  subroutine avg_elevation_classes_to_phys_state(ec_phys_state)
    use ppgrid,         only: begchunk, endchunk, pcols, pver

    ! Dummy arguments
    type(physics_state), intent(in) :: ec_phys_state(begchunk:endchunk)

    ! Local variables
    integer                         :: set, k, ic
    integer                         :: chnk, lcnhk, col
    real(r8)                        :: accum(pver)
    real(r8)                        :: area, atemp

    do chnk = begchunk, endchunk
      do set = 1, size(ec_sets, 1)
        accum(:) = 0._r8
        area = 0._r8
        do ic = 1, ec_sets(set, chnk)%num_elevation_classes
          lcnhk = ec_sets(set, chnk)%ec_loc(ic, 1)
          col = ec_sets(set, chnk)%ec_loc(ic, 2)
          atemp = ec_sets(set, chnk)%area(ic)
          ! Average physics state U (weighted by area?)
          accum(:) = accum(:) + (ec_phys_state(lcnhk)%u(col,:) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, pver
          ec_sets(set, chnk)%grid_phys_state%u(1,k) = accum(k) / area
        end do
      end do
    end do
  end subroutine avg_elevation_classes_to_phys_state

  !---------------------------------------------------------------------------
  !
  !  dyn_state_to_elevation_classes
  !
  !  Create a dynamics tendency and apply it to phys_state
  !
  !---------------------------------------------------------------------------
  subroutine dyn_state_to_elevation_classes(curr_dyn_state, prev_dyn_state, phys_state, phys_columns, dt)
    use physconst,  only: rair, gravit
    use hycoef,     only: hyam, hybm, hyai, hybi, ps0
    ! Dummy arguments
    type(ec_state_t),             intent(in)    :: curr_dyn_state
    type(ec_state_t),             intent(in)    :: prev_dyn_state
    type(ec_state_t),             intent(inout) :: phys_state
    type(phys_column_t), pointer, intent(in)    :: phys_columns(:,:)
    real(r8),                     intent(in)    :: dt

    ! Local variables
    integer                            :: ie, i, m, k, ic, nlev
    integer                            :: index, index2
    integer                            :: pt_beg, pt_end, blk_beg, blk_end
    real(r8), allocatable              :: tend(:)
    real(r8)                           :: work1(max_elevation_classes), work2

    nlev = size(curr_dyn_state%u, 2)
    allocate(tend(nlev))
    pt_beg  = LBOUND(phys_columns, 1)
    pt_end  = UBOUND(phys_columns, 1)
    blk_beg = LBOUND(phys_columns, 2)
    blk_end = UBOUND(phys_columns, 2)
    do ie = blk_beg, blk_end  ! Loop over blocks
      ! Apply U (equally to each elevation class)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%u(index,:,ie) = curr_dyn_state%u(i,:,ie)
        end do
      end do

      ! Apply V (equally to each elevation class)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%v(index,:,ie) = curr_dyn_state%v(i,:,ie)
        end do
      end do

      ! Apply T dynamics tendency
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        tend(:) = curr_dyn_state%t(i,:,ie) - prev_dyn_state%t(i,:,ie)
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%T(index,:,ie) = phys_state%T(index,:,ie) + tend(:)
        end do
      end do

      ! Apply omega 
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%omega(index,:,ie) = curr_dyn_state%omega(i,:,ie)
        end do
      end do

      ! Apply PHIS 
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%phis(index,ie) = gravit*phys_columns(i, ie)%elevation(ic)
        end do
      end do

      ! Apply PS 
      index = 0
      index2 = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        work2=0.0_r8
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index2 = index2 + 1
          work1(ic)=exp(-(phys_state%phis(index2,ie)-curr_dyn_state%phis(i,ie))/(rair*curr_dyn_state%t(i,nlev,ie)))
          work2=work2 + phys_columns(i, ie)%area(ic)*work1(ic)
        end do
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          phys_state%ps(index,ie) = curr_dyn_state%ps(i,ie)*work1(ic)/work2
        end do
      end do

      ! Apply Q dynamics tendency
      do m = 1, pcnst
        index = 0
        do i = pt_beg, pt_end ! Loop over physics columns in an element
          tend(:) = curr_dyn_state%q(i,:,m,ie) - prev_dyn_state%q(i,:,m,ie) ! change due to dynamics
          do ic = 1, phys_columns(i, ie)%num_elevation_classes
            index = index + 1
            do k=1,nlev
              if(tend(k).le. 0._r8)then
                if(prev_dyn_state%q(i,k,m,ie) > 0._r8)then
                  phys_state%q(index,k,m,ie) = phys_state%q(index,k,m,ie) + &
                       tend(k)*phys_state%q(index,k,m,ie)  / prev_dyn_state%q(i,k,m,ie)
                else
                  phys_state%q(index,k,m,ie) = phys_state%q(index,k,m,ie)
                end if
              else
                phys_state%q(index,k,m,ie)=phys_state%q(index,k,m,ie)+tend(k)
              end if
            end do
          end do
        end do
      end do

    end do ! End loop over blocks
    deallocate(tend)

  end subroutine dyn_state_to_elevation_classes

  !---------------------------------------------------------------------------
  !
  !  elevation_classes_to_dyn_tend
  !
  !  Average the physics tendencies in phys_state into dyn_tend
  !
  !---------------------------------------------------------------------------
  subroutine elevation_classes_to_dyn_tend(phys_tend, dyn_tend, phys_columns)
    ! Dummy arguments
    type(ec_state_t),             intent(in)    :: phys_tend
    type(ec_state_t),             intent(inout) :: dyn_tend
    type(phys_column_t), pointer, intent(in)    :: phys_columns(:,:)

    ! Local variables
    integer                            :: ie, i, index, m, k, ic, nlev
    integer                            :: pt_beg, pt_end, blk_beg, blk_end
    real(r8), allocatable              :: tend(:)
    real(r8)                           :: area, atemp

    nlev = size(dyn_tend%u, 2)
    allocate(tend(nlev))
    pt_beg  = LBOUND(phys_columns, 1)
    pt_end  = UBOUND(phys_columns, 1)
    blk_beg = LBOUND(phys_columns, 2)
    blk_end = UBOUND(phys_columns, 2)
    do ie = blk_beg, blk_end  ! Loop over blocks
      ! Average U physics tendency (weighted by area?)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        tend(:) = 0._r8
        area = 0._r8
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          atemp = phys_columns(i,ie)%area(ic)
          tend(:) = tend(:) + (phys_tend%u(index,:,ie) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, nlev
          dyn_tend%u(i,k,ie) = tend(k) / area
        end do
      end do

      ! Average V physics tendency (weighted by area?)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        tend(:) = 0._r8
        area = 0._r8
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          atemp = phys_columns(i,ie)%area(ic)
          tend(:) = tend(:) + (phys_tend%v(index,:,ie) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, nlev
          dyn_tend%v(i,k,ie) = tend(k) / area
        end do
      end do

      ! Average T physics tendency (weighted by area?)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        tend(:) = 0._r8
        area = 0._r8
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          atemp = phys_columns(i,ie)%area(ic)
          tend(:) = tend(:) + (phys_tend%t(index,:,ie) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, nlev
          dyn_tend%t(i,k,ie) = tend(k) / area
        end do
      end do

      do m = 1, pcnst
      ! Average q physics tendency (weighted by area?)
      index = 0
      do i = pt_beg, pt_end ! Loop over physics columns in an element
        tend(:) = 0._r8
        area = 0._r8
        do ic = 1, phys_columns(i, ie)%num_elevation_classes
          index = index + 1
          atemp = phys_columns(i,ie)%area(ic)
          tend(:) = tend(:) + (phys_tend%q(index,:,m,ie) * atemp)
          area = area + atemp
        end do ! End loop over elevation classes
        !$omp parallel do private(k)
        do k = 1, nlev
          dyn_tend%q(i,k,m,ie) = tend(k) / area
        end do
      end do
      end do

    end do ! End loop over blocks
    deallocate(tend)

  end subroutine elevation_classes_to_dyn_tend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This code should be called from dyn_comp
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine elevation_classes_readnl(NLFilename)
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: masterproc, masterprocid, mpicom
    use spmd_utils,     only: mpi_integer, mpi_logical, mpi_character
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog

    ! Dymmy arguments
    character(len=*),    intent(in)  :: NLFileName

    ! Local variables
    integer                          :: ierr
    integer                          :: unitn

    namelist /ctl_nl/ elevation_classes

    ! Default is no elevation classes
    elevation_classes = ''

    ! Read the namelist
    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(NLFilename), status='old')
      call find_group_name(unitn, 'ctl_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, ctl_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun('elevation_classes_readnl: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)

      if (len_trim(elevation_classes) > 0) then
        write(iulog, *) 'Elevation Classes will be read from ',trim(elevation_classes)
      else
        write(iulog, *) 'Elevation Classes are disabled'
      end if
    end if

    call mpi_bcast(elevation_classes, len(elevation_classes), mpi_character, masterprocid, mpicom, ierr)
      
  end subroutine elevation_classes_readnl

  subroutine elevation_classes_init(input_gridname, pts_per_block, numblocks, &
    num_subgrids, subgrid_area, subgrid_elev)
    use ncdio_atm,        only: infld
    use cam_pio_utils,    only: cam_pio_openfile, cam_pio_closefile
    use cam_pio_utils,    only: cam_pio_handle_error
    use pio,              only: PIO_NOWRITE, PIO_inq_dimid, PIO_inq_dimlen
    use pio,              only: file_desc_t

    !! Dummy arguments
    character(len=*),      intent(in)  :: input_gridname
    integer,               intent(in)  :: pts_per_block
    integer,               intent(in)  :: numblocks
    integer,  allocatable, intent(out) :: num_subgrids(:,:)
    real(r8), allocatable, intent(out) :: subgrid_area(:,:,:)
    real(r8), allocatable, intent(out) :: subgrid_elev(:,:,:)

    !! Local variables
    type(file_desc_t)                  :: fh_ec
    integer                            :: ierr
    integer                            :: dimid
    logical                            :: found
    character(len=*), parameter        :: subname = 'ELEVATION_CLASSES_INIT'

    call cam_pio_openfile(fh_ec, trim(elevation_classes), PIO_NOWRITE)

    ! We need to know the maximum number of elevation classes
    ierr = PIO_inq_dimid(fh_ec, 'MaxNoClass', dimid)
    call cam_pio_handle_error(ierr, subname//': Error finding dimension, MaxNoClass')
    ierr = PIO_inq_dimlen(fh_ec, dimid, max_elevation_classes)
    
    ! Read the relevant variables
    if (allocated(num_subgrids)) then
      deallocate(num_subgrids)
    end if
    allocate(num_subgrids(pts_per_block, numblocks))
    num_subgrids = 0
    call infld('NumSubgrids', fh_ec, 'ncol', 'MaxNoClass', 1, pts_per_block,  &
         1, numblocks, num_subgrids, found, gridname=trim(input_gridname))

    if (allocated(subgrid_area)) then
      deallocate(subgrid_area)
    end if
    allocate(subgrid_area(pts_per_block, max_elevation_classes, numblocks))
    subgrid_area = 0.0_r8
    call infld('SubgridAreaFrac', fh_ec, 'ncol', 'MaxNoClass',                &
         1, pts_per_block, 1, max_elevation_classes, 1, numblocks,            &
         subgrid_area, found, gridname=trim(input_gridname))

    if (allocated(subgrid_elev)) then
      deallocate(subgrid_elev)
    end if
    allocate(subgrid_elev(pts_per_block, max_elevation_classes, numblocks))
    subgrid_elev = 0.0_r8
    call infld('AveSubgridElv', fh_ec, 'ncol', 'MaxNoClass',                  &
         1, pts_per_block, 1, max_elevation_classes, 1, numblocks,            &
         subgrid_elev, found, gridname=trim(input_gridname))

    call cam_pio_closefile(fh_ec)
  end subroutine elevation_classes_init

end module ec_coupling
