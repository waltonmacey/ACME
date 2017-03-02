module clm_pflotran_interfaceMod

!----------------------------------------------------------------------------------------------
! NGEE-Arctic CLM-PFLOTRAN soil bgc coupling interface
! authors: Fengming YUAN1, Gautam Bishit1,2, and Guoping Tang1
!
!          1.Climate Change Science Institute & Environmental Science Division
!          Oak Ridge National Laboratory
!
!          2.Lawrence Berkley National Laboratory
!
! date: 2012 - 2015
!
! modified by Gangsheng Wang @ ORNL based on clm_bgc_interfaceMod.F90: 8/28/2015, 2/2/2017
!----------------------------------------------------------------------------------------------

#include "shr_assert.h"
!! FLEXIBLE_POOLS
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: clm_pflotran_interfaceMod
  !
  ! !DESCRIPTION:
  ! Performs
  !
  ! !USES:
  ! Most 'USES' are declaired in each subroutine
  !  use shr_const_mod, only : SHR_CONST_G
  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod    , only : bounds_type
  use filterMod    , only : clumpfilter
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg

  ! currently only works with soil columns, i.e. luntype of 'istsoil/istcrop'
  !  use landunit_varcon     , only : istsoil, istcrop

  ! (dummy) variable definitions
  !! (wgs,8/28/2015): ALM types/variables are replaced by clm_bgc_data
  use clm_bgc_interface_data, only : clm_bgc_interface_data_type


#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
  use pflotran_clm_main_module
  use pflotran_clm_setmapping_module
#endif

  ! !PUBLIC TYPES:
  implicit none

  save

  private    ! By default everything is private

#ifdef CLM_PFLOTRAN
  type(pflotran_model_type), pointer, public :: pflotran_m

#ifdef COLUMN_MODE
  integer, pointer, public :: filters_colcount_beg(:)   ! dim: column block's 1st col on filters(1:nc)
#else
  logical, pointer, public :: mapped_gcount_skip(:)     ! dim: inactive grid mask in bounds%begg:endg
#endif
#endif
  !
  character(len=256), private:: pflotran_prefix = ''
  character(len=32), private :: restart_stamp = ''

  real(r8), parameter :: rgas = 8.3144621d0                 ! m3 Pa K-1 mol-1

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_pf_readnl

  ! wrappers around '#ifdef CLM_PFLOTRAN .... #endif' block statements to maintain sane runtime behavior
  ! when pflotran is not available.
  public :: clm_pf_interface_init
  public :: clm_pf_set_restart_stamp
  public :: clm_pf_run
  public :: clm_pf_write_restart
  public :: clm_pf_finalize

  private :: pflotran_not_available

#ifdef CLM_PFLOTRAN
  ! private work functions that truely require '#ifdef CLM_PFLOTRAN .... #endif'
  private :: interface_init
  private :: pflotran_run_onestep
  private :: pflotran_write_checkpoint
  private :: pflotran_finalize
  private :: get_clm_soil_dimension     !!pf17
  private :: get_clm_soil_properties
  private :: get_clm_soil_th
  private :: get_clm_iceadj_porosity 
  private :: get_clm_bgc_conc
  private :: get_clm_bgc_rate
  private :: update_soil_bgc_pf2clm
  private :: update_bgc_gaslosses_pf2clm

!!delete the following subroutines currently for bgc-only
!   private :: update_soil_temperature_pf2clm
!   private :: update_soil_moisture_pf2clm
!   private :: get_clm_bcwflx
!   private :: get_clm_bceflx
!   private :: update_bcflow_pf2clm

#endif

contains

!-----------------------------------------------------------------------
!
! public interface functions allowing runtime behavior regardless of
! whether pflotran is compiled in.
!
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_readnl
  !
  ! !INTERFACE:
  subroutine clm_pf_readnl( NLFilename )
  !
  ! !DESCRIPTION:
  ! Read namelist for clm-pflotran interface
  !
  ! !USES:
    use clm_varctl    , only : iulog
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use clm_nlUtilsMod, only : find_nlgroup_name
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
  ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'clm_pf_readnl'  ! subroutine name
  !EOP
  !-----------------------------------------------------------------------
    namelist / clm_pflotran_inparm / pflotran_prefix

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm-pflotran namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_pflotran_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_pflotran_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg=subname //':: ERROR: reading clm_pflotran_inparm namelist.')
          end if
       end if
       call relavu( unitn )
       write(iulog, '(/, A)') " clm-pflotran namelist:"
       write(iulog, '(A, " : ", A,/)') "   pflotran_prefix", trim(pflotran_prefix)
    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(pflotran_prefix, mpicom)

  end subroutine clm_pf_readnl

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_set_restart_stamp
  !
  ! !INTERFACE:
  subroutine clm_pf_set_restart_stamp(clm_restart_filename)
  !
  ! !DESCRIPTION: Set the pflotran restart date stamp. Note we do NOT
  ! restart here, that gets handled by pflotran's internal
  ! initialization during interface_init_clm_pf()
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=256), intent(in) :: clm_restart_filename
  ! !LOCAL VARIABLES:
    integer :: name_length, start_pos, end_pos
    character(len=32) :: clm_stamp
  !EOP
  !-----------------------------------------------------------------------

    ! clm restart file name is of the form:
    !     ${CASE_NAME}.clm2.r.YYYY-MM-DD-SSSSS.nc
    ! we need to extract the: YYYY-MM-DD-SSSSS
    write(*, '("clm-pf : clm restart file name : ", A/)') trim(clm_restart_filename)
    name_length = len(trim(clm_restart_filename))
    start_pos = name_length - 18
    end_pos = name_length - 3
    clm_stamp = clm_restart_filename(start_pos : end_pos)
    write(*, '("clm-pf : clm date stamp : ", A/)') trim(clm_stamp)
    restart_stamp = clm_stamp

  end subroutine clm_pf_set_restart_stamp


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: pflotran_not_available
  !
  ! !INTERFACE:
  subroutine pflotran_not_available(subname)
  !
  ! !DESCRIPTION:
  ! Print an error message and abort.
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: subname
  ! !LOCAL VARIABLES:
  !EOP
  !-----------------------------------------------------------------------
    call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface has not been compiled " // &
         "into this version of ALM.")
  end subroutine pflotran_not_available


!******************************************************************************************!
!
! public interface function wrappers
!
!------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------
  subroutine clm_pf_interface_init(bounds, filters)

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds
    type(clumpfilter), intent(in) :: filters(:) ! filters

    character(len=256) :: subname = "clm_pf_interface_init()"

#ifdef CLM_PFLOTRAN
    call interface_init(bounds, filters)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_interface_init

  !--------------------------------------------------------------------------------------------

    subroutine clm_pf_run(clm_bgc_data, bounds, filters, ifilter)

    implicit none

    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds      ! bounds of current process
    type(clumpfilter) , intent(inout) :: filters(:)  ! filters on current process
    integer           , intent(in)    :: ifilter     ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "clm_pf_run"

#ifdef CLM_PFLOTRAN
    call pflotran_run_onestep(clm_bgc_data, bounds, filters, ifilter)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_run

  !-----------------------------------------------------------------------------
  subroutine clm_pf_write_restart(date_stamp)

    implicit none
    character(len=*), intent(in) :: date_stamp

    character(len=32) :: subname = "clm_pf_write_restart"

#ifdef CLM_PFLOTRAN
    call pflotran_write_checkpoint(date_stamp)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_write_restart


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: clm_pf_finalize
  !
  ! !INTERFACE:
  subroutine clm_pf_finalize()

    implicit none
    character(len=256) :: subname = "clm_pf_finalize"

#ifdef CLM_PFLOTRAN
    call pflotran_finalize()
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_finalize




#ifdef CLM_PFLOTRAN

!************************************************************************************!
! The following is private and requires explicit coupling between CLM and PFLOTRAN
!
!------------------------------------------------------------------------------------!
!

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: interface_init
  !
  ! !INTERFACE:
  subroutine interface_init(bounds, filters)
    !
    ! !DESCRIPTION:
    ! initialize the pflotran iterface
    !
    ! !USES:
    use clm_varctl      , only : iulog
    use GridcellType    , only : grc
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use landunit_varcon , only : istsoil, istcrop
    use decompMod       , only : get_proc_global, get_proc_clumps, ldecomp
    use spmdMod         , only : mpicom, masterproc, iam, npes
    use domainMod       , only : ldomain, lon1d, lat1d  !!lon0,lat0

    use clm_time_manager, only : get_nstep !!nsstep, nestep
    use clm_varcon      , only : dzsoi, zisoi
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevdecomp, ndecomp_pools
    use clm_varctl      , only : pf_hmode, pf_tmode, pf_cmode, pf_frzmode,  &
                                 initth_pf2clm, pf_clmnstep0,               &
                                 pf_surfaceflow  !! = pflotran_surfaceflow
    use filterMod       , only : allocFilters


    use CNDecompCascadeConType , only : decomp_cascade_con

!     use abortutils      , only : endrun


    ! pflotran
    use Option_module,                only : printErrMsg
    use Simulation_Base_class,        only : simulation_base_type
    use Simulation_Subsurface_class,  only : simulation_subsurface_type
    use Realization_Base_class,       only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use PFLOTRAN_Constants_module
    use pflotran_clm_setmapping_module
    use Mapping_module
    ! !ARGUMENTS:

    implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
! #include "finclude/petscviewer.h"

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    ! Revised by Fengming Yuan, CCSI-ORNL
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    type(bounds_type), intent(in) :: bounds     ! bounds
    type(clumpfilter), intent(in) :: filters(:) ! filters

    integer  :: global_numg             ! total number of gridcells across all processors (active)
    integer  :: global_numc             ! total number of columns across all processors (active)
    integer  :: global_nump             ! total number of pfts across all processors (active)
    integer  :: g,l,c, pid              ! indices
    integer  :: nc, nclumps, total_soilc, start_mappedgid, nx, ny, npes_pf
    integer, pointer :: total_soilc_all(:)   ! dim: npes
    integer, pointer :: mapped_gid(:)        ! dim: num_soilc
    integer  :: gid, gcount, colcount, cellcount
    integer  :: gcolumns(1:bounds%endg-bounds%begg+1)
    integer  :: ierr

    character(len= 32) :: subname = 'interface_init' ! subroutine name

    integer, pointer :: clm_all_cell_ids_nindex(:)
    integer, pointer :: clm_top_cell_ids_nindex(:)
    integer, pointer :: clm_bot_cell_ids_nindex(:)
    integer :: clm_all_npts
    integer :: clm_top_npts
    integer :: clm_bot_npts
    real(r8):: x0, x1, y0, y1, dx_global(1:ldomain%ni), dy_global(1:ldomain%nj)
    integer :: i, j
    real(r8):: lon0, lat0 !!origin of longitude/latitude
    integer :: nv  !! number of vertices
    class(realization_subsurface_type), pointer    :: realization

    associate( &
         ! Assign local pointers to derived subtypes components (landunit-level)
         ltype      =>  lun%itype      , & !  [integer (:)]  landunit type index
         lgridcell  =>  lun%gridcell   , & !  [integer (:)]  gridcell index of landunit
         ! Assign local pointer to derived subtypes components (column-level)
         cgridcell  =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
         clandunit  =>  col%landunit   , & !  [integer (:)]  landunit index of column
         cactive    =>  col%active       & !  [logic (:)]  column active or not
         )


    !!wgs:beg----------------------------------------------------------------
    !! wgs: ACME-v2: add lon0/lat0 to ldomain
    lon0 = ldomain%lon0
    lat0 = ldomain%lat0
    !!wgs:end----------------------------------------------------------------

    ! (0) determines Total grids/columns to be mapped between CLM and PFLOTRAN

#ifdef COLUMN_MODE
    ! counting 'soilc' in 'filters'
    ! NOTE: only works for 'soilc', which actually includes natural soil and crop land units
    !
    total_soilc = 0
    nclumps = get_proc_clumps()
    do nc=1,nclumps ! current process
       total_soilc = total_soilc + filters(nc)%num_soilc
    end do

    call mpi_barrier(mpicom, ierr)      ! needs all processes done first
    ! sum of all active soil columns across all processors
    call mpi_allreduce(total_soilc, global_numc, 1, MPI_INTEGER,MPI_SUM,mpicom,ierr)

    ! (active) 'soilc' global indexing across all processes
    ! and must be in order, corresponding to 'bounds%begc:endc' filtered by 'filter%soilc'
    ! bounds:           begc,  begc+1, begc+2, ..., endc               ! note: must skip those not 'ltype(l)==istsoil' or 'ltype(l)==istcrop'
    ! filters:          (1)%soilc(:), (2)%soilc(:), ..., (nclumps)%soilc(:)
    ! columnid:         start_columnid(1), start_colunmid+1,  ..., start_columnid+total_soilc-1

    allocate (total_soilc_all(0:npes-1))
    call mpi_gather(total_soilc, 1, MPI_INTEGER, &
                    total_soilc_all, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast(total_soilc_all, npes, MPI_INTEGER, 0, mpicom, ierr)

    ! CLM's natural grid id, continued and ordered across processes, for mapping to PF mesh
    ! will be assigned to calculate 'clm_cell_ids_nindex' below
    allocate(mapped_gid(1:total_soilc))

    ! corresponding starting no. of re-ordered filtered-columns (zero-based)
    ! will be used to skip block(s) of columns in 'bounds%begc:endc' for 'filters(nc)%soilc' in the interface data vec.
    allocate(filters_colcount_beg(1:nclumps))

    start_mappedgid = 0
    colcount = 0
    do pid=0, npes-1
       if (pid==iam) then

           do nc=1, nclumps
             filters_colcount_beg(nc) = 0
             do c=1,filters(nc)%num_soilc
               colcount = colcount + 1
               mapped_gid(colcount) = start_mappedgid + c

             enddo
             if (nc>1) filters_colcount_beg(nc) = colcount
             start_mappedgid = start_mappedgid + filters(nc)%num_soilc    ! skip by fileter's stride, when 'pid==iam'
           enddo

           exit ! do pid=0,npes-1
       endif
       start_mappedgid = start_mappedgid + total_soilc_all(pid)           ! skip by pes's total soil column, until 'pid==iam' and exit do-loop

    end do

#else
    ! assumption that only 1 natural/cropped soil-column allowed per grid cell NOW
    gcolumns(:) = 0   !grc%ncolumns NOT assigned values at all
    do c = bounds%begc, bounds%endc
      l = clandunit(c)
      g = cgridcell(c)

      if ((.not.(ltype(l)==istsoil)) .and. (.not.(ltype(l)==istcrop)) ) then
         write (iulog,*) 'WARNING: Land Unit type of Non-SOIL/CROP... within the domain'
         write (iulog,*) 'CLM-- PFLOTRAN does not support this land unit at present, AND will skip it'

       else
         ! count active columns for a gridcell
         gcount = g - bounds%begg + 1
         if (cactive(c)) then
            gcolumns(gcount) = gcolumns(gcount)+1
         end if

       endif

    enddo ! do c = bounds%begc, bounds%endc

    ! do checking on assumption: 1 soil col. (either natveg or crop, but not both) per grid
    do g = bounds%begg, bounds%endg
        if (gcolumns(g-bounds%begg+1) /= 1) then
          write (iulog,*) 'ERROR: More than 1 or no ACTIVE column found in gridcell:', g, gcolumns(g-bounds%begg+1)
          write (iulog,*) 'CLM-PFLOTRAN does not support this at present, AND please check your surface data, then re-run'
          write (iulog,*) ' i.e., this mode is used for user-defined CLM grid, which may be generated together with PF mesh'

          call endrun(trim(subname) // ": ERROR: Currently does not support multiple or inactive soil column per grid " // &
         "in this version of CLM-PFLOTRAN.")

        endif
    enddo


    ! counting active gridcells in current pes
    total_soilc = 0
    do l=bounds%begl, bounds%endl
        if(ltype(l)==istsoil .or. ltype(l)==istcrop) then
           total_soilc = total_soilc + 1
        endif
    enddo

    call mpi_barrier(mpicom, ierr)      ! needs all processes done first
    ! sum of all active gridcells across all processors
    call mpi_allreduce(total_soilc, global_numc, 1, MPI_INTEGER,MPI_SUM,mpicom,ierr)

    ! collect info for all pes
    allocate (total_soilc_all(0:npes-1))
    call mpi_gather(total_soilc, 1, MPI_INTEGER, &
                    total_soilc_all, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast(total_soilc_all, npes, MPI_INTEGER, 0, mpicom, ierr)

    ! CLM's natural grid id, continued and ordered across processors, for mapping to PF mesh
    ! will be assigned to calculate 'clm_cell_ids_nindex' below
    allocate(mapped_gid(1:total_soilc))

    ! mark the inactive grid (non-natveg/crop landunits)
    ! will be used to skip inactive grids in 'bounds%begg:endg'
    allocate(mapped_gcount_skip(1:bounds%endg-bounds%begg+1))

    gcount = 0
    mapped_gcount_skip(:) = .true.

!!wgs: endl = 16 > 1 in ACME
!!wgs: pflotran requires mapped_gcount_skip(g-bounds%begg+1) = .false. to pass data between ALM & pflotran
    do l=bounds%begl, bounds%endl
       g = lgridcell(l)
       if(ltype(l)==istsoil .or. ltype(l)==istcrop) then
           gcount = gcount + 1
           mapped_gid(gcount) = grc%gindex(g)      ! this is the globally grid-index, i.e. 'an' as shown below

           mapped_gcount_skip(g-bounds%begg+1) = .false.

       endif
    end do

#endif  
!! #ifdef COLUMN_MODE

    if (masterproc) then
      write(iulog,*) '%%--------------------------------------------------------%%'
      write(iulog,*) '%%                                                        %%'
      write(iulog,*) '%%               clm_pf_interface_init                    %%'
      write(iulog,*) '%%                                                        %%'
#ifdef COLUMN_MODE
      write(iulog,*) '%%              ---  1D COLUMN-MODE  ---                  %%'
      write(iulog,*) '%%  Total soil columns (natveg+crop): ',global_numc,'     %%'
#else
      write(iulog,*) '%%           ---   FULLY-3D COUPLED-MODE  ---             %%'
      write(iulog,*) '%%  Total grids with active soil columns: ',global_numc,' %%'
#endif
      write(iulog,*) '%%                                                        %%'
      write(iulog,*) '%%--------------------------------------------------------%%'
      write(iulog,*) ' '
    endif


    pf_clmnstep0 = get_nstep()

    !----------------------------------------------------------------------------------------
    ! (1) Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    !----------------------------------------------------------------------------------------
    ! (2) passing grid/mesh info to interface_data so that PF mesh can be established/mapped

    !(2a) domain/decompose
    clm_pf_idata%nzclm_mapped = nlevgrnd      ! the soil layer no. mapped btw CLM and PF for data-passing

    if (masterproc) then
      write(iulog,*) '%%                                                     %%'
      write(iulog,*) '%% CLM-Layer-No.     PF-Layer-NO.    Thickness (m)     %%'
      do j=clm_pf_idata%nzclm_mapped, 1, -1
         write(iulog,*) j, clm_pf_idata%nzclm_mapped-j+1, dzsoi(j)
      enddo
      write(iulog,*) '%%                                                     %%'
      write(iulog,*) ' '
    endif


#ifdef COLUMN_MODE
    clm_pf_idata%nxclm_mapped = global_numc  ! 1-D format, along X direction
    clm_pf_idata%nyclm_mapped = 1            ! 1-D format

    clm_pf_idata%npx = npes
    clm_pf_idata%npy = 1
    clm_pf_idata%npz = 1

    if(.not.associated(clm_pf_idata%clm_lx)) &
    allocate(clm_pf_idata%clm_lx(1:clm_pf_idata%npx))
    if(.not.associated(clm_pf_idata%clm_ly)) &
    allocate(clm_pf_idata%clm_ly(1:clm_pf_idata%npy))
    if(.not.associated(clm_pf_idata%clm_lz)) &
    allocate(clm_pf_idata%clm_lz(1:clm_pf_idata%npz))

    !
    do pid=0, npes-1
       clm_pf_idata%clm_lx(pid+1) = total_soilc_all(pid)
    end do
    clm_pf_idata%clm_ly = 1
    clm_pf_idata%clm_lz = clm_pf_idata%nzclm_mapped

    deallocate(total_soilc_all)

#else
    clm_pf_idata%nxclm_mapped = ldomain%ni  ! longitudial
    clm_pf_idata%nyclm_mapped = ldomain%nj  ! latidudial

    ! Currently, the following IS only good for user-defined non-global soil domain
    ! AND, the CLM grids ONLY over-rides PF mesh, when 'mapping_files' not provided
    ! i.e. only used for structured-grid.

    ! due to virtually 2-D surface-grid, along which PF structured-grids are decomposed,
    !the 'npes' used by PF must be some specific number
    nx = ldomain%ni
    ny = ldomain%nj
    if(npes<max(nx,ny)) then
        npes_pf = npes
    else
        npes_pf = npes - min(mod(npes,nx), mod(npes,ny))
    endif
    ! no. of processors along X-/Y-direction
    if(mod(npes_pf,nx)==0) then
        clm_pf_idata%npx = nx
        clm_pf_idata%npy = npes_pf/nx
    elseif(mod(npes_pf,ny)==0) then
        clm_pf_idata%npx = npes_pf/ny
        clm_pf_idata%npy = ny
    else
        if(nx<ny) then
            clm_pf_idata%npx = 1
            clm_pf_idata%npy = npes_pf
        else
            clm_pf_idata%npx = npes_pf
            clm_pf_idata%npy = 1
        endif
    endif
    clm_pf_idata%npz = 1                               ! No decompose along Z-direction currently

    ! calculate node no. for each processor along X-/Y-/Z-direction
    if(.not.associated(clm_pf_idata%clm_lx)) &
    allocate(clm_pf_idata%clm_lx(1:clm_pf_idata%npx))
    if(.not.associated(clm_pf_idata%clm_ly)) &
    allocate(clm_pf_idata%clm_ly(1:clm_pf_idata%npy))
    if(.not.associated(clm_pf_idata%clm_lz)) &
    allocate(clm_pf_idata%clm_lz(1:clm_pf_idata%npz))

    clm_pf_idata%clm_lx(:) = (nx-mod(nx,clm_pf_idata%npx))/clm_pf_idata%npx
    do i=1, mod(nx,clm_pf_idata%npx)
       clm_pf_idata%clm_lx(i) = clm_pf_idata%clm_lx(i)+1
    end do
    clm_pf_idata%clm_ly(:) = (ny-mod(ny,clm_pf_idata%npy))/clm_pf_idata%npy
    do j=1, mod(ny,clm_pf_idata%npy)
       clm_pf_idata%clm_ly(j) = clm_pf_idata%clm_ly(j)+1
    end do
    clm_pf_idata%clm_lz = clm_pf_idata%nzclm_mapped

    ! clean-ups
    deallocate(total_soilc_all)

#endif

    if (masterproc) then
       write(iulog,*) ''
       write(iulog,*) '----- actual PFLOTRAN mesh decompose ----- '
       write(iulog,*) 'no. of processors on X: ',clm_pf_idata%npx
       write(iulog,*) 'no. of processors on Y: ',clm_pf_idata%npy
       write(iulog,*) 'no. of processors on Z: ',clm_pf_idata%npz
    endif


    ! (2b) clm land-surface grid length/width, if provided for over-ride PF's mesh
    ! unit IS upon input (now: degree)

#ifdef COLUMN_MODE
    ! dummy values here, will be modified later

    if(.not.associated(clm_pf_idata%dxclm_global)) &
    allocate(clm_pf_idata%dxclm_global(1:clm_pf_idata%nxclm_mapped))
    if(.not.associated(clm_pf_idata%dyclm_global)) &
    allocate(clm_pf_idata%dyclm_global(1:clm_pf_idata%nyclm_mapped))
    clm_pf_idata%x0clm_global = 0.d0
    clm_pf_idata%y0clm_global = 0.d0
    clm_pf_idata%dxclm_global = 1.d0
    clm_pf_idata%dyclm_global = 1.d0

#else
    ! if given vertices,
    ! then have enough information to configure out spatial connections of CLM grids
    ! ldomain%nv missing in ACME1
    
    if (ldomain%nv == 4) then

      if(.not.associated(clm_pf_idata%dxclm_global)) &
      allocate(clm_pf_idata%dxclm_global(1:clm_pf_idata%nxclm_mapped))
      if(.not.associated(clm_pf_idata%dyclm_global)) &
      allocate(clm_pf_idata%dyclm_global(1:clm_pf_idata%nyclm_mapped))

      ! globally grids
      !NOTE: lon1d/lat1d are the centroid of a grid along longitude-axis/latitude-axis
      clm_pf_idata%x0clm_global = lon0
      x0 = lon0
      x1 = lon1d(1) + (lon1d(1) - lon0)

      do i=1,ldomain%ni-1
          dx_global(i) = abs(x1 - x0)
          x0 = x1
          x1 = lon1d(i+1) + (lon1d(i+1)-x0)
      end do
      dx_global(i) = abs(x1 - x0)

      clm_pf_idata%y0clm_global = lat0
      y0 = lat0
      y1 = lat1d(1) + (lat1d(1) - lat0)
      do j=1,ldomain%nj-1
          dy_global(j) = abs(y1 - y0)
          y0 = y1
          y1 = lat1d(j+1) + (lat1d(j+1)-y0)
      end do
      dy_global(j) = abs(y1 - y0)

      ! passing info to clm-pf-interface
      do i=1, ldomain%ni
         clm_pf_idata%dxclm_global(i)=dx_global(i)
      end do
      do j=1, ldomain%nj
         clm_pf_idata%dyclm_global(j)=dy_global(j)
      end do

    end if

#endif

    ! the soil layer thickness mapped btw CLM and PF for data-passing consistently
    ! unit IS upon CLM (now: meter)
    if(.not.associated(clm_pf_idata%dzclm_global)) &
    allocate(clm_pf_idata%dzclm_global(1:clm_pf_idata%nzclm_mapped))
    do j=1,clm_pf_idata%nzclm_mapped
       clm_pf_idata%dzclm_global(j) = dzsoi(clm_pf_idata%nzclm_mapped-j+1)
    end do
    clm_pf_idata%z0clm_global = 0._r8

    !----------------------------------------------------------------------------------------
    ! (3) passing BGC species no./constants to interface_data

    ! the CLM-CN/BGC decomposing pool size and element number, and constants
    clm_pf_idata%ndecomp_pools    = ndecomp_pools
    clm_pf_idata%ndecomp_elements = 2    ! now: C and N only

    if (.not. associated(clm_pf_idata%floating_cn_ratio)) &
    allocate(clm_pf_idata%floating_cn_ratio(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%decomp_pool_name)) &
    allocate(clm_pf_idata%decomp_pool_name(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%decomp_element_ratios)) &
    allocate(clm_pf_idata%decomp_element_ratios(1:clm_pf_idata%ndecomp_pools,1:clm_pf_idata%ndecomp_elements))

    if (.not. associated(clm_pf_idata%ispec_decomp_c)) &
    allocate(clm_pf_idata%ispec_decomp_c(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ispec_decomp_n)) &
    allocate(clm_pf_idata%ispec_decomp_n(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ispec_decomp_hr)) &
    allocate(clm_pf_idata%ispec_decomp_hr(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ispec_decomp_nmin)) &
    allocate(clm_pf_idata%ispec_decomp_nmin(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ispec_decomp_nimm)) &
    allocate(clm_pf_idata%ispec_decomp_nimm(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ispec_decomp_nimp)) &
    allocate(clm_pf_idata%ispec_decomp_nimp(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%ck_decomp_c)) &
    allocate(clm_pf_idata%ck_decomp_c(1:clm_pf_idata%ndecomp_pools))

    if (.not. associated(clm_pf_idata%fr_decomp_c)) &
    allocate(clm_pf_idata%fr_decomp_c(1:clm_pf_idata%ndecomp_pools,1:clm_pf_idata%ndecomp_pools))

    ! -----------------------------------------------------------------

    ! (4) Create PFLOTRAN model
    call pflotranModelCreate(mpicom, pflotran_prefix, pflotran_m)

    call pflotranModelSetupMappingFiles(pflotran_m)

    ! PFLOTRAN ck file NOT works well when coupling with CLM. So it's off and restart PF from CLM.
    !restart_stamp = ""
    !call pflotranModelSetupRestart(pflotran_m, restart_stamp)
    initth_pf2clm = .false.
    !if (restart_stamp .ne. "") then
    !   initth_pf2clm = .true.
    !endif

    select type (simulation => pflotran_m%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         pflotran_m%option%io_buffer = "This version of clm-pflotran only works with subsurface simulations."
         write(*, '(/A/)') pflotran_m%option%io_buffer
         call printErrMsg(pflotran_m%option)
    end select

    if(pflotran_m%option%nsurfflowdof > 0) then
       pflotran_m%option%io_buffer = "This version of clm-pflotran DOES NOT work with PF Surface simulation."
       write(*, '(/A/)') pflotran_m%option%io_buffer
       call printErrMsg(pflotran_m%option)
    endif
    pf_surfaceflow = .false.

    !------------------------------------------------

    ! Number of cells and Indexing in CLM domain's clumps on current process ('bounds')

#ifdef COLUMN_MODE
    ! soil column-wised for mapping.
    clm_all_npts = total_soilc*clm_pf_idata%nzclm_mapped
    clm_top_npts = total_soilc
    clm_bot_npts = total_soilc
    allocate(clm_all_cell_ids_nindex(1:clm_all_npts))
    allocate(clm_top_cell_ids_nindex(1:clm_top_npts))
    allocate(clm_bot_cell_ids_nindex(1:clm_bot_npts))

    cellcount = 0
    do colcount = 1, total_soilc
       gid = mapped_gid(colcount)

       ! Save cell IDs (0-based) of CLM columns in 1D array on current process ('bounds')
       clm_top_cell_ids_nindex(colcount) = (gid-1)*clm_pf_idata%nzclm_mapped
       do j = 1,clm_pf_idata%nzclm_mapped
          cellcount = cellcount + 1
          clm_all_cell_ids_nindex(cellcount) = (gid-1)*clm_pf_idata%nzclm_mapped + j - 1
       enddo
       clm_bot_cell_ids_nindex(colcount) = gid*clm_pf_idata%nzclm_mapped - 1

    end do
    deallocate(mapped_gid)

#else
    !grid-wised for mapping.
    clm_all_npts = total_soilc*clm_pf_idata%nzclm_mapped      ! 'total_soilc' was only counted with active-column grid
    clm_top_npts = total_soilc
    clm_bot_npts = total_soilc
    allocate(clm_all_cell_ids_nindex(1:clm_all_npts))
    allocate(clm_top_cell_ids_nindex(1:clm_top_npts))
    allocate(clm_bot_cell_ids_nindex(1:clm_bot_npts))

    cellcount = 0
    do gcount = 1, total_soilc

       gid = mapped_gid(gcount)

       ! Save cell IDs of CLM grid
       do j = 1,clm_pf_idata%nzclm_mapped
          cellcount = cellcount + 1
          clm_all_cell_ids_nindex(cellcount) = (gid-1)*clm_pf_idata%nzclm_mapped + j-1   ! zero-based
       enddo
       clm_top_cell_ids_nindex(gcount) = (gid-1)*clm_pf_idata%nzclm_mapped               ! zero-based
       clm_bot_cell_ids_nindex(gcount) = gid*clm_pf_idata%nzclm_mapped-1                 ! zero-based

    enddo
    deallocate(mapped_gid)

#endif


    ! CLM: 3-D Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_sub = clm_all_npts
    clm_pf_idata%ngclm_sub = clm_all_npts

    ! CLM: Surface/Bottom cells of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_2dtop = clm_top_npts
    clm_pf_idata%ngclm_2dtop = clm_top_npts

    ! CLM: bottom face of subsurface domain
    clm_pf_idata%nlclm_2dbot = clm_bot_npts
    clm_pf_idata%ngclm_2dbot = clm_bot_npts

    ! PFLOTRAN: 3-D Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
    clm_pf_idata%ngpf_sub = realization%patch%grid%ngmax

    ! For CLM/PF: ground surface NOT defined, so need to set the following to zero.
    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf  = 0
    clm_pf_idata%ngpf_srf  = 0

    ! Initialize maps for transferring data between CLM and PFLOTRAN.
    if(pflotran_m%map_clm_sub_to_pf_sub%id == CLM_3DSUB_TO_PF_3DSUB) then
       call pflotranModelInitMapping(pflotran_m, clm_all_cell_ids_nindex, &
                                  clm_all_npts, CLM_3DSUB_TO_PF_3DSUB)
    endif
    if(pflotran_m%map_pf_sub_to_clm_sub%id == PF_3DSUB_TO_CLM_3DSUB) then
       call pflotranModelInitMapping(pflotran_m, clm_all_cell_ids_nindex, &
                                  clm_all_npts, PF_3DSUB_TO_CLM_3DSUB)
    endif
    !
    if(pflotran_m%map_clm_2dtop_to_pf_2dtop%id == CLM_2DTOP_TO_PF_2DTOP) then
       call pflotranModelInitMapping(pflotran_m, clm_top_cell_ids_nindex,   &
                                      clm_top_npts, CLM_2DTOP_TO_PF_2DTOP)
    endif
    if(pflotran_m%map_pf_2dtop_to_clm_2dtop%id == PF_2DTOP_TO_CLM_2DTOP) then
       call pflotranModelInitMapping(pflotran_m, clm_top_cell_ids_nindex,   &
                                      clm_top_npts, PF_2DTOP_TO_CLM_2DTOP)
    endif
    !
    if(pflotran_m%map_clm_2dbot_to_pf_2dbot%id == CLM_2DBOT_TO_PF_2DBOT) then
       call pflotranModelInitMapping(pflotran_m, clm_bot_cell_ids_nindex, &
                                     clm_bot_npts, CLM_2DBOT_TO_PF_2DBOT)
    endif
    if(pflotran_m%map_pf_2dbot_to_clm_2dbot%id == PF_2DBOT_TO_CLM_2DBOT) then
       call pflotranModelInitMapping(pflotran_m, clm_bot_cell_ids_nindex, &
                                     clm_bot_npts, PF_2DBOT_TO_CLM_2DBOT)
    endif

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

#ifdef COLUMN_MODE
    ! if 'column-wised' mapping, vertical-flow/transport only mode IS ON by default
    if(pflotran_m%option%nflowdof > 0) then
      pflotran_m%option%flow%only_vertical_flow = PETSC_TRUE
    endif
    if(pflotran_m%option%ntrandof > 0) then
      pflotran_m%option%transport%only_vertical_tran = PETSC_TRUE
    endif

    ! checking if 'option%mapping_files' turned off by default
    if(pflotran_m%option%mapping_files) then
       pflotran_m%option%io_buffer = " COLUMN_MODE coupled clm-pflotran DOES NOT need MAPPING_FILES ON."
       write(*, '(/A/)') pflotran_m%option%io_buffer
       call printErrMsg(pflotran_m%option)
    endif

#endif

    ! if BGC is on
    if(pflotran_m%option%ntrandof > 0) then

      ! the CLM-CN/BGC decomposing pools
      clm_pf_idata%decomp_pool_name  = decomp_cascade_con%decomp_pool_name_history(1:ndecomp_pools)
      clm_pf_idata%floating_cn_ratio = decomp_cascade_con%floating_cn_ratio_decomp_pools(1:ndecomp_pools)
      ! PF bgc species names/IDs
      call pflotranModelGetRTspecies(pflotran_m)
    endif

    deallocate(clm_all_cell_ids_nindex)
    deallocate(clm_top_cell_ids_nindex)
    deallocate(clm_bot_cell_ids_nindex)

!-------------------------------------------------------------------------------------
    ! coupled module controls betweeen PFLOTRAN and CLM45 (F.-M. Yuan, Aug. 2013)
    if(pflotran_m%option%iflowmode==RICHARDS_MODE) then
      pf_hmode = .true.
      pf_tmode = .false.
      pf_frzmode = .false.

    elseif(pflotran_m%option%iflowmode==TH_MODE) then
      pf_hmode = .true.
      pf_tmode = .true.
      if (pflotran_m%option%use_th_freezing) then
         pf_frzmode = .true.
      else
         pf_frzmode = .false.
      endif
    endif

    if(pflotran_m%option%ntrandof.gt.0) then
      pf_cmode = .true.        ! initialized as '.false.' in clm initialization
    endif

    !! wgs:beg-----------------------------------------------------------
    !! wgs:  avoid to use clm_bgc_data during initialization.
    !! wgs: the following 3 subroutine-calls are moved to 'pflotran_run_onestep', where clm_bgc_data can be used

    ! force CLM soil domain into PFLOTRAN subsurface grids
!     call get_clm_soil_dimension(clm_bgc_data, bounds)

    ! Currently always set soil hydraulic/BGC properties from CLM to PF
!     call get_clm_soil_properties(clm_bgc_data, bounds, filters)
    
    ! Get top surface area of 3-D pflotran subsurface domain
!     call pflotranModelGetTopFaceArea(pflotran_m)
    !! wgs:end-----------------------------------------------------------

    ! Initialize PFLOTRAN states
    call pflotranModelStepperRunInit(pflotran_m)

    end associate
  end subroutine interface_init

  !-----------------------------------------------------------------------------
  !
  ! !SUBROUTINE: pflotran_run_onestep
  !
  ! !INTERFACE:

  subroutine pflotran_run_onestep(clm_bgc_data, bounds, filters, ifilter)
  !
  ! !DESCRIPTION:
  !
  !  F.-M. YUAN: based on Gautam's 'step_th_clm_pf',
  !              'chemistry' (PF_CMODE) added (Sept. 6, 2013)
  !
  ! !USES:


!    use clm_pflotran_interface_data
!     use PFLOTRAN_Constants_module
!     use pflotran_clm_setmapping_module
    use spmdMod           , only : mpicom, masterproc, iam, npes
    use clm_time_manager  , only : get_step_size, get_nstep, nsstep, nestep,       &
                                   is_first_step, is_first_restart_step, calc_nestep

    use clm_varctl        , only : pf_tmode, pf_hmode, pf_cmode, &
                                   pf_frzmode, pf_clmnstep0, initth_pf2clm


    implicit none

    type(bounds_type) , intent(in)    :: bounds
    type(clumpfilter) , intent(inout) :: filters(:)     ! filters on current process
    integer           , intent(in)    :: ifilter        ! which filter to be operate

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    !LOCAL VARIABLES:
    real(r8) :: dtime                         ! land model time step (sec)
    integer  :: nstep                         ! time step number
    integer  :: total_clmstep                 ! total clm time step number
    logical  :: ispfprint                     ! let PF printout or not
    logical  :: isinitpf = .FALSE.            ! (re-)initialize PF from CLM or not
    integer  :: ierr

  !-----------------------------------------------------------------------

    nstep = get_nstep() - pf_clmnstep0 !!nsstep
    dtime = get_step_size()

    if (is_first_step() .or. is_first_restart_step()) then
       isinitpf = .TRUE.
    else
       isinitpf = .FALSE.
    endif


    ! (0)
    if (isinitpf) then
      call calc_nestep()  !! nestep
       total_clmstep = nestep - pf_clmnstep0!!nestep - nsstep
       ispfprint = .true.               ! turn-on or shut-off PF's *.h5 output

       call pflotranModelUpdateFinalWaypoint(pflotran_m, total_clmstep*dtime, dtime, ispfprint)

       !! wgs:beg------------------------------------------------
       !! wgs: move from 'interface_init'
       ! force CLM soil domain into PFLOTRAN subsurface grids
       call get_clm_soil_dimension(clm_bgc_data, bounds)

       ! Currently always set soil hydraulic/BGC properties from CLM to PF
       call get_clm_soil_properties(clm_bgc_data, bounds, filters)

       ! Get top surface area of 3-D pflotran subsurface domain
       call pflotranModelGetTopFaceArea(pflotran_m)

       !! wgs:end------------------------------------------------

       ! always initializing soil 'TH' states from CLM to pflotran
       call get_clm_soil_th(clm_bgc_data, .not.initth_pf2clm, .not.initth_pf2clm, bounds, filters, ifilter)

       call pflotranModelUpdateTHfromCLM(pflotran_m, .FALSE., .FALSE.)     ! pass TH to global_auxvar

    endif


    ! (1)
    ! if PF T/H mode not available, have to pass those from CLM to global variable in PF to drive BGC/H
    if (.not. isinitpf .and. (.not.pf_tmode .or. .not.pf_hmode)) then    ! always initialize from CLM to pF, if comment out this 'if'block
       call get_clm_soil_th(clm_bgc_data, .TRUE., .TRUE., bounds, filters, ifilter)

       call pflotranModelUpdateTHfromCLM(pflotran_m, .FALSE., .FALSE.)     ! pass TH to global_auxvar

    end if

    ! ice-len adjusted porostiy
    if (.not.pf_frzmode) then
        call get_clm_iceadj_porosity(clm_bgc_data, bounds, filters, ifilter)

        call pflotranModelResetSoilPorosityFromCLM(pflotran_m)

    endif

    ! (2) pass CLM water fluxes to CLM-PFLOTRAN interface

    ! (3) CLM thermal BC to PFLOTRAN

    ! (4)
    if (pf_cmode) then

    ! (4a) always (re-)initialize PFLOTRAN soil bgc state variables from CLM-CN

       !if (isinitpf) then    ! if only initialize ONCE, uncomment this 'if ... endif' block.

          call get_clm_bgc_conc(clm_bgc_data, bounds, filters, ifilter)
          call pflotranModelSetBgcConcFromCLM(pflotran_m)
          if ((.not.pf_hmode .or. .not.pf_frzmode)) then
          ! this is needed, because at step 0, PF's interface data is empty
          ! which causes Aq. conc. adjustment balacne issue
            call pflotranModelGetSaturationFromPF(pflotran_m)
          endif

       !endif    ! if only initialize ONCE, uncomment this 'if ... endif' block.

       ! MUST reset PFLOTRAN soil aq. bgc state variables from CLM-CN due to liq. water volume change
       ! when NOT coupled with PF Hydrology or NOT in freezing-mode (porosity will be forced to vary from CLM)
        if (.not.pf_hmode .or. .not.pf_frzmode) then
          call pflotranModelUpdateAqConcFromCLM(pflotran_m)
        endif

    ! (4b) bgc rate (fluxes) from CLM to PFLOTRAN
        call get_clm_bgc_rate(clm_bgc_data, bounds, filters, ifilter)
        call pflotranModelSetBgcRatesFromCLM(pflotran_m)

    endif

    ! (5) the main callings of PFLOTRAN
    call mpi_barrier(mpicom, ierr)

    if(mod(nstep+1,48) == 0) then   ! this will allow PFLOTRAN write out every 48 CLM time-step, if relevant option is ON.
       ispfprint = .TRUE.
    else
       ispfprint = .FALSE.
    endif

    call pflotranModelStepperRunTillPauseTime( pflotran_m, (nstep+1.0d0)*dtime, dtime, ispfprint )
    call mpi_barrier(mpicom, ierr)


    ! (6) update CLM variables from PFLOTRAN

    if (pf_cmode) then
        call pflotranModelGetBgcVariablesFromPF( pflotran_m)      ! bgc variables

        call update_soil_bgc_pf2clm(clm_bgc_data, bounds, filters, ifilter)

        ! need to save the current time-step PF porosity/liq. saturation for bgc species mass conservation
        ! if CLM forced changing them into PF at NEXT timestep
        if (.not.pf_hmode .or. .not.pf_frzmode) then
           call pflotranModelGetSaturationFromPF(pflotran_m)
        endif

    endif

  end subroutine pflotran_run_onestep

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: write_checkpoint
  !
  ! !INTERFACE:
  subroutine pflotran_write_checkpoint(date_stamp)
  !
  ! !DESCRIPTION:
  ! Trigger a pflotran checkpoint file to be written
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=32), intent(in) :: date_stamp ! file name date stamp

  ! !LOCAL VARIABLES:

  !EOP
  !-----------------------------------------------------------------------

    ! temporarily OFF - it's not working well for BGC
    ! So, now must initializing PF variables from CLM each start/restart.

    !call pflotranModelStepperCheckpoint(pflotran_m, date_stamp)

  end subroutine pflotran_write_checkpoint

  !-----------------------------------------------------------------------------
  !
  ! !IROUTINE: pflotran_finalize
  !
  ! !INTERFACE:
  subroutine pflotran_finalize()
  !
  ! !DESCRIPTION:
  !
  ! finalizing pflotran runs and destroying objects
  !
  ! !USES:
    use clm_varctl           , only : use_pflotran

    implicit none

  !-----------------------------------------------------------------------

    if (use_pflotran) then
       call pflotranModelDestroy(pflotran_m)
    endif

#ifdef COLUMN_MODE
    deallocate(filters_colcount_beg)
#else
    deallocate(mapped_gcount_skip)
#endif

  end subroutine pflotran_finalize


  ! ============================= GET CLM initial/src-sink/BC to PFLOTRAN ==================================
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_soil_dimension
  !
  ! !INTERFACE:
  subroutine get_clm_soil_dimension(clm_bgc_data, bounds)
    !
    ! !DESCRIPTION:
    ! get soil column dimension to PFLOTRAN
    !
    ! !USES:
    use GridcellType    , only : grc
    use LandunitType    , only : lun
    use ColumnType      , only : col

    use clm_varpar      , only : nlevgrnd
    use domainMod       , only : ldomain
    use landunit_varcon , only : istsoil, istcrop


    !
    ! !ARGUMENTS:

    implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    type(bounds_type)                   , intent(in) :: bounds
    type(clm_bgc_interface_data_type)   , intent(in) :: clm_bgc_data
    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    ! Revised by Fengming Yuan, CCSI-ORNL, May 2015
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    integer  :: g,l,c,j  ! indices
    integer  :: gcount, cellcount   ! gcount: 0-based, cellcount: 1-based
    integer  :: v, n1, n2
    real(r8) :: p1, p2
    real(r8) :: dxsoil_clm(1:bounds%endg-bounds%begg+1)
    real(r8) :: dysoil_clm(1:bounds%endg-bounds%begg+1)
#ifndef COLUMN_MODE
    real(r8) :: wtgcell_sum(1:bounds%endg-bounds%begg+1)
    integer  :: xwtgcell_c(1:bounds%endg-bounds%begg+1)
#endif

    character(len= 32) :: subname = 'get_clm_soil_dimension' ! subroutine name

    PetscScalar, pointer :: cellid_clm_loc(:)
    PetscScalar, pointer :: zisoil_clm_loc(:)     ! 3-D PF-cell's z-node coordinates (elevation-adjusted, unit: m)
    PetscScalar, pointer :: dxsoil_clm_loc(:)     ! 3-D PF-cell's soil length (unit: degrees)
    PetscScalar, pointer :: dysoil_clm_loc(:)     ! 3-D PF-cell's soil width (unit: degrees)
    PetscScalar, pointer :: dzsoil_clm_loc(:)     ! 3-D PF-cell's soil thickness (unit: m)
    PetscScalar, pointer :: xsoil_clm_loc(:)      ! 3-D PF-cell's center x coordinates (unit: m)
    PetscScalar, pointer :: ysoil_clm_loc(:)      ! 3-D PF-cell's center y coordinates (unit: m)
    PetscScalar, pointer :: zsoil_clm_loc(:)      ! 3-D PF-cell's soil depth from ground surface at the layer center (unit: m)
    PetscScalar, pointer :: toparea_clm_loc(:)    ! 3-D PF-cell (unit: m^2)

    PetscErrorCode :: ierr


    ! for calling functions in 'geodesic.for'
    double precision a, f, dummy1, dummy2
    double precision lats(4), lons(4)
    a = 6378137.0d0           ! major-axis length of Earth Ellipsoid in metres in WGS-84
    f = 1.d0/298.257223563d0  ! flatening of Earth Ellipsoid in WGS-84
    
    associate( &
         ! Assign local pointers to derived subtypes components (gridcell-level)
         latc       =>  ldomain%latc   , & !  [real(r8) (:)]
         lonc       =>  ldomain%lonc   , & !  [real(r8) (:)]
         !! latv/lonv missing in ACME1
         latv       =>  ldomain%latv   , & !  [real(r8) (:,:)]
         lonv       =>  ldomain%lonv   , & !  [real(r8) (:,:)]

         lelev      =>  ldomain%topo   , & !  [real(r8) (:)]
         larea      =>  ldomain%area   , & !  [real(r8) (:)]
         ! landunit
         ltype      =>  lun%itype      , & !  [integer (:)]  landunit type index
         ! Assign local pointer to derived subtypes components (column-level)
         clandunit  =>  col%landunit   , & !  [integer (:)]  landunit index of column
         cgridcell  =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
         cwtgcell   =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell)
         cactive    =>  col%active     , & !  [logic (:)]  column active or not
         !
         z          =>  clm_bgc_data%z , & !  [real(r8) (:,:)]  layer depth (m) (sort of centroid from surface 0 )
         zi         =>  clm_bgc_data%zi, & !  [real(r8) (:,:)]  layer interface depth (m)
         dz         =>  clm_bgc_data%dz  & !  [real(r8) (:,:)]  layer thickness (m)
         )


#ifndef COLUMN_MODE
    ! active column weight summation for 1 grid
    wtgcell_sum(:) = 0._r8
    xwtgcell_c(:) = 0

    do c = bounds%begc, bounds%endc
       gcount = cgridcell(c) - bounds%begg + 1
       if (xwtgcell_c(gcount)<=0) xwtgcell_c(gcount) = c
       if (cactive(c)) then
          wtgcell_sum(gcount) = wtgcell_sum(gcount)+cwtgcell(c)

          if( (cwtgcell(c)>=cwtgcell(xwtgcell_c(gcount))) .or. &
              (.not.cactive(xwtgcell_c(gcount))) ) then
             xwtgcell_c(gcount) = c
          end if

       end if

    enddo
#endif

    do g = bounds%begg, bounds%endg
      gcount = g - bounds%begg                               ! 0-based
      ! re-calculating 2-D grid area if vertices are known from input file
      ! NOTE: this will over-write the grid area read-in from either 'ldomain' file or 'surfdata' file
      if (ldomain%nv==4 .or. ldomain%nv==3) then
         if (ldomain%nv==4) then
           lats = latv(g,1:4)
           lons = lonv(g,1:4)
           call area(a, f, lats, lons, 4, dummy1, dummy2)
         else if (ldomain%nv==3) then
           lats(1:3) = latv(g,1:3)
           lons(1:3) = lonv(g,1:3)
           call area(a, f, lats(1:3), lons(1:3), 3, dummy1, dummy2)
         endif

         if (dummy1 < 1.d-20) then
            call endrun(trim(subname) // ": ERROR: re-calculated ldomain%area is less than 0. " // &
            "Please check the grid vertices lat/lon in ldomain file")
         else
           larea(g) = dummy1 * 1.e-6_r8
         endif

         ! for 1-D grid, either 'dx' or 'dy' may be variable and acceptable in the model
         ! (though, currently NOT YET used for PF mesh)
         ! (NOTE: for 2-D grid, dx/dy in lon/lat can NOT be variable)
         if(.not.ldomain%isgrid2d) then

            p1 = 0._r8
            p2 = 0._r8
            n1 = 0
            n2 = 0
            do v = 1, ldomain%nv  !!wgs: ldomain%nv missing in ACME1
               if (lonv(g,v)<lonc(g)) then
                  p1 = p1 + lonv(g,v)
                  n1 = n1 + 1
               elseif (lonv(g,v)>lonc(g)) then
                  p2 = p2 + lonv(g,v)
                  n2 = n2 + 1
               end if
            end do
            if (n1 > 0 .and. n2 > 0) then
               dxsoil_clm(gcount+1) = abs(p2/n2-p1/n1)  ! degs
            endif

            p1 = 0._r8
            p2 = 0._r8
            n1 = 0
            n2 = 0
            do v = 1, ldomain%nv !!wgs: ldomain%nv missing in ACME1
               if (latv(g,v)<latc(g)) then
                  p1 = p1 + latv(g,v)
                  n1 = n1 + 1
               elseif (latv(g,v)>latc(g)) then
                  p2 = p2 + latv(g,v)
                  n2 = n2 + 1
               end if
            end do
            if (n1 > 0 .and. n2 > 0) then
               dysoil_clm(gcount+1) = abs(p2/n2-p1/n1)  ! degs
            endif

         endif  !if(.not.ldomain%isgrid2d)

       endif !if (ldomain%nv==4 .or. 3)

    end do

    !!!!
    call VecGetArrayF90(clm_pf_idata%cellid_clmp,  cellid_clm_loc,  ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%zisoil_clmp,  zisoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%dxsoil_clmp,  dxsoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%dysoil_clmp,  dysoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%dzsoil_clmp,  dzsoil_clm_loc,  ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%area_top_face_clmp,  toparea_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%xsoil_clmp,  xsoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%ysoil_clmp,  ysoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%zsoil_clmp,  zsoil_clm_loc,  ierr)
    CHKERRQ(ierr)

    zisoil_clm_loc(:)   = 0._r8
    dxsoil_clm_loc(:)   = 0._r8
    dysoil_clm_loc(:)   = 0._r8
    dzsoil_clm_loc(:)   = 0._r8

    toparea_clm_loc(:)  = 0._r8
    xsoil_clm_loc(:)    = 0._r8
    ysoil_clm_loc(:)    = 0._r8
    zsoil_clm_loc(:)    = 0._r8

    do j = 1, clm_pf_idata%nzclm_mapped

       if (j <= nlevgrnd) then


#ifdef COLUMN_MODE
        gcount = 0
        do c = bounds%begc, bounds%endc
            g = cgridcell(c)
            l = clandunit(c)

            ! note that filters%soilc includes 'istsoil' and 'istcrop'
            ! (TODO: checking col%itype and lun%itype - appears not match with each other, and col%itype IS messy)
            if (ltype(l)==istsoil .or. ltype(l)==istcrop) then
               gcount = gcount + 1                                    ! actually is the active soil column count
               cellcount = (gcount-1)*clm_pf_idata%nzclm_mapped + j   ! 1-based

               xsoil_clm_loc(cellcount) = lonc(g)
               ysoil_clm_loc(cellcount) = latc(g)
               !
               dzsoil_clm_loc(cellcount) = dz(c, j)                ! cell vertical thickness (m)
               zisoil_clm_loc(cellcount) = -zi(c, j) + lelev(g)    ! cell-node elevation (m)
               zsoil_clm_loc(cellcount)  = z(c, j)                 ! cell-center vertical depth from surface (m)

               ! top face area, scaled by active column weight and land fraction
               toparea_clm_loc(cellcount) = cwtgcell(c) * ldomain%frac(g) * larea(g) * 1.e6_r8       ! m^2


               ! after knowing 'toparea', we may get a pseudo 'dx' and 'dy' so that PF will not crash
               ! (note: PF needs these information otherwise throw-out error message, even with 'vertical_only' option)

               if (ldomain%nv == 4) then   !!wgs: ldomain%nv missing in ACME1
                  ! having 4 vertices
                  lats = latv(g,1:4)
                  lons = lonv(g,1:4)
                  dxsoil_clm_loc(cellcount) = abs(lons(1)+lons(4)-lons(2)-lons(3))/2.0_r8     &
                                               * cwtgcell(c) * ldomain%frac(g)
                    ! note: since in 'column_wise' mode, the columns are in 1D array by x-axis,
                    ! only need to scale 'dx' by column area fraction
                  dysoil_clm_loc(cellcount) = abs(lats(1)+lats(2)-lats(3)-lats(4))/2.0_r8

               else
                  dxsoil_clm_loc(cellcount) = larea(g)/(re**2)  &       ! in degrees of great circle length
                                               * cwtgcell(c) * ldomain%frac(g)
                  dysoil_clm_loc(cellcount) = larea(g)/(re**2)
               endif

            endif

        enddo ! do c = bounds%begc, bounds%endc

#else

        do g = bounds%begg, bounds%endg
            gcount = g - bounds%begg                           ! 0-based
            cellcount = gcount*clm_pf_idata%nzclm_mapped + j   ! 1-based

            cellid_clm_loc(cellcount) = (grc%gindex(g)-1)*clm_pf_idata%nzclm_mapped + j  ! 1-based

            xsoil_clm_loc(cellcount) = lonc(g)
            ysoil_clm_loc(cellcount) = latc(g)
            dxsoil_clm_loc(cellcount) = -9999.d0
            dysoil_clm_loc(cellcount) = -9999.d0

            !
            dzsoil_clm_loc(cellcount) = dz(xwtgcell_c(gcount+1), j)                ! cell vertical thickness (m), by column of max. weight in a grid
            zisoil_clm_loc(cellcount) = -zi(xwtgcell_c(gcount+1), j) + lelev(g)    ! cell-node elevation (m)
            zsoil_clm_loc(cellcount)  = z(xwtgcell_c(gcount+1), j)                 ! cell-center vertical depth from surface (m)

            ! top face area, scaled by active column weight (summed) and land fraction
            toparea_clm_loc(cellcount) = wtgcell_sum(gcount+1) * ldomain%frac(g) * larea(g) * 1.e6_r8       ! m^2

        enddo ! do g = bounds%begg, bounds%endg
#endif

      else
        call endrun(trim(subname) // ": ERROR: CLM-PF mapped soil layer numbers is greater than " // &
            " 'clm_varpar%nlevgrnd'. Please check")


      endif
    enddo

    call VecRestoreArrayF90(clm_pf_idata%cellid_clmp,  cellid_clm_loc,  ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%zisoil_clmp,  zisoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%dxsoil_clmp,  dxsoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%dysoil_clmp,  dysoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%dzsoil_clmp,  dzsoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clmp,  toparea_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%xsoil_clmp,  xsoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%ysoil_clmp,  ysoil_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%zsoil_clmp,  zsoil_clm_loc,  ierr)
    CHKERRQ(ierr)

    ! Set CLM soil domain onto PFLOTRAN grid
    call pflotranModelSetSoilDimension(pflotran_m)

    end associate
  end subroutine get_clm_soil_dimension

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_soil_properties
  !
  ! !INTERFACE:
  subroutine get_clm_soil_properties(clm_bgc_data, bounds, filters)
    !
    ! !DESCRIPTION:
    ! get soil column physical properties to PFLOTRAN
    !
    ! !USES:
    use LandunitType            , only : lun
    use ColumnType              , only : col
    use landunit_varcon         , only : istsoil, istcrop

    use clm_varpar              , only : nlevgrnd, ndecomp_pools, ndecomp_cascade_transitions
    use clm_varctl              , only : use_century_decomp
    use CNDecompCascadeConType  , only : decomp_cascade_con

    ! pflotran
    !
    ! !ARGUMENTS:

    implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscviewer.h"

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    ! Revised by Fengming Yuan, CCSI-ORNL
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    type(bounds_type), intent(in) :: bounds           ! bounds
    type(clumpfilter), intent(in) :: filters(:)       ! filters on current process

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    ! LOCAL VARAIBLES:

    integer  :: fc, c, l, g, j      ! indices
    integer  :: k,ki,kj
    integer  :: gcount, cellcount   ! gcount: 0-based, cellcount: 1-based
    integer  :: soilc1, layer1
    real(r8) :: CN_ratio_mass_to_mol
    real(r8) :: wtgcount

    character(len= 32) :: subname = 'get_clm_soil_properties' ! subroutine name


    PetscScalar, pointer :: hksat_x_clm_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_clm_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_clm_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_clm_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_clm_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_clm_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: watfc_clm_loc(:)
    PetscScalar, pointer :: bulkdensity_dry_clm_loc(:)
!     PetscScalar, pointer :: zsoi_clm_loc(:)

    PetscErrorCode :: ierr

    associate( &
         ! Assign local pointer to derived subtypes components (column-level)
         ltype                    => lun%itype                                , & !  [integer (:)]  landunit type index
         ! Assign local pointer to derived subtypes components (column-level)
         clandunit                => col%landunit                             , & !  [integer (:)]  landunit index of column
         cgridcell                => col%gridcell                             , & !  [integer (:)]  gridcell index of column
         cwtgcell                 => col%wtgcell                              , & !  [real(r8) (:)]  weight (relative to gridcell

         z                        => clm_bgc_data%z                           , & !  [real(r8) (:,:)]  layer depth (m)
         dz                       => clm_bgc_data%dz                          , & !  [real(r8) (:,:)]  layer thickness depth (m)
!          zi         =>  clm_bgc_data%zi         , & !  [real(r8) (:,:)]  interface level below a "z" level (m)
         !
         bd                       => clm_bgc_data%bd_col                      , & !
         bsw                      => clm_bgc_data%bsw_col                     , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat                    => clm_bgc_data%hksat_col                   , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat                   => clm_bgc_data%sucsat_col                  , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat                   => clm_bgc_data%watsat_col                  , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc                    => clm_bgc_data%watfc_col                   , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         
         rf_decomp_cascade        => clm_bgc_data%rf_decomp_cascade_col       , &
         pathfrac_decomp_cascade  => clm_bgc_data%pathfrac_decomp_cascade_col , &
         initial_cn_ratio         => clm_bgc_data%initial_cn_ratio            , &
         kd_decomp_pools          => clm_bgc_data%decomp_k_pools_col            &
         )

!-------------------------------------------------------------------------------------
    if(pflotran_m%option%ntrandof > 0) then
      ! the following assumes 'nclumps' in current process greater than 0 (at least 1)

      !! wgs:beg---------------------------------------------------------------------------
      !! in ACME, decomp_rate_constants have been call in 'CNEcosystemDynNoLeaching1' before pflotran-run
!       if (use_century_decomp) then
!         call decomp_rate_constants_bgc(bounds, filters(1)%num_soilc, filters(1)%soilc)  !!wgs: mismatch variables
! !         subroutine decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
! !        canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars)
!       else 
!         call decomp_rate_constants_cn(bounds, filters(1)%num_soilc, filters(1)%soilc)   !!wgs: mismatch variables
! !         subroutine decomp_rate_constants_cn(bounds, &
! !         num_soilc, filter_soilc, &
! !         canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
!       end if
      !! wgs:end---------------------------------------------------------------------------

      CN_ratio_mass_to_mol = clm_pf_idata%N_molecular_weight/clm_pf_idata%C_molecular_weight
      clm_pf_idata%decomp_element_ratios(:,1) = 1.0_r8
      clm_pf_idata%decomp_element_ratios(:,2) = 1.0_r8/initial_cn_ratio(1:ndecomp_pools) &
                                                /CN_ratio_mass_to_mol                       ! ratio in moles

      ! note: the following 'kd' already corrected by ad-factors for each pool
      clm_pf_idata%ck_decomp_c = kd_decomp_pools(1:ndecomp_pools)

      ! find the first active SOIL Column to pick up the decomposition constants
      ! NOTE: this only is good for CLM-CN reaction-network;
      !       for CLM-BGC (century-type), those constants are 'cell' dependent
      !      So, here we do the data passing by cell (although only 1 now) that can be extended by adding two loops in the future
      !
      soilc1 = filters(1)%soilc(1)
      layer1 = 1

      clm_pf_idata%fr_decomp_c = 0._r8
      do k = 1, ndecomp_cascade_transitions
        ki=decomp_cascade_con%cascade_donor_pool(k)
        kj=decomp_cascade_con%cascade_receiver_pool(k)

        if (ki>0) then
          ! taking the first 'cell' as default ('pathfrac' is 'cell'-related, which will be adjusted if needed)
          if (clm_pf_idata%fr_decomp_c(ki,ki) <=0._r8) then
            ! not-yet assign 'co2' fraction for donor-pool
            clm_pf_idata%fr_decomp_c(ki,ki) = rf_decomp_cascade(soilc1,layer1,k)             ! CO2-C respiration fraction
          elseif(clm_pf_idata%fr_decomp_c(ki,ki) .ne. rf_decomp_cascade(soilc1,layer1,k)) then
            ! have assigned 'co2' fraction for same donor-pool with different receive-pool,
            ! BUT 'co2' fraction inconsistent
            call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface finds different respiration fraction for " // &
             "same decomposition pool: " //trim(decomp_cascade_con%decomp_pool_name_history(ki)) )

          endif

          if (kj>0) then
            clm_pf_idata%fr_decomp_c(ki,kj) = (1.0_r8-rf_decomp_cascade(soilc1,layer1,k)) &
                         * pathfrac_decomp_cascade(soilc1,layer1,k)
          else
            if(clm_pf_idata%fr_decomp_c(ki,ki) .ne. 1.0) then
               ! if no receivor, respiration fraction must be 1.0
               call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface finds respiration fraction not 1.0 for " // &
                  "no-down decomposition pool: " //trim(decomp_cascade_con%decomp_pool_name_history(ki)) )
            endif

          endif

        endif

      enddo

      call pflotranModelSetSOMKfromCLM(pflotran_m)
    endif


    call VecGetArrayF90(clm_pf_idata%hksat_x_clmp, hksat_x_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_clmp, hksat_y_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_clmp, hksat_z_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_clmp,  sucsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clmp,  watsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_clmp,     bsw_clm_loc,     ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watfc_clmp,  watfc_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%bulkdensity_dry_clmp,  bulkdensity_dry_clm_loc,  ierr)
    CHKERRQ(ierr)

!     call VecGetArrayF90(clm_pf_idata%zsoi_clmp,  zsoi_clm_loc,  ierr)
!     CHKERRQ(ierr)

    hksat_x_clm_loc(:) = 0._r8
    hksat_y_clm_loc(:) = 0._r8
    hksat_z_clm_loc(:) = 0._r8
    sucsat_clm_loc(:)  = 0._r8
    watsat_clm_loc(:)  = 0._r8
    bsw_clm_loc(:)     = 0._r8
    watfc_clm_loc(:)   = 0._r8
    bulkdensity_dry_clm_loc(:) =  0._r8

    gcount = 0
    do c = bounds%begc, bounds%endc
      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)

      if (ltype(l)==istsoil .or. ltype(l)==istcrop) then

#ifdef COLUMN_MODE
        gcount = gcount + 1                                    ! 1-based column (fake grid) count
        wtgcount = 1._r8
#else
        gcount = g - bounds%begg + 1                           ! 1-based
        wtgcount = cwtgcell(c)
#endif

        do j = 1, clm_pf_idata%nzclm_mapped

          if (j <= nlevgrnd) then
            cellcount = (gcount - 1) *clm_pf_idata%nzclm_mapped + j   ! 1-based

            ! CLM calculation of wet thermal-conductivity as following:
            ! dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))*tkice**((1._r8-fl)*watsat(c,j))
            ! where, fl is the liq. saturation/total saturation
            ! so, if fl=0, it's the frozen-wet thermal-conductitivity
            !     if fl=1, it's the liq.-wet thermal-conductivity, i.e. 'tksatu'

!             tkwet_clm_loc(cellcount )  = &                      !(W/m/K)
!                 tkwet_clm_loc(cellcount ) + &
!                 tkmg(c,j)*(tkwat**watsat(c,j)) * wtgcount
!             tkdry_clm_loc(cellcount )  = &
!                 tkdry_clm_loc(cellcount ) + tkdry(c,j)*wtgcount
!             tkfrz_clm_loc(cellcount )  = &
!                 tkfrz_clm_loc(cellcount ) + &
!                 tkmg(c,j)*(tkice**watsat(c,j)) * wtgcount
!             hcvsol_clm_loc(cellcount ) = &
!                 hcvsol_clm_loc(cellcount ) + csol(c,j)*wtgcount  ! (J/m3/K)

            hksat_x_clm_loc(cellcount ) = &
                hksat_x_clm_loc(cellcount ) + hksat(c,j)*wtgcount
            hksat_y_clm_loc(cellcount ) = &
                hksat_y_clm_loc(cellcount ) + hksat(c,j)*wtgcount
            hksat_z_clm_loc(cellcount ) = &
                hksat_z_clm_loc(cellcount ) + hksat(c,j)*wtgcount

            sucsat_clm_loc( cellcount ) = &
                sucsat_clm_loc( cellcount ) + sucsat(c,j)*wtgcount
            watsat_clm_loc( cellcount ) = &
                watsat_clm_loc( cellcount ) + watsat(c,j)*wtgcount
            bsw_clm_loc(    cellcount ) = &
                bsw_clm_loc(    cellcount ) + bsw(c,j)*wtgcount
            watfc_clm_loc( cellcount )  = &
                watfc_clm_loc( cellcount ) + watfc(c,j)*wtgcount
            bulkdensity_dry_clm_loc( cellcount ) = &
                bulkdensity_dry_clm_loc( cellcount ) + bd(c,j)*wtgcount

          else
            call endrun(trim(subname) // ": ERROR: CLM-PF mapped soil layer numbers is greater than " // &
              " 'clm_varpar%nlevgrnd'. Please check")

          endif
        enddo
      endif

    enddo ! do c = bounds%begc, bounds%endc

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_clmp, hksat_x_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_clmp, hksat_y_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_clmp, hksat_z_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_clmp,  sucsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clmp,  watsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_clmp,     bsw_clm_loc,     ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%watfc_clmp,  watfc_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%bulkdensity_dry_clmp,  bulkdensity_dry_clm_loc,  ierr)
    CHKERRQ(ierr)
!     call VecRestoreArrayF90(clm_pf_idata%zsoi_clmp,  zsoi_clm_loc,  ierr)
!     CHKERRQ(ierr)

    ! Set CLM soil properties onto PFLOTRAN grid
    call pflotranModelSetSoilProp(pflotran_m)

    end associate
  end subroutine get_clm_soil_properties


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_soil_th
  !
  ! !INTERFACE:
  subroutine get_clm_soil_th(clm_bgc_data,initpftmode, initpfhmode, bounds, filters, ifilter)

  !
  ! !DESCRIPTION:
  !  update soil temperature/saturation from CLM to PFLOTRAN for driving PF's BGC
  !  if either NOT available inside PFLOTRAN
  !
  ! !USES:
    use clm_time_manager    , only : get_nstep, is_first_step, is_first_restart_step
    use shr_const_mod       , only : SHR_CONST_G
    use ColumnType          , only : col
    use clm_varctl          , only : iulog
    use clm_varcon          , only : denh2o, denice, tfrz
    use clm_varpar          , only : nlevgrnd
    use shr_infnan_mod      , only : shr_infnan_isnan

    use PFLOTRAN_Constants_module
    use clm_varctl          , only : pf_frzmode

  ! !ARGUMENTS:
    implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
    logical           , intent(in) :: initpftmode, initpfhmode
    type(bounds_type) , intent(in) :: bounds         ! bounds of current process
    type(clumpfilter) , intent(in) :: filters(:)     ! filters on current process
    integer           , intent(in) :: ifilter        ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: nstep                             ! time step number
    integer  :: fc, c, g, gcount, cellcount       ! indices

    PetscScalar, pointer :: soilpress_clmp_loc(:)
    PetscScalar, pointer :: soilpsi_clmp_loc(:)
    PetscScalar, pointer :: soillsat_clmp_loc(:)
    PetscScalar, pointer :: soilisat_clmp_loc(:)
    PetscScalar, pointer :: soilvwc_clmp_loc(:)  !
    PetscScalar, pointer :: soilt_clmp_loc(:)  !
    PetscScalar, pointer :: t_scalar_clmp_loc(:)  !
    PetscScalar, pointer :: w_scalar_clmp_loc(:)  !
    PetscScalar, pointer :: o_scalar_clmp_loc(:)  !
    PetscErrorCode :: ierr
    integer :: j
    real(r8):: sattmp, psitmp, itheta, sucmin_pa, psitmp0

    character(len= 32) :: subname = 'get_clm_soil_th' ! subroutine name

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      cgridcell       => col%gridcell               , & ! column's gridcell
      wtgcell         => col%wtgcell                , & ! column's weight relative to gridcell
      cactive         => col%active                 , & ! [logical (:)]  column active or not
      dz              => clm_bgc_data%dz            , & ! layer thickness depth (m)
    !
      sucsat          => clm_bgc_data%sucsat_col    , & ! minimum soil suction (mm) (nlevgrnd)
      bsw             => clm_bgc_data%bsw_col       , & ! Clapp and Hornberger "b"
      watsat          => clm_bgc_data%watsat_col    , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
      watmin          => clm_bgc_data%watmin_col    , & ! col minimum volumetric soil water (nlevsoi)
      sucmin          => clm_bgc_data%sucmin_col    , & ! col minimum allowable soil liquid suction pressure (mm) [Note: sucmin_col is a negative value, while sucsat_col is a positive quantity]
      soilpsi         => clm_bgc_data%soilpsi_col   , & ! soil water matric potential in each soil layer (MPa)
    !
      h2osoi_liq      => clm_bgc_data%h2osoi_liq_col, & ! liquid water (kg/m2)
      h2osoi_ice      => clm_bgc_data%h2osoi_ice_col, & ! ice lens (kg/m2)
      h2osoi_vol      => clm_bgc_data%h2osoi_vol_col, & ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
!       zwt             => clm_bgc_data%zwt          , & ! water table depth (m)
    !
      t_soisno        => clm_bgc_data%t_soisno_col  , & ! snow-soil temperature (Kelvin)
      t_scalar        => clm_bgc_data%t_scalar_col  , & ! soil temperature scalar for decomp
      w_scalar        => clm_bgc_data%w_scalar_col  , & ! soil water scalar for decomp
      o_scalar        => clm_bgc_data%o_scalar_col    & ! fraction by which decomposition is limited by anoxia
    )

    !--------------------------------------------------------------------------------------
    nstep   = get_nstep()

    call VecGetArrayF90(clm_pf_idata%press_clmp, soilpress_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilpsi_clmp, soilpsi_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilt_clmp, soilt_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%h2osoi_vol_clmp, soilvwc_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%t_scalar_clmp, t_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%w_scalar_clmp, w_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%o_scalar_clmp, o_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)


!     watmin(:,:) =  0.01_r8
!     sucmin(:,:) = -1.e8_r8
!>>>DEBUG | get_clm_soil_th--------------------------------------------------
!>>>DEBUG | watmin(:,:) & sucmin(:,:) are declared and initialized in 'SoilStateType.F90'
!>>>DEBUG | watmin=  0.328991E-01  0.351745E-01  0.402160E-01  0.568353E-01  0.984890E-01  0.164315E+00  0.165404E+00  0.414067E-01  0.339666E-01  0.328906E-01  0.328906E-01  0.328906E-01  0.328906E-01  0.328906E-01  0.328906E-01
!>>>DEBUG | sucmin= -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08 -0.101325E+08

    do fc = 1,filters(ifilter)%num_soilc
      c = filters(ifilter)%soilc(fc)
      g = cgridcell(c)
!      watmin(c,:) = 0.01_r8
!      sucmin(c,:) = 1.e8_r8

#ifdef COLUMN_MODE
      gcount = c - bounds%begc                 ! 0-based
      gcount = gcount + filters_colcount_beg(ifilter)
      ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
      ! when assigning values from 'filters(nc)%soilc(fc)'.
      ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clmp' vector.
#else
      gcount = g - bounds%begg                 ! 0-based
      if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

      do j = 1, clm_pf_idata%nzclm_mapped

        if (j<=nlevgrnd) then
          cellcount = gcount*clm_pf_idata%nzclm_mapped + j    ! 1-based

          if (initpfhmode) then
             ! this adjusting should be done first, if PF-freezing-mode off,
             ! so that the following calculation can be done correctly
             itheta = h2osoi_ice(c,j) / (dz(c,j) * denice)
!             itheta = min(itheta, watsat(c,j)-watmin(c,j)/(dz(c,j)*denh2o))
             itheta = min(itheta, watsat(c,j)-watmin(c,j))
             soilisat_clmp_loc(cellcount) = itheta/watsat(c,j)

             if(.not.pf_frzmode) then
                ! porosity will be ice-adjusted for PF, if PF freezing-mode is off,
                ! so need to adjust 'psi' so that 'saturation' in PF is correct
                sattmp = h2osoi_liq(c,j) / ((watsat(c,j)-itheta)*dz(c,j)*denh2o)
                sattmp = min(max(0.01d0, sattmp/(watsat(c,j)-itheta)),1._r8)

                ! soil matric potential re-done by Clapp-Hornburger method (this is the default used by CLM)
                ! this value IS different from what CLM used (not ice-content adjusted)
                ! So that in PF, if not ice-adjusted, the PSI is very small (negative) which implies possible water movement

                !! sucsat > 0, sucmin < 0, sucsat & sucmin have units of mm
                !! psitmp = psitmp0, as denh2o*1.e-3_r8 = 1.0
                psitmp0 = sucsat(c,j) * (-SHR_CONST_G) * (sattmp**(-bsw(c,j)))  ! Pa
!                psitmp = min(max(psitmp,-sucmin(c,j)/SHR_CONST_G),0._r8)
                psitmp = denh2o*(-SHR_CONST_G)*(sucsat(c,j)*1.e-3_r8)*(sattmp**(-bsw(c,j))) ! Pa
                sucmin_pa = denh2o*SHR_CONST_G*sucmin(c,j)*1.e-3_r8                         ! Pa
                psitmp = min(max(psitmp,sucmin_pa),0._r8)                                   !!Pa

             else
                sattmp = h2osoi_liq(c,j) / (watsat(c,j)*dz(c,j)*denh2o)
                sattmp = min(max(0.01d0, sattmp/watsat(c,j)),1._r8)

                psitmp = soilpsi(c,j)*1.e6_r8  ! MPa -> Pa
                if (shr_infnan_isnan(soilpsi(c,j)) .or. nstep<=0) then ! only for initialization, in which NOT assigned a value
                    psitmp0 = sucsat(c,j) * (-SHR_CONST_G) * ((sattmp+itheta)**(-bsw(c,j)))  ! -Pa: included both ice and liq. water as CLM does
!                    psitmp = min(max(psitmp,-sucmin(c,j)/SHR_CONST_G),0._r8)
                    psitmp = denh2o*(-SHR_CONST_G)*(sucsat(c,j)*1.e-3_r8) * ((sattmp+itheta)**(-bsw(c,j)))
                    sucmin_pa = denh2o*SHR_CONST_G*sucmin(c,j)*1.e-3_r8
                    psitmp = min(max(psitmp,sucmin_pa),0._r8)
                endif

             endif !!if(.not.pf_frzmode) 
             soillsat_clmp_loc(cellcount) = sattmp

             soilpsi_clmp_loc(cellcount)   = psitmp
             soilpress_clmp_loc(cellcount) = psitmp+clm_pf_idata%pressure_reference

             w_scalar_clmp_loc(cellcount)=w_scalar(c,j)
             o_scalar_clmp_loc(cellcount)=o_scalar(c,j)

             soilvwc_clmp_loc(cellcount) =h2osoi_vol(c,j)

          endif !!if (initpfhmode) then

          if (initpftmode) then
             soilt_clmp_loc(cellcount)=t_soisno(c,j)-tfrz
             t_scalar_clmp_loc(cellcount)=t_scalar(c,j)
          endif

#ifdef CLM_PF_OFFLINE
write(101,*) c, j, soilpress_clmp_loc(cellcount), soilt_clmp_loc(cellcount), soillsat_clmp_loc(cellcount)
#endif

        else 
            call endrun(trim(subname) // ": ERROR: CLM-PF mapped soil layer numbers is greater than " // &
              " 'clm_varpar%nlevgrnd'. Please check")

        endif !!if (j<=nlevgrnd) then

      enddo  !!do j = 1, clm_pf_idata%nzclm_mapped
    enddo !!do fc = 1,filters(ifilter)%num_soilc

!write(iulog,'(A,50(1h-))')">>>DEBUG | get_clm_soil_th"
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilpress=", soilpress_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilpsi[Pa]=", soilpsi_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soillsat=", soillsat_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilisat=", soilisat_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilt[oC]=", soilt_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilvwc=", soilvwc_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | t_scalar=", t_scalar_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | w_scalar=", w_scalar_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | o_scalar=", o_scalar_clmp_loc(1:10)

    call VecRestoreArrayF90(clm_pf_idata%press_clmp, soilpress_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilpsi_clmp, soilpsi_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilt_clmp, soilt_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%h2osoi_vol_clmp, soilvwc_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%t_scalar_clmp, t_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%w_scalar_clmp, w_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%o_scalar_clmp, o_scalar_clmp_loc, ierr)
    CHKERRQ(ierr)

   end associate


  end subroutine get_clm_soil_th


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_iceadj_porosity
  !
  ! !INTERFACE:
  subroutine get_clm_iceadj_porosity(clm_bgc_data, bounds, filters, ifilter)
  !
  ! !DESCRIPTION:
  !  update soil effective porosity from CLM to PFLOTRAN if PF freezing mode is off
  !
  ! !USES:
    use ColumnType          , only : col
    use clm_varctl          , only : iulog
    use clm_varcon          , only : denice
    use clm_varpar          , only : nlevgrnd

    use PFLOTRAN_Constants_module
    use clm_varctl          , only : pf_frzmode

  ! !ARGUMENTS:
    implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    type(bounds_type)         , intent(in) :: bounds         ! bounds
    type(clumpfilter)         , intent(in) :: filters(:)     ! filters on current process
    integer                   , intent(in) :: ifilter        ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, gcount, cellcount       ! indices
    real(r8) :: itheta

    PetscScalar, pointer :: adjporosity_clmp_loc(:)  !
    PetscScalar, pointer :: soilisat_clmp_loc(:)  !
    PetscErrorCode :: ierr

    character(len= 32) :: subname = 'get_clm_iceadj_porosity' ! subroutine name

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
    cgridcell       => col%gridcell                 , & ! column's gridcell
    dz              => clm_bgc_data%dz              , & ! layer thickness depth (m)
    watsat          => clm_bgc_data%watsat_col      , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
    h2osoi_ice      => clm_bgc_data%h2osoi_ice_col    & ! ice lens (kg/m2)
    )

    ! if 'pf_tmode' is NOT using freezing option, the phase-change of soil water done in 'SoilTemperatureMod.F90' in 'bgp2'
    ! must be included to adjust porosity (effective porosity) in pflotran
    ! This is doing prior to the real liquid water source/sink, because 'h2osoi_liq' will be updated during those calls after 'bgp2'.
    if (.not. pf_frzmode) then

        ! re-calculate the effective porosity (CLM ice-len adjusted), which should be pass to pflotran
        call VecGetArrayF90(clm_pf_idata%effporosity_clmp, adjporosity_clmp_loc,  ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
        CHKERRQ(ierr)

        do fc = 1,filters(ifilter)%num_soilc
           c = filters(ifilter)%soilc(fc)
           g = cgridcell(c)

#ifdef COLUMN_MODE
           gcount = c - bounds%begc                 ! 0-based
           gcount = gcount + filters_colcount_beg(ifilter)
           ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
           ! when assigning values from 'filters(nc)%soilc(fc)'.
           ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clmp' vector.
#else
           gcount = g - bounds%begg                 ! 0-based
           if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

           do j = 1, clm_pf_idata%nzclm_mapped
             if (j<=nlevgrnd) then
               cellcount = gcount*clm_pf_idata%nzclm_mapped + j    ! 1-based

               itheta = h2osoi_ice(c,j) / (dz(c,j) * denice)
               itheta = min(itheta, 0.99_r8*watsat(c,j))
               adjporosity_clmp_loc(cellcount ) = watsat(c,j) - itheta
               soilisat_clmp_loc(cellcount ) = itheta/watsat(c,j)
             else
               call endrun(trim(subname) // ": ERROR: CLM-PF mapped soil layer number is greater than " // &
                 " 'clm_varpar%nlevgrnd'. Please check")

             endif

           end do
        end do

!write(iulog,'(A,50(1h-))')">>>DEBUG | get_clm_iceadj_porosity"
!write(iulog,'(A30,12E14.6)')">>>DEBUG | adjporosity=", adjporosity_clmp_loc(1:10)
!write(iulog,'(A30,12E14.6)')">>>DEBUG | soilisat=", soilisat_clmp_loc(1:10)

        call VecRestoreArrayF90(clm_pf_idata%effporosity_clmp,  adjporosity_clmp_loc,  ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
        CHKERRQ(ierr)

    end if

    end associate
  end subroutine get_clm_iceadj_porosity



  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_bgc_conc(bounds)
  !
  ! !INTERFACE:

  subroutine get_clm_bgc_conc(clm_bgc_data, bounds, filters, ifilter)
!! TODO: add phosphorus vars
! #ifndef FLEXIBLE_POOLS
!     use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
! #endif
    use ColumnType          , only : col
    use clm_varctl          , only : iulog
    use clm_varpar          , only : ndecomp_pools, nlevdecomp

    implicit none

    type(bounds_type) , intent(in) :: bounds          ! bounds of current process
    type(clumpfilter) , intent(in) :: filters(:)      ! filters on current process
    integer           , intent(in) :: ifilter         ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_concentration"

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    ! Local variables
    integer  :: fc, c, g, j, k        ! do loop indices
    integer  :: gcount, cellcount
    real(r8) :: CN_ratio_mass_to_mol

! #ifdef FLEXIBLE_POOLS
    integer  :: vec_offset
    PetscScalar, pointer :: decomp_cpools_vr_clm_loc(:)      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_clm_loc(:)      ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
!     PetscScalar, pointer :: decomp_ppools_vr_clm_loc(:)      ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) P pools

    PetscScalar, pointer :: smin_no3_vr_clm_loc(:)           ! (gN/m3) vertically-resolved soil mineral NO3
    PetscScalar, pointer :: smin_nh4_vr_clm_loc(:)           ! (gN/m3) vertically-resolved soil mineral NH4
    PetscScalar, pointer :: smin_nh4sorb_vr_clm_loc(:)       ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

    PetscErrorCode :: ierr
    !
    !------------------------------------------------------------------------------------------
    !
    associate ( &
       initial_cn_ratio => clm_bgc_data%initial_cn_ratio        , &
       initial_cp_ratio => clm_bgc_data%initial_cp_ratio        , &

       decomp_cpools_vr => clm_bgc_data%decomp_cpools_vr_col    , & ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr => clm_bgc_data%decomp_npools_vr_col    , & ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       smin_no3_vr      => clm_bgc_data%smin_no3_vr_col         , & ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr      => clm_bgc_data%smin_nh4_vr_col         , & ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr  => clm_bgc_data%smin_nh4sorb_vr_col     , & ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

       decomp_ppools_vr => clm_bgc_data%decomp_ppools_vr_col    , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
       solutionp_vr     => clm_bgc_data%solutionp_vr_col        , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr       => clm_bgc_data%labilep_vr_col          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr       => clm_bgc_data%secondp_vr_col          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       occlp_vr         => clm_bgc_data%occlp_vr_col            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr         => clm_bgc_data%primp_vr_col            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
       sminp_vr         => clm_bgc_data%sminp_vr_col              & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
    )

    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_clmp, decomp_cpools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_clmp, decomp_npools_vr_clm_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%smin_no3_vr_clmp, smin_no3_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_clmp, smin_nh4_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4sorb_vr_clmp, smin_nh4sorb_vr_clm_loc, ierr)
    CHKERRQ(ierr)

    CN_ratio_mass_to_mol = clm_pf_idata%N_molecular_weight/clm_pf_idata%C_molecular_weight

     do fc = 1, filters(ifilter)%num_soilc
      c = filters(ifilter)%soilc(fc)
      g = col%gridcell(c)

#ifdef COLUMN_MODE
      gcount = c - bounds%begc             ! 0-based
      gcount = gcount + filters_colcount_beg(ifilter)
      ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
      ! when assigning values from 'filters(nc)%soilc(fc)'.
      ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clmp' vector.
#else
      gcount = g - bounds%begg                 ! 0-based
      if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

      decomp_cpools_vr_clm_loc(:) = 0._r8
      decomp_npools_vr_clm_loc(:) = 0._r8

      smin_no3_vr_clm_loc(:)      = 0._r8
      smin_nh4_vr_clm_loc(:)      = 0._r8
      smin_nh4sorb_vr_clm_loc(:)  = 0._r8

      do j = 1, clm_pf_idata%nzclm_mapped
          cellcount = gcount*clm_pf_idata%nzclm_mapped + j    ! 1-based

          ! note: all clm-pf soil layers are 'clm_pf_idata%nzclm_mapped' for both TH/BGC,
          !       but in CLM, T is within 'nlevgrnd', H is within 'nlevsoi', bgc within 'nlevdecomp'

          if(j <= nlevdecomp) then

             do k = 1, ndecomp_pools
                vec_offset = (k-1)*clm_pf_idata%nlclm_sub        ! 0-based
                ! decomp_pool vec: 'cell' first, then 'species' (i.e. cell by cell for 1 species, then species by species)
                ! Tips: then when doing 3-D data-mapping, no need to stride the vecs BUT to do segmentation.

                decomp_cpools_vr_clm_loc(vec_offset+cellcount) = decomp_cpools_vr(c,j,k)   &
                     /clm_pf_idata%C_molecular_weight

                if (clm_pf_idata%floating_cn_ratio(k)) then
                  decomp_npools_vr_clm_loc(vec_offset+cellcount) = decomp_npools_vr(c,j,k) &
                     /clm_pf_idata%N_molecular_weight
                else
                  decomp_npools_vr_clm_loc(vec_offset+cellcount) =   &
                      decomp_cpools_vr_clm_loc(vec_offset+cellcount) &
                     /(initial_cn_ratio(k)*CN_ratio_mass_to_mol)      ! initial_cn_ratio: in unit of mass
                endif

             enddo ! do k=1, ndecomp_pools

             smin_no3_vr_clm_loc(cellcount)      = smin_no3_vr(c,j) &
                        /clm_pf_idata%N_molecular_weight
             smin_nh4_vr_clm_loc(cellcount)      = smin_nh4_vr(c,j) &
                        /clm_pf_idata%N_molecular_weight
             smin_nh4sorb_vr_clm_loc(cellcount)  = smin_nh4sorb_vr(c,j) &
                        /clm_pf_idata%N_molecular_weight

           endif

       enddo ! do j = 1, clm_pf_idata%nzclm_mapped

    enddo ! do fc = 1, num_soilc

!-----------------------------------------------------------------------------
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_clmp, decomp_cpools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_clmp, decomp_npools_vr_clm_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_clmp, smin_no3_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_clmp, smin_nh4_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4sorb_vr_clmp, smin_nh4sorb_vr_clm_loc, ierr)
    CHKERRQ(ierr)

  end associate
  end subroutine get_clm_bgc_conc

  !-----------------------------------------------------------------------------
  !
  ! !IROUTINE: get_clm_bgc_rate()
  !
  ! !INTERFACE:
  subroutine get_clm_bgc_rate(clm_bgc_data,  bounds, filters, ifilter)
!! TODO: add phosphorus vars
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
    use ColumnType          , only : col
    use clm_time_manager    , only : get_step_size, get_nstep,  is_first_step, is_first_restart_step
    use clm_varpar          , only : ndecomp_pools, nlevdecomp
    use clm_varctl          , only : iulog, pf_hmode

  ! !ARGUMENTS:
    implicit none

    type(bounds_type) , intent(in) :: bounds         ! bounds of current process
    type(clumpfilter) , intent(in) :: filters(:)     ! filters on current process
    integer           , intent(in) :: ifilter        ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_rate"

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices
    integer  :: gcount, cellcount
!     real(r8) :: wtgcell, realc_gcell, realn_gcell

    real(r8) :: dtime                               ! land model time step (sec)

    ! ratios of NH4:NO3 in N deposition and fertilization (temporarily set here, will be as inputs)
   real(r8) :: r_nh4_no3_dep(bounds%begc:bounds%endc)
   real(r8) :: r_nh4_no3_fert(bounds%begc:bounds%endc)
   real(r8) :: fnh4_dep, fnh4_fert

    ! C/N source/sink rates as inputs for pflotran: Units - moles/m3/s (note: do unit conversion here for input rates)
! #ifdef FLEXIBLE_POOLS
    integer  :: vec_offset
    PetscScalar, pointer :: rate_decomp_c_clm_loc(:)       !
    PetscScalar, pointer :: rate_decomp_n_clm_loc(:)       !
!     PetscScalar, pointer :: rate_decomp_p_clm_loc(:)     !


    PetscScalar, pointer :: rate_plantndemand_clm_loc(:)   !
    PetscScalar, pointer :: rate_smin_no3_clm_loc(:)       !
    PetscScalar, pointer :: rate_smin_nh4_clm_loc(:)       !

    PetscErrorCode :: ierr

    !
    !---------------------------------------------------------------------------
    !
    associate ( &
      decomp_cpools_vr                  => clm_bgc_data%decomp_cpools_vr_col            , &      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
      decomp_npools_vr                  => clm_bgc_data%decomp_npools_vr_col            , &      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
      smin_no3_vr                       => clm_bgc_data%smin_no3_vr_col                 , &      ! (gN/m3) vertically-resolved soil mineral NO3
      smin_nh4_vr                       => clm_bgc_data%smin_nh4_vr_col                 , &      ! (gN/m3) vertically-resolved soil mineral NH4
!       smin_nh4sorb_vr => clm_bgc_data%smin_nh4sorb_vr_col      , &       ! (gN/m3) vertically-resolved soil mineral NH4 absorbed
      ! plant litering and removal + SOM/LIT vertical transport
      col_net_to_decomp_cpools_vr       => clm_bgc_data%externalc_to_decomp_cpools_col  , &
      col_net_to_decomp_npools_vr       => clm_bgc_data%externaln_to_decomp_npools_col  , &
      ! inorg. nitrogen source
      !! Inputs to NO3/NH4 is calculated in clm_bgc_interfaceMod.F90: subroutine get_clm_bgc_flux
!      ndep_to_sminn                  => nitrogenflux_vars%ndep_to_sminn_col                 , &
!      nfix_to_sminn                  => nitrogenflux_vars%nfix_to_sminn_col                 , &
!      fert_to_sminn                  => nitrogenflux_vars%fert_to_sminn_col                 , &
!      soyfixn_to_sminn               => nitrogenflux_vars%soyfixn_to_sminn_col              , &
!      supplement_to_sminn_vr         => nitrogenflux_vars%supplement_to_sminn_vr_col        , &
!      !
!      nfixation_prof                 => cnstate_vars%nfixation_prof_col                     , &
!      ndep_prof                      => cnstate_vars%ndep_prof_col                          , &
!      activeroot_prof                => cnstate_vars%activeroot_prof_col                    , &
      ! inorg. nitrogen sink (if not going to be done in PF)
      no3_net_transport_vr              => clm_bgc_data%no3_net_transport_vr_col        , &
      ! inorg. nitrogen sink potential
      col_plant_ndemand_vr              => clm_bgc_data%plant_ndemand_vr_col            , &

      externaln_to_nh4_vr               => clm_bgc_data%externaln_to_nh4_col            , &
      externaln_to_no3_vr               => clm_bgc_data%externaln_to_no3_col            , &

      col_net_to_decomp_ppools_vr       => clm_bgc_data%externalp_to_decomp_ppools_col  , &
      externalp_to_primp_vr             => clm_bgc_data%externalp_to_primp_col          , &
      externalp_to_labilep_vr           => clm_bgc_data%externalp_to_labilep_col        , &
      externalp_to_solutionp            => clm_bgc_data%externalp_to_solutionp_col      , &
      sminp_net_transport_vr            => clm_bgc_data%sminp_net_transport_vr_col      , &
      col_plant_pdemand_vr              => clm_bgc_data%plant_pdemand_vr_col              &
    )

    dtime = get_step_size()


    call VecGetArrayF90(clm_pf_idata%rate_decomp_c_clmp, rate_decomp_c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_decomp_n_clmp, rate_decomp_n_clm_loc, ierr)
    CHKERRQ(ierr)


    call VecGetArrayF90(clm_pf_idata%rate_plantndemand_clmp, rate_plantndemand_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_smin_no3_clmp, rate_smin_no3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_smin_nh4_clmp, rate_smin_nh4_clm_loc, ierr)
    CHKERRQ(ierr)

    ! Initialize to ZERO

    rate_decomp_c_clm_loc(:) = 0.0_r8
    rate_decomp_n_clm_loc(:) = 0.0_r8

    rate_smin_no3_clm_loc(:) = 0.0_r8
    rate_smin_nh4_clm_loc(:) = 0.0_r8
    rate_plantndemand_clm_loc(:) = 0.0_r8

!    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
!    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3

    do fc = 1,filters(ifilter)%num_soilc             ! only operating on soil column, although all columns (grids) are initialized above
      c = filters(ifilter)%soilc(fc)
      g = col%gridcell(c)

#ifdef COLUMN_MODE
      gcount = c - bounds%begc             ! 0-based
      gcount = gcount + filters_colcount_beg(ifilter)
      ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
      ! when assigning values from 'filters(nc)%soilc(fc)'.
      ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clmp' vector.
#else
      gcount = g - bounds%begg                 ! 0-based
      if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

      do j = 1, clm_pf_idata%nzclm_mapped
          cellcount = gcount*clm_pf_idata%nzclm_mapped + j    ! 1-based

          ! note: all clm-pf soil layers are 'clm_pf_idata%nzclm_mapped' for both TH/BGC,
          !       but in CLM, T is within 'nlevgrnd', H is within 'nlevsoi', bgc within 'nlevdecomp'

          if(j <= nlevdecomp) then
                !! just in case, we need to do some checking first (actually already done before)
              do k = 1,ndecomp_pools
                if (col_net_to_decomp_cpools_vr(c,j,k) < 0._r8) then
                  col_net_to_decomp_cpools_vr(c,j,k) =                  &
                    max(col_net_to_decomp_cpools_vr(c,j,k),             &
                        -max(decomp_cpools_vr(c,j,k)/dtime, 0._r8))
                endif

                if (col_net_to_decomp_npools_vr(c,j,k) < 0._r8) then
                  col_net_to_decomp_npools_vr(c,j,k) =                  &
                    max(col_net_to_decomp_npools_vr(c,j,k),             &
                        -max(decomp_npools_vr(c,j,k)/dtime, 0._r8))
                endif
              end do

              do k = 1, ndecomp_pools
                 vec_offset = (k-1)*clm_pf_idata%nlclm_sub        ! 0-based
                 ! decomp_pool vec: 'cell' first, then 'species' (i.e. cell by cell for 1 species, then species by species)
                 ! Tips: then when doing 3-D data-mapping, no need to stride the vecs BUT to do segmentation.

                 rate_decomp_c_clm_loc(vec_offset+cellcount)   =  &
                            col_net_to_decomp_cpools_vr(c,j,k)    &
                           /clm_pf_idata%C_molecular_weight

                 if (rate_decomp_c_clm_loc(vec_offset+cellcount)<0._r8) then
                   rate_decomp_c_clm_loc(vec_offset+cellcount) = max(     &
                     rate_decomp_c_clm_loc(vec_offset+cellcount),         &
                     -max(decomp_cpools_vr(c,j,k)/clm_pf_idata%C_molecular_weight/dtime, 0._r8))
                 endif

                 if (clm_pf_idata%floating_cn_ratio(k)) then
                   rate_decomp_n_clm_loc(vec_offset+cellcount) =  &
                           col_net_to_decomp_npools_vr(c,j,k)     &
                           /clm_pf_idata%N_molecular_weight

                   if (rate_decomp_n_clm_loc(vec_offset+cellcount)<0._r8) then
                     rate_decomp_n_clm_loc(vec_offset+cellcount) = max(     &
                       rate_decomp_n_clm_loc(vec_offset+cellcount),         &
                       -max(decomp_npools_vr(c,j,k)/clm_pf_idata%N_molecular_weight/dtime, 0._r8))
                   endif

                 endif

              enddo ! do k=1, ndecomp_pools

!! wgs:beg:-----------------------------------------------------------------
!! Inputs (externaln) to NO3/NH4 is calculated in clm_bgc_interfaceMod.F90: subroutine get_clm_bgc_flux
              rate_smin_nh4_clm_loc(cellcount) = externaln_to_nh4_vr(c,j)/ clm_pf_idata%N_molecular_weight
              if (rate_smin_nh4_clm_loc(cellcount)<0._r8) then
                 rate_smin_nh4_clm_loc(cellcount) = max(     &
                   rate_smin_nh4_clm_loc(cellcount),         &
                   -max(smin_nh4_vr(c,j)/clm_pf_idata%N_molecular_weight/dtime, 0._r8))
              endif
!              if (rate_smin_nh4_clm_loc(vec_offset+cellcount)<0._r8) then
!                 rate_smin_nh4_clm_loc(vec_offset+cellcount) = max(     &
!                   rate_smin_nh4_clm_loc(vec_offset+cellcount),         &
!                   -max(smin_nh4_vr(c,j)/clm_pf_idata%N_molecular_weight/dtime, 0._r8))
!              endif

              ! PF hydrological mode is OFF, then NO3 transport NOT to calculate in PF
              ! then it's done in CLM, so need to pass those to PF as source/sink term (RT mass transfer)
              rate_smin_no3_clm_loc(cellcount) = externaln_to_no3_vr(c,j)/ clm_pf_idata%N_molecular_weight 

!! wgs:2/27/2017: comment out:beg: do NOT consider NO3 leaching, which will be calculated after PFLOTRAN run
!              if(.not.pf_hmode) then
!                rate_smin_no3_clm_loc(cellcount) =              &
!                          rate_smin_no3_clm_loc(cellcount) -    &
!                          no3_net_transport_vr(c,j)             &
!                          / clm_pf_idata%N_molecular_weight
!              endif
              ! something unknown NOT right for NO3 leaching correction here as 'rt_mass_transfer' when 'use_century_decomp' is on

              if (rate_smin_no3_clm_loc(cellcount)<0._r8) then
                 rate_smin_no3_clm_loc(cellcount) = max(     &
                   rate_smin_no3_clm_loc(cellcount),         &
                   -max(smin_no3_vr(c,j)/clm_pf_idata%N_molecular_weight/dtime, 0._r8))
              endif

!              if (rate_smin_no3_clm_loc(vec_offset+cellcount)<0._r8) then
!                 rate_smin_no3_clm_loc(vec_offset+cellcount) = max(     &
!                   rate_smin_no3_clm_loc(vec_offset+cellcount),         &
!                   -max(smin_no3_vr(c,j)/clm_pf_idata%N_molecular_weight/dtime, 0._r8))
!              endif
!! wgs:2/27/2017: comment out:end
!! wgs:end-----------------------------------------------------------------
              ! plant N uptake rate here IS the N demand (potential uptake)
              rate_plantndemand_clm_loc(cellcount) = &
                         col_plant_ndemand_vr(c,j)/clm_pf_idata%N_molecular_weight

! #ifdef CLM_PF_DEBUG
!       write(pflotran_m%option%myrank+200,*) 'checking bgc-mass-rate - clm: ', &
!         'rank=',pflotran_m%option%myrank, 'column=',c, 'layer_id=',j, &
!         'rate_nh4_clm(layer_id)=',rate_smin_nh4_clm_loc(cellcount)
! #endif

          
          endif ! if (j<=nlevdecomp)

      enddo ! do j=1, clm_pf_idata%nzclm_mapped

    enddo ! do fc=1,numsoic

!-----------------------------------------------------------------------------
    call VecRestoreArrayF90(clm_pf_idata%rate_decomp_c_clmp, rate_decomp_c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_decomp_n_clmp, rate_decomp_n_clm_loc, ierr)
    CHKERRQ(ierr)


    call VecRestoreArrayF90(clm_pf_idata%rate_plantndemand_clmp, rate_plantndemand_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_smin_no3_clmp, rate_smin_no3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_smin_nh4_clmp, rate_smin_nh4_clm_loc, ierr)
    CHKERRQ(ierr)

    end associate
  end subroutine get_clm_bgc_rate


  ! =============================UPDATE PFLOTRAN evolving variables to CLM ==========================================




  !-----------------------------------------------------------------------------
  !
  ! !ROUTINE: update_soil_bgc_pf2clm()
  !
  ! !INTERFACE:
  !
  ! ! calculating BGC state variable changes over one time-step (rates)
  !   NOTE: Don't update the organic C/N state variables, which will be updated in those 'update' subroutines
  !           and the 'CNSoilLittVertTranspMod.F90' after 'update1'.
  !
  subroutine update_soil_bgc_pf2clm(clm_bgc_data, bounds, filters, ifilter)
!! TODO: add phosphorus vars
    use ColumnType              , only : col
    use clm_varctl              , only : iulog, use_ed
    use CNDecompCascadeConType  , only : decomp_cascade_con
    use clm_varpar              , only : ndecomp_pools, nlevdecomp, nlevdecomp_full
    use clm_varctl              , only : pf_hmode
    use clm_time_manager        , only : get_step_size,get_nstep, is_first_step, is_first_restart_step

    use clm_varcon              , only : dzsoi_decomp

    implicit none

    type(bounds_type) , intent(in) :: bounds         ! bounds of current process
    type(clumpfilter) , intent(in) :: filters(:)     ! filters on current process
    integer           , intent(in) :: ifilter        ! which filter to be operated

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_pf2clm"

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

    integer  :: fc,c,g,j,k,l
    integer  :: gcount, cellcount

    real(r8) :: dtime            ! land model time step (sec)
!!pflotran mass balance check:beg--------------------------------------------------
    integer  :: err_index                                                           ! indices
    logical  :: err_found                                                           ! error flag
    real(r8) :: pf_cinputs, pf_coutputs, pf_cdelta,pf_errcb, pf_cbeg, pf_cend
    real(r8) :: pf_ninputs, pf_noutputs, pf_ndelta,pf_errnb
    real(r8) :: pf_noutputs_gas, pf_noutputs_veg                                    !! _gas:nitrogen gases; _veg:plant uptake of NO3/NH4
    real(r8) :: pf_noutputs_nit, pf_noutputs_denit                                  !! _gas = _nit + _denit
    real(r8) :: pf_ninputs_org, pf_ninputs_min, pf_ndelta_org,pf_ndelta_min         !! _org:organic; _min:mineral nitrogen
    real(r8) :: pf_nbeg, pf_nbeg_org, pf_nbeg_min,pf_nend, pf_nend_org, pf_nend_min
    real(r8) :: pf_externaln_to_no3_col(bounds%begc:bounds%endc, 1:nlevdecomp_full)
    real(r8) :: pf_nbeg_nh4sorb,pf_nend_nh4sorb
    real(r8) :: pf_nbeg_no3,pf_nbeg_nh4
!!pflotran mass balance check:end--------------------------------------------------

    integer  :: vec_offset
    PetscScalar, pointer :: decomp_cpools_vr_clm_loc(:)      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_clm_loc(:)      ! (moleN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools

    PetscScalar, pointer :: smin_no3_vr_clm_loc(:)           ! (moleN/m3) vertically-resolved soil mineral NO3
    PetscScalar, pointer :: smin_nh4_vr_clm_loc(:)           ! (moleN/m3) vertically-resolved total soil mineral NH4
    PetscScalar, pointer :: smin_nh4sorb_vr_clm_loc(:)       ! (moleN/m3) vertically-resolved absorbed soil mineral NH4

    ! 'accextrn_vr' - accumulative (root) extracted N, i.e., actual plant N uptake from each soil layer, within a CLM timestep
    PetscScalar, pointer :: accextrnh4_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil mineral N root-extraction (accumulated)
    PetscScalar, pointer :: accextrno3_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil mineral N root-extraction (accumulated)

    ! 'acchr_vr' - accumulative CO2 proudction from decompositon (for tracking HR, not involving mass-balance)
!     PetscScalar, pointer :: acchr_vr_clm_loc(:)               ! (moleC/m3/timestep) vertically-resolved soil CO2 production (accumulated)

    ! 'accnmin_vr' - accumulative gross N mineralization within a CLM timestep
    PetscScalar, pointer :: accnmin_vr_clm_loc(:)               ! (moleN/m3/timestep) vertically-resolved soil N mineralization (accumulated)

    ! 'accnimm_vr' - accumulative N immobilization within a CLM timestep
    PetscScalar, pointer :: accnimmp_vr_clm_loc(:)              ! (moleN/m3/timestep) vertically-resolved soil N potential immoblilization (accumulated)
    PetscScalar, pointer :: accnimm_vr_clm_loc(:)               ! (moleN/m3/timestep) vertically-resolved soil N immoblilization (accumulated)

    PetscErrorCode :: ierr

!------------------------------------------------------------------------------------
     !
     associate ( &
     initial_cn_ratio             => clm_bgc_data%initial_cn_ratio             , &
     decomp_cpools_vr             => clm_bgc_data%decomp_cpools_vr_col         , &
     decomp_npools_vr             => clm_bgc_data%decomp_npools_vr_col         , &
     sminn_vr                     => clm_bgc_data%sminn_vr_col                 , &
     smin_no3_vr                  => clm_bgc_data%smin_no3_vr_col              , &
     smin_nh4_vr                  => clm_bgc_data%smin_nh4_vr_col              , &
     smin_nh4sorb_vr              => clm_bgc_data%smin_nh4sorb_vr_col          , &

     decomp_cpools_delta_vr       => clm_bgc_data%decomp_cpools_sourcesink_col , &
     decomp_npools_delta_vr       => clm_bgc_data%decomp_npools_sourcesink_col  , &

!      sminn_delta_vr               => cnf%sminn_delta_vr                    , &
!      smin_no3_delta_vr            => cnf%smin_no3_delta_vr                 , &
!      smin_nh4_delta_vr            => cnf%smin_nh4_delta_vr                 , &
!      smin_nh4sorb_delta_vr        => cnf%smin_nh4sorb_delta_vr             , &
!      hr_vr                        => ccf%hr_vr                             , &
     
     sminn_to_plant_vr            => clm_bgc_data%sminn_to_plant_vr_col         , &
     smin_no3_to_plant_vr         => clm_bgc_data%smin_no3_to_plant_vr_col      , &
     smin_nh4_to_plant_vr         => clm_bgc_data%smin_nh4_to_plant_vr_col      , &
     potential_immob_vr           => clm_bgc_data%potential_immob_vr_col        , &
     actual_immob_vr              => clm_bgc_data%actual_immob_vr_col           , &
     gross_nmin_vr                => clm_bgc_data%gross_nmin_vr_col               &
     )
! ------------------------------------------------------------------------
     dtime = get_step_size()
!!wgs:beg-----------------------------------------------------------------
!! pflotran mass blance check for Carbon & Nitrogen:
!! summary for last time-step

    do fc = 1,filters(ifilter)%num_soilc
        c = filters(ifilter)%soilc(fc)

        pf_cbeg         = 0._r8
        pf_nbeg_org     = 0._r8
        pf_nbeg_min     = 0._r8
        pf_nbeg         = 0._r8
        pf_nbeg_nh4sorb = 0._r8
        pf_nbeg_no3     = 0._r8
        pf_nbeg_nh4     = 0._r8

        do j = 1, nlevdecomp
            pf_nbeg_no3 = pf_nbeg_no3 + smin_no3_vr(c,j)*dzsoi_decomp(j)
            pf_nbeg_nh4 = pf_nbeg_nh4 + smin_nh4_vr(c,j)*dzsoi_decomp(j)
            pf_nbeg_nh4sorb = pf_nbeg_nh4sorb + smin_nh4sorb_vr(c,j)*dzsoi_decomp(j)
!!do NOT directly use sminn_vr(c,j), it does NOT always equal to (no3+nh4+nh4sorb) herein due to unknown reason
!!            pf_nbeg_min = pf_nbeg_min + sminn_vr(c,j)*dzsoi_decomp(j)
            do l = 1, ndecomp_pools
                pf_nbeg_org = pf_nbeg_org + decomp_npools_vr(c,j,l)*dzsoi_decomp(j)
                pf_cbeg     = pf_cbeg     + decomp_cpools_vr(c,j,l)*dzsoi_decomp(j)
            end do
        end do

        pf_nbeg_min = pf_nbeg_no3 + pf_nbeg_nh4 + pf_nbeg_nh4sorb
        pf_nbeg = pf_nbeg_org + pf_nbeg_min
    end do

!! do NOT consider NO3 leaching/runoff
    do fc = 1,filters(ifilter)%num_soilc
        c = filters(ifilter)%soilc(fc)
        do j = 1, nlevdecomp
            pf_externaln_to_no3_col(c,j) = clm_bgc_data%externaln_to_no3_col(c,j)
!            if(.not.pf_hmode) then
!                pf_externaln_to_no3_col(c,j) = pf_externaln_to_no3_col(c,j)  &
!                                             - clm_bgc_data%no3_net_transport_vr_col(c,j)
!                if (pf_externaln_to_no3_col(c,j)<0._r8) then
!                     pf_externaln_to_no3_col(c,j) = max(pf_externaln_to_no3_col(c,j), &
!                                                        -max(smin_no3_vr(c,j)/dtime, 0._r8))
!                endif
!            end if !!(.not.pf_hmode)

        end do
    end do
!!wgs:end-----------------------------------------------------------------

     ! soil C/N pool increments set to the previous timestep (i.e., not yet updated)
     decomp_cpools_delta_vr  = 0._r8-decomp_cpools_vr
     decomp_npools_delta_vr  = 0._r8-decomp_npools_vr

     ! clm-pf interface data updated
! FLEXIBLE_POOLS
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_clms, decomp_cpools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_clms, decomp_npools_vr_clm_loc, ierr)
     CHKERRQ(ierr)


     call VecGetArrayReadF90(clm_pf_idata%smin_no3_vr_clms, smin_no3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%smin_nh4_vr_clms, smin_nh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%smin_nh4sorb_vr_clms, smin_nh4sorb_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     call VecGetArrayReadF90(clm_pf_idata%accextrnh4_vr_clms, accextrnh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accextrno3_vr_clms, accextrno3_vr_clm_loc, ierr)
     CHKERRQ(ierr)

!      if(clm_pf_idata%ispec_hrimm>0) then
!         call VecGetArrayReadF90(clm_pf_idata%acctothr_vr_clms, acchr_vr_clm_loc, ierr)
!         CHKERRQ(ierr)
!      else
!         call VecGetArrayReadF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
!         CHKERRQ(ierr)
!      endif

     if(clm_pf_idata%ispec_nmin>0) then
        call VecGetArrayReadF90(clm_pf_idata%acctotnmin_vr_clms, accnmin_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecGetArrayReadF90(clm_pf_idata%accnmin_vr_clms, accnmin_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     if(clm_pf_idata%ispec_nimp>0) then
        call VecGetArrayReadF90(clm_pf_idata%acctotnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecGetArrayReadF90(clm_pf_idata%accnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     if(clm_pf_idata%ispec_nimm>0) then
        call VecGetArrayReadF90(clm_pf_idata%acctotnimm_vr_clms, accnimm_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecGetArrayReadF90(clm_pf_idata%accnimm_vr_clms, accnimm_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

    ! soil C/N pool increments, and actual plant N uptake, gross N mineralization and immobilization
    do fc = 1,filters(ifilter)%num_soilc          ! only operating on soil column, which then back to CLM-CN
      c = filters(ifilter)%soilc(fc)
      g = col%gridcell(c)

#ifdef COLUMN_MODE
      gcount = c - bounds%begc             ! 0-based
      gcount = gcount + filters_colcount_beg(ifilter)
      ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
      ! when assigning values from 'filters(nc)%soilc(fc)'.
      ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clms' vector.
#else
      gcount = g - bounds%begg             ! 0-based
      if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

      do j = 1, nlevdecomp

!           hr_vr(c,j)              = 0._r8
          gross_nmin_vr(c,j)      = 0._r8
          actual_immob_vr(c,j)    = 0._r8
          potential_immob_vr(c,j) = 0._r8

          if(j <= clm_pf_idata%nzclm_mapped) then

              cellcount = gcount*clm_pf_idata%nzclm_mapped + j   ! 1-based

              do k=1, ndecomp_pools

                 vec_offset = (k-1)*clm_pf_idata%ngclm_sub        ! 0-based
                 ! decomp_pool vec: 'cell' first, then 'species' (i.e. cell by cell for 1 species, then species by species)
                 ! Tips: then when doing 3-D data-mapping, no need to stride the vecs BUT to do segmentation.

                 decomp_cpools_delta_vr(c,j,k) = ( decomp_cpools_delta_vr(c,j,k)  &
                                + decomp_cpools_vr_clm_loc(vec_offset+cellcount)  &
                                * clm_pf_idata%C_molecular_weight ) !!/dtime  !!wgs: decomp_cpools_delta_vr=> clm_bgc_data%decomp_cpools_sourcesink_col

                 if (clm_pf_idata%floating_cn_ratio(k)) then
                     decomp_npools_delta_vr(c,j,k) = ( decomp_npools_delta_vr(c,j,k)  &
                                + decomp_npools_vr_clm_loc(vec_offset+cellcount)      &
                                * clm_pf_idata%N_molecular_weight ) !!/dtime
                 else
                     decomp_npools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k)/ &
                        initial_cn_ratio(k)      ! initial_cn_ratio: already in unit of mass
                 endif

                 if (abs(decomp_cpools_delta_vr(c,j,k))<=1.d-20) decomp_cpools_delta_vr(c,j,k)=0._r8
                 if (abs(decomp_npools_delta_vr(c,j,k))<=1.d-21) decomp_npools_delta_vr(c,j,k)=0._r8

!!wgs:beg:--------------------------------------------------------------------------------
!! hr_vr is moved to update_bgc_gaslosses_pf2clm()
                 ! tracking HR/NMIN/NIMM/NIMMP from SOM-C reaction network
                 ! if total 'hr/nmin/nimm/nimp' unknown, have to add up
!                  if (clm_pf_idata%ispec_decomp_hr(k)>0 .and. clm_pf_idata%ispec_hrimm<=0) then
!                     hr_vr(c,j) = hr_vr(c,j)                            &
!                           + (acchr_vr_clm_loc(vec_offset+cellcount)    &
!                           * clm_pf_idata%C_molecular_weight)/dtime
!                  endif
!!wgs:end:--------------------------------------------------------------------------------
                 !
                 if (clm_pf_idata%ispec_decomp_nmin(k)>0 .and. clm_pf_idata%ispec_nmin<=0) then
                    gross_nmin_vr(c,j)  = gross_nmin_vr(c,j)            &
                          + (accnmin_vr_clm_loc(vec_offset+cellcount)   &
                          * clm_pf_idata%N_molecular_weight)/dtime
                 endif

                 !
                 if (clm_pf_idata%ispec_decomp_nimm(k)>0 .and. clm_pf_idata%ispec_nimm<=0) then
                    actual_immob_vr(c,j) = actual_immob_vr(c,j)         &
                          + (accnimm_vr_clm_loc(vec_offset+cellcount)   &
                          * clm_pf_idata%N_molecular_weight)/dtime
                 endif

                 !
                 if (clm_pf_idata%ispec_decomp_nimp(k)>0 .and. clm_pf_idata%ispec_nimp<=0) then
                    potential_immob_vr(c,j) = potential_immob_vr(c,j)   &
                          + (accnimmp_vr_clm_loc(vec_offset+cellcount)  &
                          * clm_pf_idata%N_molecular_weight)/dtime
                 endif

              enddo ! do k=1, ndecomp_pools

              ! if total 'hr/nmin/nimm/nimp' known
!               if (clm_pf_idata%ispec_hrimm>0) then
!                  hr_vr(c,j) = hr_vr(c,j)                            &
!                           + (acchr_vr_clm_loc(cellcount)            &
!                           * clm_pf_idata%C_molecular_weight)/dtime
!               endif

              if (clm_pf_idata%ispec_nmin>0) then
                 gross_nmin_vr(c,j)  = gross_nmin_vr(c,j)           &
                          + (accnmin_vr_clm_loc(cellcount)          &
                          * clm_pf_idata%N_molecular_weight)/dtime
              endif
              if (clm_pf_idata%ispec_nimm>0) then
                 actual_immob_vr(c,j) = actual_immob_vr(c,j)        &
                          + (accnimm_vr_clm_loc(cellcount)          &
                          * clm_pf_idata%N_molecular_weight)/dtime
              endif
              if (clm_pf_idata%ispec_nimp>0) then
                 potential_immob_vr(c,j) = potential_immob_vr(c,j)  &
                          + (accnimmp_vr_clm_loc(cellcount)         &
                          * clm_pf_idata%N_molecular_weight)/dtime
              endif

!! wgs-beg:--------------------------------------------------------------------------------
!! directly update the 'smin' N pools (must bypass the 'CNNStateUpdate1,2,3' relevant to soil N)
              smin_no3_vr(c,j) = &
                           smin_no3_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              smin_nh4_vr(c,j) = &
                           smin_nh4_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              smin_nh4sorb_vr(c,j) = &
                           smin_nh4sorb_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              sminn_vr(c,j) = smin_no3_vr(c,j) + smin_nh4_vr(c,j) + smin_nh4sorb_vr(c,j)
!! wgs-end:--------------------------------------------------------------------------------

              !! flows or changes
              smin_nh4_to_plant_vr(c,j) = (accextrnh4_vr_clm_loc(cellcount)  &
                          * clm_pf_idata%N_molecular_weight)/dtime
              smin_no3_to_plant_vr(c,j) = (accextrno3_vr_clm_loc(cellcount)  &
                          * clm_pf_idata%N_molecular_weight)/dtime
              sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)

              gross_nmin_vr(c,j) = (accnmin_vr_clm_loc(cellcount)        &
                          * clm_pf_idata%N_molecular_weight)/dtime

              potential_immob_vr(c,j) = (accnimmp_vr_clm_loc(cellcount)  &
                          * clm_pf_idata%N_molecular_weight)/dtime
              actual_immob_vr(c,j)    = (accnimm_vr_clm_loc(cellcount)   &
                          * clm_pf_idata%N_molecular_weight)/dtime

          else    ! just in case 'clm_pf_idata%nzclm_mapped<nlevdcomp', all set to ZERO (different from TH)

              do k=1, ndecomp_pools
                 decomp_cpools_delta_vr(c,j,k) = 0._r8
                 decomp_npools_delta_vr(c,j,k) = 0._r8
              enddo

!               hr_vr(c,j)                 = 0._r8

!               sminn_delta_vr(c,j)        = 0._r8
!               smin_no3_delta_vr(c,j)     = 0._r8
!               smin_nh4_delta_vr(c,j)     = 0._r8
!               smin_nh4sorb_delta_vr(c,j) = 0._r8

              smin_nh4_to_plant_vr(c,j)  = 0._r8
              smin_no3_to_plant_vr(c,j)  = 0._r8
              sminn_to_plant_vr(c,j)     = 0._r8
              gross_nmin_vr(c,j)         = 0._r8
              potential_immob_vr(c,j)    = 0._r8
              actual_immob_vr(c,j)       = 0._r8

          endif

       enddo

     enddo ! do fc = 1, filters(ifilter)%num_soilc


     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_clms, decomp_cpools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_clms, decomp_npools_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     call VecRestoreArrayReadF90(clm_pf_idata%smin_no3_vr_clms, smin_no3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%smin_nh4_vr_clms, smin_nh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%smin_nh4sorb_vr_clms, smin_nh4sorb_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     call VecRestoreArrayReadF90(clm_pf_idata%accextrnh4_vr_clms, accextrnh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%accextrno3_vr_clms, accextrno3_vr_clm_loc, ierr)
     CHKERRQ(ierr)

!      if(clm_pf_idata%ispec_hrimm>0) then
!         call VecRestoreArrayReadF90(clm_pf_idata%acctothr_vr_clms, acchr_vr_clm_loc, ierr)
!         CHKERRQ(ierr)
!      else
!         call VecRestoreArrayReadF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
!         CHKERRQ(ierr)
!      endif

     if(clm_pf_idata%ispec_nmin>0) then
        call VecRestoreArrayReadF90(clm_pf_idata%acctotnmin_vr_clms, accnmin_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecRestoreArrayReadF90(clm_pf_idata%accnmin_vr_clms, accnmin_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     if(clm_pf_idata%ispec_nimp>0) then
        call VecRestoreArrayReadF90(clm_pf_idata%acctotnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecRestoreArrayReadF90(clm_pf_idata%accnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     if(clm_pf_idata%ispec_nimm>0) then
        call VecRestoreArrayReadF90(clm_pf_idata%acctotnimm_vr_clms, accnimm_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecRestoreArrayReadF90(clm_pf_idata%accnimm_vr_clms, accnimm_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     ! update bgc gas losses-------------------------------------------------
     call update_bgc_gaslosses_pf2clm(clm_bgc_data, bounds, filters, ifilter)

!!wgs:beg-----------------------------------------------------------------
!! pflotran mass blance check for Carbon & Nitrogen:
!! summary for current time-step & calculate mass balance
!! pflotran mass blance check-Carbon
    err_found = .false.
    do fc = 1,filters(ifilter)%num_soilc
        c = filters(ifilter)%soilc(fc)
        pf_cinputs = 0._r8
        pf_coutputs = 0._r8
        pf_cdelta = 0._r8

        do j = 1, nlevdecomp
            pf_coutputs = pf_coutputs + clm_bgc_data%hr_vr_col(c,j)*dzsoi_decomp(j)
            do l = 1, ndecomp_pools
                pf_cinputs = pf_cinputs + clm_bgc_data%externalc_to_decomp_cpools_col(c,j,l)*dzsoi_decomp(j)
                pf_cdelta  = pf_cdelta  + decomp_cpools_delta_vr(c,j,l)*dzsoi_decomp(j)
            end do
        end do

        pf_cend = pf_cbeg + pf_cdelta
        pf_errcb = (pf_cinputs - pf_coutputs)*dtime - pf_cdelta

        ! check for significant errors
        if (abs(pf_errcb) > 1e-8_r8) then
            err_found = .true.
            err_index = c
        end if
    end do

    if (.not. use_ed) then
         if (err_found) then
            write(iulog,'(A,100(1h-))')">>>--------PFLOTRAN Mass Balance Check:beg"
            write(iulog,'(A30,I5)')"C Balance Error in Column = ",err_index
            write(iulog,'(16A15)')"errcb", "C_in-out", "Cdelta","Cinputs","Coutputs","Cbeg","Cend"
            write(iulog,'(16E15.6)')pf_errcb, (pf_cinputs - pf_coutputs)*dtime, pf_cdelta, &
                                    pf_cinputs*dtime,pf_coutputs*dtime,pf_cbeg,pf_cend
            write(iulog,'(A,100(1h-))')">>>--------PFLOTRAN Mass Balance Check:end"
!           'endrun' is called in CNBalanceCheckMod.F90, as PFLOTRAN error will propagate into ALM
!            if((.not.is_first_step()) .and. (.not.is_first_restart_step())) then
!                call endrun(msg=errMsg(__FILE__, __LINE__))
!            end if
         end if
    end if !!(.not. use_ed)

!!-----------------------------------
!! pflotran mass blance check-Nitrogen
    do fc = 1,filters(ifilter)%num_soilc
        c = filters(ifilter)%soilc(fc)

        pf_nend_org         = 0._r8
        pf_nend_min         = 0._r8
        pf_nend             = 0._r8
        pf_nend_nh4sorb     = 0._r8

        pf_ninputs_org      = 0._r8
        pf_ninputs_min      = 0._r8
        pf_ninputs          = 0._r8

        pf_noutputs_nit     = 0._r8
        pf_noutputs_denit   = 0._r8
        pf_noutputs_gas     = 0._r8
        pf_noutputs_veg     = 0._r8
        pf_noutputs         = 0._r8

        pf_ndelta_org       = 0._r8
        pf_ndelta_min       = 0._r8
        pf_ndelta           = 0._r8

        do j = 1, nlevdecomp
            !! sminn_ver(c,j) has been calculated above
            pf_nend_min     = pf_nend_min     + sminn_vr(c,j)*dzsoi_decomp(j)
            pf_nend_nh4sorb = pf_nend_nh4sorb + smin_nh4sorb_vr(c,j)*dzsoi_decomp(j)
            pf_ninputs_min  = pf_ninputs_min  + clm_bgc_data%externaln_to_nh4_col(c,j)*dzsoi_decomp(j) &
                                              + pf_externaln_to_no3_col(c,j)*dzsoi_decomp(j)

            pf_noutputs_nit = pf_noutputs_nit + clm_bgc_data%f_ngas_decomp_vr_col(c,j)*dzsoi_decomp(j) &
                                              + clm_bgc_data%f_ngas_nitri_vr_col(c,j)*dzsoi_decomp(j)
            pf_noutputs_denit = pf_noutputs_denit + clm_bgc_data%f_ngas_denit_vr_col(c,j)*dzsoi_decomp(j)
            pf_noutputs_veg = pf_noutputs_veg + sminn_to_plant_vr(c,j)*dzsoi_decomp(j)
            do l = 1, ndecomp_pools
                pf_ndelta_org  = pf_ndelta_org  + decomp_npools_delta_vr(c,j,l)*dzsoi_decomp(j)
                pf_ninputs_org = pf_ninputs_org + clm_bgc_data%externaln_to_decomp_npools_col(c,j,l)*dzsoi_decomp(j)
            end do
        end do

        pf_nend_org     = pf_nbeg_org       + pf_ndelta_org   !!pf_ndelta_org has been calculated
        pf_nend         = pf_nend_org       + pf_nend_min
        pf_ndelta_min   = pf_nend_min       - pf_nbeg_min
        pf_ndelta       = pf_ndelta_org     + pf_ndelta_min
        pf_ninputs      = pf_ninputs_org    + pf_ninputs_min
        pf_noutputs_gas = pf_noutputs_nit   + pf_noutputs_denit
        pf_noutputs     = pf_noutputs_gas   + pf_noutputs_veg
        pf_errnb        = (pf_ninputs - pf_noutputs)*dtime - pf_ndelta
!        write(iulog,*)'>>>DEBUG | pflotran nbalance error = ', pf_errnb, c, get_nstep()
        ! check for significant errors
        if (abs(pf_errnb) > 1e-8_r8) then
            err_found = .true.
            err_index = c
        end if
    end do

    if (.not. use_ed) then
         if (err_found) then
            write(iulog,'(A,100(1h-))')">>>--------PFLOTRAN Mass Balance Check:beg"
            write(iulog,'(A30,I5)')"N Balance Error in Column = ",err_index
            write(iulog,'(20A15)')  "errnb", "N_in-out", "Ndelta",                          &
                                    "Ninputs","Noutputs", "Nbeg","Nend"
            write(iulog,'(20E15.6)')pf_errnb, (pf_ninputs - pf_noutputs)*dtime, pf_ndelta,  &
                                    pf_ninputs*dtime,pf_noutputs*dtime,pf_nbeg,pf_nend

            write(iulog,'(20A15)')  "Ninputs_org","Ninputs_min",                            &
                                    "Noutputs_nit","Noutputs_denit",                        &
                                    "Noutputs_gas","Noutputs_veg"

            write(iulog,'(20E15.6)')pf_ninputs_org*dtime,pf_ninputs_min*dtime,              &
                                    pf_noutputs_nit*dtime,pf_noutputs_denit*dtime,          &
                                    pf_noutputs_gas*dtime,pf_noutputs_veg*dtime

            write(iulog,'(20A15)')  "Ndelta_org","Nbeg_org","Nend_org",                     &
                                    "Ndelta_min","Nbeg_min","Nend_min",                     &
                                    "Ndelta_nh4sorb","Nbeg_nh4sorb","Nend_nh4sorb"
            write(iulog,'(20E15.6)')pf_ndelta_org,pf_nbeg_org,pf_nend_org,                  &
                                    pf_ndelta_min,pf_nbeg_min,pf_nend_min,                  &
                                    pf_nend_nh4sorb-pf_nbeg_nh4sorb,pf_nbeg_nh4sorb,pf_nend_nh4sorb
            write(iulog,'(A,100(1h-))')">>>--------PFLOTRAN Mass Balance Check:end"
!           'endrun' is called in CNBalanceCheckMod.F90, as PFLOTRAN error will propagate into ALM
!            if((.not.is_first_step()) .and. (.not.is_first_restart_step())) then
!                call endrun(msg=errMsg(__FILE__, __LINE__))
!            end if
        end if
    end if !!(.not. use_ed)
!!wgs:end-----------------------------------------------------------------
    end associate
  end subroutine update_soil_bgc_pf2clm


  !-----------------------------------------------------------------------------
  !
  ! !ROUTINE: update_bgc_gaslosses_pf2clm()
  !
  ! !INTERFACE:
  !
  ! This is a temporary solution to estimate pflotran bgc gaseous emission and transport loss
  ! from their aq. phase states
  ! (due to not yet available in pflotran bgc)
  !
  subroutine update_bgc_gaslosses_pf2clm(clm_bgc_data, bounds, filters, ifilter)

     use ColumnType         , only : col
     use clm_time_manager   , only : get_step_size, get_nstep
     use clm_varpar         , only : nlevdecomp
     use clm_varcon         , only : tfrz

     use clm_varctl         , only : pf_tmode, pf_hmode, pf_frzmode

     !
     implicit none

     type(bounds_type), intent(in) :: bounds         ! bounds of current process
     type(clumpfilter), intent(in) :: filters(:)     ! filters on current process
     integer          , intent(in) :: ifilter        ! which filter to be operated

     type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

     character(len=256) :: subname = "get_pf_bgc_gaslosses"

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

     integer  :: fc, c, g, j, k
     integer  :: gcount, cellcount
     real(r8) :: dtime            ! land model time step (sec)
     integer  :: nstep

     ! for gas species
     real(r8) :: tc, tk, total_p       ! temperature (oC, K), total air pressure (Pa)
     real(r8) :: co2_p, n2_p, n2o_p    ! partial pressure (Pa) of CO2, N2, N2O
     real(r8) :: cgas, cgas_p          ! mole-C(N)/m3(bulk soil)
     real(r8) :: air_vol, air_molar, wfps
     integer  :: lair_barrier(1:filters(ifilter)%num_soilc)      ! toppest soil layer that little air space for air flow into deep soil (-1: no, 0: ground, >0: soil layer)

     ! gases from PFLOTRAN are timely accumulated, so gas fluxes are calculated here if over atm. partial pressure (no explicit transport available from PF now)
     PetscScalar, pointer :: gco2_vr_clms_loc(:)               ! (M: molC/m3 bulk soil) vertically-resolved soil gas CO2 from PF's evolution
     PetscScalar, pointer :: gn2_vr_clms_loc(:)                ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2 from PF's evolution
     PetscScalar, pointer :: gn2o_vr_clms_loc(:)               ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2O from PF's evolution
     PetscScalar, pointer :: gco2_vr_clmp_loc(:)               ! (M: molC/m3 bulk soil) vertically-resolved soil gas CO2 to reset PF's CO2g
     PetscScalar, pointer :: gn2_vr_clmp_loc(:)                ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2 to reset PF's N2g
     PetscScalar, pointer :: gn2o_vr_clmp_loc(:)               ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2O to reset PF's N2Og

     ! 'acchr_vr' - accumulative CO2 proudction from decompositon (for tracking HR, not involving mass-balance)
     PetscScalar, pointer :: acchr_vr_clm_loc(:)               ! (moleC/m3/timestep) vertically-resolved soil CO2 production (accumulated)
     ! 'accngasmin_vr' - accumulative N gas proudction from mineralization (for tracking, not involving mass-balance)
     PetscScalar, pointer :: accngasmin_vr_clm_loc(:)          ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)
     ! 'accngasnitr_vr' - accumulative N gas proudction from nitrification (for tracking)
     PetscScalar, pointer :: accngasnitr_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)
     ! 'accngasdeni_vr' - accumulative N gas proudction from denitrification (for tracking)
     PetscScalar, pointer :: accngasdeni_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)

     !
     PetscScalar, pointer :: soillsat_clm_loc(:)
     PetscScalar, pointer :: soilisat_clm_loc(:)
     PetscScalar, pointer :: soilpor_clm_loc(:)
     PetscScalar, pointer :: soilt_clm_loc(:)
     PetscScalar, pointer :: soilpress_clm_loc(:)
     PetscErrorCode :: ierr

     real(r8), parameter :: rgas = 8.3144621       ! m3 Pa K-1 mol-1

!------------------------------------------------------------------------------------
    associate ( &
     forc_pco2                    => clm_bgc_data%forc_pco2_grc                 , & ! partial pressure co2 (Pa)
!      forc_pch4                    => clm_bgc_data%forc_pch4_grc                  , & ! partial pressure ch4 (Pa)
     forc_pbot                    => clm_bgc_data%forc_pbot_not_downscaled_grc  , & ! atmospheric pressure (Pa)
     frac_sno                     => clm_bgc_data%frac_sno_eff_col              , & ! fraction of ground covered by snow (0 to 1)
     frac_h2osfc                  => clm_bgc_data%frac_h2osfc_col               , & ! fraction of ground covered by surface water (0 to 1)
     dz                           => clm_bgc_data%dz                            , & ! soil layer thickness depth (m)
     hr_vr                        => clm_bgc_data%hr_vr_col                     , &
     f_co2_soil_vr                => clm_bgc_data%f_co2_soil_vr_col             , &
     f_n2o_soil_vr                => clm_bgc_data%f_n2o_soil_vr_col             , &
     f_n2_soil_vr                 => clm_bgc_data%f_n2_soil_vr_col              , &
     f_ngas_decomp_vr             => clm_bgc_data%f_ngas_decomp_vr_col          , &
     f_ngas_nitri_vr              => clm_bgc_data%f_ngas_nitri_vr_col           , &
     f_ngas_denit_vr              => clm_bgc_data%f_ngas_denit_vr_col             &
     )
! ------------------------------------------------------------------------
     dtime = get_step_size()
     nstep = get_nstep()

     ! get the current time-step state variables of aq. phase of interested species
     call VecGetArrayReadF90(clm_pf_idata%gco2_vr_clms, gco2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%gn2_vr_clms, gn2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%gn2o_vr_clms, gn2o_vr_clms_loc, ierr)
     CHKERRQ(ierr)

     call VecGetArrayF90(clm_pf_idata%gco2_vr_clmp, gco2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2_vr_clmp, gn2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2o_vr_clmp, gn2o_vr_clmp_loc, ierr)
     CHKERRQ(ierr)

     
     call VecGetArrayReadF90(clm_pf_idata%accngasmin_vr_clms, accngasmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accngasnitr_vr_clms, accngasnitr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accngasdeni_vr_clms, accngasdeni_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     if(clm_pf_idata%ispec_hrimm>0) then
        call VecGetArrayReadF90(clm_pf_idata%acctothr_vr_clms, acchr_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecGetArrayReadF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     ! env. variables to properties of gases
     if (pf_tmode) then
         call VecGetArrayReadF90(clm_pf_idata%soilt_clms, soilt_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soilt'
     else
         call VecGetArrayReadF90(clm_pf_idata%soilt_clmp, soilt_clm_loc, ierr)
         CHKERRQ(ierr)   ! CLM evolved 'soilt' - for CLM, MPI vecs and Seq. vecs should be same
     end if

     if (pf_frzmode) then
         call VecGetArrayReadF90(clm_pf_idata%soilisat_clms, soilisat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil ice saturation'
     end if

     if (pf_hmode) then
         call VecGetArrayReadF90(clm_pf_idata%soillsat_clms, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
         call VecGetArrayReadF90(clm_pf_idata%press_clms, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
     else
         call VecGetArrayReadF90(clm_pf_idata%soillsat_clmp, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
         call VecGetArrayReadF90(clm_pf_idata%press_clmp, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
     endif
     call VecGetArrayReadF90(clm_pf_idata%effporosity_clms, soilpor_clm_loc, ierr)  ! PF evolved 'soil porosity'
     CHKERRQ(ierr)

     ! find the toppest air barrier layer
     lair_barrier(:) = -1            ! (-1: no barrier, 0: ground snow/ice/water-layer barrier, >=1: barrier in soil column)
     do fc = 1,filters(ifilter)%num_soilc
       c = filters(ifilter)%soilc(fc)
       g = col%gridcell(c)

#ifdef COLUMN_MODE
       gcount = c - bounds%begc             ! 0-based
       gcount = gcount + filters_colcount_beg(ifilter)
      ! this means that skip prior clumps in 'clm_pf_idata%..._clms' vectors,
      ! when assigning values from 'filters(nc)%soilc(fc)'.
      ! NOTE: if 'nc' greater than 1, 'soilc(fc)' is a segement of '_clms' vector.
#else
       gcount = g - bounds%begg                ! 0-based
      if (mapped_gcount_skip(gcount+1)) cycle  ! skip inactive grid
#endif

      ! find the toppest air barrier layer for the current column at first

       if((frac_sno(c)+ frac_h2osfc(c))>=0.90_r8) then
         lair_barrier(c) = 0
       endif

       do j = 1, clm_pf_idata%nzclm_mapped

          cellcount = gcount*clm_pf_idata%nzclm_mapped+j

          wfps = 0._r8
          if (lair_barrier(c) >= 0) exit

          if (pf_frzmode) then
             wfps =soillsat_clm_loc(cellcount)  &
                    + soilisat_clm_loc(cellcount)
          else
             wfps =soillsat_clm_loc(cellcount)   ! note: 'lsat' from PF has been adjusted by 'isat' reduced porosity
          endif
          if (wfps > 0.95_r8) then             ! 95% total saturation as a critical-point for air-flow into deep soil
              lair_barrier(c) = j
          endif
       enddo
    

    ! gas exchanges btw atm. and non-barrierred soil layer
    ! only operating on soil column one by one, which then back to CLM-CN

       total_p = forc_pbot(g)

       do j = 1, nlevdecomp
          
          if(j <= clm_pf_idata%nzclm_mapped) then
              cellcount = gcount*clm_pf_idata%nzclm_mapped+j

              tc = soilt_clm_loc(cellcount)    ! soil layer tc (oC)
              tk = tc+tfrz

              ! total_p is soil air pressure (soil water pressure if not less than atm. pressure)
              total_p = max(total_p, soilpress_clm_loc(cellcount))

              ! the following is for adjusting air space in soil (seems not right ?? -- off)
              !air_vol = (1.0_r8 - soillsat_clm_loc(cellcount)) * &
              !                    soilpor_clm_loc(cellcount)                 ! m3 air/m3 soil
              !air_vol = max(air_vol, 0.01d0)            ! min. 0.01 to avoid math. issue

              air_vol = 1.0_r8                          ! atm used if not commented out
              air_molar = total_p*air_vol/rgas/tk       ! moles of air in a cell

              ! gas fluxes from  immobile PFLOTRAN evolving CO2imm, N2Oimm and N2imm, which are cumulative
              ! CO2 -
              cgas = gco2_vr_clms_loc(cellcount)       ! mol/m3 soil (evolving in PF, but not yet transport)
              co2_p = forc_pco2(g)                              ! assuming atm. pco2 (pa) as directly equilibrated with soil CO2(g)
              cgas_p = co2_p/forc_pbot(g) * air_molar

              f_co2_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_co2_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gco2_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gco2_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_co2_soil_vr(c,j) = f_co2_soil_vr(c,j)*clm_pf_idata%C_molecular_weight                  ! moleCO2/m3 --> gC/m3 soil
              f_co2_soil_vr(c,j) = f_co2_soil_vr(c,j)/dtime                               ! gC/m3/s

              ! N2
              cgas = gn2_vr_clms_loc(cellcount)       ! mol/m3 soil (evolving in PF, but not yet transport)
              n2_p = 0.78084_r8                                ! assuming atm. pn2 as directly equilibrated with soil n2(g)
              cgas_p = n2_p * air_molar           ! moleN2/m3

              f_n2_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_n2_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gn2_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gn2_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_n2_soil_vr(c,j) = f_n2_soil_vr(c,j)*clm_pf_idata%N_molecular_weight*2._r8    ! mole-N2/m3 --> g-N/m3 soil
              f_n2_soil_vr(c,j) = f_n2_soil_vr(c,j)/dtime                 ! gN/m3/s

              ! N2O
              cgas = gn2o_vr_clms_loc(cellcount)
              n2o_p = 310e-9_r8                   ! assuming general atm. pN2O (310ppbv in 1990) as directly equilibrated with soil N2(aq)
              cgas_p = n2o_p * air_molar          ! moleN2O/m3

              f_n2o_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_n2o_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gn2o_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gn2o_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_n2o_soil_vr(c,j) = f_n2o_soil_vr(c,j)*clm_pf_idata%N_molecular_weight*2._r8   ! mole-N2O/m3 --> g-N/m3 soil
              f_n2o_soil_vr(c,j) = f_n2o_soil_vr(c,j)/dtime                ! gN/m3/s

              ! tracking HR from SOM-C reaction network
              hr_vr(c,j)           = (acchr_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%C_molecular_weight)/dtime

              ! tracking gaseous N production from N reaction network
              f_ngas_decomp_vr(c,j)= (accngasmin_vr_clm_loc (cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

              f_ngas_nitri_vr(c,j) = (accngasnitr_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

              f_ngas_denit_vr(c,j) = (accngasdeni_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

          else    ! just in case 'clm_pf_idata%nzclm_mapped<nlevdecomp'

              f_co2_soil_vr(c,j)         = f_co2_soil_vr(c,clm_pf_idata%nzclm_mapped)
              f_n2_soil_vr(c,j)          = f_n2_soil_vr(c,clm_pf_idata%nzclm_mapped)
              f_n2o_soil_vr(c,j)         = f_n2o_soil_vr(c,clm_pf_idata%nzclm_mapped)

              hr_vr(c,j)                 = 0._r8
              f_ngas_decomp_vr(c,j)      = 0._r8
              f_ngas_nitri_vr(c,j)       = 0._r8
              f_ngas_denit_vr(c,j)       = 0._r8

          endif
       enddo !! do j = 1, nlevdecomp
     enddo !! do fc = 1,filters(ifilter)%num_soilc

     call VecRestoreArrayReadF90(clm_pf_idata%gco2_vr_clms, gco2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%gn2_vr_clms, gn2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%gn2o_vr_clms, gn2o_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gco2_vr_clmp, gco2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2_vr_clmp, gn2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_clmp, gn2o_vr_clmp_loc, ierr)
     CHKERRQ(ierr)

     if(clm_pf_idata%ispec_hrimm>0) then
        call VecRestoreArrayReadF90(clm_pf_idata%acctothr_vr_clms, acchr_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     else
        call VecRestoreArrayReadF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
        CHKERRQ(ierr)
     endif

     call VecRestoreArrayReadF90(clm_pf_idata%accngasmin_vr_clms, accngasmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%accngasnitr_vr_clms, accngasnitr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%accngasdeni_vr_clms, accngasdeni_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     if (pf_tmode) then
         call VecRestoreArrayReadF90(clm_pf_idata%soilt_clms, soilt_clm_loc, ierr)
         CHKERRQ(ierr)
     else
         call VecRestoreArrayReadF90(clm_pf_idata%soilt_clmp, soilt_clm_loc, ierr)
         CHKERRQ(ierr)  ! for CLM, MPI vecs and Seq. vecs should be same
     end if
     if (pf_frzmode) then
         call VecRestoreArrayReadF90(clm_pf_idata%soilisat_clms, soilisat_clm_loc, ierr)
         CHKERRQ(ierr)
     end if

     if (pf_hmode) then
         call VecRestoreArrayReadF90(clm_pf_idata%soillsat_clms, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
         call VecRestoreArrayReadF90(clm_pf_idata%press_clms, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
     else
         call VecRestoreArrayReadF90(clm_pf_idata%soillsat_clmp, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
         call VecRestoreArrayReadF90(clm_pf_idata%press_clmp, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
     endif
     call VecRestoreArrayReadF90(clm_pf_idata%effporosity_clms, soilpor_clm_loc, ierr)
     CHKERRQ(ierr)

     ! need to reset the PF's internal gas concentration (CLM ==> PF)
     call pflotranModelUpdateAqGasesfromCLM(pflotran_m)

    end associate
  end subroutine update_bgc_gaslosses_pf2clm


  !-----------------------------------------------------------------------------

#endif

!----------------------------------------------------------------------------------------------
 !!  END of CLM-PFLOTRAN bgc coupling interface
!----------------------------------------------------------------------------------------------

end module clm_pflotran_interfaceMod

