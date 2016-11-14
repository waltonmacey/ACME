!
! List of valid output variables and methods to access them
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define use_netcdf_interp_mod

module nc_fields

  use dimensions_mod, only: nelemd, np, ne, nc, nlev, qsize_d
  use element_mod,    only: element_t
  use kinds,          only: rl => real_kind
  use pio_io_mod,     only: pio_double
  use shr_const_mod,  only: unset => shr_const_spval
  use vertical_se,    only: v_interpolate

  implicit none

  integer, parameter :: n_dims = 4                                      ! number of dimensions
  character*(*),parameter :: dim_names(n_dims)=(/'lon ','lat ','lev ','time'/)

  !_____________________________________________________________________
  ! scalar field descriptor

  type nc_sfield_t
    integer         :: var_type           ! pio type of variable
    integer         :: dim_ids(n_dims)    ! dimension ids
    integer         :: tensor_dim         ! tensor dim:0=scalar,1=vector,2=tensor
    logical         :: is_required        ! variable required flag
    character*(5)   :: short_name         ! short name of variable
    character*(256) :: units              ! physical units
    character*(256) :: long_name          ! full name of variable
  endtype

  !_____________________________________________________________________
  ! vector field descriptor

  type nc_vfield_t
    character*(4)   :: vec_name           ! vector name
    integer         :: n_comp             ! number of components
    character*(4)   :: short_name(2)      ! name of each component
  endtype

  !_____________________________________________________________________
  !  possible output variables

  integer, parameter :: n_vars= 40 ! number of possible output fields

  type(nc_sfield_t), parameter :: nc_vars(n_vars) =&
  (/&
    nc_sfield_t(pio_double, (/1,0,0,0/),0, .true. , 'lon  ', 'degrees_east',    'column longitude'           ),&
    nc_sfield_t(pio_double, (/2,0,0,0/),0, .true. , 'lat  ', 'degrees_north',   'column latitude'            ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'lev  ', '',                'vertical coordinate'        ),&
    nc_sfield_t(pio_double, (/4,0,0,0/),0, .true. , 'time ', 'days',            'model elapsed time'         ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'hyam ', '',                'hybrid A at midpoints'      ),&
    nc_sfield_t(pio_double, (/3,0,0,0/),0, .true. , 'hybm ', '',                'hybrid B at midpoints'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'ps   ', 'Pa',              'hydrostatic surf pressure'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'z    ', 'meters',          'altitiude'                  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'geo  ', 'meters^2/sec^2',  'geopotential height'        ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'phi  ', 'meters^2/sec^2',  'geopotential height'        ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'phis ', 'meters^2/sec^2',  'surface geopotential'       ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'p    ', 'Pa',              'hydrostatic atmo pressure'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'pt   ', 'Pa',              'total atmospheric pressure' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'T    ', 'degrees Kelvin',  'temperature'                ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Tp   ', 'degrees Kelvin',  'potential temperature'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),1, .false., 'u    ', 'meters/second',   'longitudinal wind component'),&
    nc_sfield_t(pio_double, (/1,2,3,4/),1, .false., 'v    ', 'meters/second',   'latitudinal wind component' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'w    ', 'meters/second',   'vertical wind component'    ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'pp   ', 'Pa',              'pressure deviation'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q    ', 'kg/kg',           'tracer 1 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q2   ', 'kg/kg',           'tracer 2 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q3   ', 'kg/kg',           'tracer 3 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'Q4   ', 'kg/kg',           'tracer 4 mixing ratio'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'omega', '',                ''                           ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'D3   ', '',                '3d divergence of velocity'  ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'speed', 'meters/second',   'wind speed'                 ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'spd_ ', 'meters/second',   'speed - speed0'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'KE   ', 'Joules',          'kinetic energy'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'KE_  ', 'Joules',          'KE - KE0'                   ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'alpha', '',                'specific volume'            ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'rho  ', '',                'density'                    ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'dpdn ', '',                'pseudo-density'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'm    ', '',                'pseudo-density'             ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'cflux', '',                'column int horiz mass flux' ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'ndot ', '',                'vertical eta velocity'      ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'nflux', '',                'vertical eta flux'          ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O1   ', '',                'debugging output 1'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O2   ', '',                'debugging output 2'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O3   ', '',                'debugging output 3'         ),&
    nc_sfield_t(pio_double, (/1,2,3,4/),0, .false., 'O4   ', '',                'debugging output 4'         ) &
    !nc_sfield_t(type,       dim_ids,    n_spatial_dims, required?, short_name, units, long_name)
  /)

  type(nc_vfield_t),parameter::nc_vfields(1)=&                          ! list vector fields
  (/&
    !nc_vfield_t(vec_name, n_comps, (/component_names/))
    nc_vfield_t ('v'     , 2      , (/'u   ','v   '/)) &
  /)

contains

  !_____________________________________________________________________________
  function get_scalar_field(element, short_name, n0, ni) result(sfield)

    ! get 3d scalar-field data by name

    character*(*),    intent(in)  :: short_name
    type(element_t),  intent(in)  :: element
    integer,          intent(in)  :: n0       ! time level index
    integer,          intent(in)  :: ni       ! num vertical interp pts

    real(rl) :: sfield(np,np,ni)
    real(rl) :: var(np,np,nlev),var2d(np,np)
    integer i,j,qi

    if(ni==1) then

      ! get 2d scalar field by name
      select case(short_name)
        case('ps'  ); var2d = element%state%ps_v(:,:,n0)
        case default; var2d = unset
      end select
      sfield(:,:,1) = var2d
    endif

    if(ni>1) then

      ! get 3d scalar field by name
      select case(short_name)

        ! get prognostic variables
        case('T'  ); var = element%state%T(:,:,:,n0)
        case('u'  ); var = element%state%v(:,:,1,:,n0)
        case('v'  ); var = element%state%v(:,:,2,:,n0)
        case('p'  ); var = element%state%dp3d(:,:,:,n0)
        case('Q'  ); var = element%state%Q(:,:,:,1)
        case('Q2' ); qi=2; if(qsize_d>1) var = element%state%Q(:,:,:,qi)
        case('Q3' ); qi=3; if(qsize_d>2) var = element%state%Q(:,:,:,qi)
        case('Q4' ); qi=4; if(qsize_d>3) var = element%state%Q(:,:,:,qi)
        case('Q5' ); qi=5; if(qsize_d>4) var = element%state%Q(:,:,:,qi)

        case('omega');var = element%state%omega

        ! get diagnostic variables
        case('geo'); var = element%derived%phi(:,:,:)

        case default; var = unset                                      ! assign special "missing" value
      endselect

      ! interpolate each column
      do i=1,np
        do j=1,np
          sfield(i,j,:) = v_interpolate(var(i,j,:),ni)
        enddo
      enddo
    endif

  end function

  !_____________________________________________________________________________
  function get_vector_field(element, short_name, n0, ni) result(vfield)

    ! Get vector-field data by name

    character*(*),    intent(in)  :: short_name
    type(element_t),  intent(in)  :: element
    integer,          intent(in)  :: n0       ! time level index
    integer,          intent(in)  :: ni       ! num vertical interp pts

    real(rl) :: vfield(np,np,2,ni)
    real(rl) :: var(np,np,2,nlev)
    integer  :: i,j
    select case(short_name)                                             ! switch on short variable name

      case('v');
        var = element%state%v(:,:,:,:,n0)

      case default; var = unset                                         ! assign special "missing" value
    endselect

    ! interpolate each column
    do i=1,np
      do j=1,np
        vfield(i,j,1,:) = v_interpolate(var(i,j,1,:),ni)
        vfield(i,j,2,:) = v_interpolate(var(i,j,2,:),ni)
      enddo
    enddo

  end function

  end module nc_fields
