!
! namelist for dcmip2016 test 2 tropical cyclone
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test2"         ! test identifier
  ne                = 7                        ! number of elements per cube face (2dg)
  qsize             = 3                         ! num tracer fields
  ndays             = 10                        ! num simulation days: 0 = use nmax steps
  statefreq         = 5 !10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 = new run
  tstep             = 10                        ! largest timestep
  tstep_type        = 3                         ! 1 => default method
  integration       = 'explicit'                ! explicit time integration
  smooth            = 0                         ! timestep smooting (nonzero smoothing breaks this test)
  nu                = 1.0e16                    ! reduced earth hyperviz
  nu_s              = 1.0e16
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rsplit            = 1
/
&filter_nl/
&solver_nl
  precon_method     = "identity"
  maxits            = 50
  tol               = 1.e-7
/
&vert_nl
  vform             = "ccm"                     ! vertical coordinate type "ccm"=hybrid pressure/terrain
  vanalytic         = 1                         ! set vcoords in initialization routine
  vtop              = 2.73919e-1                ! vertical coordinate at top of atm (z=10000m)
/
&analysis_nl
  output_dir        = "../movies/"              ! destination dir for netcdf file
  output_timeunits  = 0 !   !2,                   ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 50 !  !3,                   ! 100 sec / 0.5 sec per step
  output_varnames1  ='ps','v','omega','Q','Q2','Q3'  !'T' ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
  interp_nlat       = 128
  interp_nlon       = 256
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
