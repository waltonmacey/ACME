!===============================================================================
! SVN $Id: shr_const_mod.F90 6354 2007-09-11 22:49:33Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk/shr/shr_const_mod.F90 $
!===============================================================================

MODULE shr_mem_mod
    
   use shr_kind_mod, only : shr_kind_r8
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none
   private
    
! PUBLIC: Public interfaces

   public ::  shr_mem_getusage, &
	      shr_mem_init, print_memory_usage
    
! PUBLIC: Public interfaces

   real(shr_kind_r8) :: mb_blk = 0.0_shr_kind_r8

!===============================================================================
CONTAINS
!===============================================================================

subroutine shr_mem_init(prt)

   implicit none

   !----- arguments -----

   logical, optional :: prt
     
   !----- local -----

   ! --- Memory stats --- 
   integer :: msize                   ! memory size (high water)
   integer :: mrss                    ! resident size (current memory use)
   integer :: msize0,msize1           ! temporary size
   integer :: mrss0,mrss1,mrss2       ! temporary rss
   integer :: mshare,mtext,mdatastack
   logical :: lprt
   integer :: ierr
 
   integer :: GPTLget_memusage

   real(shr_kind_r8),allocatable :: mem_tmp(:)
    
   !---------------------------------------------------

   lprt = .true. ! ndk
   if (present(prt)) then
      lprt = prt
   endif

   ierr = GPTLget_memusage (msize, mrss0, mshare, mtext, mdatastack)
   allocate(mem_tmp(1024*1024))    ! 1 MWord, 8 MB
   mem_tmp = -1.0
   ierr = GPTLget_memusage (msize, mrss1, mshare, mtext, mdatastack)
   deallocate(mem_tmp)
   ierr = GPTLget_memusage (msize, mrss2, mshare, mtext, mdatastack)
   mb_blk = 0.0_shr_kind_r8
   if (mrss1 - mrss0 > 0) then
      mb_blk = (8.0_shr_kind_r8)/((mrss1-mrss0)*1.0_shr_kind_r8)
   endif

   if (lprt) then
      write(s_logunit,'(A,f16.2)') '8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
      write(s_logunit,'(A,f16.2)') '8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
      write(s_logunit,'(A,f16.2)') 'Memory block size conversion in bytes is ',mb_blk*1024_shr_kind_r8*1024.0_shr_kind_r8
   endif

end subroutine shr_mem_init

!===============================================================================

subroutine shr_mem_getusage(r_msize,r_mrss)

   implicit none

   !----- arguments ---
   real(shr_kind_r8) :: r_msize,r_mrss

   !----- local ---
   integer :: msize,mrss
   integer :: mshare,mtext,mdatastack
   integer :: ierr
   integer :: GPTLget_memusage

   !---------------------------------------------------

   ierr = GPTLget_memusage (msize, mrss, mshare, mtext, mdatastack)
   r_msize = msize*mb_blk
   r_mrss  = mrss*mb_blk

end subroutine shr_mem_getusage

!===============================================================================


subroutine print_memory_usage(cname,iname)  !ndk

  implicit none
  include "mpif.h" 
  character(len=*),intent(in) :: cname
  integer,         intent(in) :: iname

  integer::rank,nprocs,ierr
  real*8 :: valueRSS, valuePeak, tvalueRSS, tvaluePeak
  real*8 :: minvalueRSS,  maxvalueRSS,  sumvalueRSS
  real*8 :: minvaluePeak, maxvaluePeak, sumvaluePeak
  real*8 :: k0=1000.0d0 
  real*8 :: k1=1024.0d0 

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  call  system_mem_usage(tvalueRSS, tvaluePeak)
  !write(*,'(a,i4, a,i4,a,f10.4,a,f10.4)') "  i=", i, " rank ", rank, &
  !     " tvalueRSS=", tvalueRSS, " tPeak=", tvaluePeak

  call MPI_reduce(tvalueRSS,  maxvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_MAX, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvalueRSS,  minvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_MIN, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvalueRSS,  sumvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, maxvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, minvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, sumvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierr)

  if(rank==0) then
     write(*,'(a,a25,i8,a,f8.2,f8.2,f8.2,f8.2,a,f8.2,f8.2,f8.2,f8.2)') &
          "PMU ", cname, iname, &
          " RSS=", sumvalueRSS/real(nprocs), minvalueRSS,maxvalueRSS, sumvalueRSS, &
          " peak=", sumvaluePeak/real(nprocs), minvaluePeak, maxvaluePeak, sumvaluePeak
  endif

  !-----------------------------------------------------------------
  call shr_mem_getusage(tvaluePeak, tvalueRSS) ! also call the ACME version
  tvaluePeak=tvaluePeak/(k0)
  tvalueRSS=tvalueRSS/(k0)

  call MPI_reduce(tvalueRSS,  maxvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_MAX, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvalueRSS,  minvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_MIN, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvalueRSS,  sumvalueRSS,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, maxvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, minvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0,MPI_COMM_WORLD, ierr)
  call MPI_reduce(tvaluePeak, sumvaluePeak, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierr)

  if(rank==0) then
     write(*,'(a,a25,i8,a,f8.2,f8.2,f8.2,f8.2,a,f8.2,f8.2,f8.2,f8.2)') &
          "PMA ", cname, iname, &
          " RSS=", sumvalueRSS/real(nprocs), minvalueRSS,maxvalueRSS, sumvalueRSS, &
          " peak=", sumvaluePeak/real(nprocs), minvaluePeak, maxvaluePeak, sumvaluePeak
  endif

end subroutine print_memory_usage

subroutine system_mem_usage(valueRSS, valuePeak)

  !use ifport !if on intel compiler
  implicit none

  real*8, intent(inout) :: valueRSS, valuePeak
  
  real*8 :: k0=1000.0d0 
  real*8 :: k1=1024.0d0 

  character(len=200):: filename=' '
  character(len=80) :: line
  character(len=8)  :: pid_char=' '
  integer :: pid, ivalueRSS, ivaluePeak
  logical :: ifxst
  
  ivalueRSS=-1    ! return negative number if not found
  ivaluePeak=-1    ! return negative number if not found
  
  !--- get process ID
  !pid=getpid()
  !write(pid_char,'(I8)') pid
  !filename='/proc/'//trim(adjustl(pid_char))//'/status'
  filename='/proc/self/status'
  
  inquire (file=filename,exist=ifxst)
  if (.not.ifxst) then
     write (*,*) 'system file does not exist'
     return
  endif
  
  open(unit=100, file=filename, action='read')
  do
     read (100,'(a)',end=120) line
     if (line(1:6).eq.'VmRSS:') then
        read (line(7:),*) ivalueRSS
        !exit
     endif
     if (line(1:7).eq.'VmPeak:') then
        read (line(8:),*) ivaluePeak
        !exit
     endif
  enddo
120 continue

  close(100)

  valueRSS = ivalueRSS/(k0*k0)
  valuePeak = ivaluePeak/(k0*k0)

  return
end subroutine system_mem_usage

END MODULE shr_mem_mod
