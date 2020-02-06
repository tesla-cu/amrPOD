! =============================================================================
! This Fortran code is to time the computation of POD on a set of binary files 
! containing the flowfield snapshots. The workflow of the program is:
!
!    - Allocate variables
!    - Initialize computation
!    - Error checking
!    - Read data
!    - Perform POD:                         Matrix syntax:
!       - Compute covariance matrix            - R     = X^T*X
!       - Perform eigenvalue decomposition     - R*Psi = Psi*Lambda
!       - Compute spatial/POD modes            - Phi   = X*Psi*Lambda^(-1/2) 
!       - Compute temporal coefficients        - A     = X^T * Phi
!    - Finalize computation
!
! This program user defined inputs from "POD.inputs" to specify the domain 
! dimensions, location of data, format of data, number of samples of CPU time, 
! variable to compute POD on, start and end of time indices, and the finest 
! resolution. This program computes the CPU time for each of the POD operations
! and the reshaping time in the AMR algorithm. 
!
! To create an executable, at the terminal, run:
!
!    gfortran -O0 -g -fcheck=all -fbacktrace -o POD.ex POD.f90 -L.. -llapack -lrefblas
!    gfortran -O3 -o POD.ex POD.f90 -L.. -llapack -lrefblas 
!
! For more timing info:
!    gfortran -pg -O3 -o POD.ex POD.f90 -L. -llapack -lrefblas 
!    gprof POD.ex > POD.stats
! 
! =============================================================================
! Begin 
! =============================================================================

include 'rshp_AMR.f90'
include 'comp_R.f90'
include 'comp_Phi.f90'
include 'comp_A.f90'

program POD

use comp_R_omp
use comp_A_omp
use comp_Phi_omp
use rshp_AMR

implicit none

! =============================================================================
! Allocate variables for POD 
! =============================================================================

! Parameters ------------------------------------------------------------------
real,    parameter :: PI = 4.0*atan(1.0)          ! pi=3.1415 ...
integer, parameter :: maxlen=256                  ! max length of string
character(len=20), parameter :: CFMT = "(F12.8)"  ! format of CPU time
character(len=20), parameter :: PFMT = "(F12.10)" ! format of p and p^m

! General ---------------------------------------------------------------------
integer               :: i, j, k, l ! spatial counting variables
integer               :: m, n       ! time counting variables
logical               :: chk        ! check variable 
integer               :: fid        ! file ID
character(len=maxlen) :: filename   ! file name
character(len=5)      :: nt_id      ! string version of nt

! Data ------------------------------------------------------------------------
integer               :: nx, ny, nz ! domain dimensions
integer               :: nspat      ! number of spatial points
integer               :: nt         ! number of time steps
character(len=maxlen) :: datadir    ! directory where data is stored
character(len=maxlen) :: datafmt    ! format of data
character(len=maxlen) :: CPUdir     ! directory to store CPU times

! Standard POD ----------------------------------------------------------------
character(len=maxlen)                         :: var   ! solution value variable
integer,                       dimension(2)   :: itime ! temporal indices (start, end)
double precision, allocatable, dimension(:,:) :: Xpod  ! snapshot matrix
double precision, allocatable, dimension(:,:) :: Rpod  ! covariance matrix
double precision, allocatable, dimension(:,:) :: Psi   ! eigenvectors
double precision, allocatable, dimension(:,:) :: Phi   ! POD modes
double precision, allocatable, dimension(:,:) :: Apod  ! temporal coefficients
double precision, allocatable, dimension(:)   :: Lam   ! eigenvalues
double precision, allocatable, dimension(:)   :: opt_work ! for dgesy
integer                                       :: eig_info ! for dgesy
integer                                       :: L_work   ! for dgesy

! AMR POD ---------------------------------------------------------------------
logical                                       :: do_amr  ! do AMR algorithm?
integer                                       :: finest  ! finest level of AMR
integer                                       :: ndim    ! dimension for POD
integer,          allocatable, dimension(:,:) :: Xgrid   ! companion matrix
double precision, allocatable, dimension(:)   :: Xcol    ! column of X
integer,          allocatable, dimension(:)   :: Xgcol   ! column of Xgrid
double precision, allocatable, dimension(:)   :: p_comp  ! computed p
double precision, allocatable, dimension(:)   :: pm_comp ! computed p^m

! Timing ----------------------------------------------------------------------
integer                            :: nsamp      ! number of samples of CPU
real                               :: Rshp_CPU   ! CPU time for reshaping
real,    allocatable, dimension(:) :: R_CPU      ! CPU time for R
real,    allocatable, dimension(:) :: Phi_CPU    ! CPU time for Phi
real,    allocatable, dimension(:) :: A_CPU      ! CPU time for A
integer                            :: cr         ! clock rate
double precision                   :: crr        ! real clock rate
integer                            :: c0, c1, c2 ! timers
integer, allocatable, dimension(:) :: cpu0, cpuf ! init and final CPU times
integer                            :: hrs, mins  ! hours, minutes
double precision                   :: secs       ! seconds
character(len=10)                  :: clock_time ! time of day

! MPI 

! Namelists - groups of variables to be read from inputs
namelist /POD_inputs/ nx, ny, nz, nsamp, itime, datafmt, datadir, var
namelist /AMR_POD_inputs/ finest

! =============================================================================
! Initialize computation
! =============================================================================

! Set up timing
call date_and_time(time=clock_time)
call system_clock(count_rate=cr)
crr = real(cr)
call system_clock(c0)

write(*, '(a)') '------------------------------------------------'//           &
                '------------------------------------------------'
write(*, '(a)') 'POD computation started at '//                                &
               clock_time(1:2)//':'//clock_time(3:4)//':'//                    &
               clock_time(5:10)//' local time.'
write(*, '(a)') '------------------------------------------------'//           &
                '------------------------------------------------'

! Set defaults for inputs
nx      = 1
ny      = 1
nz      = 1
nsamp   = 1
itime   = [0, 99999]
datafmt = '(a,a,I5.5,a)'
datadir = '/folder/'
var     = 'z_velocity'

! Read inputs
open(9, file='POD.inputs', form='formatted', status='old')
read(9, nml=POD_inputs)
close(9)

! Compute important quantities
ndim = 0
if (nx > 1) ndim = ndim + 1
if (ny > 1) ndim = ndim + 1
if (nz > 1) ndim = ndim + 1
nspat = nx*ny*nz
nt    = itime(2) - itime(1) + 1
write(nt_id, '(I5)') nt

! Allocate variables for standard POD
allocate(Xpod(nspat,nt))
allocate(Psi(nt,nt))
allocate(Lam(nt))
allocate(Rpod(nt,nt))
allocate(Phi(nspat,nt))
allocate(Apod(nt,nt))

! Allocate and assign variables for fast POD
do_amr = .true.
if (do_amr) then
   finest = 0

   open(9, file='POD.inputs', form='formatted', status='old')
   read(9, nml=AMR_POD_inputs)
   close(9)

   allocate(Xgrid(nspat,nt))
   allocate(Xcol(nspat))
   allocate(Xgcol(nspat))
end if

! Set up CPU timing
allocate(cpu0(nsamp), cpuf(nsamp))
allocate(R_cpu(nsamp), Phi_cpu(nsamp), A_cpu(nsamp))
CPUdir = trim(datadir)//'CPU_timing_nt'//trim(adjustl(nt_id))//'/'
call execute_command_line('mkdir -p '//CPUdir, wait=.true.)

! =============================================================================
! Error checking
! =============================================================================

! Check if all variable data exists in the directory
do j=1,nt
   write(filename,datafmt) trim(datadir), trim(var), itime(1)+j-1, '.bin'
   inquire(file=filename, exist=chk)
   if (.not.chk) then
      write(*,'(a)') "the file "//trim(filename)//" does not exist"
      stop
   endif
enddo

! Check if grid level data exists in the directory
if (do_amr) then
   do j=1,nt
      write(filename,datafmt) trim(datadir), 'grid_level', itime(1)+j-1, '.bin'
      inquire(file=filename, exist=chk)
      if (.not.chk) then
         write(*,'(a)') "the file "//trim(filename)//" does not exist"
         stop
      endif
   enddo
end if

! =============================================================================
! Read Data 
! =============================================================================

! Initialize CPU for reshaping
if (do_amr) then
   Rshp_CPU = 0.
endif

! Load data for each snapshot
do n=1,nt

   ! Read variable data
   write(filename,datafmt) trim(datadir), trim(var), itime(1)+n-1, '.bin'
   open(newunit=fid, file=trim(filename), action='read', access='stream',      &
      form='unformatted', status='old')
   read(fid) Xcol
   Xpod(:,n) = Xcol
   close(fid)

   ! If we are utilizing the AMR, we need to load grid data and reshape data
   if (do_amr) then
      ! Reshape original data
      Xcol = Xpod(:,n)
      call system_clock(cpu0(1))
      call reshape_AMR(nx, ny, nz, finest, Xcol, 'forward')
      call system_clock(cpuf(1))
      Rshp_CPU = Rshp_CPU + real(cpuf(1) - cpu0(1))/crr
      Xpod(:,n) = Xcol

      ! Load grid level data
      write(filename,datafmt) trim(datadir), 'grid_level', itime(1)+n-1, '.bin'
      open(newunit=fid, file=trim(filename), action='read', access='stream',   &
         form='unformatted', status='old')
      read(fid) Xcol
      close(fid)

      ! Reshape grid data
      call system_clock(cpu0(1))
      call reshape_AMR(nx, ny, nz, finest, Xcol, 'forward')
      call system_clock(cpuf(1))
      Rshp_CPU = Rshp_CPU + real(cpuf(1) - cpu0(1))/crr
      Xgrid(:,n) = int(Xcol)
   end if

enddo

! Write out CPU for reshaping
open(newunit=fid, file=trim(CPUdir)//'Reshape_CPU.txt', form='formatted')
write(*, "(a,"//CFMT//",a)") "reshaping cpu time ", Rshp_CPU, " seconds"
write(fid, CFMT) Rshp_CPU
close(fid)

! Compute p^m and p
if (do_amr) then

   allocate(p_comp(0:finest), pm_comp(0:finest))
   p_comp = 0.
   pm_comp = 0.

   ! Compute p
   do n=1,nt
      do l=0,finest
         p_comp(l) = p_comp(l) + dble(count(mask=Xgrid(:,n)==l))/dble(nt*nspat)
      enddo
   enddo
   write(*,*) "p_comp = ",p_comp
   write(*,"(a,"//PFMT//")") "sum(p_comp) = ",sum(p_comp)

   ! Compute p^m
   do i=1,nspat
      pm_comp(maxval(Xgrid(i,:))) = pm_comp(maxval(Xgrid(i,:))) + 1./dble(nspat)
   end do
   write(*,*) "pm_comp = ",pm_comp
   write(*,"(a,"//PFMT//")") "sum(pm_comp) = ",sum(pm_comp)

   ! Write out p to txt file
   open(newunit=fid, file=trim(CPUdir)//'p.txt', form='formatted')
   do i=0,finest
      write(fid, PFMT) p_comp(i)
   enddo
   close(fid)

   ! Write out p^m to txt file
   open(newunit=fid, file=trim(CPUdir)//'pm.txt', form='formatted')
   do i=0,finest
      write(fid, PFMT) pm_comp(i)
   enddo
   close(fid)

   deallocate(p_comp, pm_comp)
endif

write(*,*)

! =============================================================================
! Compute POD
! =============================================================================

! Compute R -------------------------------------------------------------------

! R with standard algorithm
write(*,*) "computing R using standard operations ..."
call system_clock(cpu0(1))
call compute_R_omp(Xpod, nspat, nt, Rpod, 0)
call system_clock(cpuf(1))

! Write out CPU time for R
! open(newunit=fid, file=trim(CPUdir)//'R_CPU_standard.txt', form='formatted')
! R_CPU = real(cpuf - cpu0)/crr
!    write(*,*) "    cpu time ", R_CPU(i), " seconds"
!    write(fid, CFMT) R_CPU(i)
! enddo
! close(fid)


! R with AMR algorithm
write(*,*) "computing R utilizing AMR ..."
call system_clock(cpu0(1))
call compute_R_omp(Xpod, nspat, nt, Rpod, 1, Xgrid, finest, ndim)
call system_clock(cpuf(1))

! Write out CPU time for R
! open(newunit=fid, file=trim(CPUdir)//'R_CPU_AMR.txt', form='formatted')

! do i=1,nsamp
!    write(*,*) "    cpu time ", R_CPU(i), " seconds"
!    write(fid, CFMT) R_CPU(i)
! enddo
! close(fid)

write(*,*)

! Compute Psi and Lambda ------------------------------------------------------

! Arguments for dsyev (double symmetric eigenvalue):
!    'V'      - compute evals and evecs
!    'U'      - use upper triangular mat
!    nt       - order of matrix
!    Rpod     - matrix to compute decomposition, outputs evecs
!    nt       - leading order of R
!    Lam      - eigenvalues
!    opt_work - 
!    L_work   - length of opt_work
!    eig_info - information if errors

L_work = 3*nt
allocate(opt_work(L_work))
call dsyev('V', 'U', nt, Rpod, nt, Lam, opt_work, L_work, eig_info)
Psi = Rpod

! Compute Phi -----------------------------------------------------------------

! Phi with standard algorithm
write(*,*) "computing Phi using standard algorithm ..."
call system_clock(cpu0(1))
call compute_Phi_omp(Xpod, Psi, Lam, nspat, nt, Phi, 0)
call system_clock(cpuf(1))

! Write out CPU time for Phi
! open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_standard.txt', form='formatted')
! Phi_CPU = real(cpuf - cpu0)/crr
! do i=1,nsamp
!    write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
!    write(fid, CFMT) Phi_CPU(i)
! enddo
! close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo

! Phi with AMR algorithm, method 1
write(*,*) "computing Phi utilizing AMR, method 1 ..."
call system_clock(cpu0(1))
call compute_Phi_omp(Xpod, Psi, Lam, nspat, nt, Phi, 1, Xgrid, finest, ndim)
call system_clock(cpuf(1))

! Write out CPU time for Phi
! open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_AMR1.txt', form='formatted')
! Phi_CPU = real(cpuf - cpu0)/crr
! do i=1,nsamp
!    write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
!    write(fid, CFMT) Phi_CPU(i)
! enddo
! close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo

! Phi with AMR algorithm, method 1
write(*,*) "computing Phi utilizing AMR, method 2 ..."
call system_clock(cpu0(1))
call compute_Phi_omp(Xpod, Psi, Lam, nspat, nt, Phi, 2, Xgrid, finest, ndim)
call system_clock(cpuf(1))

! Write out CPU time for Phi
! open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_AMR2.txt', form='formatted')
! Phi_CPU = real(cpuf - cpu0)/crr
! do i=1,nsamp
!    write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
!    write(fid, CFMT) Phi_CPU(i)
! enddo
! close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo

write(*,*)

! Compute A -------------------------------------------------------------------

! A with standard algorithm
write(*,*) "computing A using standard algorithm ..."
call system_clock(cpu0(1))
call compute_A_omp(Xpod, Phi, nspat, nt, Apod, 0)
call system_clock(cpuf(1))

! Write out CPU time for A
! open(newunit=fid, file=trim(CPUdir)//'A_CPU_standard.txt', form='formatted')
! A_CPU = real(cpuf - cpu0)/crr
! do i=1,nsamp
!    write(*,*) "    cpu time ", A_CPU(i), " seconds"
!    write(fid, CFMT) A_CPU(i)
! enddo
! close(fid)

! do i=1,4
!    write(*,*) Apod(i,1:4)
! enddo

! A with AMR algorithm
write(*,*) "computing A utilizing AMR ..."
call system_clock(cpu0(1))
call compute_A_omp(Xpod, Phi, nspat, nt, Apod, 1, Xgrid, finest, ndim)
call system_clock(cpuf(1))

! Write out CPU time for A
! open(newunit=fid, file=trim(CPUdir)//'A_CPU_AMR.txt', form='formatted')
! A_CPU = real(cpuf - cpu0)/crr
! do i=1,nsamp
!    write(*,*) "    cpu time ", A_CPU(i), " seconds"
!    write(fid, CFMT) A_CPU(i)
! enddo
! close(fid)

! do i=1,4
!    write(*,*) Apod(i,1:4)
! enddo

! =============================================================================
! Finalize computation
! =============================================================================

call date_and_time(time=clock_time)
call system_clock(c2)
c1 = (c2 - c0)/cr
hrs = int(c1/3600)
mins = int(mod(c1, 3600)/60)
secs = real(c1 - hrs*3600 - mins*60) &
         + real(mod(c2-c0, cr))/crr

write(*, '(a)') '------------------------------------------------'//           &
                '------------------------------------------------'

write(*, '(a)') 'POD computation finished at '//                               &
               clock_time(1:2)//':'//clock_time(3:4)//':'//                    &
               clock_time(5:10)//' local time.'//new_line('a')

write(*, '(A,I3.2,A,I2.2,A,I2.2,F0.3,A)')                                      &
   'Total program walltime = ', hrs, ':', mins, ':', int(secs),                &
   secs-int(secs), ' (HHH:MM:SS.SSS)'

write(*, '(a)') '------------------------------------------------'//           &
                '------------------------------------------------'

end program POD
! =============================================================================