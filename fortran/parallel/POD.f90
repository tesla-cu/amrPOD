! =============================================================================
! This Fortran code is to time the computation of POD on a set of binary files 
! containing the flowfield snapshots. The workflow of the program is:
!
!    - Allocate variables
!    - Initialize computation
!    - Error checking
!    - Read data
!    - Perform POD:                         Matlab syntax:
!       - Compute covariance matrix            - R     = X^T*X
!       - Perform eigenvalue decomposition     - R*Psi = Psi*Lambda
!       - Compute spatial/POD modes            - Phi   = X*Psi*Lambda^(-1/2) 
!       - Compute temporal coefficients        - A     = X^T * Phi
!    - Write data
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
!    mpifort -O0 -g -fopenmp -fcheck=all -fbacktrace -o POD.ex POD.f90 -L.. -llapack -lrefblas
!    mpifort -O3 -fopenmp -o POD.ex POD.f90 -L.. -llapack -lrefblas 
! 
! =============================================================================
! Begin 
! =============================================================================

include 'rshp_AMR.f90'
include 'comp_R_omp.f90'
include 'comp_Phi_omp.f90'
include 'comp_A_omp.f90'

program POD

use comp_R_omp
use comp_A_omp
use comp_Phi_omp
use rshp_AMR

use mpi
use omp_lib

implicit none

! =============================================================================
! Allocate variables for POD 
! =============================================================================

! Parameters ------------------------------------------------------------------
real,              parameter :: PI = 4.0*atan(1.0)    ! pi = 3.1415 ...
integer,           parameter :: maxlen=256            ! max length of string
character(len=20), parameter :: CFMT = "(F12.8)"      ! format of CPU time
character(len=20), parameter :: PFMT = "(F12.10)"     ! format of p and p^m
integer,           parameter :: LI  = MPI_OFFSET_KIND ! 64-bit long signed integer (presumably)

! General ---------------------------------------------------------------------
integer               :: i, j, k    ! spatial counting variables
integer               :: m, n       ! time counting variables
integer               :: l          ! level counting variable
integer               :: v          ! variable counting
logical               :: chk        ! general check variable 
integer               :: fid        ! general file ID
integer               :: fcpu       ! file ID for cpu
character(len=maxlen) :: filename   ! file name
character(len=maxlen) :: filefmt    ! format of file name
character(len=maxlen) :: var   ! solution value variable
character(len=5)      :: nt_id      ! string version of nt
integer               :: istrt, iend ! start and end indices 

! Data ------------------------------------------------------------------------
integer               :: nx, ny, nz ! domain dimensions
integer               :: nspat      ! number of spatial points
integer               :: nspp       ! number of spatial points per processor
integer               :: nstot      ! nspp * nvar = total num of spat points
integer               :: nt         ! number of time steps
character(len=maxlen) :: datadir    ! directory where data is stored

! Standard POD ----------------------------------------------------------------
character,        allocatable, dimension(:)   :: vars*(maxlen) ! variables in computation
character(len=maxlen)                         :: PODdir! directory to store POD data
integer                                       :: nvar  ! number of variable
logical                                       :: rm_mean ! remove mean?
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
integer :: Ralg
integer :: Phialg
integer :: Aalg

! AMR POD ---------------------------------------------------------------------
logical                                       :: do_amr  ! do AMR algorithm?
integer                                       :: finest  ! finest level of AMR
integer                                       :: ndim    ! dimension for POD
integer,          allocatable, dimension(:,:) :: Xgrid   ! companion matrix
double precision, allocatable, dimension(:)   :: Xcol    ! column of X
double precision, allocatable, dimension(:)   :: p_comp  ! computed p
double precision, allocatable, dimension(:)   :: pm_comp ! computed p^m

! Timing ----------------------------------------------------------------------
double precision :: Rshp_CPU   ! CPU time for reshaping
double precision :: read_CPU
double precision :: write_CPU
double precision :: R_CPU      ! CPU time for R
double precision :: Psi_CPU      ! CPU time for R
double precision :: Phi_CPU    ! CPU time for Phi
double precision :: A_CPU      ! CPU time for A
! integer                            :: cr         ! clock rate
! double precision                   :: crr        ! real clock rate
! integer                            :: c0, c1, c2 ! timers
! integer :: cpu0, cpuf ! init and final CPU times
double precision :: cpu0, cpuf
integer                            :: hrs, mins  ! hours, minutes
double precision                   :: secs       ! seconds
character(len=10)                  :: clock_time ! time of day

! MPI -------------------------------------------------------------------------
integer :: rank                         ! processor identification
integer :: nproc                        ! number of processors
integer :: mpierr                         ! MPI error flag
integer :: mpistat(MPI_STATUS_SIZE)                        ! MPI status flag
integer :: mpifid
! integer :: 
integer(LI) :: mpidisp, mpioff
integer :: nread
integer :: mpiinfo
double precision :: wtime
! type(MPI_Info) :: mpiinfo

! Namelists - groups of variables to be read from inputs
namelist /POD_inputs/ nx, ny, nz, itime, filefmt, datadir, nvar, Ralg, Phialg, &
                      Aalg
namelist /variables/ vars
namelist /AMR_POD_inputs/ finest

! =============================================================================
! Initialize computation
! =============================================================================

! Intialize MPI
call MPI_Init(mpierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)
call MPI_Comm_size(MPI_COMM_WORLD, nproc, mpierr)

! Set up timing
call date_and_time(time=clock_time)
wtime = MPI_Wtime()

if (rank == 0) then
   write(*,'(a)')'------------------------------------------------'//          &
                 '------------------------------------------------'
   write(*,'(a)')'POD computation started at '//                               &
                 clock_time(1:2)//':'//clock_time(3:4)//':'//clock_time(5:10)//&
                 ' local time.'
   write(*,'(a)') '------------------------------------------------'//         &
                  '------------------------------------------------'
endif

! Set defaults for inputs
nx      = 1
ny      = 1
nz      = 1
itime   = [0, 99999]
filefmt = '(a,a,I5.5,a)'
datadir = './'

! Read inputs
open(newunit=fid, file='POD.inputs', form='formatted', status='old')
read(fid, nml=POD_inputs)
close(fid)

! Now we know how many variables, read all variables to compute POD on
allocate(vars(nvar))
open(newunit=fid, file='POD.inputs', form='formatted', status='old')
read(fid, nml=variables)
close(fid)

! Compute important quantities
ndim = 0
if (nx > 1) ndim = ndim + 1
if (ny > 1) ndim = ndim + 1
if (nz > 1) ndim = ndim + 1
nspat = nx*ny*nz
nspp  = nspat/nproc
nstot = nspp*nvar
nt    = itime(2) - itime(1) + 1
write(nt_id, '(I5)') nt

! Allocate variables for standard POD
allocate(Xpod(nstot,nt))
allocate(Psi(nt,nt))
allocate(Lam(nt))
allocate(Rpod(nt,nt))
allocate(Phi(nstot,nt))
allocate(Apod(nt,nt))
allocate(Xcol(nspp))

! Allocate and assign variables for fast POD

! NEED TO WRITE CODE TO DETERMINE WHICH ALGORITHM TO USE
! But for now...
do_amr = .false.
if (Ralg   == -1 .or. Ralg   == 1)                  do_amr = .true.
if (Phialg == -1 .or. Phialg == 1 .or. Phialg == 2) do_amr = .true.
if (Aalg   == -1 .or. Aalg   == 1)                  do_amr = .true.

! Read in POD AMR data
if (do_amr) then
   finest = 0

   open(newunit=fid, file='POD.inputs', form='formatted', status='old')
   read(fid, nml=AMR_POD_inputs)
   close(fid)

   allocate(Xgrid(nstot,nt))
end if

! Set up output directory information and CPU timing
! allocate(cpu0(nsamp), cpuf(nsamp))
! allocate(R_cpu(nsamp), Phi_cpu(nsamp), A_cpu(nsamp))
PODdir = trim(datadir)//'POD_data_nt'//trim(adjustl(nt_id))//'/'
if (rank == 0) then
   call execute_command_line('mkdir -p '//PODdir, wait=.true.)
   open(newunit=fcpu, file=trim(PODdir)//'cpu_timing.txt', form='formatted')
endif
! =============================================================================
! Error checking
! =============================================================================

! Check the algorithms are valid
if (Ralg < -1 .or. Ralg > 1) then
   if (rank == 0) then
      write(*,'(a)') "Ralg must be -1, 0, or 1!"
   endif
   stop
endif
if (Phialg < -1 .or. Phialg > 2) then
   if (rank == 0) then
      write(*,'(a)') "Phialg must be -1, 0, 1, or 2!"
   endif
   stop
endif
if (Aalg < -1 .or. Aalg > 1) then
   if (rank == 0) then
      write(*,'(a)') "Aalg must be -1, 0, or 1!"
   endif
   stop
endif

! Check that nspat is evenly divisible by the number of MPI processors
if (abs(mod(dble(nz), dble(nproc))) > 1.e-8) then
   if (rank == 0) then
      write(*,'(a)') "the number of z spatial points is not evenly divisible"//&
                     " by the number of MPI processors"
      write(*,'(a,I10,a,I4,a,F10.3)') "nspat = ",nspat," / nproc = ",nproc,    &
                                      "  ---->  ",dble(nspat)/dble(nproc)
   endif
   stop

   ! If we are doing AMR algorithm, nspat must be evenly divisible by nproc*2^finest
   ! if (do_amr) then
   !    if (abs(mod(dble(nz), dble(nproc*(2**finest)))) > 1.e-8) then
   !       if (rank == 0) then
   !          write(*,'(a)') "the number of z spatial points is not evenly divisible "// &
   !                         "by the number of MPI processors"
   !          write(*,'(a,I10,a,I4,a,F10.3)') "nspat = ",nspat," / nproc = ",nproc,    &
   !                                          "  ---->  ",dble(nspat)/dble(nproc)
   !          endif
   !       stop
   !    endif
   ! endif
endif

! Check if all variable data exists in the directory
do v=1,nvar
   var = vars(v)
   do j=1,nt
      write(filename,filefmt) trim(datadir), trim(var), itime(1)+j-1, '.bin'
      inquire(file=filename, exist=chk)

      if (.not.chk) then
         if (rank == 0) then
            write(*,'(a)') "the file "//trim(filename)//" does not exist"
         endif
         stop
      endif
   enddo
enddo

! Check if grid level data exists in the directory
if (do_amr) then
   do j=1,nt
      write(filename,filefmt) trim(datadir), 'grid_level', itime(1)+j-1, '.bin'
      inquire(file=filename, exist=chk)

      if (.not.chk) then
         if (rank == 0) then
            write(*,'(a)') "the file "//trim(filename)//" does not exist"
         endif
         stop
      endif
   enddo
end if


call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! =============================================================================
! Read Data 
! =============================================================================

! Initialize CPU timing
read_CPU = 0.
if (do_amr) then
   Rshp_CPU = 0.
endif

! Load data for each snapshot
do v=1,nvar
   var   = vars(v)           ! variable for this loop
   istrt = (v-1)*nspp + 1    ! starting index location with Xpod
   iend  = v*nspp            ! ending index location with Xpod
   mpioff = int(rank*nspp,LI)
   mpidisp = 0

   ! Load all snapshots of that variable
   do n=1,nt
      ! Read variable data
      write(filename,filefmt) trim(datadir), trim(var), itime(1)+n-1, '.bin'

      cpu0 = MPI_Wtime()
      call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, mpifid, mpierr)
      call MPI_File_set_view(mpifid, mpidisp, MPI_DOUBLE, MPI_DOUBLE, 'native', MPI_INFO_NULL, mpierr)
      call MPI_File_read_at_all(mpifid, mpioff, Xcol, nspp, MPI_DOUBLE, mpistat, mpierr)
      call MPI_File_close(mpifid, mpierr)
      cpuf = MPI_Wtime()

      read_CPU = read_CPU + cpuf - cpu0
      Xpod(istrt:iend,n) = Xcol
   enddo

   ! If we are utilizing the AMR, we need to load grid data and reshape data
   if (do_amr) then
      ! Load grid level data
      do n=1,nt
         write(filename,filefmt) trim(datadir), 'grid_level', itime(1)+n-1, '.bin'

         cpu0 = MPI_Wtime()
         call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, mpifid, mpierr)
         call MPI_File_set_view(mpifid, mpidisp, MPI_DOUBLE, MPI_DOUBLE, 'native', MPI_INFO_NULL, mpierr)
         call MPI_File_read_at_all(mpifid, mpioff, Xcol, nspp, MPI_DOUBLE, mpistat, mpierr)
         call MPI_File_close(mpifid, mpierr)
         cpuf = MPI_Wtime()

         read_CPU = read_CPU + cpuf - cpu0
         Xgrid(istrt:iend,n) = int(Xcol)      
      enddo

      finest = maxval(Xgrid)
      call MPI_Allreduce(MPI_IN_PLACE, finest, 1, MPI_INT, MPI_MAX, &
         MPI_COMM_WORLD, mpierr)

      ! Reshape data
      cpu0 = MPI_Wtime()
      !$omp do private(Xcol)
      do n=1,nt
         ! Reshape original data
         Xcol = Xpod(istrt:iend,n)
         call reshape_AMR(nx, ny, int(nz/nproc), finest, Xcol, 'forward')
         Xpod(istrt:iend,n) = Xcol

         ! Reshape grid data
         Xcol = dble(Xgrid(istrt:iend,n))
         call reshape_AMR(nx, ny, int(nz/nproc), finest, Xcol, 'forward')
         Xgrid(istrt:iend,n) = int(Xcol)
      enddo
      !$omp end do
      cpuf = MPI_Wtime()
      Rshp_CPU = Rshp_CPU + cpuf - cpu0
   endif
enddo

! Write out CPU for reading and reshaping
if (rank == 0) then
   write(fcpu, '(a,'//CFMT//')') "read:           ", read_CPU
   if (do_amr) then
      write(fcpu, '(a,'//CFMT//')') "forward reshape:", Rshp_CPU
   endif
endif

! Compute p^m and p
if (do_amr) then

   allocate(p_comp(0:finest), pm_comp(0:finest))
   p_comp = 0.
   pm_comp = 0.

   ! Compute p
   !$omp parallel shared(p_comp)
   !$omp do reduction(+:p_comp) 
   do n=1,nt
      do l=0,finest
         p_comp(l) = p_comp(l) + dble(count(mask=Xgrid(:,n)==l))/dble(nt*nstot)
      enddo
   enddo
   !$omp end do
   !$omp end parallel
   call MPI_Allreduce(MPI_IN_PLACE, p_comp, finest+1, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)

   ! write(*,*) "p_comp = ",p_comp
   ! write(*,"(a,"//PFMT//")") "sum(p_comp) = ",sum(p_comp)

   ! Compute p^m
   !$omp parallel shared(pm_comp)
   !$omp do reduction(+:pm_comp)
   do i=1,nstot
      pm_comp(maxval(Xgrid(i,:))) = pm_comp(maxval(Xgrid(i,:))) + 1./dble(nstot)
   end do
   !$omp end do
   !$omp end parallel
   call MPI_Allreduce(MPI_IN_PLACE, pm_comp, finest+1, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)

   p_comp = p_comp/dble(nproc)
   pm_comp = pm_comp/dble(nproc)

   ! Write out relevant AMR data
   if (rank == 0) then
      open(newunit=fid, file=trim(PODdir)//'p_data.txt', form='formatted')
      write(fid,'(a,I1)') "finest = ",finest 
      write(fid,'(a)') "" 
      write(fid,'(a)') "mean(p_l):" 
      do l=0,finest
         write(fid,'(a,I1,a,'//PFMT//')') "    p_",l," = ",p_comp(l)
      enddo
      write(fid,'(a)') "" 
      write(fid,'(a)') "p_l^m:" 
      do l=0,finest
         write(fid,'(a,I1,a,'//PFMT//')') "    p_",l,"^m = ",pm_comp(l)
      enddo
      close(fid)
   endif

   deallocate(p_comp, pm_comp)
endif

if (rank == 0) write(*,*)

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! =============================================================================
! Compute POD
! =============================================================================

! Compute R -------------------------------------------------------------------

! R with standard algorithm
if (Ralg == -1 .or. Ralg == 0) then
   if (rank == 0) write(*,*) "computing R using standard operations ..."
   cpu0 = MPI_Wtime()
   call compute_R_omp(Xpod, nstot, nt, Rpod, 0)
   call MPI_Allreduce(MPI_IN_PLACE, Rpod, nt*nt, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)
   cpuf = MPI_Wtime()

   ! Write out CPU time for R
   if (rank == 0) then
      R_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "R standard:     ", R_CPU
      do i=1,4
         write(*,*) Rpod(i,1:4)
      enddo
   endif
endif

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! R with AMR algorithm
if (Ralg == -1 .or. Ralg == 1) then
   if (rank == 0) write(*,*) "computing R utilizing AMR ..."
   cpu0 = MPI_Wtime()
   call compute_R_omp(Xpod, nstot, nt, Rpod, 1, Xgrid, finest, ndim)
   call MPI_Allreduce(MPI_IN_PLACE, Rpod, nt*nt, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)
   cpuf = MPI_Wtime()

   ! Write out CPU time for R
   if (rank == 0) then
      R_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "R AMR:          ", R_CPU
      do i=1,4
         write(*,*) Rpod(i,1:4)
      enddo
   endif
endif

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

if (rank == 0) write(*,*)

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
cpu0 = MPI_Wtime()
call dsyev('V', 'U', nt, Rpod, nt, Lam, opt_work, L_work, eig_info)
cpuf = MPI_Wtime()
deallocate(opt_work)

if (rank == 0) then
   Psi_CPU = cpuf - cpu0
   write(fcpu, '(a,'//CFMT//')') "Psi:            ", Psi_CPU
endif
Psi = Rpod

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! Compute Phi -----------------------------------------------------------------

! Phi with standard algorithm
if (Phialg == -1 .or. Phialg == 0) then
   if (rank == 0) write(*,*) "computing Phi using standard algorithm ..."
   cpu0 = MPI_Wtime()
   call compute_Phi_omp(Xpod, Psi, Lam, nstot, nt, Phi, 0)
   cpuf = MPI_Wtime()

   if (rank == 0) then
      Phi_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "Phi standard:   ", Phi_CPU
      do i=1,4
         write(*,*) Phi(i,1:4)
      enddo
   endif
endif

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! Phi with AMR algorithm, method 1
if (Phialg == -1 .or. Phialg == 1) then
   if (rank == 0) write(*,*) "computing Phi utilizing AMR, method 1 ..."
   cpu0 = MPI_Wtime()
   call compute_Phi_omp(Xpod, Psi, Lam, nstot, nt, Phi, 1, Xgrid, finest, ndim)
   cpuf = MPI_Wtime()

   if (rank == 0) then
      Phi_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "Phi AMR 1:      ", Phi_CPU
      do i=1,4
         write(*,*) Phi(i,1:4)
      enddo
   endif
endif

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

! Phi with AMR algorithm, method 2
if (Phialg == -1 .or. Phialg == 2) then
   if (rank == 0) write(*,*) "computing Phi utilizing AMR, method 2 ..."
   cpu0 = MPI_Wtime()
   call compute_Phi_omp(Xpod, Psi, Lam, nstot, nt, Phi, 2, Xgrid, finest, ndim)
   cpuf = MPI_Wtime()

   if (rank == 0) then
      Phi_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "Phi AMR 2:      ", Phi_CPU
      do i=1,4
         write(*,*) Phi(i,1:4)
      enddo
   endif
endif

call MPI_Barrier(MPI_COMM_WORLD, mpierr)

if (rank == 0) write(*,*)

! Compute A -------------------------------------------------------------------

! A with standard algorithm
if (Aalg == -1 .or. Aalg == 0) then
   if (rank == 0) write(*,*) "computing A using standard algorithm ..."
   cpu0 = MPI_Wtime()
   call compute_A_omp(Xpod, Phi, nstot, nt, Apod, 0)
   call MPI_Allreduce(MPI_IN_PLACE, Apod, nt*nt, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)
   cpuf = MPI_Wtime()

   if (rank == 0) then
      A_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "A standard:     ", A_CPU
      do i=1,4
         write(*,*) Apod(i,1:4)
      enddo
   endif
endif

! A with AMR algorithm
if (Aalg == -1 .or. Aalg == 1) then
   if (rank == 0) write(*,*) "computing A utilizing AMR ..."
   cpu0 = MPI_Wtime()
   call compute_A_omp(Xpod, Phi, nstot, nt, Apod, 1, Xgrid, finest, ndim)
   call MPI_Allreduce(MPI_IN_PLACE, Apod, nt*nt, MPI_DOUBLE, MPI_SUM, &
      MPI_COMM_WORLD, mpierr)
   cpuf = MPI_Wtime()

   if (rank == 0) then
      A_CPU = cpuf - cpu0
      write(fcpu, '(a,'//CFMT//')') "A AMR:          ", A_CPU
      do i=1,4
         write(*,*) Apod(i,1:4)
      enddo
   endif
endif

! =============================================================================
! Write Data 
! =============================================================================

! Initialize CPU timing
write_CPU = 0.
if (do_amr) then
   Rshp_CPU = 0.
endif

! Write data for each POD mode
do v=1,nvar
   var   = vars(v)           ! variable for this loop
   istrt = (v-1)*nspp + 1    ! starting index location with Xpod
   iend  = v*nspp            ! ending index location with Xpod
   mpioff = int(rank*nspp,LI)
   mpidisp = 0

   ! If we were utilizing the AMR, we need to reshape data to original form
   if (do_amr) then
      cpu0 = MPI_Wtime()
      !$omp do private(Xcol)
      do n=1,nt
         ! Reshape original data
         Xcol = Phi(istrt:iend,n)
         call reshape_AMR(nx, ny, int(nz/nproc), finest, Xcol, 'reverse')
         Phi(istrt:iend,n) = Xcol
      enddo
      !$omp end do
      cpuf = MPI_Wtime()
      Rshp_CPU = Rshp_CPU + cpuf - cpu0
   end if


   ! Write POD spatial modes
   do n=1,nt
      write(filename,filefmt) trim(PODdir), 'POD_'//trim(var), itime(1)+n-1,   &
         '.bin'

      cpu0 = MPI_Wtime()
      call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpifid, mpierr)
      call MPI_File_set_view(mpifid, mpidisp, MPI_DOUBLE, MPI_DOUBLE, 'native', MPI_INFO_NULL, mpierr)
      call MPI_File_write_at_all(mpifid, mpioff, Xcol, nspp, MPI_DOUBLE, mpistat, mpierr)
      call MPI_File_close(mpifid, mpierr)
      cpuf = MPI_Wtime()

      write_CPU = write_CPU + cpuf - cpu0
   enddo
enddo

! Write data for temporal coefficients
if (rank == 0) then
   write(filename,filefmt) trim(PODdir), 'temporal_coefficients.bin'
   open(newunit=fid, file=trim(filename), action='write', access='stream',     &
      form='unformatted', status='replace')
   cpu0 = MPI_Wtime()
   write(fid) Apod
   cpuf = MPI_Wtime()
   close(fid)

   write_CPU = write_CPU + cpuf - cpu0
endif

! Write out CPU for writing and reshaping
if (rank == 0) then
   write(fcpu, '(a,'//CFMT//')') "write:          ", write_CPU
   if (do_amr) then
      write(fcpu, '(a,'//CFMT//')') "reverse reshape:", Rshp_CPU
   endif
endif

! =============================================================================
! Finalize computation
! =============================================================================

if (rank == 0) close(fcpu)

call date_and_time(time=clock_time)

wtime = MPI_Wtime() - wtime
hrs   = int(wtime/3600.)
mins  = int(mod(wtime, 3600.)/60.)
secs  = wtime - hrs*3600. - mins*60.

if (rank == 0) then
   write(*,'(a)') '------------------------------------------------'//         &
                  '------------------------------------------------'

   write(*,'(a)')'POD computation finished at '//                              &
                 clock_time(1:2)//':'//clock_time(3:4)//':'//clock_time(5:10)//&
                 ' local time.'//new_line('a')

   write(*,'(A,I3.2,A,I2.2,A,I2.2,F0.3,A)') 'Total program walltime = ',       &
      hrs,':', mins,':',int(secs),secs-int(secs),' (HHH:MM:SS.SSS)'

   write(*,'(a)') '------------------------------------------------'//         &
                  '------------------------------------------------'
endif

call MPI_Finalize(mpierr)

end program POD
! =============================================================================