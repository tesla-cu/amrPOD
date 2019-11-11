!
! This Fortran code is to compute POD on a set of binary files 
! containing the flowfield snapshots. This program is designed 
! to be extremely flexible to allow the user to choose the 
! variables, the region, the algorithm, etc. via the inputs-POD
! file. It *will* also capable of computing POD in a scalable manner and
! utilizing the AMR structure of the grid.
!
! The workflow of the program is:
!
!    - allocate variables
!    - read inputs file 
!    - read data
!    - perform POD:                         Matlab syntax:
!       - compute covariance matrix            - R=X'*X
!       - perform eigenvalue decomposition     -
!       - compute spatial modes
!       - compute temporal coefficients
!    - write POD data
! 
!
! To create an executable, at the terminal, run:
!
!    gfortran -O0 -g -fcheck=all -fbacktrace -o POD.ex POD.f90
!    gfortran -O3 -o POD.ex POD.f90
!    gfortran -O3 -o POD.ex POD.f08 -L/usr/lib -llapack -L/usr/lib -lblas
!
! The meaning of each variable name is:
!
!    algorithm - either 'normal' or 'fast', where normal is the 
!                standard method and 'fast' utilizes the AMR
!    dim - spatial dimension of the fields to compute POD
!    datafmt - format of the files to be read and used in POD
!    datadir - directory of the files to be read and used in POD
!    datatyp - type of files to be read and used in POD
!    nspat - number of spatial points (nx*ny*nz)
!    nt - number of time steps
!    nx, ny, nz - number of grid cells in the x, y, and z directions

!    
!    PODdir - directory to write the POD information
!    
!    
!


include 'reshape_AMR.f90'
include 'compute_R.f90'
include 'compute_eig.f90'
include 'compute_Phi.f90'
include 'compute_A.f90'

program POD

use comp_R
use comp_A
use comp_Phi
use rshp_AMR
! use mpi
! use omp_lib

implicit none
! include 'mpif.h'



!

! =============== Allocate variables universal for POD ==============

! ---------- General variables
integer         :: i, j, k ! spatial counting variables
integer         :: m, n ! time counting variables

integer         :: ierror, untin, iostat
real, parameter :: PI = 4.0*atan(1.0) ! pi=3.1415 ...
integer, parameter :: maxlen=256 ! max length of string
logical               :: chk ! check variable 

! ---------- Data parameters
integer :: nx, ny, nz ! number of points in each spatial direction
real :: dt ! time step



character(len=maxlen) :: datafmt, datadir, datatyp
character(len=maxlen) :: filename
integer               :: file_unit, length

real :: cpu0, cpuf ! initial and final CPU times


! ---------- POD variables ------------------------------------------

! cannot use allocatable in reading from namelist so we are setting
! the maximum number of variables to include be 6 and truncate from
! there
character, dimension(6) :: variables*(maxlen) 
integer,   dimension(6) :: var_nam_len
character(len=maxlen) :: var

integer :: nvar
integer :: nxp, nyp, nzp ! number of points in each spatial direction for POD
integer, dimension(2) :: ix, iy, iz ! spatial indices for POD (start, end)
integer :: nt ! number of time steps
integer, dimension(2) :: itime ! temporal indices for POD (start, end)

integer :: nspat ! number of spatial points
logical :: rm_mean ! should we remove the mean?

double precision, allocatable, dimension(:,:) :: Xpod, Rpod, Phi, Psi, Apod ! POD matrices
double precision, allocatable, dimension(:)   :: Lam, opt_work
double precision :: Rsum, Phisum, Asum ! single elements of R, Phi, and A
integer :: L_work, eig_info
real :: Rcpu1, Phicpu1, Acpu1 ! cpu time of each operation using normal algorithm

character(len=maxlen) :: PODdir
integer :: Ralg, Phialg, Aalg ! algorithm to use for each of the 3 steps

! ---------- Fast POD variables ------------------------------------------
integer, allocatable, dimension(:,:) :: Xgrid ! companion matrix with grid levels
double precision, allocatable, dimension(:)   :: Xcol ! column of X
integer, allocatable, dimension(:)   :: Xgcol ! column of Xgrid
real :: Rcpu2, Phicpu2, Acpu2 ! cpu time of each operation using normal algorithm
logical :: do_amr ! are we doing the fast algorithm leveraging amr anywhere?
integer :: finest, nlev ! finest level of AMR
integer :: ndim ! dimension for POD (1D, 2D, or 3D)
integer :: cval, dval
integer, allocatable, dimension(:) :: c_l, d_l
double precision, allocatable, dimension(:,:) :: test_2D



! ---------- Setup MPI

! integer :: error, rank, nproc
! real :: wtime

! call MPI_Init     (error)                        ! initialize MPI
! call MPI_Comm_size(MPI_COMM_WORLD, nproc, error) ! get the nproc
! call MPI_Comm_rank(MPI_COMM_WORLD, rank,  error) ! get the rank

! wtime = MPI_Wtime() ! wall time




! =========================== Read Inputs ===========================


namelist /POD_inputs/ nx, ny, nz, itime, ix, iy, iz, &
datafmt, datadir, datatyp , PODdir, Ralg, Phialg, Aalg, variables, &
var_nam_len, rm_mean

namelist /AMR_POD_inputs/ finest


! ---------- Set defaults
itime = [0, 99999]
! dim   = 2
nx    = 1
ny    = 1
nz    = 1
ix    = [1, 1]
iy    = [1, 1]
iz    = [1, 1]
rm_mean = .true.
Ralg = 1
Phialg = 1
Aalg = 1
datafmt = '(a,a,I5.5,a)'
datadir = '/folder/'
datatyp = 'float64'

PODdir = '/folder/'

variables = ['x_velocity', 'y_velocity', 'z_velocity', &
             '          ', '          ', '          ']
var_nam_len = (/10, 10, 10, 0, 0, 0/)




finest = 0



! ---------- Read inputs

! untin = 1000
! open(untin, file='inputs_POD2', delim='apostrophe')
! write(untin, nml=inputs_POD) ! appears to wr
! close(untin)

! This will override the defaults if they are in 'POD.inputs'
open(9, file='POD.inputs', form='formatted', status='old')
read(9, nml=POD_inputs)
close(9)

! When close to done, FIX THIS!
! nxp = ix(2) - ix(1) + 1
! nyp = iy(2) - iy(1) + 1
! nzp = iz(2) - iz(1) + 1

nxp = nx
nyp = ny
nzp = nz


ndim = 0
if (nxp > 1) ndim = ndim + 1
if (nyp > 1) ndim = ndim + 1
if (nzp > 1) ndim = ndim + 1




! Remove random variable name if not using 6 variables, set nvar
nvar = 0
do i=1,6
   if (var_nam_len(i) == 0) then
      variables(i) = ' '
   else
      nvar = nvar + 1
   endif 
end do

! ---------- Allocate and assign variables for standard POD
nt    = itime(2) - itime(1) + 1
nspat = nx*ny*nz



allocate(Xpod(nvar*nspat,nt))
allocate(Phi (nvar*nspat,nt))
allocate(Rpod(nt,        nt))
allocate(Psi(nt,nt))
allocate(Apod(nt,nt))
allocate(Lam(nt))
allocate(opt_work(nt))

! ---------- Allocate and assign variables for fast POD
if (Ralg>0 .or. Phialg>0 .or. Aalg>0) then
   do_amr = .true.
else
   do_amr = .false.
end if

if (do_amr) then
   open(9, file='POD.inputs', form='formatted', status='old')
   read(9, nml=AMR_POD_inputs)
   close(9)

   nlev = finest+1

   allocate(Xgrid(nvar*nspat,nt))
   allocate(c_l(0:finest), d_l(0:finest))
   allocate(Xcol(nvar*nspat))
   allocate(Xgcol(nvar*nspat))



   do i=0,finest
      c_l(i) = 2**(finest-i)
      d_l(i) = (2**ndim)**(finest-i)
   end do
   write(*,*) c_l
   write(*,*) d_l
   
end if

! ---------- CHECK IF WE SHOULD DO OPTIMZED






! Check (a few) that the files can be accessed


! write(*,*) "nx = ", nx
! write(*,*) "ny = ", ny
! write(*,*) "alg = ", algorithm
! write(*,*) "vars = ", variables
! write(*,*) "vars(1) = ", variables(1)





! ---------- Error checking

! Check if all variables exist in the directory
do i=1,nvar

   var = variables(i)

   do j=1,nt
      ! write file name with full directory
      write(filename,datafmt) &
         trim(datadir), trim(var), itime(1)+j-1, '.bin'

      inquire(file=filename, exist=chk)
      if (.not.chk) then
         write(*,*) "the file "//trim(filename)//" does not exist"
         stop
      endif

   enddo
enddo
! write(*,*) 'file = ', file

! Check if grid level data exist in the directory
if (do_amr) then
   do j=1,nt
      ! write file name with full directory
      write(filename,datafmt) &
         trim(datadir), 'grid_level', itime(1)+j-1, '.bin'

      inquire(file=filename, exist=chk)
      if (.not.chk) then
         write(*,*) "the file "//trim(filename)//" does not exist"
         stop
      endif
   enddo
end if

! if (rank == 0) then
write(*,*) "we will be performing POD on the following variables:"
do i=1,nvar
   write(*,*) "    ",trim(variables(i))
end do
write(*,*)
! endif



! ============================ Read Data ============================

! Note may want to read, depending on the size that X becomes,
! in the snapshots one at a time

do i=1,nvar
   var = variables(i)

   do n=1,nt
      ! Read variable data
      ! NEEDS UPDATED FOR nxp < nx, etc.
      write(filename,datafmt) &
         trim(datadir), trim(var), itime(1)+n-1, '.bin'

      open(1, file=trim(filename), action='read', &
         access='stream', form='unformatted', status='old')
      read(1) Xcol
      Xpod(:,n) = Xcol
      close(1)
      ! write(*,*) Xpod(nspat,j)

      ! If we are utilizing the AMR, we need to load grid data
      ! and reshape data
      if (do_amr) then

         Xcol = Xpod(:,n)
         call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'forward')
         Xpod(:,n) = Xcol

         write(filename,datafmt) &
            trim(datadir), 'grid_level', itime(1)+n-1, '.bin'

         open(1, file=trim(filename), action='read', &
            access='stream', form='unformatted', status='old')
         read(1) Xcol
         close(1)

         

         ! Xgcol = Xgrid(:,n)
         call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'forward')
         Xgrid(:,n) = int(Xcol)

         ! NEEDS DONE FOR XGRID, using ints
         ! Xcol = Xpod(:,n)
         ! call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'forward')
         ! Xpod = Xcol
         
         ! if (n==1) write(*,*) "NEED: RESHAPE AND XGRID"
         ! Xgrid = 0
      end if

   enddo
enddo


! write(*,*) Xpod(:,1)
write(*,*) "shape(Xpod) = ", shape(Xpod)


! =============================== POD ===============================

! For parallelization, we are most likely going to need to come up
! with a way to parallelize the outer loop when computing R, Phi, 
! and A. 

! For R, we only need to compute:
!
!  |\           i.e. do i=1,nt
!  | \                  do j=1,i
!  |  \                    compute R(i,j)
!  |   \                end do
!  |    \            end do
!  |     \
!  .      .
!  .       .
!  .        .
!  |         \
!  ------------
!
! So we need to figure out how many of these computations there are
! in total, then somehow divvy these up the to various processors.
!
! Once we do that, we send those values to the various procs, sync,
! and impose symmetry

! Serial method
! call cpu_time(stime_R1)
! do j=1,nt
!    do i=1,j
!       Rpod(i,j) = sum(Xpod(1:nspat,i)*Xpod(1:nspat,j))
!       Rpod(j,i) = Rpod(i,j)
!    enddo
! enddo
! call cpu_time(etime_R1)

! write(*,*) "cpu time, method 1: ", etime_R1-stime_R1

! write(*,*) "computing R ..."

! ---------- Compute R
write(*,*)
call cpu_time(cpu0)
call compute_R(Xpod, nspat, nt, Rpod, 0)
call cpu_time(cpuf)
Rcpu1 = cpuf - cpu0
write(*,*) "    cpu time ", Rcpu1, " seconds"

do i=1,4
   write(*,*) Rpod(i,1:4)
enddo

call cpu_time(cpu0)
call compute_R(Xpod, nspat, nt, Rpod, 1, Xgrid, finest, ndim)
call cpu_time(cpuf)
Rcpu1 = cpuf - cpu0
write(*,*) "    cpu time ", Rcpu1, " seconds"

do i=1,4
   write(*,*) Rpod(i,1:4)
enddo

! ---------- Compute Psi and Lambda


! ---------- Compute Phi


Phi = 0.1

! ---------- Compute A
write(*,*)
call cpu_time(cpu0)
call compute_A(Xpod, Phi, nspat, nt, Apod, 0)
call cpu_time(cpuf)
Acpu1 = cpuf - cpu0
write(*,*) "    cpu time ", Acpu1, " seconds"

call cpu_time(cpu0)
call compute_A(Xpod, Phi, nspat, nt, Apod, 1, Xgrid, finest, ndim)
call cpu_time(cpuf)
Acpu1 = cpuf - cpu0
write(*,*) "    cpu time ", Acpu1, " seconds"


write(*,*) "Rpod(1,1) = ", Rpod(1,1)

! Argument meanings:
!    'V' - compute evals and evecs
!    'U' - use upper triangular mat, irrelevant if computing evecs
!    nt  - order of matrix
!    Psi - evecs, if using argument 'V' above
!    nt  - leading order of A
!    Lam - eigenvalues

L_work = -1
Lam = 0.0
! write(*,*) shape(Lam)
! write(*,*) Lam
do i=1,4
   write(*,*) Rpod(i,1:4)
enddo
! call dsyev('V','L',nt,Rpod,nt,Lam,opt_work,L_work,eig_info)
! ! call ssyev('V','L',nt,Rpod,nt,Lam,opt_work,L_work,eig_info)
! write(*,*) "eig_info", eig_info
! write(*,*) "opt_work", opt_work(1)
! write(*,*) "max lam: ", Lam(nt)
! write(*,*) Lam





! First I am disappointed that method 2 is substantially slower than 
! method 1 (about 2-3 times). Surely I thought that using built in 
! functions would not speed up that much, but they did. This then begs
! the question, why even go through the hassle? Is there really that
! much to gain by doing iterations for each spatial point? By 
! implementing a fast method, how much faster do we need to be in 
! theory (i.e., # of iters) to actually see a speed up. The answer is
! probably around an order of magnitude, quite substantial.




!
! Phi should be a bit easier (wrt parallel) since it is a rectangular operation, as
! with A.




! call MPI_Finalize(error) ! shut down MPI.



! ! Set parameters (see Matlab code for more details)
! xlo =-0.1
! xhi = 0.1
! ylo =-0.1
! yhi = 0.1
! nx  = 64
! ny  = 64
! dt = 0.0001
! l  = 0.008
! U_j = 15.3
! Rii = 7.0    ! can get away using (U_j*pertmag)**2 in actual simulations

! ! Calculate basics
! xran = xhi-xlo
! yran = yhi-ylo
! dx = xran/nx
! dy = yran/ny
! nf = INT((l/dx)+0.5)
! nf2 = 2*nf
! tau = l/U_j

! ! Fill R_ij with HIT
! DO i=1,3
!     DO j=1,3
!         IF (i==j) THEN
!             R(i,j) = Rii
!         ELSE
!             R(i,j) = 0.0
!         END IF
!     END DO
! END DO

! ! Cholesky factorization to form amplitude tensor
! a(1,1) = SQRT(R(1,1))
! a(1,2) = 0.0
! a(1,3) = 0.0
! a(2,1) = R(2,1)/a(1,1)
! a(2,2) = SQRT(R(2,2) - a(2,1)**2)
! a(2,3) = 0.0
! a(3,1) = R(3,1)/a(1,1)
! a(3,2) = (R(3,2) - a(2,1)*a(3,1))/a(2,2)
! a(3,3) = SQRT(R(3,3) - a(3,1)**2 - a(3,2)**2)

! ! Calculate k_filter matrix
! ALLOCATE(k_filter(4*nf+1,4*nf+1))
! DO i=1,4*nf+1
!     DO j=1,4*nf+1
!         k_filter(i,j) = SQRT(REAL((i-1-nf2)**2 + (j-1-nf2)**2))
!     END DO
! END DO

! ! Calculate b_filter matrix
! ALLOCATE(b_filter(4*nf+1,4*nf+1))
! bj_sq = SQRT(SUM(SUM(EXP(-PI*k_filter/nf)**2, 1), 1))
! b_filter = EXP(-PI*k_filter/nf)/bj_sq

! ! Construct 2D sheets of white noise
! ALLOCATE(rx(4*nf+nx,4*nf+ny), ry(4*nf+nx,4*nf+ny), rz(4*nf+nx,4*nf+ny))
! CALL generate_white_noise_planes(4*nf+nx, 4*nf+ny, rx, ry, rz) 

! ! Calculate array from -2*nf to 2*nf
! ALLOCATE(NtoN((4*nf)+1))
! DO i=1,(4*nf)+1
!     NtoN(i) = i-nf2-1
! END DO

! ! Calculate psi fields using filter function
! ALLOCATE(psix(nx,ny), psiy(nx,ny), psiz(nx,ny))
! DO i=1,nx
!     DO j=1,ny
!         ! summing all dims of b.*r 
!         psix(i,j) = SUM(SUM(b_filter*rx(i+nf2+NtoN,j+nf2+NtoN),1),1)
!         psiy(i,j) = SUM(SUM(b_filter*ry(i+nf2+NtoN,j+nf2+NtoN),1),1)
!         psiz(i,j) = SUM(SUM(b_filter*rz(i+nf2+NtoN,j+nf2+NtoN),1),1)
!     END DO
! END DO

! ! Calculate u fields
! ALLOCATE(ux(nx,ny), uy(nx,ny), uz(nx,ny))
! ux = a(1,1)*psix
! uy = a(2,1)*psix + a(2,2)*psiy
! uz = a(3,1)*psix + a(3,2)*psiy + a(3,3)*psiz

! ! Write out u, v, w file to analyze in Matlab
! OPEN(UNIT=8, FILE='uvw.dat', STATUS='replace', ACTION='write', &
!     IOSTAT=ierror)

! do i=1,nx
!     do j=1,ny
!         WRITE(8,100) ux(i,j), uy(i,j), uz(i,j)
!         100 FORMAT (F10.6, F10.6, F10.6)
!     end do
! end do

! CLOSE(8)
! WRITE(*,*) "done computing u, v, w for inflow"


end program POD


