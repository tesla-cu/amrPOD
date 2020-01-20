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
!    gfortran -O0 -g -fcheck=all -fbacktrace -o POD.ex POD.f90 -L.. -llapack -lrefblas
!    gfortran -O3 -o POD.ex POD.f90 -L.. -llapack -lrefblas 
!
! For more timing info:
!    gfortran -pg -O3 -o POD.ex POD.f90 -L. -llapack -lrefblas 
!    gprof POD.ex > POD.stats
!
!    
!    
!


include 'rshp_AMR.f90'
include 'comp_R.f90'
include 'comp_Phi.f90'
include 'comp_A.f90'

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
integer         :: i, j, k, l ! spatial counting variables
integer         :: m, n ! time counting variables

integer         :: ierror, untin, iostat
real, parameter :: PI = 4.0*atan(1.0) ! pi=3.1415 ...
integer, parameter :: maxlen=256 ! max length of string
logical               :: chk ! check variable 

character(len=20), parameter :: CFMT = "(F12.8)"
character(len=20), parameter :: PFMT = "(F12.10)"

! ---------- Data parameters
integer :: nx, ny, nz ! number of points in each spatial direction
real :: dt ! time step



character(len=maxlen) :: datafmt, datadir, datatyp
character(len=1)      :: dataorder
character(len=maxlen) :: CPUdir
character(len=maxlen) :: filename
integer               :: file_unit, length
integer :: fid




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
character(len=5) :: nt_id
integer, dimension(2) :: itime ! temporal indices for POD (start, end)

integer :: nspat ! number of spatial points
logical :: rm_mean ! should we remove the mean?

double precision, allocatable, dimension(:,:) :: Xpod, Rpod, Phi, Psi, Apod ! POD matrices
double precision, allocatable, dimension(:)   :: Lam, opt_work
double precision :: Rsum, Phisum, Asum ! single elements of R, Phi, and A
integer :: L_work, eig_info

character(len=maxlen) :: PODdir
integer :: Ralg, Phialg, Aalg ! algorithm to use for each of the 3 steps

! ---------- Fast POD variables ------------------------------------------
integer, allocatable, dimension(:,:) :: Xgrid ! companion matrix with grid levels
double precision, allocatable, dimension(:)   :: Xcol ! column of X
integer, allocatable, dimension(:)   :: Xgcol ! column of Xgrid
logical :: do_amr ! are we doing the fast algorithm leveraging amr anywhere?
integer :: finest, nlev ! finest level of AMR
integer :: ndim ! dimension for POD (1D, 2D, or 3D)
integer :: cval, dval
integer, allocatable, dimension(:) :: c_l, d_l
double precision, allocatable, dimension(:,:) :: data_2D, data_2D_rev
double precision, allocatable, dimension(:,:,:) :: data_3D, data_3D_rev

! Timing variables -------------------------------------------------
integer :: nsamp ! number of samples for ensembles of CPU
integer :: tot_CPU ! total CPU time
real :: Rshp_CPU
real, allocatable, dimension(:) :: R_CPU
real, allocatable, dimension(:) :: Phi_CPU
real, allocatable, dimension(:) :: A_CPU ! cpu time of each operation
double precision, allocatable, dimension(:) :: p_comp, pm_comp ! computed p and p^m
integer      :: c0, c1, c2, cr
integer, allocatable, dimension(:) :: cpu0, cpuf ! init and final CPU times for each op
integer           :: hrs, mins

double precision :: clock_rate, seconds

character(len=10) :: clock_time

! real :: cpu0, cpuf ! initial and final CPU times
! real, allocatable, dimension(:) :: cpu0_arr, cpuf_arr ! initial and final CPU times

namelist /POD_inputs/ nx, ny, nz, nsamp, itime, ix, iy, iz, &
datafmt, datadir, datatyp, dataorder, PODdir, Ralg, Phialg, Aalg, variables, &
var_nam_len, rm_mean

namelist /AMR_POD_inputs/ finest


! Set up timing
call date_and_time(time=clock_time)
call system_clock(count_rate=cr)
clock_rate = real(cr)
call system_clock(c0)

write(*, '(A)') '------------------------------------------------'//   &
                '------------------------------------------------'
write(*, '(A)') 'POD computation started at '//             &
               clock_time(1:2)//':'//clock_time(3:4)//':'//           &
               clock_time(5:10)//' local time.'
write(*, '(A)') '------------------------------------------------'//   &
                '------------------------------------------------'


! ---------- Setup MPI

! integer :: error, rank, nproc
! real :: wtime

! call MPI_Init     (error)                        ! initialize MPI
! call MPI_Comm_size(MPI_COMM_WORLD, nproc, error) ! get the nproc
! call MPI_Comm_rank(MPI_COMM_WORLD, rank,  error) ! get the rank

! wtime = MPI_Wtime() ! wall time




! =========================== Read Inputs ===========================





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





! call system_clock(count_rate=rate)
! rate = real(cr, DP)
! call system_clock(c0)



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
write(nt_id, '(I5)') nt
nspat = nx*ny*nz





allocate(Xpod(nvar*nspat,nt))

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
   ! write(*,*) c_l
   ! write(*,*) d_l
   
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

! Set up CPU timing
allocate(cpu0(nsamp), cpuf(nsamp))
allocate(R_cpu(nsamp), Phi_cpu(nsamp), A_cpu(nsamp))
CPUdir = trim(datadir)//'CPU_timing_nt'//trim(adjustl(nt_id))//'/'
call execute_command_line('mkdir -p '//CPUdir, wait=.true.)

! Output important information about the POD (basically just inputs)

! ============================ Read Data ============================

! Note may want to read, depending on the size that X becomes,
! in the snapshots one at a time


Rshp_CPU = 0.

do i=1,nvar
   var = variables(i)

   do n=1,nt
      ! Read variable data
      ! NEEDS UPDATED FOR nxp < nx, etc.
      write(filename,datafmt) &
         trim(datadir), trim(var), itime(1)+n-1, '.bin'

      open(newunit=fid, file=trim(filename), action='read', &
         access='stream', form='unformatted', status='old')
      read(fid) Xcol
      ! if (dataorder=='c' .and. ndim==3) then
      !    allocate(data_3D_rev(nz,ny,nx), data_3D(nx,ny,nz))
      !    data_3D_rev = reshape(Xcol, [nz,ny,nx])
      !    data_3D = reshape(data_3D_rev, [nz,ny,nx], order=[3,2,1])
      !    Xcol = reshape(data_3D, [nspat])
      !    deallocate(data_3D_rev, data_3D)
      ! endif
      Xpod(:,n) = Xcol
      close(fid)
      ! write(*,*) Xpod(nspat,j)

      ! If we are utilizing the AMR, we need to load grid data
      ! and reshape data
      if (do_amr) then

         Xcol = Xpod(:,n)
         call system_clock(cpu0(1))
         call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'forward')
         call system_clock(cpuf(1))
         Rshp_CPU = Rshp_CPU + real(cpuf(1) - cpu0(1))/clock_rate
         Xpod(:,n) = Xcol

         write(filename,datafmt) &
            trim(datadir), 'grid_level', itime(1)+n-1, '.bin'

         open(newunit=fid, file=trim(filename), action='read', &
            access='stream', form='unformatted', status='old')
         read(fid) Xcol
         ! if (dataorder=='c' .and. ndim==3) then
         !    allocate(data_3D_rev(nz,ny,nx), data_3D(nx,ny,nz))
         !    data_3D_rev = reshape(Xcol, [nz,ny,nx])
         !    data_3D = reshape(data_3D_rev, [nz,ny,nx], order=[3,2,1])
         !    Xcol = reshape(data_3D, [nspat])
         !    deallocate(data_3D_rev, data_3D)
         ! endif
         close(fid)

         ! Xgcol = Xgrid(:,n)
         call system_clock(cpu0(1))
         call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'forward')
         call system_clock(cpuf(1))
         Rshp_CPU = Rshp_CPU + real(cpuf(1) - cpu0(1))/clock_rate
         ! call reshape_AMR(nxp, nyp, nzp*nvar, finest, Xcol, 'reverse')
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

open(newunit=fid, file=trim(CPUdir)//'Reshape_CPU.txt', &
   form='formatted')
write(*,*) "reshaping cpu time ", Rshp_CPU, " seconds"
write(fid, CFMT) Rshp_CPU
close(fid)

! ---------- Compute p^m and p
if (do_amr) then
   allocate(p_comp(0:finest), pm_comp(0:finest))
   p_comp = 0.
   pm_comp = 0.

   ! Compute p
   do n=1,nt
      do l=0,finest
         ! p_comp(l) = p_comp(l) + dble(sum(Xgrid(:,n), mask=Xgrid(:,n)==l)/nspat)
         p_comp(l) = p_comp(l) + dble(count(mask=Xgrid(:,n)==l))/dble(nt*nspat)
      enddo
   enddo
   write(*,*) "p_comp = ",p_comp
   write(*,*) "sum(p_comp) = ",sum(p_comp)

   ! Compute p^m
   do i=1,nspat
      pm_comp(maxval(Xgrid(i,:))) = pm_comp(maxval(Xgrid(i,:))) + 1./dble(nspat)
   end do
   write(*,*) "pm_comp = ",pm_comp
   write(*,*) "sum(pm_comp) = ",sum(pm_comp)

   ! Write out p to txt file
   open(newunit=fid, file=trim(CPUdir)//'p.txt', &
      form='formatted')
   do i=0,finest
      write(fid, PFMT) p_comp(i)
   enddo
   close(fid)

   ! Write out p^m to txt file
   open(newunit=fid, file=trim(CPUdir)//'pm.txt', &
      form='formatted')
   do i=0,finest
      write(fid, PFMT) pm_comp(i)
   enddo
   close(fid)

   deallocate(p_comp, pm_comp)
endif



! write(*,*) Xpod(:,1)
! write(*,*) "shape(Xpod) = ", shape(Xpod)


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


! ---------- Compute R


! call omp_set_num_threads(4)

allocate(Rpod(nt,nt))


! !$OMP PARALLEL DO PRIVATE(j)
!    do i=1,nsamp,omp_get_num_threads()
!       j = i + omp_get_thread_num()
!       call cpu_time(cpu0_arr(j))
!       call compute_R(Xpod, nspat, nt, Rpod, 0)
!       call cpu_time(cpuf_arr(j))
!       Rcpu_arr(j) = cpuf_arr(j) - cpu0_arr(j)
!       write(*,*) "    cpu time ", Rcpu_arr(j), " seconds"
!    enddo
! !$OMP END PARALLEL DO

write(*,*)
write(*,*) "computing R using standard operations ..."
!!$OMP PARALLEL DO
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_R(Xpod, nspat, nt, Rpod, 0)
   call system_clock(cpuf(i))
enddo
!!$OMP END PARALLEL DO


open(newunit=fid, file=trim(CPUdir)//'R_CPU_standard.txt', &
   form='formatted')
R_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", R_CPU(i), " seconds"
   write(fid, CFMT) R_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Rpod(i,1:4)
! enddo

write(*,*) "computing R utilizing AMR ..."
!!$OMP PARALLEL DO
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_R(Xpod, nspat, nt, Rpod, 1, Xgrid, finest, ndim)
   call system_clock(cpuf(i))
enddo
!!$OMP END PARALLEL DO
! R_CPU = cpuf - cpu0
! write(*,*) "    cpu time ", R_CPU, " seconds"

open(newunit=fid, file=trim(CPUdir)//'R_CPU_AMR.txt', &
   form='formatted')
R_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", R_CPU(i), " seconds"
   write(fid, CFMT) R_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Rpod(i,1:4)
! enddo

! ---------- Compute Psi and Lambda ---------------------------------

allocate(Psi(nt,nt))
allocate(Lam(nt))

! Arguments for dsyev (double symmetric eigenvalue):
! http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
!    'V'      - compute evals and evecs
!    'U'      - use upper triangular mat
!    nt       - order of matrix
!    Rpod     - matrix to compute decomposition, outputs evecs
!    nt       - leading order of A
!    Lam      - eigenvalues
!    opt_work - 
!    L_work   - length of opt_work
!    eig_info - information if errors

L_work = 3*nt
allocate(opt_work(L_work))
call dsyev('V','U',nt,Rpod,nt,Lam,opt_work,L_work,eig_info)
Psi = Rpod
deallocate(Rpod)


! ---------- Compute Phi --------------------------------------------
allocate(Phi(nvar*nspat,nt))


write(*,*)
write(*,*) "computing Phi using standard operations ..."
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_Phi(Xpod, Psi, Lam, nspat, nt, Phi, 0)
   call system_clock(cpuf(i))
enddo 

open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_standard.txt', &
   form='formatted')
Phi_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
   write(fid, CFMT) Phi_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo

write(*,*) "computing Phi utilizing AMR, method 1 ..."
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_Phi(Xpod, Psi, Lam, nspat, nt, Phi, 1, Xgrid, finest, ndim)
   call system_clock(cpuf(i))
enddo

open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_AMR1.txt', &
   form='formatted')
Phi_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
   write(fid, CFMT) Phi_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo

write(*,*) "computing Phi utilizing AMR, method 2 ..."
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_Phi(Xpod, Psi, Lam, nspat, nt, Phi, 2, Xgrid, finest, ndim)
   call system_clock(cpuf(i))
enddo

open(newunit=fid, file=trim(CPUdir)//'Phi_CPU_AMR2.txt', &
   form='formatted')
Phi_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", Phi_CPU(i), " seconds"
   write(fid, CFMT) Phi_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Phi(i,1:4)
! enddo


! ---------- Compute A

allocate(Apod(nt,nt))

write(*,*)
write(*,*) "computing A using standard operations ..."
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_A(Xpod, Phi, nspat, nt, Apod, 0)
   call system_clock(cpuf(i))
enddo

open(newunit=fid, file=trim(CPUdir)//'A_CPU_standard.txt', &
   form='formatted')
A_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", A_CPU(i), " seconds"
   write(fid, CFMT) A_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Apod(i,1:4)
! enddo

write(*,*) "computing A utilizing AMR ..."
do i=1,nsamp
   call system_clock(cpu0(i))
   call compute_A(Xpod, Phi, nspat, nt, Apod, 1, Xgrid, finest, ndim)
   call system_clock(cpuf(i))
enddo
open(newunit=fid, file=trim(CPUdir)//'A_CPU_AMR.txt', &
   form='formatted')
A_CPU = real(cpuf - cpu0)/clock_rate
do i=1,nsamp
   write(*,*) "    cpu time ", A_CPU(i), " seconds"
   write(fid, CFMT) A_CPU(i)
enddo
close(fid)

! do i=1,4
!    write(*,*) Apod(i,1:4)
! enddo

call date_and_time(time=clock_time)
call system_clock(c2)
c1 = (c2 - c0)/cr
hrs = int(c1/3600)
mins = int(mod(c1, 3600)/60)
seconds = real(c1 - hrs*3600 - mins*60) &
         + real(mod(c2-c0, cr))/clock_rate

write(*, '(A)') '------------------------------------------------'//   &
                '------------------------------------------------'

write(*, '(A)') 'POD computation finished at '//            &
               clock_time(1:2)//':'//clock_time(3:4)//':'//           &
               clock_time(5:10)//' local time.'//new_line('a')

write(*, '(A,I3.2,A,I2.2,A,I2.2,F0.3,A)')                              &
   'Total program walltime = ', hrs, ':', mins, ':', int(seconds),    &
   seconds-int(seconds), ' (HHH:MM:SS.SSS)'

write(*, '(A)') '------------------------------------------------'//   &
                '------------------------------------------------'


end program POD


