module rshp_AMR
contains

subroutine reshape_AMR(nx, ny, nz, finest, data, direction)
implicit none

! ---------- Arguments ----------------------------------------------
integer                              , intent(in)     :: nx, ny, nz
integer                              , intent(in)     :: finest
double precision, dimension(nx*ny*nz), intent(inout)  :: data
character(len=*)                     , intent(in)  :: direction
! ---------- General variables --------------------------------------
integer :: i, c
integer :: nspat
integer :: ndim
character(len=1)                          :: dir ! shortened direction
integer,              dimension(0:finest) :: c_l
double precision, allocatable, dimension(:)        :: data_1D
double precision, allocatable, dimension(:,:)      :: data_2Da
double precision, allocatable, dimension(:,:)      :: data_2Db
double precision, allocatable, dimension(:,:,:)    :: data_3D
double precision, allocatable, dimension(:,:,:,:)  :: data_4D


! =============== Quantities derived from user inputs ===============
nspat = nx*ny*nz
ndim = 0
if (nx > 1) ndim = ndim + 1
if (ny > 1) ndim = ndim + 1
if (nz > 1) ndim = ndim + 1
do i=0,finest
   c_l(i) = 2**(finest-i)
end do

if (trim(direction(1:1)) >= "A" .and. trim(direction(1:1)) <= "Z") then
   dir = achar(iachar(direction(1:1)) + 32)
else
   dir = direction(1:1)
endif

! ======================== Forward reshaping ========================
if (dir == "f") then

   write(*,*) "doing forward"

   ! if (ndim == 1) then
   !    write(*,*) 'data is 1D'
   if (ndim == 2) then
      write(*,*) 'data is 2D'
      ! allocate(data_2Da(nx,ny))
      ! data_2Da = reshape(data, [nx,ny])
      ! do c=c_l

      ! enddo

   elseif (ndim == 3) then
      write(*,*) 'data is 3D'
   endif


elseif (dir== "r") then

   write(*,*) "doing reverse"

else

   write(*,*) "must do forward or backward reshape!"
   stop

endif


end subroutine reshape_AMR

end module rshp_AMR