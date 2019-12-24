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
integer :: nxr, nyr, nzr
integer :: nrshp1, nrshp2, nrshp3 ! integers used for reshape
character(len=1)                          :: dir ! shortened direction
integer,              dimension(0:finest) :: c_l
double precision, allocatable, dimension(:)        :: data_1D
double precision, allocatable, dimension(:,:)      :: data_2Da
double precision, allocatable, dimension(:,:)      :: data_2Db
double precision, allocatable, dimension(:,:,:)    :: data_3Da
double precision, allocatable, dimension(:,:,:)    :: data_3Db
double precision, allocatable, dimension(:,:,:,:)  :: data_4Da
double precision, allocatable, dimension(:,:,:,:)  :: data_4Db


! =============== Quantities derived from user inputs ===============
nspat = nx*ny*nz
ndim = 0
if (nx > 1) ndim = ndim + 1
if (ny > 1) ndim = ndim + 1
if (nz > 1) ndim = ndim + 1
do i=0,finest
   c_l(i) = 2**(finest-i)
end do

! Determine direction based on lower case first letter
if (trim(direction(1:1)) >= "A" .and. trim(direction(1:1)) <= "Z") then
   dir = achar(iachar(direction(1:1)) + 32)
else
   dir = direction(1:1)
endif

! ======================== Forward reshaping ========================
if (dir == "f") then

   ! ---------- Two dimensions
   if (ndim == 2) then
      allocate(data_2Da(nx,ny))
      data_2Da = reshape(data, [nx,ny])
      do i=0,finest
         c = c_l(i)
         nxr = size(data_2Da,1)
         nyr = size(data_2Da,2)
         nrshp1 = nspat/(c*nyr)
         nrshp2 = nspat/c

         ! Permute [2,1]
         allocate(data_2Db(nyr,nxr))
         data_2Db = reshape(data_2Da, [nyr,nxr], order=[2,1])
         deallocate(data_2Da)

         ! Reshape into higher dimension
         allocate(data_3Da(nyr,c,nrshp1))
         data_3Da = reshape(data_2Db, [nyr,c,nrshp1], order=[1,2,3])
         deallocate(data_2Db)

         ! Permute [1,3,2]
         allocate(data_3Db(nyr,nrshp1,c))
         data_3Db = reshape(data_3Da, [nyr,nrshp1,c], order=[1,3,2])
         deallocate(data_3Da)

         ! Reshape back to 2D
         allocate(data_2Da(nrshp2,c))
         data_2Da = reshape(data_3Db, [nrshp2,c], order=[1,2])
         deallocate(data_3Db)

      enddo
      data = data_2Da(:,1)
      deallocate(data_2Da)

   ! ---------- Three dimensions
   elseif (ndim == 3) then
      allocate(data_3Da(nx,ny,nz))
      data_3Da = reshape(data, [nx,ny,nz])
      do i=0,finest
         c = c_l(i)
         nxr = size(data_3Da,1)
         nyr = size(data_3Da,2)
         nzr = size(data_3Da,3)
         nrshp1 = nspat/(nzr*nyr*c)
         nrshp2 = nspat/(nzr*c*c)
         nrshp3 = nspat/(c*c)

         ! Permute [3,2,1]
         allocate(data_3Db(nzr,nyr,nxr))
         data_3Db = reshape(data_3Da, [nzr,nyr,nxr], order=[3,2,1])
         deallocate(data_3Da)

         ! Reshape into higher dimension
         allocate(data_4Da(nzr,nyr,c,nrshp1))
         data_4Da = reshape(data_3Db, [nzr,nyr,c,nrshp1], order=[1,2,3,4])
         deallocate(data_3Db)

         ! Permute [1,2,4,3]
         allocate(data_4Db(nzr,nyr,nrshp1,c))
         data_4Db = reshape(data_4Da, [nzr,nyr,nrshp1,c], order=[1,2,4,3])
         deallocate(data_4Da)

         ! Reshape next dimension
         allocate(data_4Da(nzr,c,nrshp2,c))
         data_4Da = reshape(data_4Db, [nzr,c,nrshp2,c], order=[1,2,3,4])
         deallocate(data_4Db)

         ! Permute [1,3,2,4]
         allocate(data_4Db(nzr,nrshp2,c,c))
         data_4Db = reshape(data_4Da, [nzr,nrshp2,c,c], order=[1,3,2,4])
         deallocate(data_4Da)

         ! Reshape back to 3D
         allocate(data_3Da(nrshp3,c,c))
         data_3Da = reshape(data_4Db, [nrshp3,c,c], order=[1,2,3])
         deallocate(data_4Db)

      enddo
      data = data_3Da(:,1,1)
      deallocate(data_3Da)
   endif

! ======================== Reverse reshaping ========================
elseif (dir== "r") then

   write(*,*) "doing reverse"
   if (ndim == 2) then
      allocate(data_2Da(nspat,1))
      data_2Da(:,1) = data
      do i=finest-1,0,-1
         c = c_l(i)
         nxr = size(data_2Da,1)
         nyr = size(data_2Da,2)
         nrshp1 = nspat/(c*nyr)
         nrshp2 = nspat/c

         ! Reshape into higher dimension
         allocate(data_3Da(c,nrshp1,nyr))
         data_3Da = reshape(data_2Da, [c,nrshp1,nyr], order=[1,2,3])
         deallocate(data_2Da)

         ! Permute [1,3,2]
         allocate(data_3Db(c,nyr,nrshp1))
         data_3Db = reshape(data_3Da, [c,nyr,nrshp1], order=[1,3,2])
         deallocate(data_3Da)

         ! Reshape back to 2D
         allocate(data_2Db(c,nrshp2))
         data_2Db = reshape(data_3Db, [c,nrshp2], order=[1,2])
         deallocate(data_3Db)

         ! Permute [2,1]
         allocate(data_2Da(nrshp2,c))
         data_2Da = reshape(data_2Db, [nrshp2,c], order=[2,1])
         deallocate(data_2Db)
      enddo

      nxr = size(data_2Da,1)
      nyr = size(data_2Da,2)
      nrshp1 = nspat/(ny*nyr)
      ! nrshp2 = nspat/c

      ! Reshape into higher dimension
      allocate(data_3Da(ny,nrshp1,nyr))
      data_3Da = reshape(data_2Da, [ny,nrshp1,nyr], order=[1,2,3])
      deallocate(data_2Da)

      ! Permute [1,3,2]
      allocate(data_3Db(ny,nyr,nrshp1))
      data_3Db = reshape(data_3Da, [ny,nyr,nrshp1], order=[1,3,2])
      deallocate(data_3Da)

      ! Reshape back to 2D
      allocate(data_2Db(ny,nx))
      data_2Db = reshape(data_3Db, [ny,nx], order=[1,2])
      deallocate(data_3Db)

      ! Permute [2,1]
      allocate(data_2Da(nx,ny))
      data_2Da = reshape(data_2Db, [nx,ny], order=[2,1])

      data = reshape(data_2Da, [nspat])
      deallocate(data_2Da)
   elseif (ndim == 3) then
      write(*,*) "NOT DONE!!!!!!!!"
      stop

      allocate(data_3Da(nspat,1,1))
      data_3Da(:,1,1) = data
      do i=finest-1,0,-1
         c = c_l(i)
         nxr = size(data_3Da,1)
         nyr = size(data_3Da,2)
         nzr = size(data_3Da,3)
         nrshp1 = nspat/(c*nyr*nzr)
         nrshp2 = nspat/(c*c*nzr)
         nrshp3 = nspat/(c*c)

         ! Reshape into higher dimension
         allocate(data_4Da(c,nrshp1,nyr,nzr))
         data_4Da = reshape(data_3Da, [c,nrshp1,nyr,nzr], order=[1,2,3,4])
         deallocate(data_3Da)

         ! Permute [1,3,2,4]
         allocate(data_4Db(c,nyr,nrshp1,nzr))
         data_4Db = reshape(data_4Da, [c,nyr,nrshp1,nzr], order=[1,3,2,4])
         deallocate(data_4Da)

         ! Reshape other dimension
         allocate(data_4Da(c,c,nrshp2,nzr))
         data_4Da = reshape(data_4Db, [c,c,nrshp2,nzr], order=[1,2,3,4])
         deallocate(data_4Db)

         ! Permute [1,2,4,3]
         allocate(data_4Db(c,c,nzr,nrshp2))
         data_4Db = reshape(data_4Da, [c,c,nzr,nrshp2], order=[1,2,4,3])
         deallocate(data_4Da)

         ! Reshape to three dimensions
         allocate(data_3Db(c,c,nrshp3))
         data_3Db = reshape(data_4Db, [c,c,nrshp3], order=[1,2,3])
         deallocate(data_4Db)

         ! Permute [3,2,1]
         allocate(data_3Da(nrshp3,c,c))
         data_3Da = reshape(data_3Db, [nrshp3,c,c], order=[3,2,1])
         deallocate(data_3Db)
      enddo

      nxr = size(data_3Da,1)
      nyr = size(data_3Da,2)
      nzr = size(data_3Da,3)
      nrshp1 = nspat/(nz*nyr*nzr)
      nrshp2 = nspat/(nz*ny*nzr)

      ! Reshape into higher dimension
      allocate(data_4Da(nz,nrshp1,nyr,nzr))
      data_4Da = reshape(data_3Da, [nz,nrshp1,nyr,nzr], order=[1,2,3,4])
      deallocate(data_3Da)

      ! Permute [1,3,2,4]
      allocate(data_4Db(nz,nyr,nrshp1,nzr))
      data_4Db = reshape(data_4Da, [nz,nyr,nrshp1,nzr], order=[1,3,2,4])
      deallocate(data_4Da)

      ! Reshape other dimension
      allocate(data_4Da(nz,ny,nrshp2,nzr))
      data_4Da = reshape(data_4Db, [nz,ny,nrshp2,nzr], order=[1,2,3,4])
      deallocate(data_4Db)

      ! Permute [1,2,4,3]
      allocate(data_4Db(nz,ny,nzr,nrshp2))
      data_4Db = reshape(data_4Da, [nz,ny,nzr,nrshp2], order=[1,2,4,3])
      deallocate(data_4Da)

      ! Reshape to three dimensions
      allocate(data_3Db(nz,ny,nx))
      data_3Db = reshape(data_4Db, [nz,ny,nx], order=[1,2,3])
      deallocate(data_4Db)

      ! Permute [3,2,1]
      allocate(data_3Da(nx,ny,nz))
      data_3Da = reshape(data_3Db, [nx,ny,nz], order=[3,2,1])
      deallocate(data_3Db)

      data = reshape(data_3Da, [nspat])
      deallocate(data_3Da)
   endif

else

   write(*,*) "must do forward or backward reshape!"
   stop

endif


end subroutine reshape_AMR

end module rshp_AMR