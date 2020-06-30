module rshp_AMR
contains

subroutine reshape_AMR(nx, ny, nz, finest, data, direction)
implicit none

! Arguments -------------------------------------------------------------------
integer,                               intent(in)     :: nx, ny, nz
integer,                               intent(in)     :: finest
double precision, dimension(nx*ny*nz), intent(inout)  :: data
character(len=*),                      intent(in)     :: direction

! General variables -----------------------------------------------------------
character(len=1) :: dir ! shortened direction
integer          :: c
integer          :: nspat, ndim
integer          :: nx1, ny1, nz1, nx2, ny2, nz2
integer          :: i, ii, i1, j1, k1, i2, j2, k2

integer,                       dimension(0:finest) :: c_l
double precision, allocatable, dimension(:)        :: data_1D
double precision, allocatable, dimension(:,:)      :: data_2Da
double precision, allocatable, dimension(:,:)      :: data_2Db
double precision, allocatable, dimension(:,:,:)    :: data_3Da
double precision, allocatable, dimension(:,:,:)    :: data_3Db

! =============================================================================
! Quantities derived from user inputs
! =============================================================================

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

! =============================================================================
! Forward reshaping
! =============================================================================
if (dir == "f") then

   ! Two dimensions -----------------------------------------------------------
   if (ndim == 2) then

      allocate(data_2Da(nx,ny))
      data_2Da = reshape(data, [nx,ny])
      do i=0,finest
         c = c_l(i)
         nx1 = size(data_2Da,1)
         ny1 = size(data_2Da,2)
         nx2 = nspat/c
         ny2 = c
         allocate(data_2Db(nx2,ny2))

         do j2=1,c
            i2 = 1
            do j1=j2,ny1,c
               do ii=1,nx1,c
                  do i1=ii,ii+c-1
                     data_2Db(i2,j2) = data_2Da(i1,j1)
                     i2 = i2 + 1
                  enddo
                  i2 = i2 + ny1 - c
               enddo
               i2 = i2 - nx2 + c
            enddo
         enddo

         deallocate(data_2Da)
         allocate(data_2Da(nx2,ny2))
         data_2Da = data_2Db
         deallocate(data_2Db)
      enddo
      data = data_2Da(:,1)
      deallocate(data_2Da)

   ! Three dimensions ---------------------------------------------------------
   elseif (ndim == 3) then

      allocate(data_3Da(nx,ny,nz))
      data_3Da = reshape(data, [nx,ny,nz])
      do i=0,finest
         c = c_l(i)
         nx1 = size(data_3Da,1)
         ny1 = size(data_3Da,2)
         nz1 = size(data_3Da,3)
         nx2 = nspat/(c*c)
         ny2 = c
         nz2 = c
         allocate(data_3Db(nx2,ny2,nz2))

         do k2=1,c
            do j2=1,c
               i2 = 1
               do k1=k2,nz1,c
                  do j1=j2,ny1,c
                     do ii=1,nx1,c
                        do i1=ii,ii+c-1
                           data_3Db(i2,j2,k2) = data_3Da(i1,j1,k1)
                           i2 = i2 + 1
                        enddo
                        i2 = i2 + ny1*nz1/c - c
                     enddo
                     i2 = i2 - nx2 + c
                  enddo
               enddo
            enddo
         enddo

         deallocate(data_3Da)
         allocate(data_3Da(nx2,ny2,nz2))
         data_3Da = data_3Db
         deallocate(data_3Db)

      enddo
      data = data_3Da(:,1,1)
      deallocate(data_3Da)

   endif

! =============================================================================
! Reverse reshaping
! =============================================================================
elseif (dir== "r") then

   ! Two dimensions -----------------------------------------------------------
   if (ndim == 2) then

      allocate(data_2Da(nspat,1))
      data_2Da(:,1) = data
      do i=finest-1,-1,-1
         nx1 = size(data_2Da,1)
         ny1 = size(data_2Da,2)
         if (i >= 0) then
            c = c_l(i)
            nx2 = nspat/c
            ny2 = c
         else
            nx2 = nx
            ny2 = ny
         endif
         allocate(data_2Db(nx2,ny2))

         do j1=1,ny1
            i1 = 1
            do j2=j1,ny2,ny1
               do ii=1,nx2,ny1
                  do i2=ii,ii+ny1-1
                     data_2Db(i2,j2) = data_2Da(i1,j1)
                     i1 = i1 + 1
                  enddo
                  i1 = i1 + ny2 - ny1
               enddo
               i1 = i1 - nx1 + ny1
            enddo
         enddo

         deallocate(data_2Da)
         allocate(data_2Da(nx2,ny2))
         data_2Da = data_2Db
         deallocate(data_2Db)
      enddo

      data = reshape(data_2Da, [nspat])
      deallocate(data_2Da)

   ! Three dimensions -----------------------------------------------------------
   elseif (ndim == 3) then

      allocate(data_3Da(nspat,1,1))
      data_3Da(:,1,1) = data
      do i=finest-1,-1,-1
         nx1 = size(data_3Da,1)
         ny1 = size(data_3Da,2)
         nz1 = size(data_3Da,3)
         if (i >= 0) then
            c = c_l(i)
            nx2 = nspat/(c*c)
            ny2 = c
            nz2 = c
         else
            nx2 = nx
            ny2 = ny
            nz2 = nz
         endif
         allocate(data_3Db(nx2, ny2, nz2))

         do k1=1,nz1
            do j1=1,ny1
               i1 = 1
               do k2=k1,nz2,nz1
                  do j2=j1,ny2,ny1
                     do ii=1,nx2,ny1
                        do i2=ii,ii+ny1-1
                           data_3Db(i2,j2,k2) = data_3Da(i1,j1,k1)
                           i1 = i1 + 1
                        enddo
                        i1 = i1 + ny2*nz2/ny1 - ny1
                     enddo
                     i1 = i1 - nx1 + ny1
                  enddo
               enddo
            enddo
         enddo

         deallocate(data_3Da)
         allocate(data_3Da(nx2,ny2,nz2))
         data_3Da = data_3Db
         deallocate(data_3Db)
      enddo

      data = reshape(data_3Da, [nspat])
      deallocate(data_3Da)
   endif

else

   write(*,*) "must do forward or reverse reshape!"
   stop

endif

end subroutine reshape_AMR

end module rshp_AMR