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

integer :: nx1, ny1, nz1, nx2, ny2, nz2
integer :: i1, j1, k1, i2, j2, k2, ii

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
      ! allocate(data_2Da(nx,ny))
      ! data_2Da = reshape(data, [nx,ny])
      ! do i=0,finest
      !    c = c_l(i)
      !    nxr = size(data_2Da,1)
      !    nyr = size(data_2Da,2)
      !    nrshp1 = nspat/(c*nyr)
      !    nrshp2 = nspat/c

      !    ! Permute [2,1]
      !    allocate(data_2Db(nyr,nxr))
      !    data_2Db = reshape(data_2Da, [nyr,nxr], order=[2,1])
      !    deallocate(data_2Da)

      !    ! Reshape into higher dimension
      !    allocate(data_3Da(nyr,c,nrshp1))
      !    data_3Da = reshape(data_2Db, [nyr,c,nrshp1], order=[1,2,3])
      !    deallocate(data_2Db)

      !    ! Permute [1,3,2]
      !    allocate(data_3Db(nyr,nrshp1,c))
      !    data_3Db = reshape(data_3Da, [nyr,nrshp1,c], order=[1,3,2])
      !    deallocate(data_3Da)

      !    ! Reshape back to 2D
      !    allocate(data_2Da(nrshp2,c))
      !    data_2Da = reshape(data_3Db, [nrshp2,c], order=[1,2])
      !    deallocate(data_3Db)

      ! enddo
      ! data = data_2Da(:,1)
      ! deallocate(data_2Da)

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

   ! ---------- Three dimensions
   elseif (ndim == 3) then
      ! allocate(data_3Da(nx,ny,nz))
      ! data_3Da = reshape(data, [nx,ny,nz])
      ! do i=0,finest
      !    c = c_l(i)
      !    nxr = size(data_3Da,1)
      !    nyr = size(data_3Da,2)
      !    nzr = size(data_3Da,3)
      !    nrshp1 = nspat/(nzr*nyr*c)
      !    nrshp2 = nspat/(nzr*c*c)
      !    nrshp3 = nspat/(c*c)

      !    ! Permute [3,2,1]
      !    allocate(data_3Db(nzr,nyr,nxr))
      !    data_3Db = reshape(data_3Da, [nzr,nyr,nxr], order=[3,2,1])
      !    deallocate(data_3Da)

      !    ! Reshape into higher dimension
      !    allocate(data_4Da(nzr,nyr,c,nrshp1))
      !    data_4Da = reshape(data_3Db, [nzr,nyr,c,nrshp1], order=[1,2,3,4])
      !    deallocate(data_3Db)

      !    ! Permute [1,2,4,3]
      !    allocate(data_4Db(nzr,nyr,nrshp1,c))
      !    data_4Db = reshape(data_4Da, [nzr,nyr,nrshp1,c], order=[1,2,4,3])
      !    deallocate(data_4Da)

      !    ! Reshape next dimension
      !    allocate(data_4Da(nzr,c,nrshp2,c))
      !    data_4Da = reshape(data_4Db, [nzr,c,nrshp2,c], order=[1,2,3,4])
      !    deallocate(data_4Db)

      !    ! Permute [1,3,2,4]
      !    allocate(data_4Db(nzr,nrshp2,c,c))
      !    data_4Db = reshape(data_4Da, [nzr,nrshp2,c,c], order=[1,3,2,4])
      !    deallocate(data_4Da)

      !    ! Reshape back to 3D
      !    allocate(data_3Da(nrshp3,c,c))
      !    data_3Da = reshape(data_4Db, [nrshp3,c,c], order=[1,2,3])
      !    deallocate(data_4Db)

      ! enddo
      ! data = data_3Da(:,1,1)
      ! deallocate(data_3Da)

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

      ! This is a slight improvement over the old method but only about 10%
      !   - still needs to be verified as well for all compositions

      ! allocate(data_3Da(nx,ny,nz))
      ! data_3Da = reshape(data, [nx,ny,nz])
      ! do i=0,finest
      !    if (i > 0) then
      !       deallocate(data_3Da)
      !       allocate(data_3Da(nxr,nyr,nyr))
      !       data_3Da = data_3Db
      !       deallocate(data_3Db)
      !    endif
      !    c = c_l(i)
      !    nxr = size(data_3Da,1)
      !    nyr = size(data_3Da,2)
      !    nzr = size(data_3Da,3)
      !    nrshp1 = nspat/(nzr*nyr*c)
      !    nrshp2 = nspat/(nzr*c*c)
      !    nrshp3 = nspat/(c*c)

      !    allocate(data_3Db(nrshp3,c,c))

      !    data_3Db = reshape( & ! Reshape back to 3D                              
      !               reshape( & ! Permute [1,3,2,4]
      !               reshape( & ! Reshape next dimension
      !               reshape( & ! Permute [1,2,4,3]
      !               reshape( & ! Reshape into higher dimension
      !               reshape( & ! Permute [3,2,1]
      !                       data_3Da, [nzr,nyr,nxr],      order=[3,2,1]),      &
      !                                 [nzr,nyr,c,nrshp1], order=[1,2,3,4]),    &
      !                                 [nzr,nyr,nrshp1,c], order=[1,2,4,3]),    &
      !                                 [nzr,c,nrshp2,c],   order=[1,2,3,4]),    &
      !                                 [nzr,nrshp2,c,c],   order=[1,3,2,4]),    &
      !                                 [nrshp3,c,c],       order=[1,2,3])
      ! enddo
      ! data = data_3Db(:,1,1)
      ! deallocate(data_3Da, data_3Db)
   endif

! ======================== Reverse reshaping ========================
elseif (dir== "r") then

   write(*,*) "doing reverse. NEED TO VALIDATE."
   if (ndim == 2) then
   !    allocate(data_2Da(nspat,1))
   !    data_2Da(:,1) = data
   !    do i=finest-1,0,-1
   !       c = c_l(i)
   !       nxr = size(data_2Da,1)
   !       nyr = size(data_2Da,2)
   !       nrshp1 = nspat/(c*nyr)
   !       nrshp2 = nspat/c

   !       ! Reshape into higher dimension
   !       allocate(data_3Da(c,nrshp1,nyr))
   !       data_3Da = reshape(data_2Da, [c,nrshp1,nyr], order=[1,2,3])
   !       deallocate(data_2Da)

   !       ! Permute [1,3,2]
   !       allocate(data_3Db(c,nyr,nrshp1))
   !       data_3Db = reshape(data_3Da, [c,nyr,nrshp1], order=[1,3,2])
   !       deallocate(data_3Da)

   !       ! Reshape back to 2D
   !       allocate(data_2Db(c,nrshp2))
   !       data_2Db = reshape(data_3Db, [c,nrshp2], order=[1,2])
   !       deallocate(data_3Db)

   !       ! Permute [2,1]
   !       allocate(data_2Da(nrshp2,c))
   !       data_2Da = reshape(data_2Db, [nrshp2,c], order=[2,1])
   !       deallocate(data_2Db)
   !    enddo

   !    nxr = size(data_2Da,1)
   !    nyr = size(data_2Da,2)
   !    nrshp1 = nspat/(ny*nyr)
   !    ! nrshp2 = nspat/c

   !    ! Reshape into higher dimension
   !    allocate(data_3Da(ny,nrshp1,nyr))
   !    data_3Da = reshape(data_2Da, [ny,nrshp1,nyr], order=[1,2,3])
   !    deallocate(data_2Da)

   !    ! Permute [1,3,2]
   !    allocate(data_3Db(ny,nyr,nrshp1))
   !    data_3Db = reshape(data_3Da, [ny,nyr,nrshp1], order=[1,3,2])
   !    deallocate(data_3Da)

   !    ! Reshape back to 2D
   !    allocate(data_2Db(ny,nx))
   !    data_2Db = reshape(data_3Db, [ny,nx], order=[1,2])
   !    deallocate(data_3Db)

   !    ! Permute [2,1]
   !    allocate(data_2Da(nx,ny))
   !    data_2Da = reshape(data_2Db, [nx,ny], order=[2,1])

   !    data = reshape(data_2Da, [nspat])
   !    deallocate(data_2Da)

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

   elseif (ndim == 3) then

      allocate(data_3Da(nspat,1,1))
      data_3Da(:,1,1) = data
      do i=finest-1,-1,-1
         nxr = size(data_3Da,1)
         nyr = size(data_3Da,2)
         nzr = size(data_3Da,3)
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