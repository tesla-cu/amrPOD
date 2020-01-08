module comp_Phi
contains

subroutine compute_Phi(Xpod, Psi, Lambda, nspat, nt, Phi, method, Xgrid, finest, ndim)
implicit none

! ---------- General variables --------------------------------------
integer :: i, j, jj, k, l, m, n, p
double precision :: Phisum, Lsum
! ---------- Standard POD variables ---------------------------------
integer,                               intent(in)  :: nspat
integer,                               intent(in)  :: nt
double precision, dimension(nspat,nt), intent(in)  :: Xpod
double precision, dimension(nt,nt),    intent(in)  :: Psi
double precision, dimension(nt),       intent(in)  :: Lambda
double precision, dimension(nspat,nt), intent(out) :: Phi
integer,                               intent(in)  :: method
double precision                                   :: temp
! ---------- AMR POD variables --------------------------------------
integer, optional, dimension(nspat,nt), intent(in) :: Xgrid
integer, optional,                      intent(in) :: finest
integer, optional,                      intent(in) :: ndim
integer                                            :: dval
integer                                            :: Gval
integer                                            :: d_0, d_1
integer                                            :: Xgrid_max
integer                                            :: lvl, nlev
integer,          allocatable, dimension(:)        :: Gmat1
integer,          allocatable, dimension(:)        :: d_l
integer,          allocatable, dimension(:,:)      :: nl2
double precision, allocatable, dimension(:,:)      :: Hmat2
integer,          allocatable, dimension(:,:,:)    :: Gmat2
double precision                                   :: Hsum
integer,          allocatable, dimension(:)        :: Xgrid_max2
double precision, allocatable, dimension(:)        :: Phisum2

! ========================== Standard POD ===========================
if (method == 0) then

   ! Old method
   ! do i=1,nspat
   !    do m=1,nt
   !       Phisum= 0.
   !       do n=1,nt
   !          Phisum = Phisum + Xpod(i,n)*Psi(n,m)
   !       enddo
   !       Phi(i,m) = Phisum/sqrt(Lambda(m))
   !    enddo
   ! enddo

   ! New method based on LAPACK
   Phi = 0.
   do m=1,nt
      do n=1,nt
         temp = Psi(n,m)
         do i=1,nspat
            Phi(i,m) = Phi(i,m) + temp*Xpod(i,n)
         enddo
      enddo
      ! Phi(:,m) = Phi(:,m)/sqrt(Lambda(m))
      temp = sqrt(Lambda(m))
      do i=1,nspat
         Phi(i,m) = Phi(i,m)/temp
      enddo
   enddo

! ============================= AMR POD =============================
elseif ((method==1) .or. (method==2)) then

   ! Check if all optional arguments are inputed
   if (present(Xgrid) .and. present(finest) .and. present(ndim)) then
      allocate(d_l(0:finest))
      do i=0,finest
         d_l(i) = (2**ndim)**(finest-i)
      end do
   else
      write(*,*) "Not all optional arguments are present for AMR!"
      stop
   endif

   d_0 = d_l(0)
   d_1 = d_l(1)

   ! ---------- Method 1 --------------------------------------------
   if (method == 1) then

      ! Old method
      ! allocate(Gmat1(d_0))
      ! do i=1,nspat,d_0
      !    ! Precompute maximum grid level for each spatial location
      !    Gmat1 = 0
      !    jj = 1
      !    j = i
      !    do while(jj <= d_0)
      !       Xgrid_max = Xgrid(j,1)
      !       do m=2,nt
      !          if (Xgrid(j,m) > Xgrid_max) then
      !             Xgrid_max = Xgrid(j,m)
      !             if (Xgrid_max == finest) then
      !                exit
      !             endif
      !          endif
      !       enddo
      !       Gval = d_l(Xgrid_max)
      !       Gmat1(jj) = Gval
      !       jj = jj + Gval
      !       j = j + Gval
      !    enddo

      !    ! Compute elements of Phi within coarse cell
      !    do m=1,nt
      !       jj = 1
      !       j=i
      !       do while(jj <= d_0)
      !          Phisum = 0.
      !          do n=1,nt
      !             Phisum = Phisum + Xpod(j,n)*Psi(n,m)
      !          enddo
      !          Gval = Gmat1(jj)

      !          ! Take this out
      !          Phisum = Phisum/sqrt(Lambda(m))
      !          do k=j,j+Gval-1
      !             Phi(k,m) = Phisum
      !          enddo
      !          ! and do after completion
      !          jj = jj + Gval
      !          j = j + Gval
      !       enddo
      !    enddo
      ! enddo
      ! deallocate(Gmat1)

      ! New method based on LAPACK
      allocate(Gmat1(d_0))
      allocate(Xgrid_max2(d_0))
      allocate(Phisum2(d_0))
      Phi = 0.
      do i=1,nspat,d_0
         ! Precompute maximum grid level for each spatial location
         Gmat1 = 0
         jj = 1
         j = i
         Xgrid_max2 = 0
         do m=1,nt
            do while(jj <= d_0)
               if (Xgrid(j,m) > Xgrid_max2(jj)) then
                  Xgrid_max2(jj) = Xgrid(j,m)
               endif
               Gval = d_l(Xgrid_max2(jj))
               Gmat1(jj) = Gval
               jj = jj + Gval
               j = j + Gval
            enddo
         enddo

         ! Compute elements of Phi within coarse cell
         do m=1,nt
            Phisum2 = 0.
            do n=1,nt
               jj = 1
               j = i
               temp = Psi(n,m)
               do while(jj <= d_0)
                  Phisum2(jj) = Phisum2(jj) + temp*Xpod(j,n)
                  Gval = Gmat1(jj)
                  jj = jj + Gval
                  j = j + Gval
               enddo
            enddo
            jj = 1
            j = i
            do while(jj <= d_0)
               temp = Phisum2(jj)/sqrt(Lambda(m))
               Gval = Gmat1(jj)
               do k=j,j+Gval-1
                  Phi(k,m) = temp
               enddo
               jj = jj + Gval
               j = j + Gval
            enddo
         enddo
      enddo
      deallocate(Gmat1, Xgrid_max2, Phisum2)

   ! ---------- Method 2 --------------------------------------------
   elseif (method == 2) then

      nlev = finest + 1

      allocate(Gmat2(0:finest, d_1, nt))
      allocate(nl2(0:finest, d_1))
      allocate(Hmat2(d_0,0:finest))

      do i=1,nspat,d_0
         Gmat2 = 0
         nl2 = 0

         ! Tabulate where grid levels are located 
         do n=1,nt
            lvl = Xgrid(i,n)
            if (lvl == 0) then
               nl2(0,1) = nl2(0,1) + 1
               Gmat2(0,1,nl2(0,1)) = n
            else
               if (finest > 1) then
                  jj = 1
                  j = i
                  do while(jj <= d_1)
                     lvl = Xgrid(j,n)
                     if (lvl == finest) then
                        nl2(finest,jj) = nl2(finest,jj) + 1
                        Gmat2(finest,jj,nl2(finest,jj)) = n
                        jj = jj + 1
                        j = j + d_l(finest-1)
                     else
                        do l=1,finest-1
                           if (lvl==l) then
                              nl2(l,jj) = nl2(l,jj) + 1
                              Gmat2(l,jj,nl2(l,jj)) = n
                              jj = jj + d_l(l+1)
                              j = j + d_l(l)
                           endif
                        enddo
                     endif
                  enddo
               else
                  nl2(1,1) = nl2(1,1) + 1
                  Gmat2(1,1,nl2(1,1)) = n
               endif
            endif
         enddo

         ! Compute elements of H
         do n=1,nt
            Hmat2 = 0.
            if (nl2(0,1) > 0) then
               Lsum = 0.
               do m=1,nl2(0,1)
                  k = Gmat2(0,1,m)
                  Lsum = Lsum + Xpod(i,k)*Psi(k,n)
               enddo
               do m=1,d_0
                  Hmat2(m,0) = Lsum
               enddo
            endif
            if (nl2(0,1) < nt) then
               if (finest > 2) then
                  do l=1,finest-2
                     jj = 1
                     do j=i,i+d_0-1,d_l(l)
                        if (nl2(l,jj) > 0) then
                           Lsum = 0.
                           do m=1,nl2(l,jj)
                              k = Gmat2(l,jj,m)
                              Lsum = Lsum + Xpod(j,k)*Psi(k,n)
                           enddo
                           do m=(jj-1)*d_l(finest-1)+1,(jj-1)*d_l(finest-1)+d_l(l)
                              Hmat2(m,l) = Lsum
                           enddo
                        endif
                        jj = jj + d_l(l+1)
                     enddo
                  enddo
               endif
               if (finest > 1) then
                  jj = 1
                  do j=i,i+d_0-1,d_l(finest-1)
                     if (nl2(finest-1,jj) > 0) then
                        l = finest-1
                        Lsum = 0.
                        do m=1,nl2(l,jj)
                           k = Gmat2(l,jj,m)
                           Lsum = Lsum + Xpod(j,k)*Psi(k,n)
                        enddo
                        do m=(jj-1)*d_l(l)+1,jj*d_l(l)
                           Hmat2(m,l) = Lsum
                        enddo
                     endif
                     if (nl2(finest,jj) > 0) then
                        do k=j,j+d_l(finest-1)-1
                           Lsum = 0.
                           do m=1,nl2(finest,jj)
                              p = Gmat2(finest,jj,m)
                              Lsum = Lsum + Xpod(k,p)*Psi(p,n)
                           enddo
                           Hmat2(k-j+d_l(finest-1)*(jj-1)+1,finest) = Lsum
                        enddo
                     endif
                     jj = jj + 1
                  enddo
               else
                  do k=i,i+d_0-1
                     Lsum = 0.
                     do m=1,nl2(1,1)
                        p = Gmat2(1,1,m)
                        Lsum = Lsum + Xpod(k,p)*Psi(p,n)
                     enddo
                     Hmat2(k-i+1,1) = Lsum
                  enddo
               endif
            endif

            ! Compute elements of Phi
            do m=1,d_0
               Hsum = 0.
               do l=0,finest
                  Hsum = Hsum + Hmat2(m,l)
               enddo
               Phi(m+i-1,n) = Hsum/sqrt(Lambda(n))
            enddo
         enddo
      enddo
      deallocate(Gmat2, nl2, Hmat2)
   endif

   deallocate(d_l)
endif

end subroutine compute_Phi

end module comp_Phi