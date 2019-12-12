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
double precision, dimension(nt,nspat)              :: XpodT
double precision, dimension(nt,nt),    intent(in)  :: Psi
double precision, dimension(nt),       intent(in)  :: Lambda
double precision, dimension(nspat,nt), intent(out) :: Phi
integer,                               intent(in)  :: method
! ---------- AMR POD variables --------------------------------------
integer, optional, dimension(nspat,nt), intent(in) :: Xgrid
integer, optional,                      intent(in) :: finest
integer, optional,                      intent(in) :: ndim
integer                                            :: dval
integer                                            :: Gval
integer                                            :: d_0, d_1
integer                                            :: Xgrid_max
integer                                            :: lvl, nlev
integer, allocatable, dimension(:)                 :: d_l
integer, allocatable, dimension(:)                 :: Gmat1
integer, allocatable, dimension(:,:)               :: nl2
double precision, allocatable, dimension(:,:)      :: Hmat2
integer, allocatable, dimension(:,:,:)             :: Gmat2
double precision                                   :: Hsum

! ========================== Standard POD ===========================
if (method == 0) then
   XpodT = reshape(Xpod, [nt,nspat], order=[2,1])
   do i=1,nspat
      do m=1,nt
         Phisum= 0.
         do n=1,nt
            ! Phisum = Phisum + Xpod(i,n)*Psi(n,m)
            Phisum = Phisum + XpodT(n,i)*Psi(n,m)
         enddo
         Phi(i,m) = Phisum/sqrt(Lambda(m))
      enddo
   enddo

! ============================= AMR POD =============================
elseif ((method == 1) .or. (method==2)) then

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

      allocate(Gmat1(nspat))

      do i=1,nspat,d_0
         ! Precompute maximum grid level for each spatial location
         Gmat1 = 0
         jj = 1
         j = i
         do while(jj <= d_0)
            Xgrid_max = Xgrid(j,1)
            do m=2,nt
               if (Xgrid(j,m) > Xgrid_max) then
                  Xgrid_max = Xgrid(j,m)
                  if (Xgrid_max == finest) then
                     exit
                  endif
               endif
            enddo
            Gval = d_l(Xgrid_max)
            Gmat1(jj) = Gval
            jj = jj + Gval
            j = j + Gval
         enddo

         ! Compute elements of Phi within coarse cell
         do m=1,nt
            jj = 1
            j=i
            do while(jj <= d_0)
               Phisum = 0.
               do n=1,nt
                  Phisum = Phisum + Xpod(j,n)*Psi(n,m)
               enddo
               Gval = Gmat1(jj)

               Phisum = Phisum/sqrt(Lambda(m))
               do k=j,j+Gval-1
                  Phi(k,m) = Phisum
               enddo
               jj = jj + Gval
               j = j + Gval
            enddo
         enddo
      enddo
      deallocate(Gmat1)

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
                     Hmat2(k-i,1) = Lsum
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