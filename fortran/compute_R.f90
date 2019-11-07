module comp_R
contains

subroutine compute_R(Xpod, nspat, nt, Rpod, method, Xgrid, finest, ndim)
implicit none

! ---------- General variables --------------------------------------
integer :: i, m, n
double precision :: Rsum
! ---------- Standard POD variables ---------------------------------
integer                              , intent(in)  :: nspat
integer                              , intent(in)  :: nt
double precision, dimension(nspat,nt), intent(in)  :: Xpod
double precision, dimension(nt,nt)   , intent(out) :: Rpod
integer                              , intent(in)  :: method
! ---------- AMR POD variables --------------------------------------
! integer, optional, dimension(nspat,nt) :: Xgrid
! integer, optional                      :: finest
! integer, optional, dimension(0:finest) :: d_l
! integer :: dval
integer, optional, dimension(nspat,nt), intent(in) :: Xgrid
integer, optional,                      intent(in) :: finest
integer, optional,                      intent(in) :: ndim
integer                                            :: dval
integer, allocatable, dimension(:)                 :: d_l

! ========================== Standard POD ===========================
if (method == 0) then
   write(*,*) "computing R using standard operations ..."
   do n=1,nt
      do m=1,n
         Rsum= 0.
         do i=1,nspat
            Rsum = Rsum + Xpod(i,m)*Xpod(i,n)
         enddo
         Rpod(m,n) = Rsum
         Rpod(n,m) = Rsum
      enddo
   enddo

! ============================= AMR POD =============================
elseif (method == 1) then
   write(*,*) "computing R utilizing AMR ..."

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

   do n=1,nt
      do m=1,n
         Rsum= 0.
         i = 1
         if (m==n) then
            do while(i <= nspat)
               dval = d_l(Xgrid(i,m))
               Rsum = Rsum + dble(dval)*Xpod(i,m)*Xpod(i,n)
               i    = i + dval
            end do
         else
            do while(i <= nspat)
               if (Xgrid(i,m) > Xgrid(i,n)) then
                  dval = d_l(Xgrid(i,m))
               else
                  dval = d_l(Xgrid(i,n))
               end if
               Rsum = Rsum + dble(dval)*Xpod(i,m)*Xpod(i,n)
               i    = i + dval
            end do
         end if     
         Rpod(m,n) = Rsum
         Rpod(n,m) = Rsum
      enddo
   enddo

   deallocate(d_l)
endif

end subroutine compute_R

end module comp_R