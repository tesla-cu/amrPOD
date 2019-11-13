module comp_A
contains

subroutine compute_A(Xpod, Phi, nspat, nt, Apod, method, Xgrid, finest, ndim)
implicit none

! ---------- General variables --------------------------------------
integer :: i, m, n
double precision :: Asum
! ---------- Standard POD variables ---------------------------------
integer                              , intent(in)  :: nspat
integer                              , intent(in)  :: nt
double precision, dimension(nspat,nt), intent(in)  :: Xpod, Phi
double precision, dimension(nt,nt)   , intent(out) :: Apod
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
integer                                            :: Xgrid_max
integer, allocatable, dimension(:)                 :: d_l
integer, allocatable, dimension(:)                 :: Gmat

! ========================== Standard POD ===========================
if (method == 0) then
   write(*,*) "computing A using standard operations ..."
   do m=1,nt
      do n=1,nt
         Asum= 0.
         do i=1,nspat
            Asum = Asum + Xpod(i,m)*Phi(i,n)
         enddo
         Apod(m,n) = Asum
      enddo
   enddo

! ============================= AMR POD =============================
elseif (method == 1) then
   write(*,*) "computing A utilizing AMR ..."

   ! Check if all optional arguments are inputed
   if (present(Xgrid) .and. present(finest) .and. present(ndim)) then
      allocate(Gmat(nspat))
      allocate(d_l(0:finest))
      do i=0,finest
         d_l(i) = (2**ndim)**(finest-i)
      end do
   else
      write(*,*) "Not all optional arguments are present for AMR!"
      stop
   endif

   ! Precompute maximum grid level for each spatial location
   i = 1
   do while(i <= nspat)
      Xgrid_max = Xgrid(i,1)
      do m=2,nt
         if (Xgrid(i,m) > Xgrid_max) then
            Xgrid_max = Xgrid(i,m)
            if (Xgrid_max == finest) then
               exit
            endif
         endif
      enddo
      Gmat(i) = d_l(Xgrid_max)
      i = i + d_l(Xgrid_max)
   enddo

   ! Compute A with proper weighting
   do m=1,nt
      do n=1,nt
         Asum= 0.
         i = 1
         do while(i <= nspat)
            Asum = Asum + dble(Gmat(i))*Xpod(i,m)*Phi(i,n)
            i = i + Gmat(i)
         end do 
         Apod(m,n) = Asum
      enddo
   enddo

   deallocate(Gmat, d_l)
endif

end subroutine compute_A

end module comp_A