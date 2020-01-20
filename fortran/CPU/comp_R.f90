! =============================================================================
! Description:
!
!     This module computes:
!
!         A = X^T * X
!     
!     where X is of dimension (nspat x nt). This module supports skipping 
!     operations with supplementary variables (AMR POD).
! =============================================================================

module comp_R
contains

subroutine compute_R(Xpod, nspat, nt, Rpod, method, Xgrid, finest, ndim)
implicit none

! =============================================================================
! Allocate variables
! =============================================================================

! General variables -----------------------------------------------------------
integer, intent(in) :: method     ! standard - 0, AMR - 1
integer             :: i, j, m, n ! indices
double precision    :: Rsum

! Standard POD ----------------------------------------------------------------
integer,                               intent(in)  :: nspat
integer,                               intent(in)  :: nt
double precision, dimension(nspat,nt), intent(in)  :: Xpod
double precision, dimension(nt,nt),    intent(out) :: Rpod

! AMR POD variables -----------------------------------------------------------
integer, optional, dimension(nspat,nt), intent(in) :: Xgrid
integer, optional,                      intent(in) :: finest
integer, optional,                      intent(in) :: ndim
integer                                            :: dval
integer                                            :: d_f1
integer, allocatable, dimension(:)                 :: d_l

! =============================================================================
! Standard POD 
! =============================================================================
if (method == 0) then
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

! =============================================================================
! AMR POD 
! =============================================================================
elseif (method == 1) then

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

   d_f1 = d_l(finest-1)

   do n=1,nt
      do m=1,n
         Rsum= 0.
         i = 1

         ! Compute diagonal elements of R
         if (m==n) then
            do while(i <= nspat)
               if (Xgrid(i,m) == finest) then
                  do j=i,i+d_f1-1
                     Rsum = Rsum + Xpod(j,m)*Xpod(j,m)
                  enddo
                  i = i + d_f1
               else
                  dval = d_l(Xgrid(i,m))
                  Rsum = Rsum + dble(dval)*Xpod(i,m)*Xpod(i,m)
                  i    = i + dval
               endif
            end do

         ! Compute off-diagonal elements of R
         else
            do while(i <= nspat)
               if ((Xgrid(i,m) == finest) .or. (Xgrid(i,n) == finest)) then
                  do j=i,i+d_f1-1
                     Rsum = Rsum + Xpod(j,m)*Xpod(j,n)
                  enddo
                  i = i + d_f1
               else
                  if (Xgrid(i,m) > Xgrid(i,n)) then
                     dval = d_l(Xgrid(i,m))
                  else
                     dval = d_l(Xgrid(i,n))
                  end if
                  Rsum = Rsum + dble(dval)*Xpod(i,m)*Xpod(i,n)
                  i    = i + dval
               endif
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
! =============================================================================