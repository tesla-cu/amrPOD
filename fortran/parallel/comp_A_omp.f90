! =============================================================================
! Description:
!
!     This module computes:
!
!         A = X^T * Phi
!     
!     where X and Phi are of dimension (nspat x nt). This module supports 
!     skipping operations with supplementary variables (AMR POD).
! =============================================================================

module comp_A_omp
contains

subroutine compute_A_omp(Xpod, Phi, nspat, nt, Apod, method, Xgrid, finest, ndim)
implicit none

! =============================================================================
! Allocate variables
! =============================================================================

! General ---------------------------------------------------------------------
integer, intent(in) :: method     ! standard - 0, AMR - 1
integer             :: i, j, m, n ! indices

! Standard POD ----------------------------------------------------------------
integer                              , intent(in)  :: nspat
integer                              , intent(in)  :: nt
double precision, dimension(nspat,nt), intent(in)  :: Xpod, Phi ! X and Phi
double precision, dimension(nt,nt)   , intent(out) :: Apod      ! A 
double precision                                   :: Asum      ! temporary sum

! AMR POD ---------------------------------------------------------------------
integer, optional, dimension(nspat,nt), intent(in) :: Xgrid
integer, optional,                      intent(in) :: finest
integer, optional,                      intent(in) :: ndim
integer                                            :: d_f1
integer                                            :: Xgrid_max
integer, allocatable, dimension(:)                 :: d_l
integer, allocatable, dimension(:)                 :: Gmat

! =============================================================================
! Standard POD 
! =============================================================================
if (method == 0) then
   !$omp do collapse(2)
   do m=1,nt
      do n=1,nt
         Apod(m,n) = dot_product(Xpod(:,m), Phi(:,n))
      enddo
   enddo
   !$omp end do

! =============================================================================
! AMR POD 
! =============================================================================
elseif (method == 1) then

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

   d_f1 = d_l(finest-1)

   ! Precompute maximum grid level for each spatial location

   ! i = 1
   ! do while(i <= nspat)
   !    Xgrid_max = Xgrid(i,1)
   !    do m=2,nt
   !       if (Xgrid(i,m) > Xgrid_max) then
   !          Xgrid_max = Xgrid(i,m)
   !          if (Xgrid_max == finest) then
   !             exit
   !          endif
   !       endif
   !    enddo
   !    Gmat(i) = d_l(Xgrid_max)
   !    i = i + d_l(Xgrid_max)
   ! enddo

   i = 1
   do while(i <= nspat)
      Xgrid_max = Xgrid(i,1)
      !$omp do
      do m=2,nt
         if (Xgrid(i,m) > Xgrid_max) then
            Xgrid_max = Xgrid(i,m)
            ! if (Xgrid_max == finest) then
            !    exit
            ! endif
         endif
      enddo
      !$omp end do
      if (Xgrid_max == finest) then
         Gmat(i) = 1
         i = i + d_f1
      else
         Gmat(i) = d_l(Xgrid_max)
         i = i + d_l(Xgrid_max)
      endif
   enddo

   ! Gmat(1:nspat:d_f1) = maxval(Xgrid(1:nspat:d_f1,:), dim=2)
   ! do i=1,nspat,d_f1
   !    Gmat(i) = d_l(Gmat(i))
   ! enddo
   
   ! Gmat = maxval(Xgrid, dim=2)
   ! do i=1,nspat
   !    Gmat(i) = d_l(Gmat(i))
   ! enddo

   ! Compute A with proper weighting
   !$omp do collapse(2) private(i, Asum)
   do m=1,nt
      do n=1,nt
         Asum = 0.
         i = 1
         do while(i <= nspat)
            if (Gmat(i) == 1) then
               Asum = Asum + dot_product(Xpod(i:i+d_f1-1,m), Phi(i:i+d_f1-1,n))
               i = i + d_f1
            else
               Asum = Asum + dble(Gmat(i))*Xpod(i,m)*Phi(i,n)
               i = i + Gmat(i)
            endif
            
         end do 
         Apod(m,n) = Asum
      enddo
   enddo
   !$omp end do

   deallocate(Gmat, d_l)
endif

end subroutine compute_A_omp
end module comp_A_omp
! =============================================================================