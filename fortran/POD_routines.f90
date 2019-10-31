MODULE POD_routines
IMPLICIT NONE

CONTAINS

    SUBROUTINE generate_white_noise_planes(nx, ny, rx, ry, rz)
    !
    ! Purpose:
    !     To generate two dimensional planes of white noise with 
    !     dimensions nx by ny
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx, ny
    REAL, PARAMETER :: PI = 4.0*ATAN(1.0)
    REAL, ALLOCATABLE, DIMENSION(:,:) :: rx1, rx2, ry1, ry2, rz1, rz2
    REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: rx, ry, rz

    ALLOCATE(rx1(nx,ny), ry1(nx,ny), rz1(nx,ny))
    ALLOCATE(rx2(nx,ny), ry2(nx,ny), rz2(nx,ny))
    ALLOCATE(rx (nx,ny), ry (nx,ny), rz (nx,ny))


    CALL RANDOM_NUMBER(rx1)
    CALL RANDOM_NUMBER(ry1)
    CALL RANDOM_NUMBER(rz1)
    CALL RANDOM_NUMBER(rx2)
    CALL RANDOM_NUMBER(ry2)
    CALL RANDOM_NUMBER(rz2)

    rx = SQRT(-2.0*LOG(rx1)) * COS(2.0*PI*rx2)
    ry = SQRT(-2.0*LOG(ry1)) * COS(2.0*PI*ry2)
    rz = SQRT(-2.0*LOG(rz1)) * COS(2.0*PI*rz2)

    ! x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

    END SUBROUTINE generate_white_noise_planes

END MODULE POD_routines



! SUBROUTINE generate_white_noise_planes(nx, ny, rx, ry, rz)
! !
! ! Purpose:
! !     To generate two dimensional planes of white noise with 
! !     dimensions nx by ny
! !
! IMPLICIT NONE

! ! Inputs
! INTEGER, INTENT(IN) :: nx, ny

! ! Outputs
! REAL, DIMENSION(5,5), INTENT(OUT) :: rx, ry, rz

! ! ! Subroutine variables
! ! INTEGER :: i, j


! CALL RANDOM_NUMBER(rx)
! CALL RANDOM_NUMBER(ry)
! CALL RANDOM_NUMBER(rz)

! END SUBROUTINE generate_white_noise_planes

