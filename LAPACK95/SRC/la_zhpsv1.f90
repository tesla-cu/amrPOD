 SUBROUTINE ZHPSV1_F95( AP, B, UPLO, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: HPSV_F77 => LA_HPSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    COMPLEX(WP), INTENT(INOUT) :: AP(:), B(:)
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_HPSV'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, N, NN, SIPIV, ISTAT, ISTAT1
    COMPLEX(WP) :: WW
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LPIV(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, REAL, INT, AIMAG
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; NN = SIZE(AP)
    WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
    IF( PRESENT(UPLO) )THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
    IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
!   .. TEST THE ARGUMENTS
    IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
    ELSE IF( SIZE( B ) /= N ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF( SIPIV /= N )THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
      IF( PRESENT(IPIV) )THEN; LPIV => IPIV
      ELSE; ALLOCATE( LPIV(N), STAT = ISTAT ); END IF
      IF( ISTAT == 0 ) THEN
         CALL HPSV_F77( LUPLO, N, 1, AP, LPIV, B, N, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(IPIV) )DEALLOCATE(LPIV, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE ZHPSV1_F95
