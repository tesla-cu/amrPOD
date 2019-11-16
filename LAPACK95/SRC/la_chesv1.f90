 SUBROUTINE CHESV1_F95( A, B, UPLO, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => SP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: HESV_F77 => LA_HESV, ILAENV_F77 => ILAENV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:)
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_HESV'
    CHARACTER(LEN=6), PARAMETER :: BSNAME = 'CHETRF'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV, N, LWORK, NB
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LPIV(:)
    COMPLEX(WP), POINTER :: WORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; N = SIZE(A,1)
    IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
    IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = SIZE(A,1); END IF
!   .. TEST THE ARGUMENTS
    IF( SIZE( A, 2 ) /= N .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B ) /= N ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF( SIPIV /= N )THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
!  .. DETERMINE THE WORKSPACE
      IF( PRESENT(IPIV) )THEN; LPIV => IPIV
      ELSE; ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT ); END IF
      IF( ISTAT == 0 )THEN
         NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
         IF( NB <= 1 .OR. NB >= N ) NB = 1; LWORK = N*NB
         ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN
            DEALLOCATE(WORK, STAT=ISTAT1); LWORK = 3*N-1
            ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF( ISTAT /= 0 ) THEN; LINFO = - 100
            ELSE; CALL ERINFO( -200, SRNAME, LINFO ); ENDIF
         ENDIF
         IF ( ISTAT == 0 ) &
!           .. CALL LAPACK77 ROUTINE
            CALL HESV_F77( LUPLO, N, 1, A, N, LPIV, B, N, WORK, LWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(IPIV) )DEALLOCATE(LPIV, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE CHESV1_F95
