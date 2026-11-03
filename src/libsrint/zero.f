      SUBROUTINE ZERO(A,LEN)
C
C ZEROS OUT THE FIRST LEN ELEMENTS OF DOUBLE PRECISION VECTOR A.
C
CEND
      implicit none 
      integer, intent(in) :: LEN
      REAL*8, intent(inout) :: A(LEN)
      INTEGER M, I, MP1
C
C     CALL DSCAL(LEN, 0.d0, A, 1)
      M=MOD(LEN,5)
      IF(M.NE.0) THEN
       DO I=1,M
        A(I)=0.D0
       ENDDO
      ENDIF
      IF(LEN.GE.5) THEN
       MP1=M+1
!$omp parallel do 
       DO I=MP1,LEN,5
        A(I)=0.D0
        A(I+1)=0.D0
        A(I+2)=0.D0
        A(I+3)=0.D0
        A(I+4)=0.D0
       ENDDO
!$omp end parallel do 
      ENDIF
       DO I=1,LEN
        A(I)=0.D0 
       ENDDO
      RETURN
      END
