      SUBROUTINE ZEROSG(A,LEN)
C
C ZEROS OUT THE FIRST LEN ELEMENTS OF DOUBLE PRECISION VECTOR A.
C
CEND
      use cublas
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*4, device :: A(LEN)
      integer*4 one, len4
      len4 = len; one = 1 
C
      CALL CUBLASSSCAL(LEN4, 0.0, A, one)
      RETURN
      END
