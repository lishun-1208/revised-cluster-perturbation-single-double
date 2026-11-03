      SUBROUTINE ZERODG(A,LEN)
C
C ZEROS OUT THE FIRST LEN ELEMENTS OF DOUBLE PRECISION VECTOR A.
C
CEND
      use cublas
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision, device :: A(LEN)
      integer*4 one, len4
      len4 = len; one = 1 
C
      CALL CUBLASDSCAL(LEN4, 0.d0, A, one)
      RETURN
      END
