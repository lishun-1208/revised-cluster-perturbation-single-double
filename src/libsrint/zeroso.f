      SUBROUTINE ZEROSO(A,LEN)
C
C ZEROS OUT THE FIRST LEN ELEMENTS OF DOUBLE PRECISION VECTOR A.
C
CEND  
      use mod_iop
      IMPLICIT NONE
      REAL*4,managed :: A(*)
      INTEGER LEN, i
C
      IF(PREC==4.and.igpu==0.and.len/=0) THEN  
!$omp parallel do 
          do i = 1, len
             a(i) = 0.0 
          enddo 
!$omp end parallel do 
C          CALL SCALSO(LEN, 0.d0, A, 1)
      ELSEIF(PREC==8.and.igpu==0.and.len/=0) THEN 
C         CALL DSCAL(LEN, 0.d0, A, 1)
         CALL ZEROSOD(A, LEN)
      ELSEIF(PREC==4.and.igpu==1.and.len/=0) THEN 
          CALL SSCALG(LEN, 0.0, A, 1)
      ELSEIF(PREC==8.and.igpu==1.and.len/=0) THEN 
         CALL DSCALG(LEN, 0.d0, A, 1)
      ENDIF
      RETURN
      END
C
      SUBROUTINE ZEROSOD(A,LEN)
      IMPLICIT NONE
      REAL*8,managed :: A(*)
      INTEGER LEN, i
C
!$omp parallel do 
          do i = 1, len
             a(i) = 0.d0
          enddo
!$omp end parallel do 
      return
      END
