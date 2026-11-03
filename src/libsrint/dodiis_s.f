

      SUBROUTINE DODIIS_S(ERR,AUGERR,TMP,NPHYS,NDIM)
C
C GENERATE DIIS EXPANSION COEFFICIENTS, GIVEN THE ERROR MATRIX
C
      IMPLICIT NONE 
      integer, intent(in) :: NPHYS, NDIM
      real*4, intent(inout) :: ERR(NPHYS,NPHYS)
      real*4, intent(inout) :: AUGERR(NPHYS+1,NPHYS+2)
C
      real*4 DET(2)
      real*4 TMP(*)
      integer I 
C
      CALL ZERO_S(AUGERR,(NPHYS+1)*(NPHYS+2))
      CALL BLKCPYR(ERR,NPHYS,NPHYS,AUGERR,NPHYS+1,NPHYS+1,2,2)
CFL
C   Apply a diagonal preconditioner (1.d-4,1,...,1) in order to slightly improve
C   numerical stability:
C
      DO 10 I=2,NDIM+1
       AUGERR(I,1)=-1.e-4
       AUGERR(1,I)=-1.e-4
10    CONTINUE
      AUGERR(1,NDIM+2)=-1.e-4
      CALL ZERO_S(AUGERR(2,NDIM+2),NDIM)
!      write(6,*)'test in dodiis'
      CALL MINVR(AUGERR,NDIM+1,NPHYS+1,TMP,DET,0.00,1,0)
!      write(6,*)'test in dodiis'
      CALL SCOPY(NDIM,AUGERR(2,NDIM+2),1,TMP,1)
      RETURN
      END
