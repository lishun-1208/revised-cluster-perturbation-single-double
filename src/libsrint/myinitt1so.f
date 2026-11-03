      SUBROUTINE MYINITT1SO(RHS, FAI, IRREPX, NV, NO, VRT, POP)
C
      use mod_iop
      IMPLICIT NONE
C
      INTEGER, intent(in) :: IRREPX, NV, NO, POP(8), VRT(8)
      DOUBLE PRECISION, intent(in) :: FAI(*)
C
      DOUBLE PRECISION, intent(out) :: RHS(*)
C
      IF(PREC==8) THEN 
        CALL MYINITT1_D(RHS, FAI, IRREPX, NV, NO, VRT, POP)
      ELSEIF(PREC==4) THEN 
        CALL MYINITT1_S(RHS, FAI, IRREPX, NV, NO, VRT, POP)
      ENDIF
C
      RETURN
      END
      SUBROUTINE MYINITT1_D(RHS, FAI, IRREPX, NV, NO, VRT, POP)
C
      use mod_symm
      IMPLICIT NONE
C
      INTEGER, intent(in) :: IRREPX, NV, NO, POP(8), VRT(8)
      DOUBLE PRECISION, intent(in) :: FAI(NV, NO)
C
      DOUBLE PRECISION, intent(out) :: RHS(*)
C
      INTEGER :: IRREPL, IRREPR, I, A, IDX, IDXA, IDXI
      INTEGER :: IOFFO(8), IOFFV(8)
C
      IOFFO(1) = 0
      DO IRREPR = 1, NIRREP-1
         IOFFO(IRREPR+1) = IOFFO(IRREPR) + POP(IRREPR)
      ENDDO
      IOFFV(1) = 0
      DO IRREPR = 1, NIRREP-1
         IOFFV(IRREPR+1) = IOFFV(IRREPR) + VRT(IRREPR)
      ENDDO
C
      IDX = 0
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPX,IRREPR)
         DO I = 1, POP(IRREPR)
         IDXI = IOFFO(IRREPR) + I
         DO A = 1, VRT(IRREPL)
            IDXA = IOFFV(IRREPL) + A
            IDX = IDX + 1
            RHS(IDX) = FAI(IDXA, IDXI)
         ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE MYINITT1_S(RHS, FAI, IRREPX, NV, NO, VRT, POP)
C
      use mod_symm
      IMPLICIT NONE
C
      INTEGER, intent(in) :: IRREPX, NV, NO, POP(8), VRT(8)
      REAL*4, intent(in) :: FAI(NV, NO)
C
      REAL*4, intent(out) :: RHS(*)
C
      INTEGER :: IRREPL, IRREPR, I, A, IDX, IDXA, IDXI
      INTEGER :: IOFFO(8), IOFFV(8)
C
      IOFFO(1) = 0
      DO IRREPR = 1, NIRREP-1
         IOFFO(IRREPR+1) = IOFFO(IRREPR) + POP(IRREPR)
      ENDDO
      IOFFV(1) = 0
      DO IRREPR = 1, NIRREP-1
         IOFFV(IRREPR+1) = IOFFV(IRREPR) + VRT(IRREPR)
      ENDDO
C
      IDX = 0
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPX,IRREPR)
         DO I = 1, POP(IRREPR)
         IDXI = IOFFO(IRREPR) + I
         DO A = 1, VRT(IRREPL)
            IDXA = IOFFV(IRREPL) + A
            IDX = IDX + 1
            RHS(IDX) = FAI(IDXA, IDXI)
         ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
