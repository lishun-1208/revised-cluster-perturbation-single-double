       SUBROUTINE SSTGENTMP(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
C
C THIS ROUTINE RESORTS A W(PQ,RS) SYMMETRY-PACKED QUANTITY TO
C W(PS,RQ) [INTERCHANGE OF Q AND S].  
C
C INPUT :
C         WIN - THE W(PQ,RS) QUANTITY [ALL IRREPS]
C       NSIZE - THE SIZE OF THE W QUANTITY
C        POP1 - POPULATION VECTOR FOR P
C        POP2 - POPULATION VECTOR FOR Q
C        POP3 - POPULATION VECTOR FOR R
C        POP4 - POPULATION VECTOR FOR S
C      IRREPX - THE OVERALL SYMMETRY OF THE LIST (USUALLY 1)
C
C OUTPUT:
C        WOUT - THE RESORTED P(PS,RQ) QUANTITY
C
CEND
      use omp_lib
      use mod_iop
      IMPLICIT NONE 
C
      CHARACTER*4, intent(in) :: TYPE
      real*4, intent(in)  :: WIN(*)
      real*4, intent(out) :: WOUT(*), IPOS(*)
      INTEGER, intent(in) :: POP1(8),POP2(8),POP3(8),POP4(8), IRREPX
      INTEGER, intent(in) :: NSIZE 
      if(prec==8) then
      CALL SSTGEND(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      elseif(prec==4) then
      CALL SSTGENS(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      endif
      return
      end

