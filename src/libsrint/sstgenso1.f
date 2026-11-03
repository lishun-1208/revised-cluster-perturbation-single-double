       SUBROUTINE SSTGENSO1(WIN,WOUT,POP1,POP2,POP3,
     &                   POP4,IRREPX,TYPE)
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
      use mod_iop
      use omp_lib
      IMPLICIT NONE 
C
      CHARACTER*4, intent(in) :: TYPE
      real*4, intent(in),managed  :: WIN(*)
      real*4, intent(out),managed :: WOUT(*)
      INTEGER, intent(in) :: POP1(8),POP2(8),POP3(8),POP4(8), IRREPX
      integer nsize, ipos
      if(prec==8.and.igpu==0) then
      CALL SSTGEND(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      elseif(prec==4.and.igpu==0) then
      CALL SSTGENS(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      elseif(prec==8.and.igpu==1) then
      CALL SSTGENDG(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      elseif(prec==4.and.igpu==1) then
      CALL SSTGENSG(WIN,WOUT,NSIZE,POP1,POP2,POP3,
     &                   POP4,IPOS, IRREPX,TYPE)
      endif
      return
      end

