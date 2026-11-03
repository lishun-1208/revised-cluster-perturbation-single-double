      SUBROUTINE RDDIRSO_S(IUNIT,IREC,IVEC,LENGTH)
C
C READS A RECORD FROM DIRECT ACCESS UNIT IUNIT.
C
C  INPUT:
C        IUNIT - UNIT NUMBER FOR DIRECT ACCESS FILE.
C        IREC  - NUMBER OF RECORD TO BE READ.
C        LENGTH- LENGTH OF VECTOR IN *INTEGER* WORDS.
C
C        IVEC  - CONTENTS OF RECORD IREC.
C
CEND
      IMPLICIT NONE 
      integer, intent(in) :: LENGTH, IREC, IUNIT
      real*4, intent(inout) :: IVEC(LENGTH)
C
      integer irealunit, irealrec, IER
C      nsplit=max(1,numsplit(iunit-49))
C      ioff=mod(irec-1,nsplit)
C      irealunit=iunit+ioff*5
C      irealrec =1+(irec-1)/nsplit
      irealunit = iunit 
      irealrec = irec
C      write(6,*)'begin read', irealunit, irealrec
      READ(IrealUNIT,REC=IrealREC,ERR=555,IOSTAT=IER)IVEC
      RETURN
555   WRITE(6,*)'ERROR IN RDDIRSO_S' 
      END
