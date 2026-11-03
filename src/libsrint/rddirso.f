      SUBROUTINE RDDIRSO(IUNIT,IREC,IVEC,LENGTH)
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
      double precision, intent(inout) :: IVEC(LENGTH)
C
      integer irealunit, irealrec, IER
      irealunit = iunit 
      irealrec = irec
      READ(IrealUNIT,REC=IrealREC,ERR=555,IOSTAT=IER)IVEC
      RETURN
555   WRITE(6,*)'ERROR IN RDDRISO'
      stop 
      END
