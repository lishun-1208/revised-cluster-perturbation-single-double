      SUBROUTINE CLMOIODIIS(ICORE,NUML,NUMR, nlist, filename, dsize)
C-index of the file, 412.
C-NUML distribution size of icore
C-NUMR number of distribution
C-IOPSIZE the number of iop
C dsize was suggested to be set as 4096
C     
      IMPLICIT NONE
      integer*8, intent(in) :: NUML
      integer, intent(in) :: NUMR, nlist, dsize
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(1000),MOIONUMR(1000),RECTMP(1000)
C
C----------------------------------
C
      integer iop, i, list 
      character*2 tmpc
      character*10 filename2
C      write(6,*)'num', numl, numr, nlist, dsize
      OPEN(UNIT=412, FILE=FILENAME, FORM='UNFORMATTED', access='direct',
     &            recl=dsize*8)
    
      CALL IZERO(MOIONUML, nlist)
      CALL IZERO(MOIONUMR, nlist)
      CALL IZERO(rectmp, nlist)
C      write(6,*)'clmoioso by wz'
C***********************************************************************
      do list = 1,nlist 
         CALL UPDMOIDIIS(ICORE,NUML, NUMR,list, nlist, dsize)
      enddo
C
      CLOSE(412)
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=filename2,
     &                 FORM='UNFORMATTED',status='unknown')
      rewind 413 
      write(413) nlist 
      write(413) dsize 
      write(413) numl
      write(413) numr
      write(413) (rectmp(i),i=1,nlist)
      close(413)
C***********************************************************************
      RETURN
      END
