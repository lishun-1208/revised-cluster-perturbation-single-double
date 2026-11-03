      SUBROUTINE CLMOIODIIS_S(ICORE,NUML,NUMR, IOPSIZE, filename, dsize)
C-index of the file, 412.
      IMPLICIT NONE
      integer*8, intent(in) :: NUML
      integer, intent(in) :: NUMR, IOPSIZE, dsize
      real*4, intent(in) :: icore(*)
      character*8, intent(in) :: filename
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(1000),MOIONUMR(1000),RECTMP(1000)
C
      INTEGER :: SYMRX, SYMRY, SYMRZ
      COMMON/SYMRXYZ/ SYMRX, SYMRY, SYMRZ
C
      INTEGER :: POP, VRT, NT, NF1, NF2
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      integer iop, i 
      character*2 tmpc
      character*10 filename2
C
C----------------------------------
C
      character*200 filepath

      call get_f_path(filepath)
C      write(6,*)'num', numl, numr, iopsize, dsize
      OPEN(UNIT=412, FILE=trim(filepath)//trim(FILENAME),
     &            FORM='UNFORMATTED', access='direct',
     &            recl=dsize*4)
    
      CALL IZERO(MOIONUML, iopsize)
      CALL IZERO(MOIONUMR, iopsize)
      CALL IZERO(rectmp, iopsize)
C      write(6,*)'clmoioso by wz'
C***********************************************************************
      do iop = 1, iopsize
         CALL UPDMOIDIIS_S(ICORE,NUML, NUMR,iop, iopsize, dsize)
      enddo
C
      CLOSE(412)
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=trim(filepath)//trim(filename2),
     &                 FORM='UNFORMATTED',status='unknown')
      rewind 413 
      write(413) iopsize 
      write(413) dsize 
      write(413) numl
      write(413) numr
      write(413) (rectmp(i),i=1,iopsize)
      close(413)
C***********************************************************************
      RETURN
      END
