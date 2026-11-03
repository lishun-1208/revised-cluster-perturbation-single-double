      SUBROUTINE GETALLSO_S(Z,list,filename)
C
C THIS ROUTINE RETRIEVES ALL IRREPS OF A SYMMETRY PACKED DISTRIBUTION
C  AND PLACES THEM IN THE VECTOR Z.  LENGTH IS THE TOTAL LENGTH OF THE
C  LIST, LSTNUM IS THE NUMBER OF THE LIST ON THE MOINTS FILE.
C
CEND
      use mod_symm
      use mod_orbit
      IMPLICIT NONE
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      INTEGER, external :: IRPDSO
C----------------------------------
      INTEGER, intent(in) :: list
      real*4, intent(inout) :: Z(*)
      character(len=8), intent(in) :: filename
C
      integer*8 ioff, tmpi1, tmpi2
      integer irrepx, i, j, ilist, irrepr, irrepl
      integer numr, numl
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
      character*200 filepath

      call get_f_path(filepath)
C
      tmpc = ".c"
      filename2(1:8) = filename(1:8)
      filename2(9:10) = tmpc(1:2)
C      write(6,*)'begin read',' ', trim(filename2),' in getalltmp'
      OPEN(UNIT=413, FILE=trim(filepath)//trim(FILENAME2), 
     &               FORM='UNFORMATTED')
      rewind 413
      read(413) IRREPX
      read(413) (popp(i), i=1,8)
      read(413) (popq(i), i=1,8)
      read(413) (popr(i), i=1,8)
      read(413) (pops(i), i=1,8)
      read(413) (ioppq(i), i = 1,8)
      read(413) (ioprs(i), i = 1,8)
      read(413) (irrep0(i), i = 1,8)
      read(413) (irrep0(i), i = 9,16)
      read(413) ((moionuml(i,j),i=1,10),j=1,20)
      read(413) ((moionumr(i,j),i=1,10),j=1,20)
      read(413) ((rectmp(i,j),i=1,10),j=1,20)
      close(413)
C      write(6,*)'end read',' ', trim(filename2),' in getalltmp'
C
C  LOOP OVER IRREPS.
C
      IOFF=1
      ILIST=list
      if(ilist>8) ilist = list-8
      DO 10 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREP0(list))
C
C GET NUMBER OF DISTRIBUTIONS AND DISTRIBUTION SIZE.
C
       NUMR=IRPDSO(IRREPR, POPR, POPS, IOPRS(ilist))
       NUML=IRPDSO(IRREPL, POPP, POPQ, IOPPQ(ilist))
       tmpi1 = numl; tmpi2 = numr
C       write(6,*)'irrepr', IRREPR, POPR(1), POPS(1) 
C       write(6,*)'irrepl', IRREPl, POPP(1), POPQ(1) 
C       write(6,*)'ioprs and ioppq', ilist, ioprs(ilist), ioppq(ilist)
C       write(6,*)'numdis, dissiz', numl, numr
C
C STUFF VECTOR INTO Z.
C
       IF(NUMR.EQ.0)GOTO 10
C       write(6,*)'begin getlsttmp', numl*numr
       CALL GETLSTSO_S(Z(IOFF),1,NUMR,IRREPR,list,filename)
C-----This step is because that numl and numr are both integer*4 so that 
C-----numl*numr are also a number with integer*4. If numl*numr is too
C-----large, some error will occur.
C-----tmpi1 and tmpi2 are both integer*8
       
       IOFF=IOFF+tmpi1*tmpi2  
10    CONTINUE
      RETURN
      END
