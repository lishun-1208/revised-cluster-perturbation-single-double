      SUBROUTINE CLMOIOSO0(IRREPX, ICORE, 
     &           POPP, POPQ, IOPPQIN, POPR, POPS, IOPRSIN, filename)
C-IRREPX: the symmetry of the tmpfile W(P,Q,R,S) 
C-ICORE: Memory that need to allocate the space 
C-POPP, Q, R, S: the orbital type of the four index. Like POP and VRT in
C--------------- common block
C-filename: the filename to save the W(P,Q,R,S), one should not that the
Clength of the filename must be (length=8)
C-IOPPQ and IOPRS: to define P<Q or P<=Q or P,Q
      use mod_symm
      use mod_orbit
      use mod_size
      IMPLICIT NONE
C-------------------
      INTEGER, external :: IRPDSO
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOSO/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C----------------------------------

      integer, intent(in) :: popp(8), popq(8)
      integer, intent(in) :: popr(8), pops(8)
      integer, intent(in) :: irrepx, ioppqin, ioprsin
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
      LOGICAL DOALL,ABIJ,NOABCD
C
      INTEGER IRREPR(8), IRREPI(8)
      integer ioppq(8), ioprs(8), i, j
      character*2 tmpc
      character*10 filename2
      integer ilength
      character*200 filepath

      call get_f_path(filepath)
      OPEN(UNIT=412, FILE=trim(filepath)//trim(filename), 
     &           FORM='UNFORMATTED', 
     &            access='direct', recl=iprcln)
      CALL IZERO(MOIONUML, 200)
      CALL IZERO(MOIONUMR, 200)
      CALL IZERO(rectmp, 200)
C***********************************************************************
        IRREPR(1) = IRREPX 
        ioppq(1) = ioppqin; ioprs(1) = ioprsin
        CALL INIPCKSO(IRREPX,ICORE,POPP,POPQ,IOPPQ(1), 
     &                                POPR,POPS,IOPRS(1),1)
C
      CLOSE(412)
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=trim(filepath)//trim(filename2), 
     &                 FORM='UNFORMATTED',status='unknown')
      rewind 413 
      write(413) IRREPX
      write(413) (popp(i), i=1,8)
      write(413) (popq(i), i=1,8)
      write(413) (popr(i), i=1,8)
      write(413) (pops(i), i=1,8)
      write(413) (ioppq(i), i = 1,8)
      write(413) (ioprs(i), i = 1,8)
      write(413) (irrepr(i), i = 1,8)
      write(413) (irrepi(i), i = 1,8)
      write(413) ((moionuml(i,j),i=1,10),j=1,20)
      write(413) ((moionumr(i,j),i=1,10),j=1,20)
      write(413) ((rectmp(i,j),i=1,10),j=1,20)
      close(413)
C***********************************************************************
      RETURN
      END
