      SUBROUTINE CLMOIOSO0_S(IRREPX, ICORE, 
     &           POPP, POPQ,ioppqin, POPR, POPS,ioprsin, filename)
C-index of the file, 412.
C-IRREPX: the symmetry of the tmpfile W(P,Q,R,S) 
C-ICORE: Memory that need to allocate the space 
C-POPP, Q, R, S: the orbital type of the four index. Like POP and VRT in
C--------------- common block
C-filename: the filename to save the W(P,Q,R,S), one should not that the
Clength of 
C--------------- the filename must be (length=8)
C-IOP: to define the type of Wpqrs. 
C-iop==0, then --for P=Q and R=S like Wmnij and Wabef
C----------------------------R-----------I
C IPQRSAAAA       1      1+8
C IPQRSABAB       2      2+8
C IPQRSAAAB       3      3+8
C IPQRSABAA       4      4+8
C IPQRSAABB       5      5+8
C-iop==1, then ---for P=Q and R/=S like Wmnie and Wabei
C IPQRSAABA       6      6+8
C-iop==2, then ---for P/=Q and R/=S like Wmbej
C IPQRSABBA       7      7+8
C IPQRSBAAA       8      8+8
C
C-iop negative for scalcar relativistic calculation
C-iop== -1,    ---for P=Q and R=S
C-iop== -2,    ---for P=Q and R/=S
C-iop== -3,    ---for P/=Q and R/=S
C IPQRSAAAA       1      
C IPQRSABAB       2  
C IPQRSABBA       3   
C
C-iop =-10, just for assym of I<Ab||Cd> and T2(Ab,Ij) 
C   A<=B, and I<=J
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
      real*4, intent(in) :: icore(*)
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
     &            FORM='UNFORMATTED', access='direct', recl=iprcln/2)
!     OPEN(UNIT=412, FILE=filename, FORM='UNFORMATTED', 
!    &            access='direct', recl=iprcln/2)
      CALL IZERO(MOIONUML, 200)
      CALL IZERO(MOIONUMR, 200)
      CALL IZERO(rectmp, 200)
C***********************************************************************
        IRREPR(1) = IRREPX 
        ioppq(1) = ioppqin; ioprs(1) = ioprsin
        CALL INIPCKSO_S(IRREPX,ICORE,POPP,POPQ,IOPPQ(1), 
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
