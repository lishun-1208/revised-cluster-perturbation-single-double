      SUBROUTINE CLMOIOSO(IRREPX, ICORE, 
     &           POPP, POPQ, POPR, POPS, filename, iop)
      use mod_orbit
      use mod_symm
      use mod_size
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
C-------------- --R-------I
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
C
C-iop = 11, for AAR, ABR, AAI, ABI
C IPQAA       1      1+2
C IPQAB       2      2+2
C

      IMPLICIT NONE
C-------------------
      INTEGER, external :: IDSYMSZ
      INTEGER, external :: IRPDSO
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOSO/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
!----------------------------------

      integer, intent(in) :: popp(8,2), popq(8,2)
      integer, intent(in) :: popr(8,2), pops(8,2)
      integer, intent(in) :: iop, irrepx
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
      LOGICAL DOALL,ABIJ,NOABCD
C
      INTEGER IRREPR(8), IRREPI(8)
      integer ioppq(8), ioprs(8), i, j
      character*2 tmpc
      character*10 filename2
      character*200 filepath
      integer ilength
!      character*200 filepath

!      call get_f_path(filepath)
C
!      ZORB(1:8) = 0
!      ZORB(1) = 1
      call get_f_path(filepath)
      OPEN(UNIT=412, FILE=trim(filepath)//trim(filename), 
     &            FORM='UNFORMATTED', access='direct', recl=iprcln)
      CALL IZERO(MOIONUML, 200)
      CALL IZERO(MOIONUMR, 200)
      CALL IZERO(rectmp, 200)
C***********************************************************************
C iop==0, then ---for P=Q  and R=S like Wmnij and Wabef
C 1 IPQRSAAAA     
C 2 IPQRSABAB
C 3 IPQRSAAAB     
C 4 IPQRSABAA
C 5 IPQRSAABB
      if(iop==0.or.iop==1.or.iop==2) then 
       ioppq(1) = 0; ioprs(1) = 0
       ioppq(2) = 1; ioprs(2) = 1
       ioppq(3) = 0; ioprs(3) = 1 
       ioppq(4) = 1; ioprs(4) = 0 
       ioppq(5) = 0; ioprs(5) = 0
C iop==1, then ---for P=Q  and R/=S like Wmnie and Wabei
C 6 IPQRSAABA
      if(iop==1) then  
        ioprs(1) = 1; ioprs(5) = 1; ioprs(4) = 1 
        ioppq(6) = 0; ioprs(6) = 1
C iop==2, then ---for P/=Q and R/=S like Wmbej
C 7 IPQRSABBA
C 8 IPQRSBAAA
      elseif(iop==2) then
        ioppq(1) = 1; ioppq(3) = 1; ioppq(5) = 1
        ioprs(1) = 1; ioprs(5) = 1; ioprs(4) = 1 
        ioppq(6) = 1; ioprs(6) = 1
        ioppq(7) = 1; ioprs(7) = 1
        ioppq(8) = 1; ioprs(8) = 1
      endif
C=======
      IRREPR(1)=IRREPX;               IRREPI(1)=DIRPRD(SYMRZ, IRREPR(1))
      IRREPR(2)=IRREPX;               IRREPI(2)=DIRPRD(SYMRZ, IRREPR(2))
      IRREPR(3)=DIRPRD(IRREPX, SYMRY);IRREPI(3)=DIRPRD(SYMRZ, IRREPR(3))
      IRREPR(4)=DIRPRD(IRREPX, SYMRY);IRREPI(4)=DIRPRD(SYMRZ, IRREPR(4))
      IRREPR(5)=IRREPX;               IRREPI(5)=DIRPRD(SYMRZ, IRREPR(5))
      IRREPR(6)=DIRPRD(IRREPX, SYMRY);IRREPI(6)=DIRPRD(SYMRZ, IRREPR(6))
      IRREPR(7)=IRREPX;               IRREPI(7)=DIRPRD(SYMRZ, IRREPR(7))
      IRREPR(8)=DIRPRD(IRREPX, SYMRY);IRREPI(8)=DIRPRD(SYMRZ, IRREPR(8))
C      write(6,*)'sym', symrx, symry, symrz
C      write(6,*)'irrepr in clmoioso', filename, irrepr
C      write(6,*)'irrepi in clmoioso', filename, irrepi
       CALL INIPCKSO(IRREPR(1),ICORE,POPP,POPQ,IOPPQ(1), 
     &                             POPR,POPS,IOPRS(1),1)
       CALL INIPCKSO(IRREPR(2),ICORE,POPP,POPQ,IOPPQ(2), 
     &                             POPR,POPS,IOPRS(2),2)
       CALL INIPCKSO(IRREPR(3),ICORE,POPP,POPQ,IOPPQ(3), 
     &                             POPR,POPS,IOPRS(3),3)
       CALL INIPCKSO(IRREPR(4),ICORE,POPP,POPQ,IOPPQ(4), 
     &                             POPR,POPS,IOPRS(4),4)
       CALL INIPCKSO(IRREPR(5),ICORE,POPP,POPQ,IOPPQ(5), 
     &                                POPR,POPS,IOPRS(5),5)
C     write(6,*)'rectmp in clmoiotmp', rectmp(1,6)
       if(iop==0) then 
         rectmp(1, 9)=rectmp(1,6)
         goto 100
       endif
       CALL INIPCKSO(IRREPR(6),ICORE,POPP,POPQ,IOPPQ(6), 
     &                                POPR,POPS,IOPRS(6),6)
       if(iop==1) then 
         rectmp(1,9)=rectmp(1,7) 
         goto 100 
       endif
       CALL INIPCKSO(IRREPR(7),ICORE,POPP,POPQ,IOPPQ(7), 
     &                                POPR,POPS,IOPRS(7),7)
       CALL INIPCKSO(IRREPR(8),ICORE,POPP,POPQ,IOPPQ(8), 
     &                                POPR,POPS,IOPRS(8),8)
C
100    CALL INIPCKSO(IRREPI(1),ICORE,POPP,POPQ,IOPPQ(1), 
     &                                POPR,POPS,IOPRS(1),9)
       CALL INIPCKSO(IRREPI(2),ICORE,POPP,POPQ,IOPPQ(2), 
     &                                POPR,POPS,IOPRS(2),10)
       CALL INIPCKSO(IRREPI(3),ICORE,POPP,POPQ,IOPPQ(3), 
     &                                POPR,POPS,IOPRS(3),11)
       CALL INIPCKSO(IRREPI(4),ICORE,POPP,POPQ,IOPPQ(4), 
     &                                POPR,POPS,IOPRS(4),12)
       CALL INIPCKSO(IRREPI(5),ICORE,POPP,POPQ,IOPPQ(5), 
     &                                POPR,POPS,IOPRS(5),13)
       if(iop==0) goto 200
       CALL INIPCKSO(IRREPI(6),ICORE,POPP,POPQ,IOPPQ(6), 
     &                                POPR,POPS,IOPRS(6),14)
       if(iop==1) goto 200 
       CALL INIPCKSO(IRREPI(7),ICORE,POPP,POPQ,IOPPQ(7), 
     &                                POPR,POPS,IOPRS(7),15)
       CALL INIPCKSO(IRREPI(8),ICORE,POPP,POPQ,IOPPQ(8), 
     &                                POPR,POPS,IOPRS(8),16)
200    endif
C----for AA-two-index
      if(iop==11) then
        ioppq(1) = 1; ioprs(1) = 1
        ioppq(2) = 1; ioprs(2) = 1 
        ioppq(3) = 1; ioprs(3) = 1 
        ioppq(4) = 1; ioprs(4) = 1 
      IRREPR(1)=IRREPX;               IRREPR(3)=DIRPRD(SYMRZ, IRREPR(1))
      IRREPR(2)=DIRPRD(IRREPX, SYMRY);IRREPR(4)=DIRPRD(SYMRZ, IRREPR(2))
       CALL INIPCKSO(IRREPR(1),ICORE,POPP,POPQ,IOPPQ(1), 
     &                             ZORB,ZORB,IOPRS(1),1)
       CALL INIPCKSO(IRREPR(2),ICORE,POPP,POPQ,IOPPQ(2), 
     &                             ZORB,ZORB,IOPRS(2),2)
       CALL INIPCKSO(IRREPR(3),ICORE,POPP,POPQ,IOPPQ(3), 
     &                             ZORB,ZORB,IOPRS(3),3)
       CALL INIPCKSO(IRREPR(4),ICORE,POPP,POPQ,IOPPQ(4), 
     &                             ZORB,ZORB,IOPRS(4),4)
      endif

C----for scalar relativistic   
C----iop==-1 p=q,r=s
C----iop==-2 p=q,r/=s
C----iop==-3 p/=q,r/=s
      if(iop==-1.or.iop==-2.or.iop==-3) then 
      IRREPR(1)=IRREPX      
      IRREPR(2)=IRREPX      
      IRREPR(3)=IRREPX      
      if(iop==-1) then 
        ioppq(1) = 0; ioprs(1) = 0
        ioppq(2) = 1; ioprs(2) = 1 
      elseif(iop==-2) then       
        ioppq(1) = 0; ioprs(1) = 1
        ioppq(2) = 1; ioprs(2) = 1
      elseif(iop==-3) then       
        ioppq(1) = 1; ioprs(1) = 1
        ioppq(2) = 1; ioprs(2) = 1 
        ioppq(3) = 1; ioprs(3) = 1 
      endif
       CALL INIPCKSO(IRREPX,ICORE,POPP,POPQ,IOPPQ(1), 
     &                                POPR,POPS,IOPRS(1),1)
       CALL INIPCKSO(IRREPX,ICORE,POPP,POPQ,IOPPQ(2), 
     &                                POPR,POPS,IOPRS(2),2)
       if(iop==-3) then 
       CALL INIPCKSO(IRREPX,ICORE,POPP,POPQ,IOPPQ(3), 
     &                                POPR,POPS,IOPRS(3),3)
       endif 
      endif
C----for scalar relativistic uhf
C----AAAA, ABAB, BBBB, BABA, AABB, BBAA
C----iop==-11 p=q,r=s
C----iop==-22 p=q,r/=s
C----iop==-33 p/=q,r/=s
      if(iop==-11.or.iop==-22.or.iop==-33) then
      IRREPR(1)=IRREPX
      IRREPR(2)=IRREPX
      IRREPR(3)=IRREPX
      if(iop==-11) then
        ioppq(1) = 0; ioprs(1) = 0
        ioppq(2) = 1; ioprs(2) = 1
        ioppq(3) = 0; ioprs(3) = 0
      elseif(iop==-22) then
        ioppq(1) = 0; ioprs(1) = 1
        ioppq(2) = 1; ioprs(2) = 1
        ioppq(3) = 0; ioprs(3) = 1
        ioppq(4) = 1; ioprs(4) = 1
      elseif(iop==-33) then
        ioppq(1) = 1; ioprs(1) = 1
        ioppq(2) = 1; ioprs(2) = 1
        ioppq(3) = 1; ioprs(3) = 1
        ioppq(4) = 1; ioprs(4) = 1
        ioppq(5) = 1; ioprs(5) = 1
        ioppq(6) = 1; ioprs(6) = 1
      endif
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,1),POPQ(1,1),IOPPQ(1),
     &                            POPR(1,1),POPS(1,1),IOPRS(1),1)
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,1),POPQ(1,2),IOPPQ(2),
     &                            POPR(1,1),POPS(1,2),IOPRS(2),2)
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,2),POPQ(1,2),IOPPQ(3),
     &                            POPR(1,2),POPS(1,2),IOPRS(3),3)
       if(iop==-22.or.iop==-33) then
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,2),POPQ(1,1),IOPPQ(4),
     &                            POPR(1,2),POPS(1,1),IOPRS(4),4)
       endif
       if(iop==-33) then
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,1),POPQ(1,1),IOPPQ(5),
     &                            POPR(1,2),POPS(1,2),IOPRS(5),5)
       CALL INIPCKSO(IRREPX,ICORE,POPP(1,2),POPQ(1,2),IOPPQ(6),
     &                            POPR(1,1),POPS(1,1),IOPRS(6),6)
       endif
      endif
C
      CLOSE(412)
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=trim(filepath)//trim(filename2), 
     &                 FORM='UNFORMATTED',status='unknown')
      rewind 413 
      write(413) IRREPX
      write(413) ((popp(i,j), i=1,8),j=1,2)
      write(413) ((popq(i,j), i=1,8),j=1,2)
      write(413) ((popr(i,j), i=1,8),j=1,2)
      write(413) ((pops(i,j), i=1,8),j=1,2)
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
