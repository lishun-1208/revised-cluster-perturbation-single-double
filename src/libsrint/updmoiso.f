      SUBROUTINE UPDMOISO(numl,numr,irrepx,list,ICORE)
C
C UPDATES THE MOIO AND MOIOWD VECTORS of tmpfile 
C icore: Memory that need to allocate the space.
C numl: the left side number of Wpqrs with symmetry = irrepx 
C numr: the right side number of Wpqrs with symmetry = irrepx 
C list: to define the spincase of tmpfile 
CEND
C      use commonmod
      use mod_size
      use mod_symm
      IMPLICIT NONE 
C-------------------
      INTEGER, external :: IRPDSO
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOSO/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C----------------------------------
      integer, intent(in) :: numl, numr, list, irrepx
      double precision, intent(inout) :: icore(*)
      real*8 itmp(ichcsz)
      LOGICAL YESNO,YESNO1,YESNO2
C
      integer ifile
      integer*8 nwords 
      integer nfull, npart, totrec, totwrd, listrec
      integer iwrite, iout, irecord 
C
C 
C ASSIGN MOIO(LSTSPN,LSTNUM) AND MOIOWD(LSTSPN,LSTNUM) AND THE FILE NUMBER
C
      IFILE=412
      MOIONUML(IRREPX, LIST)=NUML
      MOIONUMR(IRREPX, LIST)=NUMR
C
C COMPUTE TOTAL NUMBER OF WORDS IN LIST.
C
      NWORDS=NUML*NUMR
C
C COMPUTE THE NUMBER OF FULL RECORDS THIS WILL TAKE UP.  USUALLY THIS
C   WILL BE ZERO.  ALSO COMPUTE HOW MANY WORDS OF A PARTIAL RECORD ARE
C   REQUIRED (OFTEN JUST NWORDS).
C   
      NFULL=NWORDS/IPRCWD
      NPART=MOD(NWORDS,IPRCWD)
C
C INCREMENT TOTREC BY THE NUMBER OF FULL RECORDS THAT THIS WILL TAKE UP.
C
      TOTREC = NFULL 
      IF(NPART/=0) TOTREC = TOTREC+1
C
C COMPUTE POSITION OF LAST WORD WRITTEN IN PARTIAL RECORD.  IF THE LOGIC
C  RECORD GOES ACROSS A PHYSICAL RECORD BOUNDARY, INCREMENT TOTREC BY ON
C
      TOTWRD = NPART
       if(irrepx==nirrep) then 
       rectmp(1, list+1) = rectmp(irrepx, list) 
     &        + totrec
       else
       rectmp(irrepx+1, list) = rectmp(irrepx, list) 
     &        + totrec
       endif
       listrec=rectmp(irrepx, list)
C
C WRITE AN EMPTY BUFFER TO RECORD RECORD IF IT IS BEYOND EOF.
C
        CALL ZERO(itmp, ICHCSZ)
        IWRITE=listrec + TOTREC
        DO 1055 IOUT=1, totrec 
         IRECORD=IWRITE-IOUT+1
          CALL WRDIRSO(IFILE,IRECORD,itmp,ICHCSZ)
1055    CONTINUE
      RETURN
      END
