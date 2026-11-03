      SUBROUTINE UPDMOIDIIS_S(icore, numl, numr, list, lsize, dsize)
C
C UPDATES THE MOIO AND MOIOWD VECTORS of tmpfile 
C icore: Memory that need to allocate the space.
C numl: the left side number of Wpqrs with symmetry = irrepx 
C numr: the right side number of Wpqrs with symmetry = irrepx 
C list: to define the spincase of tmpfile 
CEND
C      use commonmod
      IMPLICIT NONE
      integer*8, intent(in) :: numl 
      integer, intent(in) :: numr, list, lsize, dsize 
      real*4, intent(inout) :: icore(*)
C-------------------
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(1000),MOIONUMR(1000),RECTMP(1000)
C----------------------------------
      integer ifile
      integer*8 nwords 
      integer nfull, npart, totrec, totwrd, listrec
      integer iwrite, iout, irecord 
C
C 
C ASSIGN MOIO(LSTSPN,LSTNUM) AND MOIOWD(LSTSPN,LSTNUM) AND THE FILE NUMBER
C
      IFILE=412
C
C COMPUTE TOTAL NUMBER OF WORDS IN LIST.
C
      NWORDS=NUML*NUMR
C      write(6,*)'num in upddiis', numl ,numr, nwords
C
C COMPUTE THE NUMBER OF FULL RECORDS THIS WILL TAKE UP.  USUALLY THIS
C   WILL BE ZERO.  ALSO COMPUTE HOW MANY WORDS OF A PARTIAL RECORD ARE
C   REQUIRED (OFTEN JUST NWORDS).
C
      NFULL=NWORDS/dsize
      NPART=MOD(NWORDS,dsize)
C
C INCREMENT TOTREC BY THE NUMBER OF FULL RECORDS THAT THIS WILL TAKE UP.
C
      TOTREC = 1+NFULL
C
C COMPUTE POSITION OF LAST WORD WRITTEN IN PARTIAL RECORD.  IF THE LOGIC
C  RECORD GOES ACROSS A PHYSICAL RECORD BOUNDARY, INCREMENT TOTREC BY ON
C
      TOTWRD = NPART
       if(list<lsize) then 
          rectmp(list+1) = rectmp(list)+totrec
       endif
       listrec=rectmp(list)
C
C WRITE AN EMPTY BUFFER TO RECORD RECORD IF IT IS BEYOND EOF.
C
       CALL ZERO_S(ICORE, dsize)
C
       IWRITE=listrec + TOTREC
C       write(6,*)'totrec', totrec
        DO 1055 IOUT=1, totrec 
           IRECORD=IWRITE-IOUT+1
           CALL WRDIRSO_S(IFILE,IRECORD,ICORE,dsize)
1055    CONTINUE
      RETURN
      END
