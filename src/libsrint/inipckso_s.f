      SUBROUTINE INIPCKSO_S(IRREPX, icore, POPP, POPQ, IOPPQ, 
     &                    POPR, POPS, IOPRS, LIST)
C-this is to DRIVES THE CREATION OF MOIO POINTERS FOR SYMMETRY PACKED
C-LISTS
C-IRREPX: the symmetry of the tmpfile Wpqrs with defined spin case
C-icore: Memory that need to allocate the space
C-POPP, Q, R, S: the orbital type of the four index. Like 
C--------------------POP and VRT in common block
C-IOPPQ: if IOPPQ=0, then POPP<POPQ, if IOPPQ=1, then POPP, POPQ
C-IOPPQ: if IOPPQ=2, then POPP<=POPQ
C-IOPRS: if IOPRS=0, then POPR<POPS, if IOPRS=1, then POPR, POPS
C-IOPRS: if IOPRS=2, then POPR<=POPS
C-LIST: to define the spincase of tmpfile

      use mod_symm
      use mod_orbit
      IMPLICIT none
C-------------------
      INTEGER, external :: IRPDSO
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOSO/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
C----------------------------------
      integer, intent(in) :: ioppq, ioprs 
      integer, intent(in) :: popp(8), popq(8)
      integer, intent(in) :: popr(8), pops(8)
      integer, intent(in) :: LIST, IRREPX
      real*4, intent(in) :: icore(*)
C
      integer irrepr, irrepl, NUML, NUMR
C
      DO 10 IRREPR=1,NIRREP
            IRREPL=DIRPRD(IRREPR,IRREPX)
            NUML=IRPDSO(IRREPL, POPP, POPQ, IOPPQ) 
            NUMR=IRPDSO(IRREPR, POPR, POPS, IOPRS)
!           call zeroso(icore, numl*numr)
            CALL UPDMOISO_S(NUML,NUMR,IRREPR,LIST,ICORE)
10    CONTINUE
      RETURN
      end
