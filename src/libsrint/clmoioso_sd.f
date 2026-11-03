      SUBROUTINE CLMOIOSO_SD(IRREPX, ICORE, 
     &           POPP, POPQ, IOPPQIN, POPR, POPS, IOPRSIN, filename)
C-IRREPX: the symmetry of the tmpfile W(P,Q,R,S) 
C-ICORE: Memory that need to allocate the space 
C-POPP, Q, R, S: the orbital type of the four index. Like POP and VRT in
C--------------- common block
C-filename: the filename to save the W(P,Q,R,S), one should not that the
Clength of the filename must be (length=8)
C-IOPPQ and IOPRS: to define P<Q or P<=Q or P,Q
      use mod_iop
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

      if(prec==4) then 
       call CLMOIOSO0_s(IRREPX, ICORE, 
     &           POPP, POPQ, IOPPQIN, POPR, POPS, IOPRSIN, filename)
      elseif(prec==8) then 
       call CLMOIOSO0(IRREPX, ICORE, 
     &           POPP, POPQ, IOPPQIN, POPR, POPS, IOPRSIN, filename)
      endif
      return
      end
