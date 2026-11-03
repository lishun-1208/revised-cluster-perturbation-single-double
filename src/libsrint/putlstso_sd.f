      SUBROUTINE PUTLSTSO_SD(ZLIST,NFIRST,NBLOCK,
     &                           ISPIN,ILIST,filename)
C
C THIS SUBROUTINE RETRIEVES NBLOCK DISTRIBUTIONS OF A PARTICULAR INTEGRA
C  LIST, BEGINNING WITH DISTRIBUTION NFIRST.  THE VALUES ARE RETURNED IN
C  VECTOR ZLIST.
C
C   OUTPUT:
C
C         ZLIST - THE TARGET VECTOR USED IN I/O OPERATIONS.
C
C   INPUT PARAMETERS:
C
C         NFIRST- THE NUMBER OF THE FIRST DISTRIBUTION WHICH WILL BE
C                 WRITTEN TO OR READ FROM DISK.
C         NBLOCK- THE NUMBER OF DISTRIBUTIONS WHICH WILL BE HANDLED.
C         USECIN- THE CACHE BUFFER WHICH IS TO BE USED FOR I/O OPERATION
C                 BEST PERFORMANCE IS ACHIEVED IF THIS IS SET TO 1 FOR
C                 T AMPLITUDE I/O AND 2 FOR INTEGRAL I/O.
C         ISPIN - THE LEFT-HAND INDEX OF THE LIST [MOIO(ISPIN,ITYPE)].
C         ILIST - THE RIGHT-HAND INDEX OF THE LIST [MOIO(ISPIN,ITYPE)].
C
CEND
      use mod_iop
      use mod_symm
      use mod_size
      IMPLICIT NONE
      integer, intent(in)::nfirst, nblock, ispin, ilist 
      character(len=8), intent(in) :: filename
C
      double precision, intent(inout) :: ZLIST(*) 
C     
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      double precision, allocatable :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
C
      integer*8 nwords, nstart, nend, nwrdf, nwrdf2
      integer irrepx, i, j, ifile, listrec
      integer recstart, recend, istat
      integer iout, irecord, istart
      integer*8 tmpi1, tmpi2
      character*200 filepath

      if(prec==8) then
       call PUTLSTSO_D(ZLIST,NFIRST,NBLOCK,
     &                           ISPIN,ILIST,filename)
      elseif(prec==4) then 
       call PUTLSTSO_S(ZLIST,NFIRST,NBLOCK,
     &                           ISPIN,ILIST,filename)
      endif
300   RETURN
      END
