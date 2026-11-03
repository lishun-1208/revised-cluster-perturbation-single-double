      SUBROUTINE GETLSTSO_SD(ZLIST,NFIRST,NBLOCK,
     &                     ISPIN,ILIST,filename)
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
C     NFIRST- THE NUMBER OF THE FIRST DISTRIBUTION WHICH WILL BE
C             WRITTEN TO OR READ FROM DISK.
C     NBLOCK- THE NUMBER OF DISTRIBUTIONS WHICH WILL BE HANDLED.
C             T AMPLITUDE I/O AND 2 FOR INTEGRAL I/O.
C     ISPIN - THE SYMMETRY OF THE RIGHT SIDE (IRREPRS) OF WPQRS 
C     ILIST - THE LIST of the saved Wpqrs with defined spincase. 
C
CEND
      use mod_iop
      use mod_symm
      use mod_orbit
      use mod_size
      IMPLICIT NONE 
      integer, intent(in)::nfirst, nblock, ispin, ilist 
      character(len=8), intent(in) :: filename
     
      double precision, intent(inout) :: ZLIST(*)
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      DOUBLE PRECISION, ALLOCATABLE :: TMP(:)
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      character(len=2) tmpc
      character(len=10) filename2
C
      integer(kind=8) :: nwords, nstart, nend, nwrdf, nwrdf2
      integer*8 test1, test2
      integer irrepx, i, j, ifile, listrec
      integer recstart,  recend, istat
      integer iout, irecord, istart
      character*200 filepath

      if(prec==8) then
        call GETLSTSO_D(ZLIST,NFIRST,NBLOCK,
     &                     ISPIN,ILIST,filename)
      elseif(prec==4) then 
        call GETLSTSO_S(ZLIST,NFIRST,NBLOCK,
     &                     ISPIN,ILIST,filename)
      endif
C write(6,*)'========end of getlst======='
      RETURN
      END
