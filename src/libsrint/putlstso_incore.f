      SUBROUTINE PUTLSTSO_incore(ZLIST,NFIRST,NBLOCK,
     &             irreprs0, popp, popq, ioppq, popr, pops, ioprs, ipqrs,
     &             isympqrs)
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
      integer, intent(in)::nfirst, nblock, isympqrs, irreprs0
      integer, intent(in)::popp(8), popq(8), popr(8), pops(8)
      integer, intent(in)::ioppq, ioprs 
     
      real*4, intent(inout) :: ZLIST(*), ipqrs(*)
C
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      DOUBLE PRECISION, ALLOCATABLE :: TMP(:)
C
      integer test1, test2, numpq, numrs 
      integer ioffpqrs(8) 
      integer irrepx, i, j, ifile, listrec
      integer recstart,  recend, istat
      integer irreppq, irreprs, ioff0, ioff, irreppq0
      integer iout, irecord, istart, numpq0
      integer, external :: irpdso
C
      ioffpqrs(1) = 1 
      do irreprs= 1, nirrep - 1 
         irreppq = dirprd(irreprs, isympqrs) 
         numpq = irpdso(irreppq, popp, popq, ioppq) 
         numrs = irpdso(irreprs, popr, pops, ioprs)
         ioffpqrs(irreprs+1) = ioffpqrs(irreprs) + numpq*numrs*isd
      enddo 
      
      ioff0 = ioffpqrs(irreprs0) 
      irreppq0 = dirprd(irreprs0, isympqrs) 
      numpq0 = irpdso(irreppq0, popp, popq, ioppq) 
      ioff = ioff0 - 1 + (nfirst-1)*numpq0*isd + 1 
      call myicopyso(zlist, ipqrs(ioff), numpq0*nblock)
 
      RETURN
      END
