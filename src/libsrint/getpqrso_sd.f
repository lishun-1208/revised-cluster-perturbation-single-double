      SUBROUTINE GETPQRSO_SD(Z, NFIRST, IRREPS, LIST, filename)
C
C THIS ROUTINE RETRIEVES ALL IRREPS OF PQR IN PQRS 
C  IN THE VECTOR Z.  LENGTH IS THE TOTAL LENGTH OF THE
C  LIST, LSTNUM IS THE NUMBER OF THE LIST ON THE MOINTS FILE.
C
CEND
      use mod_iop
      use mod_symm
      use mod_orbit
      IMPLICIT NONE
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      INTEGER, external :: IRPDSO, ISYMOFFSO
C----------------------------------
      INTEGER, intent(inout) :: list, NFIRST,  IRREPS
      DOUBLE PRECISION, intent(inout) :: Z(*)
      character(len=8), intent(in) :: filename
C
      integer ioffpq(8), ioffpr(8), ioplist
      integer*8 ioff, tmpi1, tmpi2, ioff2
      integer irrepx, i, j, ilist, irrepr, irrepl
      integer numr, numl, nums, numpq, numrs, ncount
      integer irreppr, irrepqs, numpr, numqs, maxpq
      integer irreprs, irreppq, nfirstrs, nreadrs
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      double precision, allocatable:: tmpdata(:)
      character(len=2) tmpc
      character(len=10) filename2
C
      if(prec==8) then 
         call GETPQRSO_d(Z, NFIRST, IRREPS, LIST, filename)
      elseif(prec==4) then 
         call GETPQRSO_s(Z, NFIRST, IRREPS, LIST, filename)
      endif
      return
      end
