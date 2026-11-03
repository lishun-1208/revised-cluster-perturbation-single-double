      SUBROUTINE GETPQRSO_INCORE(Z, NFIRST, IRREPS, popp, popq, ioppq,
     &               popr, pops, ioprs, ipqrs, isympqrs)
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
      INTEGER, external :: IRPDSO, ISYMOFFSO
C----------------------------------
      INTEGER, intent(inout) :: NFIRST,  IRREPS
      DOUBLE PRECISION, intent(inout) :: Z(*)
      integer, intent(inout)::popp(8), popq(8), popr(8), pops(8)
      integer, intent(inout)::ioppq, ioprs, isympqrs
      real*8, intent(inout) :: ipqrs(*)
C
      integer ioffpq(8), ioffpr(8), ioplist
      integer*8 ioff, tmpi1, tmpi2, ioff2
      integer irrepx, i, j, ilist, irrepr, irrepl
      integer numr, numl, nums, numpq, numrs, ncount
      integer irreppr, irrepqs, numpr, numqs, maxpq
      integer irreprs, irreppq, nfirstrs, nreadrs
      double precision, allocatable:: tmpdata(:)
      character(len=2) tmpc
      character(len=10) filename2

      if(prec==8) then 
         call GETPQRSOD_INCORE(Z, NFIRST, IRREPS, popp, popq, ioppq,
     &               popr, pops, ioprs, ipqrs, isympqrs)
      elseif(prec==4) then 
         call GETPQRSOS_INCORE(Z, NFIRST, IRREPS, popp, popq, ioppq,
     &               popr, pops, ioprs, ipqrs, isympqrs)
      endif
C 
      RETURN
      END
