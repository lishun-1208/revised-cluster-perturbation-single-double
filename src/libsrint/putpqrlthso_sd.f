      SUBROUTINE PUTPQRLTHSO_SD(Z, NFIRST,lth, IRREPS, LIST, filename)
C
C THIS ROUTINE RETRIEVES ALL IRREPS OF PQR IN PQRS 
C  IN THE VECTOR Z.  LENGTH IS THE TOTAL LENGTH OF THE
C  LIST, LSTNUM IS THE NUMBER OF THE LIST ON THE MOINTS FILE.
C
CEND
      use mod_symm
      use mod_orbit
      use mod_iop
      IMPLICIT NONE
C-------------------
      INTEGER :: MOIONUML, MOIONUMR, rectmp
      COMMON/MOIOTMP/MOIONUML(10,20), MOIONUMR(10,20), RECTMP(10,20)
C
      INTEGER, external :: IRPDSO, ISYMOFFSO
C----------------------------------
      INTEGER, intent(inout) :: list, lth, NFIRST,  IRREPS
      real*4, intent(inout) :: Z(*)
      character(len=8), intent(in) :: filename
C
      integer ioffpq(8), ioffpr(8), ioplist
      integer*8 ioff, tmpi1, tmpi2, ioff2
      integer irrepx, i, j, ilist, irrepr, irrepl, irreppqr
      integer numr, numl, nums, numpq, numrs, ncount
      integer irreppr, irrepqs, numpr, numqs, maxpq
      integer irreprs, irreppq, nfirstrs, nreadrs
      integer popp(8), popq(8), popr(8), pops(8)
      integer ioppq(8), ioprs(8), irrep0(16)
      double precision, allocatable:: tmpdata(:)
      character(len=2) tmpc
      character(len=10) filename2
      character*200 filepath
      integer, external :: idsymsoc

      call get_f_path(filepath)
C 
      ioplist=list
      if(list>8) then 
         ioplist = list-8 
      endif
      tmpc = ".c"
      filename2(1:8)=filename(1:8)
      filename2(9:10)=tmpc(1:2)
      OPEN(UNIT=413, FILE=trim(filepath)//trim(FILENAME2), 
     &                 FORM='UNFORMATTED')
      rewind 413
      read(413) IRREPX
      read(413) (popp(i), i=1,8)
      read(413) (popq(i), i=1,8)
      read(413) (popr(i), i=1,8)
      read(413) (pops(i), i=1,8)
      read(413) (ioppq(i), i = 1,8)
      read(413) (ioprs(i), i = 1,8)
      read(413) (irrep0(i), i = 1,8)
      read(413) (irrep0(i), i = 9,16)
      read(413) ((moionuml(i,j),i=1,10),j=1,20)
      read(413) ((moionumr(i,j),i=1,10),j=1,20)
      read(413) ((rectmp(i,j),i=1,10),j=1,20)
      close(413)
      maxpq = 0
      do irrepr = 1, nirrep 
         numpq = irpdso(irrepr, popp, popq, 1) 
         maxpq = max(maxpq, numpq) 
      enddo 
      allocate(tmpdata(maxpq))
C----I will first try to get PQRS term.
C----With symmetry of S be defined, the symmetry of PQR block is 
C----also determined. i will first resort PQR and then sstgen pqr
C----to prq
C      write(6,*)'begin get pqr'
!      irrepx = dirprd(irrepx, irrep0(list))
      irrepx = irrep0(list)
C      write(6,*)'irrepx', irrepx
      ioff = 1
      irreppqr = dirprd(irreps, irrepx)
      do i = 1, lth 
!        write(6,*)'ioff', ioff
         call putpqrso_sd(z(ioff), nfirst-1+i, irreps, list, filename)
         ioff = ioff + idsymsoc(irreppqr, popp, popq, ioppq(list),
     &                                    popr, zorb,1)*isd
      enddo 
      deallocate(tmpdata)
      RETURN
      END
