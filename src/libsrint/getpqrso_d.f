      SUBROUTINE GETPQRSO_D(Z, NFIRST, IRREPS, LIST, filename)
C
C THIS ROUTINE RETRIEVES ALL IRREPS OF PQR IN PQRS 
C  IN THE VECTOR Z.  LENGTH IS THE TOTAL LENGTH OF THE
C  LIST, LSTNUM IS THE NUMBER OF THE LIST ON THE MOINTS FILE.
C
CEND
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
      character*200 filepath

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
      if(ioprs(ioplist)==1) then 
      ioff = 1
      do irrepr = 1, nirrep 
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq(ioplist)) 
         nums = pops(irreps) 
         numr = popr(irrepr) 
          nfirstrs = isymoffso(irreps,irreprs, popr, pops, ioprs(ioplist))
     &                 -1+(nfirst-1)*numr + 1
          nreadrs = numr 
            call getlstso(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, list, filename) 
           ioff = ioff + numpq*numr
      enddo 
C
C----to r<s
      elseif(ioprs(ioplist)==2) then
      ioff = 1
      ncount = 1
      do irrepr = 1, nirrep
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq(ioplist))
         nums = pops(irreps)
         numr = popr(irrepr)
         if(irrepr == irreps) then 
           nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)-1
     &             +(nfirst-1)*nfirst/2 + 1
           nreadrs = nfirst 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
           call getlstso(z(ioff), nfirstrs, nreadrs,
     &                irreprs, list, filename)
         do i = nfirst+1, nums
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)-1 
     &                + i*(i-1)/2 + nfirst 
            nreadrs = 1
            ioff2 = ioff + (i-1)*numpq 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
            call getlstso(z(ioff2), nfirstrs, nreadrs,
     &                    irreprs, list, filename)
          enddo
          elseif(irrepr <irreps) then 
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)- 1
     &              +(nfirst-1)*numr + 1
            nreadrs = numr 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
            call getlstso(z(ioff), nfirstrs, nreadrs,
     &                irreprs, list, filename)
          elseif(irrepr > irreps) then 
            do i = 1, numr 
               nfirstrs = isymoffso(irrepr, irreprs, pops,popr,2)-1
     &           + (i-1)*nums + nfirst
               nreadrs = 1
           ncount = ncount + nreadrs
               ioff2 = ioff + (i-1)*numpq 
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
               call getlstso(z(ioff2), nfirstrs, nreadrs,
     &                    irreprs, list, filename)
            enddo
          endif  
          ioff = ioff + numpq*numr
      enddo
C
      elseif(ioprs(ioplist)==0) then
C     write(6,*)'test here', irrepx
      ioff = 1
      ncount = 1
      do irrepr = 1, nirrep
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq(ioplist))
         nums = pops(irreps)
         numr = popr(irrepr)
         if(irrepr == irreps) then
           nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)-1
     &             +(nfirst-1)*(nfirst-2)/2 + 1
           nreadrs = nfirst-1
           ncount = ncount + nreadrs
           call getlstso(z(ioff), nfirstrs, nreadrs,
     &                irreprs, list, filename)
         i = nfirst 
         ioff2 = ioff + (i-1)*numpq
         call zero(z(ioff2), numpq)
         do i = nfirst+1, nums
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)-1
     &                + (i-2)*(i-1)/2 + nfirst
            nreadrs = 1
            ioff2 = ioff + (i-1)*numpq
           ncount = ncount + nreadrs
            call getlstso(z(ioff2), nfirstrs, nreadrs,
     &                    irreprs, list, filename)
            call dscal(numpq, -1.d0, z(ioff2), 1)
          enddo
          elseif(irrepr <irreps) then
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)- 1
     &              +(nfirst-1)*numr + 1
            nreadrs = numr
           ncount = ncount + nreadrs
            call getlstso(z(ioff), nfirstrs, nreadrs,
     &                irreprs, list, filename)
          elseif(irrepr > irreps) then
            do i = 1, numr
               nfirstrs = isymoffso(irrepr, irreprs, pops,popr,0)-1
     &           + (i-1)*nums + nfirst
               nreadrs = 1
           ncount = ncount + nreadrs
               ioff2 = ioff + (i-1)*numpq
               call getlstso(z(ioff2), nfirstrs, nreadrs,
     &                    irreprs, list, filename)
               call dscal(numpq, -1.d0, z(ioff2), 1)
            enddo
          endif
          ioff = ioff + numpq*numr
      enddo

      endif
          
10    CONTINUE
      deallocate(tmpdata)
      RETURN
      END
