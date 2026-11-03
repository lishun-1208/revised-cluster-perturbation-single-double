      SUBROUTINE GETPQRSOD_INCORE(Z, NFIRST, IRREPS, popp, popq, ioppq,
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
      character(len=2) tmpc
      character(len=10) filename2
C 
      do irrepr = 1, nirrep 
         numpq = irpdso(irrepr, popp, popq, 1) 
         maxpq = max(maxpq, numpq) 
      enddo 
C----I will first try to get PQRS term.
C----With symmetry of S be defined, the symmetry of PQR block is 
C----also determined. i will first resort PQR and then sstgen pqr
C----to prq
C      write(6,*)'begin get pqr'
!      irrepx = dirprd(irrepx, irrep0(list))
      irrepx = isympqrs
C      write(6,*)'irrepx', irrepx
      if(ioprs==1) then 
      ioff = 1
      do irrepr = 1, nirrep 
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq) 
         nums = pops(irreps) 
         numr = popr(irrepr) 
          nfirstrs = isymoffso(irreps,irreprs, popr, pops, ioprs)
     &                 -1+(nfirst-1)*numr + 1
          nreadrs = numr 
            call getlstso_incore(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
           ioff = ioff + numpq*numr
      enddo 
C
C----to r<s
      elseif(ioprs==2) then
      ioff = 1
      ncount = 1
      do irrepr = 1, nirrep
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq)
         nums = pops(irreps)
         numr = popr(irrepr)
         if(irrepr == irreps) then 
           nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)-1
     &             +(nfirst-1)*nfirst/2 + 1
           nreadrs = nfirst 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
!          call getlstso(z(ioff), nfirstrs, nreadrs,
!    &                irreprs, list, filename)
            call getlstso_incore(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
         do i = nfirst+1, nums
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)-1 
     &                + i*(i-1)/2 + nfirst 
            nreadrs = 1
            ioff2 = ioff + (i-1)*numpq 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
!           call getlstso(z(ioff2), nfirstrs, nreadrs,
!    &                    irreprs, list, filename)
            call getlstso_incore(z(ioff2), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
          enddo
          elseif(irrepr <irreps) then 
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 2)- 1
     &              +(nfirst-1)*numr + 1
            nreadrs = numr 
           ncount = ncount + nreadrs
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
!           call getlstso(z(ioff), nfirstrs, nreadrs,
!    &                irreprs, list, filename)
            call getlstso_incore(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
          elseif(irrepr > irreps) then 
            do i = 1, numr 
               nfirstrs = isymoffso(irrepr, irreprs, pops,popr,2)-1
     &           + (i-1)*nums + nfirst
               nreadrs = 1
           ncount = ncount + nreadrs
               ioff2 = ioff + (i-1)*numpq 
C           write(6,'(A, I6, I6, I6, I6, I6, I6)')
C     &      'num in getpqr', nfirstrs, nreadrs, numr, nums, irrepr, irreps
!              call getlstso(z(ioff2), nfirstrs, nreadrs,
!    &                    irreprs, list, filename)
            call getlstso_incore(z(ioff2), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
            enddo
          endif  
          ioff = ioff + numpq*numr
      enddo
C
      elseif(ioprs==0) then
C     write(6,*)'test here', irrepx
      ioff = 1
      ncount = 1
      do irrepr = 1, nirrep
         irreprs = dirprd(irrepr, irreps)
         irreppq = dirprd(irreprs, irrepx)
         numpq = irpdso(irreppq, popp, popq, ioppq)
         nums = pops(irreps)
         numr = popr(irrepr)
         if(irrepr == irreps) then
           nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)-1
     &             +(nfirst-1)*(nfirst-2)/2 + 1
           nreadrs = nfirst-1
           ncount = ncount + nreadrs
!          call getlstso(z(ioff), nfirstrs, nreadrs,
!    &                irreprs, list, filename)
            call getlstso_incore(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
         i = nfirst 
         ioff2 = ioff + (i-1)*numpq
         call zeroso(z(ioff2), numpq)
         do i = nfirst+1, nums
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)-1
     &                + (i-2)*(i-1)/2 + nfirst
            nreadrs = 1
            ioff2 = ioff + (i-1)*numpq
           ncount = ncount + nreadrs
!           call getlstso(z(ioff2), nfirstrs, nreadrs,
!    &                    irreprs, list, filename)
            call getlstso_incore(z(ioff2), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
            call scalso(numpq, -1.d0, z(ioff2), 1)
          enddo
          elseif(irrepr <irreps) then
            nfirstrs = isymoffso(irreps, irreprs, popr, pops, 0)- 1
     &              +(nfirst-1)*numr + 1
            nreadrs = numr
           ncount = ncount + nreadrs
!           call getlstso(z(ioff), nfirstrs, nreadrs,
!    &                irreprs, list, filename)
            call getlstso_incore(z(ioff), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
          elseif(irrepr > irreps) then
            do i = 1, numr
               nfirstrs = isymoffso(irrepr, irreprs, pops,popr,0)-1
     &           + (i-1)*nums + nfirst
               nreadrs = 1
           ncount = ncount + nreadrs
               ioff2 = ioff + (i-1)*numpq
!              call getlstso(z(ioff2), nfirstrs, nreadrs,
!    &                    irreprs, list, filename)
            call getlstso_incore(z(ioff2), nfirstrs, nreadrs, 
     &                irreprs, popp, popq, ioppq, popr, pops, ioprs, 
     &                ipqrs, isympqrs) 
               call scalso(numpq, -1.d0, z(ioff2), 1)
            enddo
          endif
          ioff = ioff + numpq*numr
      enddo

      endif
          
10    CONTINUE
      RETURN
      END
