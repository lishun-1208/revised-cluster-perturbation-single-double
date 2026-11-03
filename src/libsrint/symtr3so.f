      subroutine symtr3so(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS, 
     &                     IIN, IOUT) 
C     this program is to do the expand of I(p,q,L) to I(q,p,L) 
C-------------------
      use mod_iop
      use omp_lib
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), IRREPPQ
      integer, intent(in) :: NUMRS, NUMPQ
      double precision, intent(in),managed :: iin(*)
      double precision, intent(out),managed :: iout(*)
C
      if(prec==8.and.igpu==0) then
         call symtr3sod(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS,
     &                     IIN, IOUT)
      elseif(prec==4.and.igpu==0) then
         call symtr3sos(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS,
     &                     IIN, IOUT)
      elseif(prec==8.and.igpu==1) then
         call symtr3sodg(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS,
     &                     IIN, IOUT)
      elseif(prec==4.and.igpu==1) then
         call symtr3sosg(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS,
     &                     IIN, IOUT)
      endif
      return
      end
      subroutine symtr3sos(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS, 
     &                     IIN, IOUT) 
C     this program is to do the expand of I(p,q,L) to I(q,p,L) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), IRREPPQ
      integer, intent(in) :: NUMRS, NUMPQ
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
      ioffo(1)=0; ioffn(1)=0
         do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            numq = popq(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popq, 1)
            vstartn =  isymoffso(irrepp, irreppq, popq, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, nump 
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq+ vstartn-1+(p-1)*numq+q
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
          enddo
C
      return 
      end  
      subroutine symtr3sod(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS, 
     &                     IIN, IOUT) 
C     this program is to do the expand of I(p,q,L) to I(q,p,L) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), IRREPPQ
      integer, intent(in) :: NUMRS, NUMPQ
      double precision, intent(in) :: iin(*)
      double precision, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
      ioffo(1)=0; ioffn(1)=0
         do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            numq = popq(irrepq)
            nump = popp(irrepp)
            if(numq==0.or.numq==0.or.numrs==0) cycle
            vstarto =  isymoffso(irrepq, irreppq, popp, popq, 1)
            vstartn =  isymoffso(irrepp, irreppq, popq, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, nump 
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq+ vstartn-1+(p-1)*numq+q
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
          enddo
C
      return 
      end  
      subroutine symtr3sosg(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS, 
     &                     IIN, IOUT) 
C     this program is to do the expand of I(p,q,L) to I(q,p,L) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), IRREPPQ
      integer, intent(in) :: NUMRS, NUMPQ
      real*4, intent(in),device :: iin(*)
      real*4, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
      ioffo(1)=0; ioffn(1)=0
         do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            numq = popq(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popq, 1)
            vstartn =  isymoffso(irrepp, irreppq, popq, popp, 1)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs 
!$acc loop independent          
            do q = 1, numq 
!$acc loop independent          
            do p = 1, nump 
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq+ vstartn-1+(p-1)*numq+q
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
          enddo
C
      return 
      end  
      subroutine symtr3sodg(IRREPPQ, POPP, POPQ, NUMPQ, NUMRS, 
     &                     IIN, IOUT) 
C     this program is to do the expand of I(p,q,L) to I(q,p,L) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), IRREPPQ
      integer, intent(in) :: NUMRS, NUMPQ
      double precision, intent(in),device :: iin(*)
      double precision, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
      ioffo(1)=0; ioffn(1)=0
         do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            numq = popq(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popq, 1)
            vstartn =  isymoffso(irrepp, irreppq, popq, popp, 1)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs 
!$acc loop independent          
            do q = 1, numq 
!$acc loop independent          
            do p = 1, nump 
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq+ vstartn-1+(p-1)*numq+q
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
          enddo
C
      return 
      end  
