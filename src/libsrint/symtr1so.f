      subroutine symtr1so(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT) 
C     this program is to transfer I(pqrs) to I(pqsr) 
C 
C-------------------
      use mod_iop
      use omp_lib
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8) 
      integer, intent(in) :: IRREPRS, NUMPQ 
      double precision, intent(in),managed :: iin(*)
      double precision, intent(out),managed :: iout(*)
C
      if(prec==8.and.igpu==0) then 
          call symtr1sod(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT)
      elseif(prec==4.and.igpu==0) then 
          call symtr1sos(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT)
      elseif(prec==8.and.igpu==1) then 
          call symtr1sodg(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT)
      elseif(prec==4.and.igpu==1) then 
          call symtr1sosg(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT)
      endif
      return
      end
      subroutine symtr1sos(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT) 
C     this program is to transfer I(pqrs) to I(pqsr) 
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), pops(8) 
      integer, intent(in) :: IRREPRS, NUMPQ 
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer ioffo(8), ioffn(8), vstarto, vstartn
      integer vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreppq, numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, pops, 1)
            vstartn = isymoffso(irrepr, irreprs, pops, popr, 1)
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
          enddo
C
      return 
      end  
      subroutine symtr1sod(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT) 
C     this program is to transfer I(pqrs) to I(pqsr) 
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), pops(8) 
      integer, intent(in) :: IRREPRS, NUMPQ 
      double precision, intent(in) :: iin(*)
      double precision, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreppq, numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, pops, 1)
            vstartn = isymoffso(irrepr, irreprs, pops, popr, 1)
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
          enddo
C
      return 
      end  
      subroutine symtr1sosg(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT) 
C     this program is to transfer I(pqrs) to I(pqsr) 
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), pops(8) 
      integer, intent(in) :: IRREPRS, NUMPQ 
      real*4, intent(in),device :: iin(*)
      real*4, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreppq, numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, pops, 1)
            vstartn = isymoffso(irrepr, irreprs, pops, popr, 1)
!$acc kernels
!$acc loop independent          
            do s = 1, nums
!$acc loop independent          
            do r = 1, numr
!$acc loop independent          
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
          enddo
C
      return 
      end  
      subroutine symtr1sodg(IRREPRS, POPR, POPS, NUMPQ, IIN, IOUT) 
C     this program is to transfer I(pqrs) to I(pqsr) 
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8) 
      integer, intent(in) :: IRREPRS, NUMPQ 
      double precision, intent(in),device :: iin(*)
      double precision, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreppq, numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = pops(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, pops, 1)
            vstartn = isymoffso(irrepr, irreprs, pops, popr, 1)
!$acc kernels
!$acc loop independent          
            do s = 1, nums
!$acc loop independent          
            do r = 1, numr
!$acc loop independent          
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
          enddo
C
      return 
      end  
