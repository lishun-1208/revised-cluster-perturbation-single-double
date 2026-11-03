      subroutine EXP3SYMSO(sym, POPR, NUMPQ, IIN, IOUT)
C     expand r<=s to r,s  
C 
C-------------------
      use mod_iop
      use omp_lib
      implicit none
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: sym, numpq
      real*4, intent(in),managed :: iin(*)
      real*4, intent(out),managed :: iout(*)
      if(prec==8.and.igpu==0) then 
         call EXP3SYMSOd(sym, POPR, NUMPQ, IIN, IOUT)
      elseif(prec==4.and.igpu==0) then 
         call exp3symsos(sym, popr, numpq, iin, iout)
      elseif(prec==8.and.igpu==1) then 
         call EXP3SYMSOdg(sym, POPR, NUMPQ, IIN, IOUT)
      elseif(prec==4.and.igpu==1) then 
         call exp3symsosg(sym, popr, numpq, iin, iout)
      endif
      return
      end
C
      subroutine EXP3SYMSOS(sym, POPR, NUMPQ, IIN, IOUT)
C     expand r<=s to r,s  
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: sym, numpq
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs, irreppq, numrs,  numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         irreprs = sym
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr>irreps) cycle
            if(irrepr==irreps) then
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, s
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq 
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO
            do s = 2, nums
            do r = 1, s-1
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq 
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            elseif(irrepr<irreps) then
C
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            endif
          enddo
C
      return 
      end  
      subroutine EXP3SYMSOD(sym, POPR, NUMPQ, IIN, IOUT)
C     expand r<=s to r,s  
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: sym, numpq
      double precision, intent(in) :: iin(*)
      double precision, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs, irreppq, numrs,  numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         irreprs = sym
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr>irreps) cycle
            if(irrepr==irreps) then
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, s
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq 
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO
            do s = 2, nums
            do r = 1, s-1
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq 
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            elseif(irrepr<irreps) then
C
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$OMP PARALLEL DO
            do s = 1, nums
            do r = 1, numr
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            endif
          enddo
C
      return 
      end  
      subroutine EXP3SYMSOSG(sym, POPR, NUMPQ, IIN, IOUT)
C     expand r<=s to r,s  
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: sym, numpq
      real*4, intent(in),device :: iin(*)
      real*4, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs, irreppq, numrs,  numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         irreprs = sym
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr>irreps) cycle
            if(irrepr==irreps) then
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$acc kernels
!$acc loop independent
            do s = 1, nums
!$acc loop independent
            do r = 1, s
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq 
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
!$acc kernels
!$acc loop independent
            do s = 2, nums
!$acc loop independent
            do r = 1, s-1
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq 
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            elseif(irrepr<irreps) then
C
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$acc kernels
!$acc loop independent
            do s = 1, nums
!$acc loop independent
            do r = 1, numr
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            endif
          enddo
C
      return 
      end  
      subroutine EXP3SYMSODG(sym, POPR, NUMPQ, IIN, IOUT)
C     expand r<=s to r,s  
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: sym, numpq
      double precision, intent(in),device :: iin(*)
      double precision, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs, irreppq, numrs,  numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         irreprs = sym
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr>irreps) cycle
            if(irrepr==irreps) then
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$acc kernels
!$acc loop independent
            do s = 1, nums
!$acc loop independent
            do r = 1, s
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq 
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
!$acc kernels
!$acc loop independent
            do s = 2, nums
!$acc loop independent
            do r = 1, s-1
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*s/2+r-1)*numpq + pq 
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq 
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            elseif(irrepr<irreps) then
C
            numr = popr(irrepr)
            nums = popr(irreps)
            vstarto = isymoffso(irreps, irreprs, popr, popr, 2)
            vstartn = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartn2= isymoffso(irrepr, irreprs, popr, popr, 1)
C
!$acc kernels
!$acc loop independent
            do s = 1, nums
!$acc loop independent
            do r = 1, numr
!$acc loop independent
            do pq = 1, numpq
               numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
               numout = (vstartn-1 + (s-1)*numr +r-1)*numpq + pq
               numout2 = (vstartn2-1 + (r-1)*nums +s-1)*numpq + pq
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            endif
          enddo
C
      return 
      end  
