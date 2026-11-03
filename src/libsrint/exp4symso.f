      subroutine exp4symso(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
C     this program is to do the expand of I(p<=q,L) to I(p,q) 
C-------------------
      use omp_lib
      use mod_iop
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
      real*4, intent(in),managed :: iin(*)
      real*4, intent(out),managed :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      if(prec==8.and.igpu==0) then
         call exp4symsod(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
      elseif(prec==4.and.igpu==0) then
         call exp4symsos(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
      elseif(prec==8.and.igpu==1) then
         call exp4symsodg(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
      elseif(prec==4.and.igpu==1) then
         call exp4symsosg(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
      endif
C
      return
      end
      subroutine exp4symsos(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
C     this program is to do the expand of I(p<=q,L) to I(p,q) 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
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
            if(irrepp>irrepq) cycle  
            if(irrepp==irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C        
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs
            do q = 1, numq
            do p = 1, nump
                numin = idxin+rs*numpq+q*nump+p
                numout= idxout+rs*numpq2+q*nump+p
                numout2 = idxout2+rs*numpq2+p*numq+q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
C               write(6,*)'numin', numin
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
      subroutine exp4symsod(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
C     this program is to do the expand of I(p<=q,L) to I(p,q) 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
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
            if(irrepp>irrepq) cycle  
            if(irrepp==irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO
            do rs = 1, numrs 
            do q = 1, numq 
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C        
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$OMP PARALLEL DO
            do rs = 1, numrs
            do q = 1, numq
            do p = 1, nump
                numin = idxin+rs*numpq+q*nump+p
                numout= idxout+rs*numpq2+q*nump+p
                numout2 = idxout2+rs*numpq2+p*numq+q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
C               write(6,*)'numin', numin
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
      subroutine exp4symsosg(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
C     this program is to do the expand of I(p<=q,L) to I(p,q) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
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
            if(irrepp>irrepq) cycle  
            if(irrepp==irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent
            do rs = 1, numrs 
!$acc loop independent
            do q = 1, numq 
!$acc loop independent
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
!$acc kernels
!$acc loop independent
            do rs = 1, numrs 
!$acc loop independent
            do q = 1, numq 
!$acc loop independent
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C        
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent
            do rs = 1, numrs
!$acc loop independent
            do q = 1, numq
!$acc loop independent
            do p = 1, nump
                numin = idxin+rs*numpq+q*nump+p
                numout= idxout+rs*numpq2+q*nump+p
                numout2 = idxout2+rs*numpq2+p*numq+q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
C               write(6,*)'numin', numin
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
      subroutine exp4symsodg(IRREPPQ, POPP, NUMPQ, 
     &                     NUMPQ2, NUMRS, IIN, IOUT) 
C     this program is to do the expand of I(p<=q,L) to I(p,q) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
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
            if(irrepp>irrepq) cycle  
            if(irrepp==irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent
            do rs = 1, numrs 
!$acc loop independent
            do q = 1, numq 
!$acc loop independent
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
!$acc kernels
!$acc loop independent
            do rs = 1, numrs 
!$acc loop independent
            do q = 1, numq 
!$acc loop independent
            do p = 1, q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*q/2+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout2)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C        
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 2)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent
            do rs = 1, numrs
!$acc loop independent
            do q = 1, numq
!$acc loop independent
            do p = 1, nump
                numin = idxin+rs*numpq+q*nump+p
                numout= idxout+rs*numpq2+q*nump+p
                numout2 = idxout2+rs*numpq2+p*numq+q
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               numout2 = (rs-1)*numpq2 + vstartn2-1+(p-1)*numq + q
               iout(numout)=iin(numin)
               iout(numout2)=iin(numin)
C               write(6,*)'numin', numin
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
