      subroutine sqsymso(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
C     this program is to do the compress of I(p,q,L) to I(p<q) 
C     numpq2, p<q
C-------------------
      use mod_iop
      use omp_lib
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
      double precision, intent(in),managed :: iin(*)
      double precision, intent(out),managed :: iout(*)
      if(prec==8.and.igpu==0) then
      call sqsymsod(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
      elseif(prec==4.and.igpu==0) then
      call sqsymsos(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
      elseif(prec==8.and.igpu==1) then
      call sqsymsodg(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
      elseif(prec==4.and.igpu==1) then
      call sqsymsosg(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
      endif
      return
      end
      subroutine sqsymsos(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
C     this program is to do the compress of I(p,q,L) to I(p<q) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ 
      integer, intent(in) :: NUMRS
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
C
      integer nsize1, nsize2
      integer ioffo(8), ioffn(8), vstarto, vstartn
      integer vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreprs
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo
      integer vstarts, vstartr, vstartp, vstartq
      integer idxin, idxout, idxout2
      integer numout, numout2, numin, numout3, numout4
      integer, external :: idsymsoc
      integer, external ::  isymoffso
C
      ioffo(1)=0; ioffn(1)=0
         do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            if(irrepp>irrepq) cycle  
            if(irrepp==irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$OMP PARALLEL DO
         
            do rs = 1, numrs 
            do q = 2, numq 
            do p = 1, q-1
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*(q-2)/2+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$OMP PARALLEL DO
            do rs = 1, numrs
            do q = 1, numq
            do p = 1, nump
                numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
                numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
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
      subroutine sqsymsod(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
C     this program is to do the compress of I(p,q,L) to I(p<q) 
C-------------------
      use mod_symm
      use omp_lib
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
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$OMP PARALLEL DO
         
            do rs = 1, numrs 
            do q = 2, numq 
            do p = 1, q-1
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*(q-2)/2+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
C
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$OMP PARALLEL DO
            do rs = 1, numrs
            do q = 1, numq
            do p = 1, nump
                numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
                numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
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
      subroutine sqsymsosg(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
C     this program is to do the compress of I(p,q,L) to I(p<q) 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
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
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs 
!$acc loop independent          
            do q = 2, numq 
!$acc loop independent          
            do p = 1, q-1
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*(q-2)/2+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs
!$acc loop independent          
            do q = 1, numq
!$acc loop independent          
            do p = 1, nump
                numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
                numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
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
      subroutine sqsymsodg(IRREPPQ, POPP, NUMPQ2, 
     &                     NUMPQ, NUMRS, IOUT, IIN) 
C     this program is to do the compress of I(p,q,L) to I(p<q) 
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
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs 
!$acc loop independent          
            do q = 2, numq 
!$acc loop independent          
            do p = 1, q-1
               numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
               numout = (rs-1)*numpq2+ vstartn-1+(q-1)*(q-2)/2+p
               iout(numout)=iin(numin)
            enddo
            enddo
            enddo
!$acc end kernels
C
            elseif(irrepp<irrepq) then 
            numq = popp(irrepq)
            nump = popp(irrepp)
            vstarto =  isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartn =  isymoffso(irrepq, irreppq, popp, popp, 0)
!$acc kernels
!$acc loop independent          
            do rs = 1, numrs
!$acc loop independent          
            do q = 1, numq
!$acc loop independent          
            do p = 1, nump
                numin = (rs-1)*numpq+ vstarto-1+(q-1)*nump+p
                numout = (rs-1)*numpq2+ vstartn-1+(q-1)*nump+p
               iout(numout)=iin(numin)
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
