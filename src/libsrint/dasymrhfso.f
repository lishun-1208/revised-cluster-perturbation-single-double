      subroutine dasymrhfso(irreppq, irreprs, popp, popr, numpq, iout, iin)
C
C-----iout(p,q,r,s) = iin(p,q,r,s)-iin(q,p,s,r)
C 
C-------------------
      use mod_iop
      use omp_lib
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), popp(8)
      integer, intent(in) :: IRREPRS, IRREPPQ, NUMPQ
      double precision, intent(in),managed :: iin(*)
      double precision, intent(out),managed :: iout(*)
      double precision, EXTERNAL :: DNRM2
      if(prec==8.and.igpu==0) then 
        call dasymrhfsod(irreppq, irreprs, popp, popr, numpq, iout, iin)
      elseif(prec==4.and.igpu==0) then 
        call dasymrhfsos(irreppq, irreprs, popp, popr, numpq, iout, iin)
      elseif(prec==8.and.igpu==1) then 
        call dasymrhfsodg(irreppq, irreprs, popp, popr, numpq, iout, iin)
      elseif(prec==4.and.igpu==1) then 
        call dasymrhfsosg(irreppq, irreprs, popp, popr, numpq, iout, iin)
      endif
C
      return 
      end  
      subroutine dasymrhfsod(irreppq, irreprs, popp, popr, numpq, iout, iin)
C
C-----iout(p,q,r,s) = iin(p,q,r,s)-iin(q,p,s,r)
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), popp(8)
      integer, intent(in) :: IRREPRS, IRREPPQ, NUMPQ
      double precision, intent(in) :: iin(*)
      double precision, intent(out) :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartrs1, vstartrs2, vstartpq1, vstartpq2
      integer*8 idxin, idxout, idxout2
      integer*8 numin2, numin
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = popr(irreps)
            vstartrs1 = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartrs2 = isymoffso(irrepr, irreprs, popr, popr, 1)
         do irrepq = 1, nirrep 
            irrepp = dirprd(irreppq, irrepq) 
            nump = popp(irrepp)
            numq = popp(irrepq)
            vstartpq1 = isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartpq2 = isymoffso(irrepp, irreppq, popp, popp, 1)
            if(numr==0.or.nums==0.or.numq==0.or.numq==0) cycle
!$OMP PARALLEL DO
            do s = 1, nums 
            do r = 1, numr 
            do q = 1, numq
            do p = 1, nump 
                  numin = (vstartrs1-1 + (s-1)*numr +r-1)*numpq 
     &                   + vstartpq1 - 1 +(q-1)*nump + p 
                  numin2 = (vstartrs2-1 + (r-1)*nums +s-1)*numpq 
     &                   + vstartpq2 - 1 +(p-1)*numq + q 
                  iout(numin)=iin(numin) - iin(numin2)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
          enddo
          enddo
C
      return 
      end  
      subroutine dasymrhfsos(irreppq, irreprs, popp, popr, numpq, iout, iin)
C
C-----iout(p,q,r,s) = iin(p,q,r,s)-iin(q,p,s,r)
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), popp(8)
      integer, intent(in) :: IRREPRS, IRREPPQ, NUMPQ
      real(kind=4), intent(in) :: iin(*)
      real(kind=4), intent(out) :: iout(*)
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartrs1, vstartrs2, vstartpq1, vstartpq2
      integer*8 idxin, idxout, idxout2
      integer*8 numin2, numin
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = popr(irreps)
            vstartrs1 = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartrs2 = isymoffso(irrepr, irreprs, popr, popr, 1)
         do irrepq = 1, nirrep 
            irrepp = dirprd(irreppq, irrepq) 
            nump = popp(irrepp)
            numq = popp(irrepq)
            vstartpq1 = isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartpq2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$OMP PARALLEL DO
            do s = 1, nums 
            do r = 1, numr 
            do q = 1, numq
            do p = 1, nump 
                  numin = (vstartrs1-1 + (s-1)*numr +r-1)*numpq 
     &                   + vstartpq1 - 1 +(q-1)*nump + p 
                  numin2 = (vstartrs2-1 + (r-1)*nums +s-1)*numpq 
     &                   + vstartpq2 - 1 +(p-1)*numq + q 
                  iout(numin)=iin(numin) - iin(numin2)
            enddo
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO
          enddo
          enddo
C
      return 
      end  
      subroutine dasymrhfsodg(irreppq, irreprs, popp, popr, numpq, iout, iin)
C
C-----iout(p,q,r,s) = iin(p,q,r,s)-iin(q,p,s,r)
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), popp(8)
      integer, intent(in) :: IRREPRS, IRREPPQ, NUMPQ
      double precision, intent(in),device :: iin(*)
      double precision, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartrs1, vstartrs2, vstartpq1, vstartpq2
      integer*8 idxin, idxout, idxout2
      integer*8 numin2, numin
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = popr(irreps)
            vstartrs1 = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartrs2 = isymoffso(irrepr, irreprs, popr, popr, 1)
         do irrepq = 1, nirrep 
            irrepp = dirprd(irreppq, irrepq) 
            nump = popp(irrepp)
            numq = popp(irrepq)
            vstartpq1 = isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartpq2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent          
            do s = 1, nums 
!$acc loop independent          
            do r = 1, numr 
!$acc loop independent          
            do q = 1, numq
!$acc loop independent          
            do p = 1, nump 
                  numin = (vstartrs1-1 + (s-1)*numr +r-1)*numpq 
     &                   + vstartpq1 - 1 +(q-1)*nump + p 
                  numin2 = (vstartrs2-1 + (r-1)*nums +s-1)*numpq 
     &                   + vstartpq2 - 1 +(p-1)*numq + q 
                  iout(numin)=iin(numin) - iin(numin2)
            enddo
            enddo
            enddo
            enddo
!$acc end kernels
          enddo
          enddo
C
      return 
      end  
      subroutine dasymrhfsosg(irreppq, irreprs, popp, popr, numpq, iout, iin)
C
C-----iout(p,q,r,s) = iin(p,q,r,s)-iin(q,p,s,r)
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8), popp(8)
      integer, intent(in) :: IRREPRS, IRREPPQ, NUMPQ
      REAL*4, intent(in),device :: iin(*)
      REAL*4, intent(out),device :: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartrs1, vstartrs2, vstartpq1, vstartpq2
      integer*8 idxin, idxout, idxout2
      integer*8 numin2, numin
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            numr = popr(irrepr)
            nums = popr(irreps)
            vstartrs1 = isymoffso(irreps, irreprs, popr, popr, 1)
            vstartrs2 = isymoffso(irrepr, irreprs, popr, popr, 1)
         do irrepq = 1, nirrep 
            irrepp = dirprd(irreppq, irrepq) 
            nump = popp(irrepp)
            numq = popp(irrepq)
            vstartpq1 = isymoffso(irrepq, irreppq, popp, popp, 1)
            vstartpq2 = isymoffso(irrepp, irreppq, popp, popp, 1)
!$acc kernels
!$acc loop independent          
            do s = 1, nums 
!$acc loop independent          
            do r = 1, numr 
!$acc loop independent          
            do q = 1, numq
!$acc loop independent          
            do p = 1, nump 
                  numin = (vstartrs1-1 + (s-1)*numr +r-1)*numpq 
     &                   + vstartpq1 - 1 +(q-1)*nump + p 
                  numin2 = (vstartrs2-1 + (r-1)*nums +s-1)*numpq 
     &                   + vstartpq2 - 1 +(p-1)*numq + q 
                  iout(numin)=iin(numin) - iin(numin2)
            enddo
            enddo
            enddo
            enddo
!$acc end kernels
          enddo
          enddo
C
      return 
      end  
