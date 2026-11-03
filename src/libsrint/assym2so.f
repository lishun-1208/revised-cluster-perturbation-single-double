      subroutine assym2so(IRREPRS, POPR, NUMPQ, iout, iin) 
C
C-----iout(pq,r,s) = iin(pq,r,s)-iin(pq,s,r)
C 
C-------------------
      use omp_lib
      use mod_iop
      use mod_symm
      implicit none 
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: IRREPRS, NUMPQ 
      real(kind=4), intent(in),managed :: iin(*)
      real(kind=4), intent(out),managed:: iout(*)
      double precision, EXTERNAL :: DNRM2
C
      if(prec==8.and.igpu==0) then 
         call assym2sod(IRREPRS, POPR, NUMPQ, iout, iin) 
      elseif(prec==4.and.igpu==0) then 
         call assym2sos(IRREPRS, POPR, NUMPQ, iout, iin) 
      elseif(prec==8.and.igpu==1) then 
         call assym2sodg(IRREPRS, POPR, NUMPQ, iout, iin) 
      elseif(prec==4.and.igpu==1) then 
         call assym2sosg(IRREPRS, POPR, NUMPQ, iout, iin) 
      endif
C
      return 
      end  

      subroutine assym2sod(IRREPRS, POPR, NUMPQ, iout, iin) 
C
C-----iout(pq,r,s) = iin(pq,r,s)-iin(pq,s,r)
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: IRREPRS, NUMPQ 
      real(kind=8), intent(in) :: iin(*)
      real(kind=8), intent(out) :: iout(*)
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
      integer vstartso, vstartqo, vstarto2, numin2
      integer vstarts, vstartr, vstartp, vstartq
      integer idxin, idxout, idxout2
      integer numout, numout2, numin, numout3, numout4
      integer, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr==irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$OMP PARALLEL DO
               do s = 2, nums
               do r = 1, s-1 
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*(s-2)/2+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$OMP END PARALLEL DO
            elseif(irrepr<irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$OMP PARALLEL DO
               do s = 1, nums
               do r = 1, numr
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*numr+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$OMP END PARALLEL DO
            endif
          enddo
C
      return 
      end  

      subroutine assym2sos(IRREPRS, POPR, NUMPQ, iout, iin) 
C
C-----iout(pq,r,s) = iin(pq,r,s)-iin(pq,s,r)
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: IRREPRS, NUMPQ 
      real(kind=4), intent(in) :: iin(*)
      real(kind=4), intent(out) :: iout(*)
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
      integer vstartso, vstartqo, vstarto2, numin2
      integer vstarts, vstartr, vstartp, vstartq
      integer idxin, idxout, idxout2
      integer numout, numout2, numin, numout3, numout4
      integer, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr==irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$OMP PARALLEL DO
               do s = 2, nums
               do r = 1, s-1 
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*(s-2)/2+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$OMP END PARALLEL DO
            elseif(irrepr<irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$OMP PARALLEL DO
               do s = 1, nums
               do r = 1, numr
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*numr+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$OMP END PARALLEL DO
            endif
          enddo
C
      return 
      end  

      subroutine assym2sodg(IRREPRS, POPR, NUMPQ, iout, iin) 
C
C-----iout(pq,r,s) = iin(pq,r,s)-iin(pq,s,r)
C 
C-------------------
      use mod_symm
      use omp_lib
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popr(8)
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
      integer vstartso, vstartqo, vstarto2, numin2
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr==irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$acc kernels
!$acc loop independent          
               do s = 2, nums
!$acc loop independent          
               do r = 1, s-1 
!$acc loop independent          
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*(s-2)/2+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$acc end kernels
            elseif(irrepr<irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$acc kernels
!$acc loop independent          
               do s = 1, nums
!$acc loop independent          
               do r = 1, numr
!$acc loop independent          
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*numr+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$acc end kernels
            endif
          enddo
C
      return 
      end  

      subroutine assym2sosg(IRREPRS, POPR, NUMPQ, iout, iin) 
C
C-----iout(pq,r,s) = iin(pq,r,s)-iin(pq,s,r)
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C
C----------------------------------

      integer, intent(in) :: popr(8)
      integer, intent(in) :: IRREPRS, NUMPQ 
      real*4, intent(in),device :: iin(*)
      real*4, intent(out),device :: iout(*)
C
      integer irrepr, irrepl, i, numl,  numl2, numr2, popl
      integer nsize1, nsize2
      integer*8 ioffo(8), ioffn(8), vstarto, vstartn
      integer*8 vstartn2
      integer irrepp, irrepq, irrepr1, irreps
      integer irreppq, numrs, numrs2, numpq2 
      integer nump, numq, numr, nums
      integer p, q, r, s, pq, rs
      integer vstartso, vstartqo, vstarto2, numin2
      integer vstarts, vstartr, vstartp, vstartq
      integer*8 idxin, idxout, idxout2
      integer*8 numout, numout2, numin, numout3, numout4
      integer*8, external :: idsymsoc
      integer, external ::  isymoffso
C
         do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if(irrepr==irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$acc kernels
!$acc loop independent          
               do s = 2, nums
!$acc loop independent          
               do r = 1, s-1 
!$acc loop independent          
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*(s-2)/2+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$acc end kernels
            elseif(irrepr<irreps) then 
              numr = popr(irrepr)
              nums = popr(irreps)
              vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
              vstarto2 = isymoffso(irrepr, irreprs, popr, popr, 1)
              vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
!$acc kernels
!$acc loop independent          
               do s = 1, nums
!$acc loop independent          
               do r = 1, numr
!$acc loop independent          
               do pq = 1, numpq
                  numin = (vstarto-1 + (s-1)*numr +r-1)*numpq + pq
                  numin2 = (vstarto2-1 + (r-1)*nums +s-1)*numpq + pq
                  numout = (vstartn-1 + (s-1)*numr+r-1)*numpq + pq
                  iout(numout)=iin(numin)-iin(numin2)
               enddo
               enddo
               enddo
!$acc end kernels
            endif
          enddo
C
      return 
      end  
