      subroutine sumpq(symw, popp, popr, pops, ioprs,iin, iout) 
C
C-----iout(rs) = sum_p[iin(p,p,r,s)]
C 
C-------------------
      use omp_lib
      use mod_iop
      implicit none 
C----------------------------------

      integer, intent(in) :: popr(8), pops(8), popp(8)
      integer, intent(in) :: symw, ioprs
      real(kind=4), intent(in),managed :: iin(*)
      real(kind=4), intent(out),managed:: iout(*)
C
      integer,external :: irpdso,isymoffso 
      integer numrs, numpq, ioffw(8), ioff 
      integer irreprs, irreppq, irrepp, irrepq, nump, numq, ioffi, i, p
C
      if(prec==8.and.igpu==0) then 
         call sumpqd(symw, popp, popr, pops, ioprs,iin, iout) 
      elseif(prec==4.and.igpu==0) then 
         call sumpqs(symw, popp, popr, pops, ioprs,iin, iout) 
      elseif(prec==8.and.igpu==1) then 
         call sumpqdg(symw, popp, popr, pops, ioprs,iin, iout) 
      elseif(prec==4.and.igpu==1) then 
         call sumpqsg(symw, popp, popr, pops, ioprs,iin, iout) 
      endif 
C
      return 
      end  

      subroutine sumpqs(symw, popp, popr, pops, ioprs,iin, iout) 
C
C-----iout(rs) = sum_p[iin(p,p,r,s)]
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8), popp(8)
      integer, intent(in) :: symw, ioprs
      real(kind=4), intent(in),managed :: iin(*)
      real(kind=4), intent(out),managed:: iout(*)
C
      integer,external :: irpdso,isymoffso 
      integer numrs, numpq, ioffw(8), ioff 
      integer irreprs, irreppq, irrepp, irrepq, nump, numq, ioffi, i, p
C
      ioffw(1) = 1
      do irreprs = 1, nirrep - 1
         irreppq = dirprd(irreprs, symw) 
         numpq = irpdso(irreppq, popp, popp, 1)
         numrs = irpdso(irreprs, popr, pops, ioprs) 
         ioffw(irreprs+1) = ioffw(irreprs) + numpq*numrs 
      enddo
C 
      irreprs = symw
      numrs = irpdso(irreprs, popr, pops, ioprs)
      numpq = irpdso(1, popp, popp, 1)
      call zeroso(iout, numrs) 
      do irrepp = 1, nirrep 
         irrepq = irrepp 
         nump = popp(irrepp)
         numq = popp(irrepq) 
         ioffi = isymoffso(irrepq, 1, popp, popp, 1)+ioffw(irreprs)-1 
!$omp parallel do 
         do i = 1, numrs 
         do p = 1, nump 
            ioff = ioffi-1 + (p-1)*nump+p +(i-1)*numpq
            iout(i) = iout(i) + iin(ioff) 
         enddo
         enddo
!$omp end parallel do 
      enddo 
C
      return 
      end  

      subroutine sumpqd(symw, popp, popr, pops, ioprs,iin, iout) 
C
C-----iout(rs) = sum_p[iin(p,p,r,s)]
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8), popp(8)
      integer, intent(in) :: symw, ioprs
      real(kind=8), intent(in),managed :: iin(*)
      real(kind=8), intent(out),managed:: iout(*)
C
      integer,external :: irpdso,isymoffso 
      integer numrs, numpq, ioffw(8), ioff 
      integer irreprs, irreppq, irrepp, irrepq, nump, numq, ioffi, i, p
C
      ioffw(1) = 1
      do irreprs = 1, nirrep - 1
         irreppq = dirprd(irreprs, symw) 
         numpq = irpdso(irreppq, popp, popp, 1)
         numrs = irpdso(irreprs, popr, pops, ioprs) 
         ioffw(irreprs+1) = ioffw(irreprs) + numpq*numrs 
      enddo
C 
      irreprs = symw
      numrs = irpdso(irreprs, popr, pops, ioprs)
      numpq = irpdso(1, popp, popp, 1)
      call zeroso(iout, numrs) 
      do irrepp = 1, nirrep 
         irrepq = irrepp 
         nump = popp(irrepp)
         numq = popp(irrepq) 
         ioffi = isymoffso(irrepq, 1, popp, popp, 1)+ioffw(irreprs)-1 
         do i = 1, numrs 
         do p = 1, nump 
            ioff = ioffi-1 + (p-1)*nump+p +(i-1)*numpq
            iout(i) = iout(i) + iin(ioff) 
         enddo
         enddo
      enddo 
C
      return 
      end  

      subroutine sumpqsg(symw, popp, popr, pops, ioprs,iin, iout) 
C
C-----iout(rs) = sum_p[iin(p,p,r,s)]
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8), popp(8)
      integer, intent(in) :: symw, ioprs
      real(kind=4), intent(in),device :: iin(*)
      real(kind=4), intent(out),device:: iout(*)
C
      integer,external :: irpdso,isymoffso 
      integer numrs, numpq, ioffw(8), ioff 
      integer irreprs, irreppq, irrepp, irrepq, nump, numq, ioffi, i, p
C
      ioffw(1) = 1
      do irreprs = 1, nirrep - 1
         irreppq = dirprd(irreprs, symw) 
         numpq = irpdso(irreppq, popp, popp, 1)
         numrs = irpdso(irreprs, popr, pops, ioprs) 
         ioffw(irreprs+1) = ioffw(irreprs) + numpq*numrs 
      enddo
C 
      irreprs = symw
      numrs = irpdso(irreprs, popr, pops, ioprs)
      numpq = irpdso(1, popp, popp, 1)
      call zeroso(iout, numrs) 
      do irrepp = 1, nirrep 
         irrepq = irrepp 
         nump = popp(irrepp)
         numq = popp(irrepq) 
         ioffi = isymoffso(irrepq, 1, popp, popp, 1)+ioffw(irreprs)-1 
!$acc kernels 
         do i = 1, numrs 
         do p = 1, nump 
            ioff = ioffi-1 + (p-1)*nump+p +(i-1)*numpq
            iout(i) = iout(i) + iin(ioff) 
         enddo
         enddo
!$acc end kernels 
      enddo 
C
      return 
      end  

      subroutine sumpqdg(symw, popp, popr, pops, ioprs,iin, iout) 
C
C-----iout(rs) = sum_p[iin(p,p,r,s)]
C 
C-------------------
      use omp_lib
      use mod_symm
      implicit none 
C
C----------------------------------

      integer, intent(in) :: popr(8), pops(8), popp(8)
      integer, intent(in) :: symw, ioprs
      real(kind=8), intent(in),device :: iin(*)
      real(kind=8), intent(out),device:: iout(*)
C
      integer,external :: irpdso,isymoffso 
      integer numrs, numpq, ioffw(8), ioff 
      integer irreprs, irreppq, irrepp, irrepq, nump, numq, ioffi, i, p
C
      ioffw(1) = 1
      do irreprs = 1, nirrep - 1
         irreppq = dirprd(irreprs, symw) 
         numpq = irpdso(irreppq, popp, popp, 1)
         numrs = irpdso(irreprs, popr, pops, ioprs) 
         ioffw(irreprs+1) = ioffw(irreprs) + numpq*numrs 
      enddo
C 
      irreprs = symw
      numrs = irpdso(irreprs, popr, pops, ioprs)
      numpq = irpdso(1, popp, popp, 1)
      call zeroso(iout, numrs) 
      do irrepp = 1, nirrep 
         irrepq = irrepp 
         nump = popp(irrepp)
         numq = popp(irrepq) 
         ioffi = isymoffso(irrepq, 1, popp, popp, 1)+ioffw(irreprs)-1 
!$acc kernels 
         do i = 1, numrs 
         do p = 1, nump 
            ioff = ioffi-1 + (p-1)*nump+p +(i-1)*numpq
            iout(i) = iout(i) + iin(ioff) 
         enddo
         enddo
!$acc end kernels 
      enddo 
C
      return 
      end  

