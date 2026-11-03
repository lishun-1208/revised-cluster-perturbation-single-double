      subroutine itranspso(iin, iout, popp,popq, ioppq, 
     &                     popr, pops, ioprs,  sym) 
C     this program is to do the transp of I(pqrs) 
C     transp Iin(pq,rs) to Iout(rs,pq)
C----ioppq=0, p<q, ioppq=1, p,q, =2, p<=q---
C----ioprs=0, r<s, ioprs=1, r,s, =2, r<=s---
C----sym: the symmetry of Ipqrs
C
      use mod_iop
      use mod_symm
      implicit none 
C-------------------
      INTEGER, external :: IRPDSO
C----------------------------------
      integer, intent(in) :: popp(8), popq(8), ioppq 
      integer, intent(in) :: popr(8), pops(8), ioprs 
      integer, intent(in) :: sym
      real*4, intent(in),managed :: iin(*)
      real*4, intent(out),managed :: iout(*)
C
      integer irrepr, irrepl, i, numl, numr
      integer nsize1, nsize2, ioff, ioff1(8)
      double precision, external :: dnrm2
C
      if(ioppq/=0.and.ioppq/=1.and.ioppq/=2) then
         write(6,*)'error of ioppq in isymtrso'
         stop
      endif
      if(ioprs/=0.and.ioprs/=1.and.ioprs/=2) then
         write(6,*)'error of ioprs in isymtrso'
         stop
      endif

      ioff1(1)=1
      do irrepr = 1, nirrep-1 
         irrepl = dirprd(irrepr, sym)
         numl = irpdso(irrepl, popr, pops, ioprs) 
         numr = irpdso(irrepr, popp, popq, ioppq)
         ioff1(irrepr+1)=ioff1(irrepr)+numl*numr*isd
      enddo
         
C
      ioff = 1
      do irrepr = 1, nirrep 
         irrepl = dirprd(irrepr, sym)
         numl = irpdso(irrepl, popp, popq, ioppq) 
         numr = irpdso(irrepr, popr, pops, ioprs)
         nsize1 = numl*numr
         call transpso(iin(ioff), iout(ioff1(irrepl)), numr, numl)
         ioff = ioff + nsize1*isd 
      enddo
      return 
      end  
