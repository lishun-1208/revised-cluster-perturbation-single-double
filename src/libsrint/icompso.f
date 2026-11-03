      subroutine icompso(iin, iout, popp, popq, ioppq, icpq,  
     &   popr, pops, ioprs, icrs, sym) 
C     this program is to do the compress of I(pqrs) with iop 
C     icpq=0, compress p,q to p<q; icpq=2, compress p,q to p<=q
C     icrs=0, compress r,s to r<s; icrs=2, compress r,s to r<=s
C-----attention! this program can only deal with on side. 
C-----for example: icpq/=1, then icrs must equal to 1 
C-------------------
      use mod_iop
      use mod_symm
      implicit none 
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp(8), popq(8), popr(8), pops(8) 
      integer, intent(in) :: ioppq, icpq, ioprs, icrs 
      integer, intent(in) :: sym
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
      real*4, allocatable :: itmp(:)
      DOUBLE PRECISION, EXTERNAL :: DNRM2
C
      integer irrepr, irrepl, i, numl, numr, numl2, numr2, popl
      integer nsize1, nsize2, ioff, ioff1, nsize
      if(ioppq/=0.and.ioppq/=1.and.ioppq/=2) then
         write(6,*)'error of ioppq in isymtrso'
         stop
      endif
      if(ioprs/=0.and.ioprs/=1.and.ioprs/=2) then
         write(6,*)'error of ioprs in isymtrso'
         stop
      endif
C
      if(icpq==0.and.icrs==1) then 
      ioff = 1 
      ioff1 = 1
      do irrepr = 1, nirrep 
         irrepl = dirprd(irrepr, sym) 
         numr = irpdso(irrepr, popr, pops, ioprs)
         numl = irpdso(irrepl, popp, popq, ioppq) 
         numl2 = irpdso(irrepl, popp, popq, icpq) 
         call sqsymso(irrepl, popp, numl2, numl, numr, 
     &              iout(ioff1),iin(ioff))  
         ioff = ioff + numl*numr*isd
         ioff1 = ioff1 + numl2*numr*isd
      enddo
C 
      elseif(icpq==1.and.icrs==0) then
      ioff = 1 
      ioff1 = 1
      do irrepr = 1, nirrep 
         irrepl = dirprd(irrepr, sym) 
         numl = irpdso(irrepl, popp, popq, ioppq) 
         numr = irpdso(irrepr, popr, pops, ioprs)
         numr2 = irpdso(irrepr, popr, pops, icrs) 
         call sq2symso(irrepr, popr, numl, iout(ioff1), iin(ioff))
         ioff = ioff + numl*numr*isd 
         ioff1 = ioff1 + numl*numr2*isd
      enddo
C
      elseif(icpq==2.and.icrs==1) then
      ioff = 1
      ioff1 = 1
      do irrepr = 1, nirrep
         irrepl = dirprd(irrepr, sym)
         numr = irpdso(irrepr, popr, pops, ioprs)
         numl = irpdso(irrepl, popp, popq, ioppq)
         numl2 = irpdso(irrepl, popp, popq, icpq)
         call sq4symso(irrepl, popp, numl2, numl, numr,
     &              iout(ioff1),iin(ioff))
         ioff = ioff + numl*numr*isd
         ioff1 = ioff1 + numl2*numr*isd
      enddo
C     
      elseif(icpq==1.and.icrs==2) then
      ioff = 1
      ioff1 = 1
      do irrepr = 1, nirrep
         irrepl = dirprd(irrepr, sym)
         numl = irpdso(irrepl, popp, popq, ioppq)
         numr = irpdso(irrepr, popr, pops, ioprs)
         numr2 = irpdso(irrepr, popr, pops, icrs)   
         call sq3symso(irrepr, popr, numl, iout(ioff1), iin(ioff))
         ioff = ioff + numl*numr*isd
         ioff1 = ioff1 + numl*numr2*isd
      enddo
C
      else
        write(6,*)'attention!!!, error in icompso'
        stop
      endif
      return 
      end  
