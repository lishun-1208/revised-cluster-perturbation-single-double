      subroutine isymtrso(iin,iout, popp,popq,ioppq,icpq,popr, pops, ioprs,
     &             icrs, sym) 
C     this program is to do the symtr of I(pqrs) 
C     if ioppq = 0, p<q; if ioprs = 0, r<s
C     if icpq = 1, change p,q to q,p
C     if icrs = 1, change r,s to s,r
C     sym, the symmetry of Ipqrs
C
      use mod_iop
      implicit none 
C----------------------------------

      integer, intent(in) :: popp, popq, popr, pops 
      integer, intent(in) :: ioppq, ioprs, icpq, icrs
      integer, intent(in) :: sym
      real*4, intent(in),managed :: iin(*)
      real*4, intent(out),managed :: iout(*)
      if(ioppq/=0.and.ioppq/=1.and.ioppq/=2) then 
         write(6,*)'error of ioppq in isymtrso'
         stop
      endif
      if(ioprs/=0.and.ioprs/=1.and.ioprs/=2) then 
         write(6,*)'error of ioprs in isymtrso'
         stop
      endif
C
      if(prec==8.and.igpu==0) then 
         call isymtrsod(iin,iout, popp,popq,ioppq,icpq,popr, pops,ioprs,
     &             icrs, sym) 
      elseif(prec==4.and.igpu==0) then 
         call isymtrsos(iin,iout, popp,popq,ioppq,icpq,popr, pops,ioprs,
     &             icrs, sym)
      elseif(prec==8.and.igpu==1) then 
         call isymtrsodg(iin,iout, popp,popq,ioppq,icpq,popr, pops,ioprs,
     &             icrs, sym) 
      elseif(prec==4.and.igpu==1) then 
         call isymtrsosg(iin,iout, popp,popq,ioppq,icpq,popr, pops,ioprs,
     &             icrs, sym)
      endif
      return 
      end  
      subroutine isymtrsos(iin,iout, popp,popq,ioppq,icpq,popr, pops, ioprs,
     &             icrs, sym) 
C     this program is to do the symtr of I(pqrs) 
C     if ioppq = 0, p<q; if ioprs = 0, r<s
C     if icpq = 1, change p,q to q,p
C     if icrs = 1, change r,s to s,r
C     sym, the symmetry of Ipqrs
C
      use mod_symm
      implicit none 
C-------------------
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp, popq, popr, pops 
      integer, intent(in) :: ioppq, ioprs, icpq, icrs
      integer, intent(in) :: sym
      real*4, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
      real*4, allocatable :: itmp(:)
C
      integer irrepr, irrepl, i
      integer nsize1, nsize2, ioff, ioff1
      integer nummax, istat, numl, numr
      nummax = 0
      if(icpq/=1.and.icrs==1) then 
       ioff = 1
       do irrepr = 1, nirrep 
          irrepl = dirprd(irrepr, sym) 
          numl = irpdso(irrepl, popp, popq, ioppq) 
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
C          write(6,*)'symtr1so 1'
          call symtr1sos(irrepr, popr, pops, numl, iin(ioff),
     &            iout(ioff)) 
C          write(6,*)'symtr1so 2'
          ioff = ioff + nsize1 
       enddo
C
      elseif(icpq==1.and.icrs/=1) then
C
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr3sos(irrepl, popp, popq, numl, numr, iin(ioff),
     &            iout(ioff))
          ioff = ioff + nsize1
       enddo
C
      elseif(icpq==1.and.icrs==1) then 
C
       nsize1= 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = max(nsize1, numl*numr) 
       enddo 
       allocate(itmp(nsize1))
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr1sos(irrepr, popr, pops, numl, iin(ioff),
     &            itmp) 
          call symtr3sos(irrepl, popp, popq, numl, numr, itmp,iout(ioff))
          ioff = ioff + nsize1
       enddo
       deallocate(itmp)
      endif
C
      return 
      end  
      subroutine isymtrsod(iin,iout, popp,popq,ioppq,icpq,popr, pops, ioprs,
     &             icrs, sym) 
C     this program is to do the symtr of I(pqrs) 
C     if ioppq = 0, p<q; if ioprs = 0, r<s
C     if icpq = 1, change p,q to q,p
C     if icrs = 1, change r,s to s,r
C     sym, the symmetry of Ipqrs
C
      use mod_symm
      implicit none 
C-------------------
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp, popq, popr, pops 
      integer, intent(in) :: ioppq, ioprs, icpq, icrs
      integer, intent(in) :: sym
      double precision, intent(in) :: iin(*)
      double precision, intent(out) :: iout(*)
      double precision, allocatable :: itmp(:)
C
      integer irrepr, irrepl, i
      integer nsize1, nsize2, ioff, ioff1
      integer nummax, istat, numl, numr
      nummax = 0
      if(icpq/=1.and.icrs==1) then 
       ioff = 1
       do irrepr = 1, nirrep 
          irrepl = dirprd(irrepr, sym) 
          numl = irpdso(irrepl, popp, popq, ioppq) 
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
C          write(6,*)'symtr1so 1'
          call symtr1sod(irrepr, popr, pops, numl, iin(ioff),
     &            iout(ioff)) 
C          write(6,*)'symtr1so 2'
          ioff = ioff + nsize1 
       enddo
C
      elseif(icpq==1.and.icrs/=1) then
C
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr3sod(irrepl, popp, popq, numl, numr, iin(ioff),
     &            iout(ioff))
          ioff = ioff + nsize1
       enddo
C
      elseif(icpq==1.and.icrs==1) then 
C
       nsize1= 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = max(nsize1, numl*numr) 
       enddo 
       allocate(itmp(nsize1))
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr1sod(irrepr, popr, pops, numl, iin(ioff),
     &            itmp) 
          call symtr3sod(irrepl, popp, popq, numl, numr, itmp,iout(ioff))
          ioff = ioff + nsize1
       enddo
       deallocate(itmp)
      endif
C
      return 
      end  
      subroutine isymtrsosg(iin,iout, popp,popq,ioppq,icpq,popr, pops, ioprs,
     &             icrs, sym) 
C     this program is to do the symtr of I(pqrs) 
C     if ioppq = 0, p<q; if ioprs = 0, r<s
C     if icpq = 1, change p,q to q,p
C     if icrs = 1, change r,s to s,r
C     sym, the symmetry of Ipqrs
C
      use mod_symm
      implicit none 
C-------------------
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp, popq, popr, pops 
      integer, intent(in) :: ioppq, ioprs, icpq, icrs
      integer, intent(in) :: sym
      real*4, intent(in),device :: iin(*)
      real*4, intent(out),device :: iout(*)
      real*4, allocatable,device :: itmp(:)
C
      integer irrepr, irrepl, i
      integer nsize1, nsize2, ioff, ioff1
      integer nummax, istat, numl, numr
      nummax = 0
      if(icpq/=1.and.icrs==1) then 
       ioff = 1
       do irrepr = 1, nirrep 
          irrepl = dirprd(irrepr, sym) 
          numl = irpdso(irrepl, popp, popq, ioppq) 
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr1sosg(irrepr, popr, pops, numl, iin(ioff),
     &            iout(ioff)) 
          ioff = ioff + nsize1 
       enddo
C
      elseif(icpq==1.and.icrs/=1) then
C
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr3sosg(irrepl, popp, popq, numl, numr, iin(ioff),
     &            iout(ioff))
          ioff = ioff + nsize1
       enddo
C
      elseif(icpq==1.and.icrs==1) then 
C
       nsize1= 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = max(nsize1, numl*numr) 
       enddo 
       allocate(itmp(nsize1))
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr1sosg(irrepr, popr, pops, numl, iin(ioff),
     &            itmp) 
          call symtr3sosg(irrepl, popp, popq, numl, numr, itmp,iout(ioff))
          ioff = ioff + nsize1
       enddo
       deallocate(itmp)
      endif
C
      return 
      end  
      subroutine isymtrsodg(iin,iout, popp,popq,ioppq,icpq,popr, pops, ioprs,
     &             icrs, sym) 
C     this program is to do the symtr of I(pqrs) 
C     if ioppq = 0, p<q; if ioprs = 0, r<s
C     if icpq = 1, change p,q to q,p
C     if icrs = 1, change r,s to s,r
C     sym, the symmetry of Ipqrs
C
      use mod_symm
      implicit none 
C-------------------
      INTEGER, external :: IRPDSO
C----------------------------------

      integer, intent(in) :: popp, popq, popr, pops 
      integer, intent(in) :: ioppq, ioprs, icpq, icrs
      integer, intent(in) :: sym
      double precision, intent(in),device :: iin(*)
      double precision, intent(out),device :: iout(*)
      double precision, allocatable,device :: itmp(:)
C
      integer irrepr, irrepl, i
      integer nsize1, nsize2, ioff, ioff1
      integer nummax, istat, numl, numr
      nummax = 0
      if(icpq/=1.and.icrs==1) then 
       ioff = 1
       do irrepr = 1, nirrep 
          irrepl = dirprd(irrepr, sym) 
          numl = irpdso(irrepl, popp, popq, ioppq) 
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
C          write(6,*)'symtr1so 1'
          call symtr1sodg(irrepr, popr, pops, numl, iin(ioff),
     &            iout(ioff)) 
C          write(6,*)'symtr1so 2'
          ioff = ioff + nsize1 
       enddo
C
      elseif(icpq==1.and.icrs/=1) then
C
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr3sodg(irrepl, popp, popq, numl, numr, iin(ioff),
     &            iout(ioff))
          ioff = ioff + nsize1
       enddo
C
      elseif(icpq==1.and.icrs==1) then 
C
       nsize1= 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = max(nsize1, numl*numr) 
       enddo 
       allocate(itmp(nsize1))
       ioff = 1
       do irrepr = 1, nirrep
          irrepl = dirprd(irrepr, sym)
          numl = irpdso(irrepl, popp, popq, ioppq)
          numr = irpdso(irrepr, popr, pops, ioprs)
          nsize1 = numl*numr
          call symtr1sodg(irrepr, popr, pops, numl, iin(ioff),
     &            itmp) 
          call symtr3sodg(irrepl, popp, popq, numl, numr, itmp,iout(ioff))
          ioff = ioff + nsize1
       enddo
       deallocate(itmp)
      endif
C
      return 
      end  
