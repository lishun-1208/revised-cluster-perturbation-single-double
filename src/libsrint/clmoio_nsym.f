      subroutine clmoio_nsym(icore, 
     &           popp, popq, ioppq, popr, pops, ioprs, filename)
      use mod_symm
      use mod_orbit
      use mod_size
      use mod_iop
      implicit none
!-------------------
      integer, external :: irpdso
!----------------------------------

      integer, intent(in) :: popp, popq
      integer, intent(in) :: popr, pops
      integer, intent(in) :: ioppq, ioprs
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
!
      integer i, j
      character*2 tmpc
      character*10 filename2
      integer ilength, numl, numr
      integer nwords, numpq, numrs 
      integer nfull, npart, totrec, totwrd, listrec
      integer iwrite, iout, irecord, ifile
      if(ioppq==1) then 
         numpq = popp*popq
      elseif(ioppq==0) then 
         numpq = popp*(popq-1)/2
      elseif(ioppq==2) then 
         numpq = popp*(popp+1)/2
      endif

      if(ioprs==1) then 
         numrs = popr*pops
      elseif(ioprs==0) then 
         numrs = popr*(popr-1)/2
      elseif(ioprs==2) then 
         numrs = popr*(popr+1)/2
      endif

      open(unit=412, file=filename, form='unformatted', 
     &            access='direct', recl=iprcwd*prec)
      numl = numpq; numr = numrs
!***********************************************************************
!
! compute total number of words in list.
!
      nwords=numl*numr
!
! compute the number of full records this will take up.  usually this
!   will be zero.  also compute how many words of a partial record are
!   required (often just nwords).
!   
      nfull=nwords/iprcwd
      npart=mod(nwords,iprcwd)
!
!!increment totrec by the number of full records that this will take up.
!
      totrec = 1+nfull
!
! compute position of last word written in partial record.  if the logic
!  record goes across a physical record boundary, increment totrec by on
!
      totwrd = npart
      listrec= 1
!
! write an empty buffer to record record if it is beyond eof.
!
      call zero(icore, iprcwd)
        iwrite=listrec + totrec
        do iout=1, totrec
           irecord=iwrite-iout+1
           call wrdir_nsym(ifile,irecord,icore,iprcwd)
        enddo

      close(412)
      return
      end

      subroutine wrdir_nsym(iunit,irec,ivec,length)
!
      implicit none
      integer, intent(in) :: length, irec, iunit
      double precision, intent(inout) :: ivec(length)
!
      integer irealunit, irealrec, ier
      write(iunit,rec=irec,IOSTAT=IER)ivec
      return
      end
!
      subroutine rddir_nsym(iunit,irec,ivec,length)
!
      implicit none
      integer, intent(in) :: length, irec, iunit
      double precision, intent(inout) :: ivec(length)
!
      integer irealunit, irealrec, ier
      read(iunit,rec=irec,IOSTAT=IER)ivec
      return
      end
!
      subroutine putlst_nsym(zlist,nfirst,nblock,filename,
     &                   popp, popq, ioppq, popr, pops, ioprs)
      use mod_size
      use mod_iop
      integer, intent(in)::ioppq, ioprs
      integer, intent(in)::popp, popq, popr, pops
      integer, intent(in)::nfirst, nblock
      character(len=8), intent(in) :: filename
!
      real*8, intent(inout) :: zlist(*)
!
      integer nwords, nstart, nend, nwrdf, nwrdf2
      integer  i, j, ifile, listrec, numl, numr
      integer recstart, recend, istat
      integer iout, irecord, istart
      real*8, allocatable :: tmp(:)
!
      open(unit=415, file=filename, form='unformatted', access='direct',
     &            recl=prec*iprcwd)
      ifile = 415
!
      if(ioppq==1) then
         numpq = popp*popq
      elseif(ioppq==0) then
         numpq = popp*(popq-1)/2
      elseif(ioppq==2) then
         numpq = popp*(popp+1)/2
      endif

      if(ioprs==1) then
         numrs = popr*pops
      elseif(ioprs==0) then
         numrs = popr*(popr-1)/2
      elseif(ioprs==2) then
         numrs = popr*(popr+1)/2
      endif
!
      allocate(tmp(iprcwd))
!
      listrec = 1 
      nwords = nblock*numpq
      nstart = (nfirst-1)*numpq
      nend = (nfirst+nblock-1)*numpq
!
      recstart = (nstart-1)/iprcwd
      nwrdf = nstart - iprcwd*recstart
!
      recstart = recstart + 1
!
      recend = (nend-1)/iprcwd 
      nwrdf2 = nend - iprcwd*recend 
!
      recend = recend + 1
! write the first line 
      if(recend > recstart) then 
         call rddir_nsym(ifile, recstart, tmp, iprcwd)
         call dcopy(iprcwd-nwrdf,  zlist, 1, tmp(nwrdf+1), 1)
         call wrdir_nsym(ifile, recstart, tmp, ichcsz)
         do iout = 1, recend - recstart - 1
            irecord = recstart + iout 
            istart = (iout-1)*iprcwd + iprcwd - nwrdf + 1
            call wrdir_nsym(ifile, irecord, zlist(istart), iprcwd)
         enddo
         call rddir_nsym(ifile, recend, tmp, iprcwd)  
         istart = nwords - nwrdf2 + 1 
         call dcopy(nwrdf2, zlist(istart), 1, tmp, 1) 
         call wrdir_nsym(ifile, recend, tmp, iprcwd) 
!
       elseif(recend == recstart) then 
         call rddir_nsym(ifile, recstart, tmp, ichcsz)
         call dcopy(nwrdf2-nwrdf, zlist,1, tmp(nwrdf+1), 1)
         call wrdir_nsym(ifile, recstart,tmp,ichcsz)
       endif
       deallocate(tmp)
       close(415) 
       return
       end
!
      subroutine getlst_nsym(zlist,nfirst,nblock,filename,
     &                   popp, popq, ioppq, popr, pops, ioprs)
      use mod_size
      use mod_iop
      integer, intent(in)::ioppq, ioprs
      integer, intent(in)::popp, popq, popr, pops
      integer, intent(in)::nfirst, nblock
      character(len=8), intent(in) :: filename
!
      real*8, intent(inout) :: zlist(*)
!
      integer nwords, nstart, nend, nwrdf, nwrdf2
      integer  i, j, ifile, listrec, numl, numr
      integer recstart, recend, istat
      integer iout, irecord, istart
      real*8, allocatable :: tmp(:)
!
      open(unit=415, file=filename, form='unformatted', access='direct',
     &            recl=prec*iprcwd)
      ifile = 415
!
      if(ioppq==1) then
         numpq = popp*popq
      elseif(ioppq==0) then
         numpq = popp*(popq-1)/2
      elseif(ioppq==2) then
         numpq = popp*(popp+1)/2
      endif

      if(ioprs==1) then
         numrs = popr*pops
      elseif(ioprs==0) then
         numrs = popr*(popr-1)/2
      elseif(ioprs==2) then
         numrs = popr*(popr+1)/2
      endif
!
      allocate(tmp(iprcwd))
!
      listrec = 1
      nwords = nblock*numpq
      nstart = (nfirst-1)*numpq
      nend = (nfirst+nblock-1)*numpq
!
      recstart = (nstart-1)/iprcwd
      nwrdf = nstart - iprcwd*recstart
!
      recstart = recstart + 1
!
      recend = (nend-1)/iprcwd
      nwrdf2 = nend - iprcwd*recend
!
      recend = recend + 1
! write the first line 
      if(recend > recstart) then
         call rddir_nsym(ifile, recstart, tmp, iprcwd)
         call dcopy(iprcwd-nwrdf,  tmp(nwrdf+1), 1, zlist, 1)
         do iout = 1, recend - recstart - 1
            irecord = recstart + iout
            istart = (iout-1)*iprcwd + iprcwd - nwrdf + 1
            call rddir_nsym(ifile, irecord, zlist(istart), iprcwd)
         enddo
         call rddir_nsym(ifile, recend, tmp, iprcwd)
         istart = nwords - nwrdf2 + 1
         call dcopy(nwrdf2, tmp, 1, zlist(istart), 1)
!
       elseif(recend == recstart) then
         call rddir_nsym(ifile, recstart, tmp, ichcsz)
         call dcopy(nwrdf2-nwrdf, tmp(nwrdf+1), 1, zlist,1)
       endif
       deallocate(tmp)
       close(415)
       return
       end

