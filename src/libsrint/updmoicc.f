      subroutine updmoicc(numl,numr,irrepx,list,icore)
c
c updates the moio and moiowd vectors of tmpfile 
c icore: memory that need to allocate the space.
c numl: the left side number of wpqrs with symmetry = irrepx 
c numr: the right side number of wpqrs with symmetry = irrepx 
c list: to define the spincase of tmpfile 
cend
c      use commonmod
      use mod_size
      use mod_symm
      implicit none 
c-------------------
      integer, external :: irpdso
      integer :: moionuml, moionumr, rec_cc
      common/moiocc/moionuml(10), moionumr(10), rec_cc(10)
c----------------------------------
      integer, intent(in) :: numl, numr, list, irrepx
      double precision, intent(inout) :: icore(*)
c
      integer ifile
      integer nwords 
      integer nfull, npart, totrec, totwrd, listrec
      integer iwrite, iout, irecord 
c
c 
c assign moio(lstspn,lstnum) and moiowd(lstspn,lstnum) and the file number
c
      ifile=412
      moionuml(irrepx)=numl
      moionumr(irrepx)=numr
c
c compute total number of words in list.
c
      nwords=numl*numr
c
c compute the number of full records this will take up.  usually this
c   will be zero.  also compute how many words of a partial record are
c   required (often just nwords).
c   
      nfull=nwords/iprcwd
      npart=mod(nwords,iprcwd)
c
c increment totrec by the number of full records that this will take up.
c
      totrec = 1+nfull
c
c compute position of last word written in partial record.  if the logic
c  record goes across a physical record boundary, increment totrec by on
c
      totwrd = npart
      rec_cc(irrepx+1) = rec_cc(irrepx) + totrec
      listrec=rec_cc(irrepx)
c
c write an empty buffer to record record if it is beyond eof.
c
      call zero(icore, iprcwd)
        iwrite=listrec + totrec
        do iout=1, totrec 
           irecord=iwrite-iout+1
           call wrdirso(ifile,irecord,icore,iprcwd)
        enddo
      return
      end
