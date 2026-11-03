      subroutine putlstcc(zlist,nfirst,nblock,ispin,irrepx,filename,
     &                   popp, popq, ioppq, popr, pops, ioprs)
!
! this subroutine retrieves nblock distributions of a particular integra
!  list, beginning with distribution nfirst.  the values are returned in
!  vector zlist.
!
!   output:
!
!         zlist - the target vector used in i/o operations.
!
!   input parameters:
!
!         nfirst- the number of the first distribution which will be
!                 written to or read from disk.
!         nblock- the number of distributions which will be handled.
!         ispin - the left-hand index of the list [moio(ispin,itype)].
!         filename  
!
!end
      use mod_symm
      use mod_size
      use mod_iop
      implicit none
      integer, intent(in)::irrepx, ioppq, ioprs
      integer, intent(in)::popp(8), popq(8), popr(8), pops(8)
      integer, intent(in)::nfirst, nblock, ispin
      character(len=8), intent(in) :: filename
!
      double precision, intent(inout) :: zlist(*) 
!     
      integer :: moionuml, moionumr, rec_cc
      common/moiotmp/moionuml(10), moionumr(10), rec_cc(10)
!
      double precision, allocatable :: tmp(:)
!
      integer nwords, nstart, nend, nwrdf, nwrdf2
      integer  i, j, ifile, listrec
      integer recstart, recend, istat
      integer iout, irecord, istart
      integer, external :: irpdso
!
! deal with ridiculous list numbers
!
      call moiocc_calc(popp, popq, ioppq, popr, pops, ioprs, irrepx) 
      open(unit=415, file=filename, form='unformatted', access='direct',
     &            recl=prec*iprcwd)
      ifile = 415 
      listrec = rec_cc(ispin)
! define the started rec and ended rec , one should note that recstart
! and recend are both integer

      nwords = nblock*moionuml(ispin) 

      nstart = (nfirst-1)*moionuml(ispin) 

      nend = (nfirst+nblock-1)*moionuml(ispin) 

      recstart=(nstart-1)/iprcwd 
      nwrdf = nstart - iprcwd*recstart
! 
      recstart=(nstart-1)/iprcwd + 1 + listrec

! 
      recend = (nend-1)/iprcwd
      nwrdf2 = nend - iprcwd*recend
!
      recend = (nend-1)/iprcwd + 1 + listrec
!
      allocate(tmp(iprcwd), stat = istat) 
! write the first line
      if(recend>recstart) then 
           call rddirso(ifile, recstart, tmp, ichcsz)
           call myicopy(zlist, tmp(nwrdf+1), ichcsz-nwrdf)
           call wrdirso(ifile, recstart, tmp, ichcsz)
! write the second to the penul-timate line
           do iout = 1, recend-recstart-1 
              irecord = recstart + iout
              istart = (iout-1)*ichcsz + ichcsz - nwrdf + 1 
              call wrdirso(ifile,irecord,zlist(istart),ichcsz)
           enddo 
! write the last line 
        call rddirso(ifile, recend, tmp, ichcsz)
        istart = nwords - nwrdf2 + 1
        call myicopy(zlist(istart), tmp, nwrdf2)
        call wrdirso(ifile,recend,tmp,ichcsz)
!
        elseif(recend==recstart) then             
          call rddirso(ifile, recstart, tmp, ichcsz)
          call myicopy(zlist, tmp(nwrdf+1), nwrdf2-nwrdf)
          call wrdirso(ifile,recstart,tmp,ichcsz)
        endif
        deallocate(tmp)
!       endif
!      write(6,*)'========end of putlst======='
      close(415)
300   return
      end
