      subroutine clmoiocc(irrepx, icore, 
     &           popp, popq, ioppq, popr, pops, ioprs, filename)
!-irrepx: the symmetry of the tmpfile w(p,q,r,s) 
!-icore: memory that need to allocate the space 
!-popp, q, r, s: the orbital type of the four index. like pop and vrt in
!--------------- common block
!-filename: the filename to save the w(p,q,r,s), one should not that the
!length of the filename must be (length=8)
!-ioppq and ioprs: to define p<q or p<=q or p,q
      use mod_symm
      use mod_orbit
      use mod_size
      use mod_iop
      implicit none
!-------------------
      integer, external :: irpdso
!
      integer :: moionuml, moionumr, rec_cc
      common/moiocc/moionuml(10), moionumr(10), rec_cc(10)
!----------------------------------

      integer, intent(in) :: popp(8), popq(8)
      integer, intent(in) :: popr(8), pops(8)
      integer, intent(in) :: irrepx, ioppq, ioprs
      double precision, intent(in) :: icore(*)
      character*8, intent(in) :: filename
!
      integer ioppq(8), ioprs(8), i, j 
      integer irrepr, irrepl, numr, numl
      character*2 tmpc
      character*10 filename2
      integer ilength
      open(unit=412, file=filename, form='unformatted', 
     &            access='direct', recl=iprcwd*prec)
      call izero(moionuml, 10)
      call izero(moionumr, 10)
      call izero(rec_cc, 10)
      rec_cc(1) = 1
!***********************************************************************
      do    irrepr=1,nirrep
            irrepl=dirprd(irrepr,irrepx)
            numl=irpdso(irrepl, popp, popq, ioppq)
            numr=irpdso(irrepr, popr, pops, ioprs)
            call updmoicc(numl,numr,irrepr,icore)
      enddo
      close(412)
!***********************************************************************
      return
      end
