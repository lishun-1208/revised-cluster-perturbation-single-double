      subroutine putallcc(z,irrepx,popp,popq, ioppq,
     &                      popr, pops, ioprs, filename)
!
! this routine writes out all irreps of a symmetry packed distribution.
!
!      z  - vector holding the list.
!      filename - which file to put vector z
!
!end
      use mod_symm
      use mod_orbit
      implicit none 
!-------------------
      integer, external :: irpdso
!
      integer :: moionuml, moionumr, rec_cc
      common/moiocc/moionuml(10), moionumr(10), rec_cc(10)
!
!----------------------------------
      integer, intent(in) :: irrepx, ioppq, ioprs
      integer, intent(in) :: popp(8), popq(8), popr(8), pops(8)
      real*8, intent(inout):: z(*)
      character(len=8), intent(in) :: filename
!
      integer ioff
      integer irrepx, i, j, irrepr, irrepl
      integer numr, numl
!
!  loop over irreps.
!
      ioff=1
      do irrepr=1,nirrep
         irrepl=dirprd(irrepr, irrepx)
!
! get number of distributions and distribution size.
!
         numr=irpdso(irrepr, popr, pops, ioprs)
         numl=irpdso(irrepl, popp, popq, ioppq)
!
! plunk vector on moints.
!
         call putlstcc(z(ioff),1,numr,irrepr,filename, 
     &                 popp, popq, ioppq, popr, pops, ioprs)
         ioff=ioff+numl*numr
      enddo
      return
      end
