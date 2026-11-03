      subroutine moiocc_calc(popp, popq, ioppq, popr, pops, ioprs, irrepx)
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
!
      integer  i, j, irrepr, irrepl
      character*2 tmpc
      character*10 filename2
      integer ilength
      integer nwords, nfull, npart, totrec, totwrd
      call izero(moionuml, 10)
      call izero(moionumr, 10)
      call izero(rec_cc, 10)
      do irrepr = 1, nirrep 
         irrepl = dirprd(irrepr, irrepx)
         moionumr(irrepr)=irpdso(irrepr, popr, pops, ioprs)
         moionuml(irrepr)=irpdso(irrepl, popp, popq, ioppq)
         nwords = moionumr(irrepr)*moionuml(irrepr) 
         nfull = nwords/iprcwd
         npart = mod(nwords,iprcwd)
         totrec = 1 + nfull
         totwrd = npart 
         rec_cc(irrepr+1)=rec_cc(irrepr) + totrec
      enddo
      return
!
      end
