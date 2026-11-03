      subroutine inipckcc(irrepx, icore, popp, popq, ioppq, 
     &                    popr, pops, ioprs)
c-this is to drives the creation of moio pointers for symmetry packed
c-lists
c-irrepx: the symmetry of the tmpfile wpqrs with defined spin case
c-icore: memory that need to allocate the space
c-popp, q, r, s: the orbital type of the four index. like 
c--------------------pop and vrt in common block
c-ioppq: if ioppq=0, then popp<popq, if ioppq=1, then popp, popq
c-ioppq: if ioppq=2, then popp<=popq
c-ioprs: if ioprs=0, then popr<pops, if ioprs=1, then popr, pops
c-ioprs: if ioprs=2, then popr<=pops
c-list: to define the spincase of tmpfile

      use mod_symm
      implicit none
c-------------------
      integer, external :: irpdso
c
      integer :: moionuml, moionumr, rec_cc
      common/moiocc/moionuml(10), moionumr(10), rec_cc(10)
c----------------------------------
      integer, intent(in) :: ioppq, ioprs 
      integer, intent(in) :: popp(8), popq(8)
      integer, intent(in) :: popr(8), pops(8)
      integer, intent(in) :: irrepx
      double precision, intent(in) :: icore(*)
c
      integer irrepr, irrepl, numl, numr
c
      do 10 irrepr=1,nirrep
            irrepl=dirprd(irrepr,irrepx)
            numl=irpdso(irrepl, popp, popq, ioppq) 
            numr=irpdso(irrepr, popr, pops, ioprs)
            call updmoicc(numl,numr,irrepr,icore)
10    continue
      return
      end
