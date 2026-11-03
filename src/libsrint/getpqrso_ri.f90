 subroutine getpqrso_ri(icore, din, irrepdin, rivv,rivv2, naux8,ioppq)
use mod_iop
use mod_symm
use mod_orbit
use mod_ri 
implicit none 
!-------------------
 
!     INTEGER :: SYMRX, SYMRY, SYMRZ
!     COMMON/SYMRXYZ/ SYMRX, SYMRY, SYMRZ
!
!     INTEGER :: POP, VRT, NT, NF1, NF2
!     COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
!
!     INTEGER :: NSTART, NIRREP, IRREPS, DIRPRD
!     COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
!
!     INTEGER :: IRPDPD, ISYTYP, NTOT
!     COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
!
!     INTEGER :: ISYMOFF
!     COMMON/SYMLOC/ISYMOFF(8,8,25)
!     INTEGER :: PRECSO, ISD, IGPU, TOTGPU
!     COMMON/PRECSOC/PRECSO, ISD, IGPU, TOTGPU
!----------------------------------
 integer, intent(in) ::  din, irrepdin,ioppq 
 real*4, intent(in) :: rivv(*), rivv2(*)
 real*4, intent(out) :: icore(*)
 integer, intent(in) :: naux8(8)

 integer irrep, n, a, b, c, d
 integer numa, numb, numc, numd
 integer numab, numcd, ioff, nabc, ioffd,ioffd2(2000, 8), ioffout
 integer irrepa, irrepb, irrepc, irrepd, numnb, i0, i1
 integer irrepab, irrepcd, nsize, ioffabc, numna
 integer ioffab(8), ioffaa(8), nreadab, ioffvv(8) 
 integer, external :: irpdso, idsymsoc, isymoffso
 real*8, external :: nrm2so

 vrta(1:8)=vrt(1:8,1)
! zorb(1:8) = 0 
! zorb(1) = 1
 ioffab(1) = 1; ioffaa(1) = 1
 ioffvv(1) = 1
 do irrep = 1, nirrep - 1
    numna = irpdso(irrep, vrta, naux8, 1) 
    numb = vrta(irrep)
    ioffvv(irrep+1) = ioffvv(irrep) + numna*numb*isd
 enddo

 do irrepcd = 1, nirrep-1 
    irrepab = dirprd(irrepcd, 1) 
    numab = irpdso(irrepab, vrta, vrta, 1) 
    numcd = irpdso(irrepcd, vrta, vrta, 1) 
    ioffab(irrepcd+1) = ioffab(irrepcd) + numab*numcd*isd
    numab = irpdso(irrepab, vrta, vrta, 0) 
    numcd = irpdso(irrepcd, vrta, vrta, 0) 
    ioffaa(irrepcd+1) = ioffaa(irrepcd) + numab*numcd*isd
 enddo
! <ab|cd> = (ac|bd) 
 nsize = idsymsoc( 1, naux8, zorb, 1, vrta, vrta, 1)
! write(6,*)'norm of rivv', nsize, nrm2so(nsize, rivv, 1), nrm2so(nsize, rivv2,1)
 do irrepd = 1, nirrep
    numnb = irpdso(irrepd, naux8, vrta, 1)
    numd = vrta(irrepd)
!   do irrepb = 1, nirrep
!      irrepac = dirprd(irrepabc, irrepb) 
!    write(6,*)'ioffvv', ioffvv(irrepd)
!    (ac|bd) with special d
!    write(6,*)'vrta', vrta, zorb
    do d = 1, numd
       ioffd = ioffvv(irrepd)-1 + (d-1)*numnb*isd + 1
       ioffd2(d,irrepd) = ioffd
    enddo 
  enddo
       nabc = idsymsoc(irrepdin, vrta, vrta, 1, vrta, zorb, 1)
       ioffd = ioffd2(din, irrepdin)
       irrepd = irrepdin
       i0 = 1
       i1 = i0 + nabc*isd
       call VmnpqVmnrs_to_Vpqrs(rivv, rivv2(ioffd), 1, irrepd, naux8, zorb, 1, vrta, vrta, 1, &
                     vrta, zorb, 1, icore(i1), irrepd, 1.d0, 0.d0)
       call sstgenso(icore(i1), icore,  nsize, vrta, vrta, vrta, zorb, icore,&
                     irrepd, '1324')
!<AB||CD> = <AB|CD> - <BA|CD>
       if(ioppq==0) then 
          call isymtrso(icore, icore(i1), vrta, vrta, 1, 1, vrta, zorb, 1, 0, irrepd)
          call axpyso(nabc, -1.d0, icore(i1), 1, icore, 1)
          call icompso(icore, icore(i1), vrta, vrta, 1, 0, vrta, zorb, 1, 1,irrepd) 
          nabc = idsymsoc(irrepd, vrta, vrta, 0, vrta, zorb, 1) 
          call myicopyso(icore(i1), icore, nabc)
       elseif(ioppq==2) then 
          call isymtrso(icore, icore(i1), vrta, vrta, 1, 1, vrta, zorb, 1, 0, irrepd)
          call axpyso(nabc, 1.d0, icore(i1), 1, icore, 1)
          call icompso(icore, icore(i1), vrta, vrta, 1, 2, vrta, zorb, 1, 1,irrepd) 
          nabc = idsymsoc(irrepd, vrta, vrta, 2, vrta, zorb, 1) 
          call myicopyso(icore(i1), icore, nabc)
       endif
  return
  end
