 subroutine combt1t2(t1r, t1i, t2r, t2i, icore) 
!! combine spin integral to nspin integral.
! iout for spin integral, iin for nspin integral
! intype for close shell
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA
!
! spin case need to be put in nspin integral time-reversal symmetry must be employed
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB! ABBB, BBAA, BBAB, BBBB 
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB, AABA, ABBA, BBBA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA, BBAA, BBAB, BBBB, BAAA, BAAB, BABB
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA, BBBB, 
!                    AABA, ABBA, BBBA, BAAA, BAAB, BABB, BABA
! example::
! example:: T2(abij) ABIJ, ABIj, AbIJ, ABij, AbIj
! example:: R2(abci) ABCI, ABCi, AbCI, ABci, AbCi, ABcI
! example:: W2(mbej) MBEJ, MBEj, MbEJ, MBej, MbEj, MBeJ, mBEJ, MbeJ
!     ioppq_spin = 1; ioprs_spin = 1; sym_spin = symry
!     ioppq_ns   = 0; ioprs_ns   = 0; sym_nspin = 1 
!     call zeroso(t2r,idsymsoc(sym_spin,nv,nv,ioppq_ns,no,no,ioprs_ns)
!     call combintegral(t2r, t2aaabr, vrta, vrtb, vrta, vrta, ioppq_spin, &
!                        popa, popb, popa, popb, ioprs_spin, symspin, &
!                        nv, nv, 0, no, no, 0, sym_nspin, 'AAAB')
! type for open shell
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, ABBB, BBAA, BBAB, BBBB, AABA, ABBA, BBBA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA, BBAA, BBAB, BBBB, BAAA, BAAB, BABB
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA, BBBB, 
!                    AABA, ABBA, BBBA, BAAA, BAAB, BABB, BABA
 use mod_iop
 use mod_size
 use mod_orbit
 use mod_ioff
 use mod_symm

 implicit none
!
 double precision, intent(inout) :: t1r(*), t1i(*), t2r(*), t2i(*) 
 double precision, intent(inout) :: icore(*)
!
 integer :: popp_s(8), popq_s(8), popr_s(8), pops_s(8)
 integer :: irreppq, irreprs, irrepp, irrepq, irrepr, irreps
 integer :: irreppqn, irreprsn, irreppn, irrepqn, irreprn, irrepsn
 integer :: numpn, numqn, numrn, numsn, numpqn, numrsn
 integer :: i, p, q, r, s, numpq, numrs, nump, numq, numr, nums
 integer :: ioffns(8), ioffnspin, ioffspin
 integer :: irrepqp, irrepsr, numqp, numsr 
 integer :: ioffin(9), ioffout(9), ioffin2(9)
 integer :: ioffp, ioffq, ioffr, ioffs, idx, idxn
 integer :: nsize, sym2, idxn2, i0, i1, i2, i3, nsize2
 integer :: numpqn2, numrsn2, ioffns2(8), ioffnspin2
 integer :: irrepna, irrepnb
 character*4 str
!
 integer, external :: irpdso, idsymsoc, isymoffso
 real*8, external :: nrm2so 
!
 do irrepna = 1, nirrep
    irrepnb = dirprd(irrepna, symry)
    nv(irrepna) = vrta(irrepna)+vrtb(irrepnb)
    no(irrepna) = popa(irrepna)+popb(irrepnb)
 enddo

! zorb(1:8) = 0 
! zorb(1) = 1
 nsize = 0
 nsize2 = 0 
 do irrepp = 1, nirrep 
    nsize = max(nsize, idsymsoc(irrepp, nv, nv, 1, no, no, 1))
    nsize2 = max(nsize, idsymsoc(irrepp, vrtb, vrtb, 1, popa, popa, 1)) 
 enddo
 i0 = 1 
 i1 = i0 + nsize 
 i2 = i1 + nsize2
 !t1r AA, AB, BA=-conj(AB), BB=conj(AA)
 call getallso(icore(i1), 1, 'T1AAR___')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, 1, &
                     nv, no, 1, zorb, zorb, 1, 1, 'AAAA')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, 1, &
                     nv, no, 1, zorb, zorb, 1, 1, 'BBAA')
 call getallso(icore(i1), 1, 'T1ABR___')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symry, &
                     nv, no, 1, zorb, zorb, 1, 1, 'ABAA')
 call scalso(irpdso(symry, nv, no, 1), -1.d0, icore(i1), 1)
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symry, &
                     nv, no, 1, zorb, zorb, 1, 1, 'BAAA')
 write(6,*)'norm of t1r 1', nrm2so(irpdso(1,nv,no,1),icore(i0),1)
 write(6,*)'norm of t1r 2', nrm2so(irpdso(1,nv,no,1),t1r      ,1)
!do i = 1, irpdso(1, nv, no, 1) 
!   if(abs(icore(i0-1+i)-t1r(i))>1.d-10) then
!      write(6,'(A, I10, 2F18.12)')'test',i,t1r(i), icore(i)
!   endif
!enddo
 ! t1i AA, AB, BA=conj(AB), BB=-conj(AA)
 call getallso(icore(i1), 1, 'T1AAI___')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrz, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'AAAA')
 call scalso(irpdso(symrz, nv, no, 1), -1.d0, icore(i1), 1)
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrz, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'BBAA')
 call getallso(icore(i1), 1, 'T1ABI___')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrx, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'ABAA')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrx, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'BAAA')
 write(6,*)'norm of t1i 1', nrm2so(irpdso(symrz,nv,no,1),icore(i0),1)
 write(6,*)'norm of t1i 2', nrm2so(irpdso(symrz,nv,no,1),t1i      ,1)

! FOR T2. ABIJ, ABIj, AbIJ, ABij, AbIj
! Abij = -conj(aBIJ) = conj(BaIJ); abIJ = conj(ABij); abIj=-conj(ABiJ)=conj(ABJi), abij=conj(ABIJ)
 call zeroso(icore(i0), nsize)
 call getallso(icore(i1), 1, 'T2AAAAR_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1    , &
                     nv, nv, 0, no  , no  , 0, 1    , 'AAAA')
! write(6,*)'norm of t2r 1', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abij=conj(ABIJ)
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1    , &
                     nv, nv, 0, no  , no  , 0, 1    , 'BBBB')
! write(6,*)'norm of t2r 2', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2AAABR_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symry, &
                     nv, nv, 0, no  , no  , 0, 1    , 'AAAB')
! write(6,*)'norm of t2r ', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abIj=-conj(ABiJ)=conj(ABJi)
 call isymtrso(icore(i1), icore(i2), vrta, vrta, 0, 0, popa, popb, 1, 1, symry) 
 call combintegral(icore(i0), icore(i2), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symry, &
                     nv, nv, 0, no  , no  , 0, 1    , 'BBAB')
! write(6,*)'norm of t2r 4', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2ABAAR_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symry, &
                     nv, nv, 0, no  , no  , 0, 1    , 'ABAA')
! write(6,*)'norm of t2r 5', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !Abij = -conj(aBIJ) = conj(BaIJ)
 call isymtrso(icore(i1), icore(i2), vrta, vrta, 1, 1, popa, popb, 0, 0, symry) 
 call combintegral(icore(i0), icore(i2), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symry, &
                     nv, nv, 0, no  , no  , 0, 1    , 'ABBB')
! write(6,*)'norm of t2r 6', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2AABBR_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1    , &
                     nv, nv, 0, no  , no  , 0, 1    , 'AABB')
! write(6,*)'norm of t2r 7', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abIJ = conj(ABij) 
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1    , &
                     nv, nv, 0, no  , no  , 0, 1    , 'BBAA')

! write(6,*)'norm of t2r 8', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 call getallso(icore(i1), 1, 'T2ABABR_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 1, 1    , &
                     nv, nv, 0, no  , no  , 0, 1    , 'ABAB')
 write(6,*)'norm of t2r 9', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 write(6,*)'norm of t2r 0', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),t2r,1)
!----------------------------------------------------------------
!----------------------------------------------------------------
 call zeroso(icore(i0), nsize)
 call getallso(icore(i1), 1, 'T2AAAAI_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no  , no  , 0, symrz, 'AAAA')
! write(6,*)'norm of t2r 1', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abij=conj(ABIJ)
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no  , no  , 0, symrz, 'BBBB')
! write(6,*)'norm of t2r 2', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2AAABI_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symrx, &
                     nv, nv, 0, no  , no  , 0, symrz, 'AAAB')
! write(6,*)'norm of t2r ', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abIj=-conj(ABiJ)=conj(ABJi)
 call isymtrso(icore(i1), icore(i2), vrta, vrta, 0, 0, popa, popb, 1, 1, symrx)
 call scalso(idsymsoc(symrx,vrta,vrta,0,popa,popb,1),-1.d0, icore(i2), 1)
 call combintegral(icore(i0), icore(i2), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symrx, &
                     nv, nv, 0, no  , no  , 0, symrz, 'BBAB')
! write(6,*)'norm of t2r 4', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2ABAAI_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symrx, &
                     nv, nv, 0, no  , no  , 0, symrz, 'ABAA')
! write(6,*)'norm of t2r 5', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !Abij = -conj(aBIJ) = conj(BaIJ)
 call isymtrso(icore(i1), icore(i2), vrta, vrta, 1, 1, popa, popb, 0, 0, symry)
 call combintegral(icore(i0), icore(i2), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symrx, &
                     nv, nv, 0, no  , no  , 0, symrz, 'ABBB')
! write(6,*)'norm of t2r 6', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)

 call getallso(icore(i1), 1, 'T2AABBI_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no  , no  , 0, symrz, 'AABB')
! write(6,*)'norm of t2r 7', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 !abIJ = conj(ABij) 
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no  , no  , 0, symrz, 'BBAA')

! write(6,*)'norm of t2r 8', nrm2so(idsymsoc(1,nv,nv,0,no,no,0),icore(i0),1)
 call getallso(icore(i1), 1, 'T2ABABI_')
 call combintegral(icore(i0), icore(i1), vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 1, symrz, &
                     nv, nv, 0, no  , no  , 0, symrz, 'ABAB')
 write(6,*)'norm of t2i 9', nrm2so(idsymsoc(symrz,nv,nv,0,no,no,0),icore(i0),1)
 write(6,*)'norm of t2i 0', nrm2so(idsymsoc(symrz,nv,nv,0,no,no,0),t2i,1)
 
 !------------------------------------------------------------------ 
 return
 end
