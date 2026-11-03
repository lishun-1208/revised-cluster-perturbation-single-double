 subroutine splitt1t2(t1r, t1i, t2r, t2i, icore) 
! split nspin integral to spin integral.
! intype for close shell
! FOR P==Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB
! FOR P==Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA
! FOR P, Q.and.R==S. AAAA, AAAB, ABAA, AABB, ABAB, BAAA
! FOR P, Q.and.R, S. AAAA, AAAB, ABAA, AABB, ABAB, AABA, BAAA, ABBA
!
! example::
! example:: T2(abij) ABIJ, ABIj, AbIJ, ABij, AbIj
! example:: R2(abci) ABCI, ABCi, AbCI, ABci, AbCi, ABcI
! example:: W2(mbej) MBEJ, MBEj, MbEJ, MBej, MbEj, MBeJ, mBEJ, MbeJ
! ioppq_spin = 1; ioprs_spin = 1; sym_spin = symry
! ioppq_ns   = 0; ioprs_ns   = 0; sym_nspin = 1 
! call splitintegral(t2r, t2aaabr, vrta, vrtb, vrta, vrta, ioppq_spin, &
!                    popa, popb, popa, popb, ioprs_spin, symspin, &
!                    nv, nv, 0, no, no, 0, sym_nspin, 'AAAB')
!
! intype for open shell
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
 double precision, intent(inout) :: icore
!
 integer :: popp_s(8), popq_s(8), popr_s(8), pops_s(8)
 integer :: irreppq, irreprs, irrepp, irrepq, irrepr, irreps
 integer :: irreppqn, irreprsn, irreppn, irrepqn, irreprn, irrepsn
 integer :: numpn, numqn, numrn, numsn, numpqn, numrsn
 integer :: p, q, r, s, numpq, numrs, nump, numq, numr, nums
 integer :: ioffns(8), ioffnspin, ioffspin
 integer :: irrepqp, irrepsr, numqp, numsr 
 integer :: ioffin(9), ioffout(9), ioffin2(9)
 integer :: ioffp, ioffq, ioffr, ioffs, idx, idxn
 integer :: nsize, sym2, idxn2
 integer :: numpqn2, numrsn2, ioffns2(8), ioffnspin2
 character*4 str
!
 integer, external :: irpdso, idsymsoc, isymoffso
 real*8, external :: nrm2so 
!
 call clmoioso_sd(1    , icore, vrta, zorb, 1, popa, zorb, 1, 'T1AAR___')
 call clmoioso_sd(symrz, icore, vrta, zorb, 1, popa, zorb, 1, 'T1AAI___')
 call clmoioso_sd(symry, icore, vrta, zorb, 1, popa, zorb, 1, 'T1ABR___')
 call clmoioso_sd(symrx, icore, vrta, zorb, 1, popa, zorb, 1, 'T1ABI___')

 call clmoioso_sd(1    , icore, vrta, vrta, 0, popa, popa, 0, 'T2AAAAR_')
 call clmoioso_sd(symrz, icore, vrta, vrta, 0, popa, popa, 0, 'T2AAAAI_')
 call clmoioso_sd(1    , icore, vrta, vrta, 0, popb, popb, 0, 'T2AABBR_')
 call clmoioso_sd(symrz, icore, vrta, vrta, 0, popb, popb, 0, 'T2AABBI_')
 call clmoioso_sd(1    , icore, vrta, vrtb, 1, popa, popb, 1, 'T2ABABR_')
 call clmoioso_sd(symrz, icore, vrta, vrtb, 1, popa, popb, 1, 'T2ABABI_')
 call clmoioso_sd(symry, icore, vrta, vrta, 0, popa, popb, 1, 'T2AAABR_')
 call clmoioso_sd(symrx, icore, vrta, vrta, 0, popa, popb, 1, 'T2AAABI_')
 call clmoioso_sd(symry, icore, vrta, vrtb, 1, popa, popb, 0, 'T2ABAAR_')
 call clmoioso_sd(symrx, icore, vrta, vrtb, 1, popa, popb, 0, 'T2ABAAI_')


 call splitintegral(t1r, icore, vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, 1, &
                     nv, no, 1, zorb, zorb, 1, 1, 'AAAA')
 call putallso(icore, 1, 'T1AAR___')
 write(6,*)'norm of t1aar', nrm2so(irpdso(1,vrta,popa,1),icore,1)
 call splitintegral(t1r, icore, vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symry, &
                     nv, no, 1, zorb, zorb, 1, 1, 'ABAA')
 call putallso(icore, 1, 'T1ABR___')
 write(6,*)'norm of t1abr', nrm2so(irpdso(symry,vrta,popa,1),icore,1)

 call splitintegral(t1i, icore, vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrz, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'AAAA')
 call putallso(icore, 1, 'T1AAI___')
 write(6,*)'norm of t1aai', nrm2so(irpdso(symrz,vrta,popa,1),icore,1)
 call splitintegral(t1i, icore, vrta, vrtb, popa, popb, 1, &
                     zorb, zorb, zorb, zorb, 1, symrx, &
                     nv, no, 1, zorb, zorb, 1, symrz, 'ABAA')
 call putallso(icore, 1, 'T1ABI___')
 write(6,*)'norm of t1abi', nrm2so(irpdso(symrx,vrta,popa,1),icore,1)
!--------------------------------------------------
 call splitintegral(t2r, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1, &
                     nv, nv, 0, no, no, 0, 1, 'AAAA')
 call putallso(icore, 1, 'T2AAAAR_')
 write(6,*)'norm of t2aaaar ', nrm2so(idsymsoc(1,vrta,vrta,0,popa,popa,0),icore,1)*2.d0
 call splitintegral(t2r, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symry, &
                     nv, nv, 0, no, no, 0, 1, 'AAAB')
 call putallso(icore, 1, 'T2AAABR_')
 write(6,*)'norm of t2aaabr ', nrm2so(idsymsoc(symry,vrta,vrta,0,popa,popa,1),icore,1)*dsqrt(2.d0)
 call splitintegral(t2r, icore, vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symry, &
                     nv, nv, 0, no, no, 0, 1, 'ABAA')
 call putallso(icore, 1, 'T2ABAAR_')
 write(6,*)'norm of t2abaar ', nrm2so(idsymsoc(symry,vrta,vrta,1,popa,popa,0),icore,1)*dsqrt(2.d0)
 call splitintegral(t2r, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, 1, &
                     nv, nv, 0, no, no, 0, 1, 'AABB')
 call putallso(icore, 1, 'T2AABBR_')
 write(6,*)'norm of t2aabbr ', nrm2so(idsymsoc(1,vrta,vrta,0,popa,popa,0),icore,1)*2.d0
 call splitintegral(t2r, icore, vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 1, 1, &
                     nv, nv, 0, no, no, 0, 1, 'ABAB')
 call putallso(icore, 1, 'T2ABABR_')
 write(6,*)'norm of t2ababr ', nrm2so(idsymsoc(1,vrta,vrta,1,popa,popa,1),icore,1)



 call splitintegral(t2i, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no, no, 0, symrz, 'AAAA')
 call putallso(icore, 1, 'T2AAAAI_')
 write(6,*)'norm of t2aaaai ', nrm2so(idsymsoc(symrz,vrta,vrta,0,popa,popa,0),icore,1)*2.d0
 call splitintegral(t2i, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 1, symrx, &
                     nv, nv, 0, no, no, 0, symrz, 'AAAB')
 call putallso(icore, 1, 'T2AAABI_')
 write(6,*)'norm of t2aaabi ', nrm2so(idsymsoc(symrx,vrta,vrta,0,popa,popa,1),icore,1)*dsqrt(2.d0)
 call splitintegral(t2i, icore, vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 0, symrx, &
                     nv, nv, 0, no, no, 0, symrz, 'ABAA')
 call putallso(icore, 1, 'T2ABAAI_')
 write(6,*)'norm of t2abaai ', nrm2so(idsymsoc(symrx,vrta,vrta,1,popa,popa,0),icore,1)*dsqrt(2.d0)
 call splitintegral(t2i, icore, vrta, vrtb, vrta, vrtb, 0, &
                     popa, popb, popa, popb, 0, symrz, &
                     nv, nv, 0, no, no, 0, symrz, 'AABB')
 call putallso(icore, 1, 'T2AABBI_')
 write(6,*)'norm of t2aabbi ', nrm2so(idsymsoc(symrz,vrta,vrta,0,popa,popa,0),icore,1)*2.d0
 call splitintegral(t2i, icore, vrta, vrtb, vrta, vrtb, 1, &
                     popa, popb, popa, popb, 1, symrz, &
                     nv, nv, 0, no, no, 0, symrz, 'ABAB')
 call putallso(icore, 1, 'T2ABABI_')
 write(6,*)'norm of t2ababi ', nrm2so(idsymsoc(symrz,vrta,vrta,1,popa,popa,1),icore,1)
 
 return
 end
