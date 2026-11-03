! write by lishun, 2024.12.31
! V(pr,qs) = fact2*V(pr,qs)+fact*V(mp,q)*V(mr,s)
! need irrepmp=irrepmr=irreppr=irrepqs=1
!---------------------------------
! modified at 2025.5.19
subroutine VmpqVmrs2Vprqs(Vmpq, Vmrs, symVmpq, symVmrs, pop_m, &
                          pop_p, pop_q, pop_r, pop_s, Vprqs, symVprqs, fact, fact2)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(in) :: Vmpq, Vmrs
    integer*8, intent(in) :: symVmpq, symVmrs
    integer*8, dimension(8), intent(in) :: pop_m, pop_p, pop_q, pop_r, pop_s
    real*8, dimension(*), intent(inout) :: Vprqs
    integer*8, intent(in) :: symVprqs
    real*8, intent(in) :: fact, fact2

    integer*8 :: irrepp, irrepr, irrepm
    integer*8 :: nump, numr, numm
    integer*8 :: irrepq, numq, q
    integer*8 :: irreps, nums, s
    integer*8 :: irrepmp, nummp
    integer*8 :: irrepmr, nummr
    integer*8 :: irreppr, numpr
    integer*8 :: irrepqs, numqs
    integer*8, dimension(8) :: ioffVmpq, ioffVmrs, ioffVprqs
    integer*8, dimension(5) :: idxVprqs ! pointer of V(pr,qs) belong to irrepqs,irreps,s,q,irrepr
    integer*8, dimension(3) :: idxVmpq, idxVmrs ! pointer of V(ab,c) belong to irrepc,c,irrepb
    integer*8, external :: irpdso, idsymsoc, isymoffso

    ioffVmpq = 1
    ioffVmrs = 1
    ioffVprqs = 1
    do irreps = 1, nirrep - 1
        ! pointer of Vmpq belong to irrepq
        irrepq = irreps
        irrepmp = dirprd(symVmpq, irrepq)
        numq = pop_q(irrepq)
        nummp = irpdso(irrepmp, pop_m, pop_p, 1)
        ioffVmpq(irrepq + 1) = ioffVmpq(irrepq) + nummp*numq
        ! pointer of Vmrs belong to irreps
        irrepmr = dirprd(symVmrs, irreps)
        nums = pop_s(irreps)
        nummr = irpdso(irrepmr, pop_m, pop_r, 1)
        ioffVmrs(irreps + 1) = ioffVmrs(irrepmr) + nummr*nums
        ! pointer of Vprqs belong to irrepqs
        irrepqs = irreps
        irreppr = dirprd(symVprqs, irrepqs)
        numqs = irpdso(irrepqs, pop_q, pop_s, 1)
        numpr = irpdso(irreppr, pop_p, pop_r, 1)
        ioffVprqs(irrepqs + 1) = ioffVprqs(irrepqs) + numpr*numqs
    end do

    do irreps = 1, nirrep
        nums = pop_s(irreps)
        irrepmr = dirprd(symVmrs, irreps)
        nummr = irpdso(irrepmr, pop_m, pop_r, 1)

        idxVmrs(1) = ioffVmrs(irreps)
        do s = 1, nums
            idxVmrs(2) = idxVmrs(1) + nummr*(s - 1)

            do irrepq = 1, nirrep
                numq = pop_q(irrepq)
                irrepmp = dirprd(symVmpq, irrepq)
                nummp = irpdso(irrepmp, pop_m, pop_p, 1)
                irrepqs = dirprd(irrepq, irreps)
                irreppr = dirprd(symVprqs, irrepqs)
                numpr = irpdso(irreppr, pop_p, pop_r, 1)

                idxVmpq(1) = ioffVmpq(irrepq)
                idxVprqs(1) = ioffVprqs(irrepqs)
                idxVprqs(2) = idxVprqs(1) + (isymoffso(irreps, irrepqs, pop_q, pop_s, 1) - 1)*numpr
                idxVprqs(3) = idxVprqs(2) + (s - 1)*numq*numpr
                do q = 1, numq
                    idxVmpq(2) = idxVmpq(1) + nummp*(q - 1)
                    idxVprqs(4) = idxVprqs(3) + numpr*(q - 1)

                    do irrepr = 1, nirrep
                        irrepp = dirprd(irreppr, irrepr)
                        irrepm = dirprd(irrepmr, irrepr)
                        numm = pop_m(irrepm)
                        nump = pop_p(irrepp)
                        numr = pop_r(irrepr)

                        idxVmpq(3) = idxVmpq(2) + isymoffso(irrepp, irrepmp, pop_m, pop_p, 1) - 1
                        idxVmrs(3) = idxVmrs(2) + isymoffso(irrepr, irrepmr, pop_m, pop_r, 1) - 1
                        idxVprqs(5) = idxVprqs(4) + isymoffso(irrepr, irreppr, pop_p, pop_r, 1) - 1

                        ! V(pr,qs) = fact2*V(pr,qs)+fact*V(mp,q)*V(mr,s)
                        call dgemm('T', 'N', nump, numr, numm, &
                                   fact, Vmpq(idxVmpq(3)), max(1, numm), Vmrs(idxVmrs(3)), max(1, numm), &
                                   fact2, Vprqs(idxVprqs(5)), max(1, nump))
                    end do
                end do
            end do
        end do
    end do

end subroutine VmpqVmrs2Vprqs
