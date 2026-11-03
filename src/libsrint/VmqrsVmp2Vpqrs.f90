! write by lishun, 2024.12.5
! V3(pq,rs) = fact2*V3(pq,rs)+fact*V1(mq,rs)*V2(m,p)
! need ioppq=iopmq=1
!--------------------------------
! modified at 2025.5.19
subroutine VmqrsVmp2Vpqrs(Vmqrs, Vmp, symVmqrs, symVmp, pop_m, &
                          pop_p, pop_q, pop_r, pop_s, ioprs, Vpqrs, symVpqrs, fact, fact2)
    use mod_symm
    implicit none

    real*8, dimension(*), intent(in) :: Vmqrs, Vmp
    integer*8, intent(in) :: symVmqrs, symVmp
    integer*8, dimension(8), intent(in) :: pop_m, pop_p, pop_q, pop_r, pop_s
    integer*8, intent(in) :: ioprs ! 0 for r<s, 1 for rs
    real*8, dimension(*), intent(inout) :: Vpqrs
    integer*8, intent(in) :: symVpqrs
    real*8, intent(in) :: fact, fact2

    integer*8 :: symm
    integer*8 :: irrepRS, numRS, rs
    integer*8 :: irrepMQ, numMQ
    integer*8 :: irrepPQ, numPQ
    integer*8 :: irrepQ, numQ
    integer*8 :: irrepM, numM
    integer*8 :: irrepP, numP
    integer*8, dimension(8) :: ioffV1, ioffV3 ! position of Vmqrs and Vpqrs belong to an irrepRS
    integer*8, dimension(3) :: idxV1, idxV3 ! position of V belong to an irrepRS,RS,irrepQ
    integer*8 :: idxV2 ! position of Vmp belong to an irrepP
    integer*8 :: lenV1, lenV2, lenV3, iend! length of V
    integer*8, external :: irpdso, idsymsoc, isymoffso

    symm = dirprd(symVmqrs, symVmp)
    if (symm /= symVpqrs) then
        write (6, '(A,3I2)') 'ERROR: symmetry error in VmqrsVmp2Vpqrs', symVmqrs, symVmp, symVpqrs
        stop
    end if
    lenV1 = idsymsoc(symVmqrs, pop_m, pop_q, 1, pop_r, pop_s, ioprs)
    lenV2 = irpdso(symVmp, pop_m, pop_p, 1)
    lenV3 = idsymsoc(symVpqrs, pop_p, pop_q, 1, pop_r, pop_s, ioprs)

    !------ give ioffV1 and ioffV3
    ioffV1(1) = 1
    ioffV3(1) = 1
    do irrepRS = 1, nirrep - 1
        numRS = irpdso(irrepRS, pop_r, pop_s, ioprs)
        ! for V1
        irrepMQ = dirprd(irrepRS, symVmqrs)
        numMQ = irpdso(irrepMQ, pop_m, pop_q, 1)
        ioffV1(irrepRS + 1) = ioffV1(irrepRS) + numMQ*numRS
        ! for V3
        irrepPQ = dirprd(irrepRS, symVpqrs)
        numPQ = irpdso(irrepPQ, pop_p, pop_q, 1)
        ioffV3(irrepRS + 1) = ioffV3(irrepRS) + numPQ*numRS
    end do

    do irrepRS = 1, nirrep
        numRS = irpdso(irrepRS, pop_r, pop_s, ioprs)
        irrepMQ = dirprd(irrepRS, symVmqrs)
        numMQ = irpdso(irrepMQ, pop_m, pop_q, 1)
        irrepPQ = dirprd(irrepRS, symVpqrs)
        numPQ = irpdso(irrepPQ, pop_p, pop_q, 1)

        idxV1(1) = ioffV1(irrepRS)
        idxV3(1) = ioffV3(irrepRS)
        do rs = 1, numRS

            idxV1(2) = idxV1(1) + numMQ*(rs - 1)
            idxV3(2) = idxV3(1) + numPQ*(rs - 1)
            do irrepQ = 1, nirrep
                numQ = pop_q(irrepQ)
                irrepM = dirprd(irrepQ, irrepMQ)
                numM = pop_m(irrepM)
                irrepP = dirprd(irrepQ, irrepPQ)
                numP = pop_p(irrepP)

                idxV1(3) = idxV1(2) + isymoffso(irrepQ, irrepMQ, pop_m, pop_q, 1) - 1
                idxV3(3) = idxV3(2) + isymoffso(irrepQ, irrepPQ, pop_p, pop_q, 1) - 1
                idxV2 = isymoffso(irrepP, symVmp, pop_m, pop_p, 1)

                ! V3(pq,rs) = fact2*V3(pq,rs)+fact*V2(m,p)*V1(x,rs)
                call dgemm('T', 'N', numP, numQ, numM, &
                           fact, Vmp(idxV2), max(1, numM), Vmqrs(idxV1(3)), max(1, numM), &
                           fact2, Vpqrs(idxV3(3)), max(1, numP))
            end do
        end do
    end do
end subroutine VmqrsVmp2Vpqrs
