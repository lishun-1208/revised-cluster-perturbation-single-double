! write by lishun, 2024.10.28
! resort 2e-intergral <pq||rs> form spin to no spin case, in hard disk
! ATTENTION: have problem for open shell abcd
!
! for <ab||ci> (also <ij||ka>), just need <ab||ci>, <AB||CI>, <aB||cI> and <Ab||Ci>
! because <aB||Ci>=-<Ba||Ci>, <Ab||cI>=-<bA||cI>
subroutine resortVabci(file_ns, file_aaaa, file_bbbb, file_abab, file_baba, &
                       popPa, popQa, popRa, popSa, &
                       popPb, popQb, popRb, popSb, &
                       popP, popQ, popR, popS, &
                       iopPQaa, iopRSaa, iopPQbb, iopRSbb, iopPQ, iopRS, &
                       iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    ! 2e-integral file
    character(len=8), intent(in) :: file_ns
    character(len=8), intent(in) :: file_aaaa, file_bbbb, file_abab, file_baba
    ! pop of P,Q,R,S
    integer*8, dimension(8), intent(in) :: popPa, popQa, popRa, popSa
    integer*8, dimension(8), intent(in) :: popPb, popQb, popRb, popSb
    integer*8, dimension(8), intent(in) :: popP, popQ, popR, popS
    ! iopPQ and iopRS for aaaa, bbbb and ns(no symmetry)
    integer*8, intent(in) :: iopPQaa, iopRSaa, iopPQbb, iopRSbb
    integer*8, intent(in) :: iopPQ, iopRS
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: irrepPQ, irrepRS
    integer*8 :: irrepPQaa, irrepPQbb, irrepPQab, irrepPQba
    integer*8 :: numPQ, numPQaa, numPQbb, numPQab, numPQba
    integer*8 :: irrepP, irrepQ, irrepR, irrepS
    integer*8 :: numP, numQ, numR, numS
    integer*8 :: P, Q, R, S
    integer*8 :: irrepPa, irrepQa, irrepRa, irrepSa
    integer*8 :: numPa, numQa, numRa, numSa
    integer*8 :: Pa, Qa, Ra, Sa
    integer*8 :: irrepPb, irrepQb, irrepRb, irrepSb
    integer*8 :: numPb, numQb, numRb, numSb
    integer*8 :: Pb, Qb, Rb, Sb
    ! pointer of <pq||R*> belong to an irrepR,irrepS
    integer*8, dimension(8, 8) :: ioffaaaa, ioffbbbb, ioffabab, ioffbaba
    integer*8, dimension(8, 8) :: ioff
    ! pointer of V(pq,R*) belong to irrepR,R,irrepQ,Q,P, with and with out spin
    integer*8, dimension(5) :: idx, idx_ns
    integer*8 :: nSize
    integer*8 :: ins, iaaaa, ibbbb, iabab, ibaba
    integer*8, external :: idsymsoc, isymoffso, irpdso

    !------- give ioff(irrepR,irrepS), for V(pq,R*)
    ioff(1, 1:8) = 0
    ioffaaaa(1, 1:8) = 0
    ioffbbbb(1, 1:8) = 0
    ioffabab(1, 1:8) = 0
    ioffbaba(1, 1:8) = 0
    do irrepS = 1, nirrep
        do irrepR = 1, nirrep - 1
            irrepRS = dirprd(irrepR, irrepS)
            irrepPQ = dirprd(irrepRS, 1)
            ! ns
            numPQ = irpdso(irrepPQ, popP, popQ, 1)
            numR = popR(irrepR)
            ioff(irrepR + 1, irrepS) = ioff(irrepR, irrepS) + numPQ*numR
            ! aaaa
            numPQaa = irpdso(irrepPQ, popPa, popQa, 1)
            numRa = popRa(irrepR)
            ioffaaaa(irrepR + 1, irrepS) = ioffaaaa(irrepR, irrepS) + numPQaa*numRa
            ! bbbb
            numPQbb = irpdso(irrepPQ, popPb, popQb, 1)
            numRb = popRb(irrepR)
            ioffbbbb(irrepR + 1, irrepS) = ioffbbbb(irrepR, irrepS) + numPQbb*numRb
            ! abab
            numPQab = irpdso(irrepPQ, popPa, popQb, 1)
            numRa = popRa(irrepR)
            ioffabab(irrepR + 1, irrepS) = ioffabab(irrepR, irrepS) + numPQab*numRa
            ! baba
            numPQba = irpdso(irrepPQ, popPb, popQa, 1)
            numRb = popRb(irrepR)
            ioffbaba(irrepR + 1, irrepS) = ioffbaba(irrepR, irrepS) + numPQba*numRb
        end do
    end do

    ins = 1
    ! cycle for irrepS, no symmetry
    do irrepS = 1, nirrep
        irrepSa = irrepS
        irrepSb = dirprd(irrepSa, symry)
        numS = popS(irrepS)
        numSa = popSa(irrepSa)
        numSb = popSb(irrepSb)
        if ((numSa + numSb) /= numS) then
            write (6, *) 'numS error'
            stop
        end if

        !------------- for alpha S
        ! need AAAA, BABA, ABBA=-BABA
        do Sa = 1, numSa
            S = Sa
            nSize = idsymsoc(irrepS, popP, popQ, 1, popR, zorb, 1)
            call dscal(nSize, 0.0d0, iCore(ins), 1)
            iaaaa = ins + nSize
            ibaba = iaaaa + idsymsoc(irrepSa, popPa, popQa, 1, popRa, zorb, 1)
            !------- read integrals for an alpha S
            ! read V(pq,r*), AAAA
            if (iopPQaa == 1) then
                call getpqrso(iCore(iaaaa), Sa, irrepSa, 1, file_aaaa)
            else if (iopPQaa == 0) then
                call getpqrso(iCore(ibaba), Sa, irrepSa, 1, file_aaaa)
                call iexpso(iCore(ibaba), iCore(iaaaa), popPa, popQa, 0, 1, popRa, zorb, 1, 0, irrepSa)
            end if
            ! read V(Pq,R*), BABA
            call getpqrso(iCore(ibaba), Sa, irrepSa, 1, file_baba)

            ! cycle for irrepRS
            do irrepRS = 1, nirrep
                irrepPQ = dirprd(irrepRS, 1)
                numPQ = irpdso(irrepPQ, popP, popQ, 1)

                irrepR = dirprd(irrepS, irrepRS)
                irrepRa = irrepR
                irrepRb = dirprd(irrepRa, symry)
                numR = popR(irrepR)
                numRa = popRa(irrepRa)
                numRb = popRb(irrepRb)
                if ((numRa + numRb) /= numR) then
                    write (6, *) 'numR error'
                    stop
                end if

                idx_ns(1) = ioff(irrepR, irrepS)
                ! cycle for irrepQ
                do irrepQ = 1, nirrep
                    irrepQa = irrepQ
                    irrepQb = dirprd(irrepQa, symry)
                    numQ = popQ(irrepQ)
                    numQa = popQa(irrepQa)
                    numQb = popQb(irrepQb)
                    if ((numQa + numQb) /= numQ) then
                        write (6, *) 'numQ error'
                        stop
                    end if

                    irrepP = dirprd(irrepQ, irrepPQ)
                    irrepPa = irrepP
                    irrepPb = dirprd(irrepPa, symry)
                    numP = popP(irrepP)
                    numPa = popPa(irrepPa)
                    numPb = popPb(irrepPb)
                    if ((numPa + numPb) /= numP) then
                        write (6, *) 'numP error'
                        stop
                    end if

                    idx_ns(2) = idx_ns(1) + isymoffso(irrepQ, irrepPQ, popP, popQ, 1) - 1
                    !------ for AAAA
                    irrepPQaa = dirprd(irrepPa, irrepQa)
                    numPQaa = irpdso(irrepPQaa, popPa, popQa, 1)
                    idx(1) = ioffaaaa(irrepRa, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepQa, irrepPQaa, popPa, popQa, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQaa*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + numPa*(Qa - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + Pa
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(iaaaa - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for BABA
                    irrepPQba = dirprd(irrepPb, irrepQa)
                    numPQba = irpdso(irrepPQba, popPb, popQa, 1)
                    idx(1) = ioffbaba(irrepRb, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepQa, irrepPQba, popPb, popQa, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQba*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + numPb*(Qa - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + Pb
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(ibaba - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for ABBA=-BABA, <pQ||Rs>=-<Qp||Rs>
                    irrepPQab = dirprd(irrepPa, irrepQb)
                    numPQab = irpdso(irrepPQab, popPa, popQb, 1)
                    idx(1) = ioffbaba(irrepRb, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepPa, irrepPQab, popQb, popPa, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQab*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + Qb
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + numQb*(Pa - 1)
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = -iCore(ibaba - 1 + idx(5))
                            end do
                        end do
                    end do
                end do
            end do
            if (iopPQ == 1) then
                ! write V(pq,r*)
                call putpqrso(iCore(ins), S, irrepS, 1, file_ns)
            else if (iopPQ == 0) then
                call icompso(iCore(ins), iCore(iaaaa), popP, popQ, 1, 0, popR, zorb, 1, 1, irrepS)
                call putpqrso(iCore(iaaaa), S, irrepS, 1, file_ns)
            end if
        end do

        !------------- for beta S
        ! need BBBB, ABAB, BAAB=-ABAB
        do Sb = 1, numSb
            S = numSa + Sb
            nSize = idsymsoc(irrepS, popP, popQ, 1, popR, zorb, 1)
            call dscal(nSize, 0.0d0, iCore(ins), 1)
            ibbbb = ins + nSize
            iabab = ibbbb + idsymsoc(irrepSb, popPb, popQb, 1, popRb, zorb, 1)
            !------- read integrals for an beta S
            ! read V(PQ,R*), BBBB
            if (iopPQbb == 1) then
                call getpqrso(iCore(ibbbb), Sb, irrepSb, 1, file_bbbb)
            else if (iopPQbb == 0) then
                call getpqrso(iCore(iabab), Sb, irrepSb, 1, file_bbbb)
                call iexpso(iCore(iabab), iCore(ibbbb), popPb, popQb, 0, 1, popRb, zorb, 1, 0, irrepSb)
            end if
            ! read V(pQ,r*), ABAB
            call getpqrso(iCore(iabab), Sb, irrepSb, 1, file_abab)

            ! cycle for irrepRS
            do irrepRS = 1, nirrep
                irrepPQ = dirprd(irrepRS, 1)
                numPQ = irpdso(irrepPQ, popP, popQ, 1)

                irrepR = dirprd(irrepS, irrepRS)
                irrepRa = irrepR
                irrepRb = dirprd(irrepRa, symry)
                numR = popR(irrepR)
                numRa = popRa(irrepRa)
                numRb = popRb(irrepRb)
                if ((numRa + numRb) /= numR) then
                    write (6, *) 'numR error'
                    stop
                end if

                idx_ns(1) = ioff(irrepR, irrepS)
                ! cycle for irrepQ
                do irrepQ = 1, nirrep
                    irrepQa = irrepQ
                    irrepQb = dirprd(irrepQa, symry)
                    numQ = popQ(irrepQ)
                    numQa = popQa(irrepQa)
                    numQb = popQb(irrepQb)
                    if ((numQa + numQb) /= numQ) then
                        write (6, *) 'numQ error'
                        stop
                    end if

                    irrepP = dirprd(irrepQ, irrepPQ)
                    irrepPa = irrepP
                    irrepPb = dirprd(irrepPa, symry)
                    numP = popP(irrepP)
                    numPa = popPa(irrepPa)
                    numPb = popPb(irrepPb)
                    if ((numPa + numPb) /= numP) then
                        write (6, *) 'numP error'
                        stop
                    end if

                    idx_ns(2) = idx_ns(1) + isymoffso(irrepQ, irrepPQ, popP, popQ, 1) - 1
                    !------ for BBBB
                    irrepPQbb = dirprd(irrepPb, irrepQb)
                    numPQbb = irpdso(irrepPQbb, popPb, popQb, 1)
                    idx(1) = ioffbbbb(irrepRb, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepQb, irrepPQbb, popPb, popQb, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQbb*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + numPb*(Qb - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + Pb
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(ibbbb - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for ABAB
                    irrepPQab = dirprd(irrepPa, irrepQb)
                    numPQab = irpdso(irrepPQab, popPa, popQb, 1)
                    idx(1) = ioffabab(irrepRa, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepQb, irrepPQab, popPa, popQb, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQab*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + numPa*(Qb - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + Pa
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for BAAB, <Pq||rS>=-<qP||rS>
                    irrepPQba = dirprd(irrepPb, irrepQa)
                    numPQba = irpdso(irrepPQba, popPb, popQa, 1)
                    idx(1) = ioffabab(irrepRa, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepPb, irrepPQba, popQa, popPb, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQba*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + Qa
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + numQa*(Pb - 1)
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = -iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                end do
            end do
            if (iopPQ == 1) then
                ! write V(pq,r*)
                call putpqrso(iCore(ins), S, irrepS, 1, file_ns)
            else if (iopPQ == 0) then
                call icompso(iCore(ins), iCore(ibbbb), popP, popQ, 1, 0, popR, zorb, 1, 1, irrepS)
                call putpqrso(iCore(ibbbb), S, irrepS, 1, file_ns)
            end if
        end do
    end do

end subroutine resortVabci

! for <ab||cd> (also <ij||kl> and <ab||ij>), just need <ab||cd>, <AB||CD>, <aB||cD>
! because <Ab||Cd>=<bA||dC>, <aB||Cd>=-<aB||dC>, <Ab||cD>=-<bA||cD>
subroutine resortVabcd(file_ns, file_aaaa, file_bbbb, file_abab, &
                       popPa, popQa, popRa, popSa, &
                       popPb, popQb, popRb, popSb, &
                       popP, popQ, popR, popS, &
                       iopPQaa, iopRSaa, iopPQbb, iopRSbb, iopPQ, iopRS, &
                       iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    ! 2e-integral file
    character(len=8), intent(in) :: file_ns
    character(len=8), intent(in) :: file_aaaa, file_bbbb, file_abab
    ! pop of P,Q,R,S
    integer*8, dimension(8), intent(in) :: popPa, popQa, popRa, popSa
    integer*8, dimension(8), intent(in) :: popPb, popQb, popRb, popSb
    integer*8, dimension(8), intent(in) :: popP, popQ, popR, popS
    ! iopPQ and iopRS for aaaa, bbbb and ns(no symmetry)
    integer*8, intent(in) :: iopPQaa, iopRSaa, iopPQbb, iopRSbb
    integer*8, intent(in) :: iopPQ, iopRS
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: irrepPQ, irrepRS
    integer*8 :: irrepPQaa, irrepPQbb, irrepPQab, irrepPQba
    integer*8 :: numPQ, numPQaa, numPQbb, numPQab, numPQba
    integer*8 :: irrepP, irrepQ, irrepR, irrepS
    integer*8 :: numP, numQ, numR, numS
    integer*8 :: P, Q, R, S
    integer*8 :: irrepPa, irrepQa, irrepRa, irrepSa
    integer*8 :: numPa, numQa, numRa, numSa
    integer*8 :: Pa, Qa, Ra, Sa
    integer*8 :: irrepPb, irrepQb, irrepRb, irrepSb
    integer*8 :: numPb, numQb, numRb, numSb
    integer*8 :: Pb, Qb, Rb, Sb
    ! pointer of <pq||R*> belong to an irrepR,irrepS
    integer*8, dimension(8, 8) :: ioffaaaa, ioffbbbb, ioffabab, ioffbaba
    integer*8, dimension(8, 8) :: ioff
    ! pointer of V(pq,R*) belong to irrepR,R,irrepQ,Q,P, with and with out spin
    integer*8, dimension(5) :: idx, idx_ns
    integer*8 :: i0, nSize
    integer*8 :: ins, iaaaa, ibbbb, iabab
    integer*8, external :: idsymsoc, isymoffso, irpdso

    !------- give ioff(irrepR,irrepS), for V(pq,r*)
    ioff(1, 1:8) = 0
    ioffaaaa(1, 1:8) = 0
    ioffbbbb(1, 1:8) = 0
    ioffabab(1, 1:8) = 0
    ioffbaba(1, 1:8) = 0
    do irrepS = 1, nirrep
        do irrepR = 1, nirrep - 1
            irrepRS = dirprd(irrepR, irrepS)
            irrepPQ = dirprd(irrepRS, 1)
            ! ns
            numPQ = irpdso(irrepPQ, popP, popQ, 1)
            numR = popR(irrepR)
            ioff(irrepR + 1, irrepS) = ioff(irrepR, irrepS) + numPQ*numR
            ! aaaa
            numPQaa = irpdso(irrepPQ, popPa, popQa, 1)
            numRa = popRa(irrepR)
            ioffaaaa(irrepR + 1, irrepS) = ioffaaaa(irrepR, irrepS) + numPQaa*numRa
            ! bbbb
            numPQbb = irpdso(irrepPQ, popPb, popQb, 1)
            numRb = popRb(irrepR)
            ioffbbbb(irrepR + 1, irrepS) = ioffbbbb(irrepR, irrepS) + numPQbb*numRb
            ! abab
            numPQab = irpdso(irrepPQ, popPa, popQb, 1)
            numRa = popRa(irrepR)
            ioffabab(irrepR + 1, irrepS) = ioffabab(irrepR, irrepS) + numPQab*numRa
            ! baba
            numPQba = irpdso(irrepPQ, popPb, popQa, 1)
            numRb = popRb(irrepR)
            ioffbaba(irrepR + 1, irrepS) = ioffbaba(irrepR, irrepS) + numPQba*numRb
        end do
    end do

    ins = 1
    ! cycle for irrepS, no symmetry
    do irrepS = 1, nirrep
        irrepSa = irrepS
        irrepSb = dirprd(irrepSa, symry)
        numS = popS(irrepS)
        numSa = popSa(irrepSa)
        numSb = popSb(irrepSb)
        if ((numSa + numSb) /= numS) then
            write (6, *) 'numS error'
            stop
        end if

        !------------- for alpha S
        ! need AAAA, BABA=ABAB, ABBA=-ABAB
        ! <Pq||Rs>=<qP||sR>, <pQ||Rs>=-<pQ||sR>
        do Sa = 1, numSa
            S = Sa
            nSize = idsymsoc(irrepS, popP, popQ, 1, popR, zorb, 1)
            call dscal(nSize, 0.0d0, iCore(ins), 1)
            iaaaa = ins + nSize
            iabab = iaaaa + idsymsoc(irrepSa, popPa, popQa, 1, popRa, zorb, 1)
            i0 = iabab + idsymsoc(irrepSa, popPb, popQa, 1, popRb, zorb, 1)
            !------- read integrals for an alpha S
            ! read V(pq,r*), AAAA
            if (iopPQaa == 1) then
                call getpqrso(iCore(iaaaa), Sa, irrepSa, 1, file_aaaa)
            else if (iopPQaa == 0) then
                call getpqrso(iCore(iabab), Sa, irrepSa, 1, file_aaaa)
                call iexpso(iCore(iabab), iCore(iaaaa), popPa, popQa, 0, 1, popRa, zorb, 1, 0, irrepSa)
            end if
            ! read V(qP,sR)=V(Pq,Rs)=-V(qP,Rs) for one s
            call getpqs(iCore(iabab), Sa, irrepSa, popQa, popPb, popSa, popRb, &
                        1, 1, file_abab, 1)

            ! cycle for irrepRS
            do irrepRS = 1, nirrep
                irrepPQ = dirprd(irrepRS, 1)
                numPQ = irpdso(irrepPQ, popP, popQ, 1)

                irrepR = dirprd(irrepS, irrepRS)
                irrepRa = irrepR
                irrepRb = dirprd(irrepRa, symry)
                numR = popR(irrepR)
                numRa = popRa(irrepRa)
                numRb = popRb(irrepRb)
                if ((numRa + numRb) /= numR) then
                    write (6, *) 'numR error'
                    stop
                end if

                idx_ns(1) = ioff(irrepR, irrepS)
                ! cycle for irrepQ
                do irrepQ = 1, nirrep
                    irrepQa = irrepQ
                    irrepQb = dirprd(irrepQa, symry)
                    numQ = popQ(irrepQ)
                    numQa = popQa(irrepQa)
                    numQb = popQb(irrepQb)
                    if ((numQa + numQb) /= numQ) then
                        write (6, *) 'numQ error'
                        stop
                    end if

                    irrepP = dirprd(irrepQ, irrepPQ)
                    irrepPa = irrepP
                    irrepPb = dirprd(irrepPa, symry)
                    numP = popP(irrepP)
                    numPa = popPa(irrepPa)
                    numPb = popPb(irrepPb)
                    if ((numPa + numPb) /= numP) then
                        write (6, *) 'numP error'
                        stop
                    end if

                    idx_ns(2) = idx_ns(1) + isymoffso(irrepQ, irrepPQ, popP, popQ, 1) - 1
                    !------ for AAAA
                    irrepPQaa = dirprd(irrepPa, irrepQa)
                    numPQaa = irpdso(irrepPQaa, popPa, popQa, 1)
                    idx(1) = ioffaaaa(irrepRa, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepQa, irrepPQaa, popPa, popQa, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQaa*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + numPa*(Qa - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + Pa
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(iaaaa - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for BABA, V(Pq,R*)=V(qP,*R)
                    irrepPQba = dirprd(irrepPb, irrepQa)
                    numPQba = irpdso(irrepPQba, popPb, popQa, 1)
                    idx(1) = ioffbaba(irrepRb, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepPb, irrepPQba, popQa, popPb, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQba*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + Qa
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + numQa*(Pb - 1)
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for ABBA, <pQ||R*>=-<pQ||*R>
                    irrepPQab = dirprd(irrepPa, irrepQb)
                    numPQab = irpdso(irrepPQab, popPa, popQb, 1)
                    idx(1) = ioffbaba(irrepRb, irrepSa)
                    idx(2) = idx(1) + isymoffso(irrepQb, irrepPQab, popPa, popQb, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQab*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + numPa*(Qb - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + Pa
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = -iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                end do
            end do
            if (iopPQ == 1) then
                ! write V(pq,r*)
                call putpqrso(iCore(ins), S, irrepS, 1, file_ns)
            else if (iopPQ == 0) then
                call icompso(iCore(ins), iCore(iaaaa), popP, popQ, 1, 0, popR, zorb, 1, 1, irrepS)
                call putpqrso(iCore(iaaaa), S, irrepS, 1, file_ns)
            end if
        end do

        !------------- for beta S
        ! need BBBB, ABAB, BAAB
        ! <Pq||rS>=-<qP||rS>
        do Sb = 1, numSb
            S = numSa + Sb
            nSize = idsymsoc(irrepS, popP, popQ, 1, popR, zorb, 1)
            call dscal(nSize, 0.0d0, iCore(ins), 1)
            ibbbb = ins + nSize
            iabab = ibbbb + idsymsoc(irrepSb, popPb, popQb, 1, popRb, zorb, 1)
            !------- read integrals for an beta S
            ! read V(pq,r*), BBBB
            if (iopPQbb == 1) then
                call getpqrso(iCore(ibbbb), Sb, irrepSb, 1, file_bbbb)
            else if (iopPQbb == 0) then
                call getpqrso(iCore(iabab), Sb, irrepSb, 1, file_bbbb)
                call iexpso(iCore(iabab), iCore(ibbbb), popPb, popQb, 0, 1, popRb, zorb, 1, 0, irrepSb)
            end if
            ! read V(pq,r*), ABAB
            call getpqrso(iCore(iabab), Sb, irrepSb, 1, file_abab)

            ! cycle for irrepRS
            do irrepRS = 1, nirrep
                irrepPQ = dirprd(irrepRS, 1)
                numPQ = irpdso(irrepPQ, popP, popQ, 1)

                irrepR = dirprd(irrepS, irrepRS)
                irrepRa = irrepR
                irrepRb = dirprd(irrepRa, symry)
                numR = popR(irrepR)
                numRa = popRa(irrepRa)
                numRb = popRb(irrepRb)
                if ((numRa + numRb) /= numR) then
                    write (6, *) 'numR error'
                    stop
                end if

                idx_ns(1) = ioff(irrepR, irrepS)
                ! cycle for irrepQ
                do irrepQ = 1, nirrep
                    irrepQa = irrepQ
                    irrepQb = dirprd(irrepQa, symry)
                    numQ = popQ(irrepQ)
                    numQa = popQa(irrepQa)
                    numQb = popQb(irrepQb)
                    if ((numQa + numQb) /= numQ) then
                        write (6, *) 'numQ error'
                        stop
                    end if

                    irrepP = dirprd(irrepQ, irrepPQ)
                    irrepPa = irrepP
                    irrepPb = dirprd(irrepPa, symry)
                    numP = popP(irrepP)
                    numPa = popPa(irrepPa)
                    numPb = popPb(irrepPb)
                    if ((numPa + numPb) /= numP) then
                        write (6, *) 'numP error'
                        stop
                    end if

                    idx_ns(2) = idx_ns(1) + isymoffso(irrepQ, irrepPQ, popP, popQ, 1) - 1
                    !------ for BBBB
                    irrepPQbb = dirprd(irrepPb, irrepQb)
                    numPQbb = irpdso(irrepPQbb, popPb, popQb, 1)
                    idx(1) = ioffbbbb(irrepRb, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepQb, irrepPQbb, popPb, popQb, 1) - 1
                    do Rb = 1, numRb
                        R = numRa + Rb
                        idx(3) = idx(2) + numPQbb*(Rb - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + numPb*(Qb - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + Pb
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(ibbbb - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for ABAB
                    irrepPQab = dirprd(irrepPa, irrepQb)
                    numPQab = irpdso(irrepPQab, popPa, popQb, 1)
                    idx(1) = ioffabab(irrepRa, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepQb, irrepPQab, popPa, popQb, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQab*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qb = 1, numQb
                            Q = numQa + Qb
                            idx(4) = idx(3) + numPa*(Qb - 1)
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pa = 1, numPa
                                P = Pa
                                idx(5) = idx(4) + Pa
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                    !------ for BAAB, <Pq||rS>=-<qP||rS>
                    irrepPQba = dirprd(irrepPb, irrepQa)
                    numPQba = irpdso(irrepPQba, popPb, popQa, 1)
                    idx(1) = ioffabab(irrepRa, irrepSb)
                    idx(2) = idx(1) + isymoffso(irrepPb, irrepPQba, popQa, popPb, 1) - 1
                    do Ra = 1, numRa
                        R = Ra
                        idx(3) = idx(2) + numPQba*(Ra - 1)
                        idx_ns(3) = idx_ns(2) + numPQ*(R - 1)
                        do Qa = 1, numQa
                            Q = Qa
                            idx(4) = idx(3) + Qa
                            idx_ns(4) = idx_ns(3) + numP*(Q - 1)
                            do Pb = 1, numPb
                                P = numPa + Pb
                                idx(5) = idx(4) + numQa*(Pb - 1)
                                idx_ns(5) = idx_ns(4) + P
                                iCore(ins - 1 + idx_ns(5)) = -iCore(iabab - 1 + idx(5))
                            end do
                        end do
                    end do
                end do
            end do
            if (iopPQ == 1) then
                ! write V(pq,r*)
                call putpqrso(iCore(ins), S, irrepS, 1, file_ns)
            else if (iopPQ == 0) then
                call icompso(iCore(ins), iCore(ibbbb), popP, popQ, 1, 0, popR, zorb, 1, 1, irrepS)
                call putpqrso(iCore(ibbbb), S, irrepS, 1, file_ns)
            end if
        end do
    end do

end subroutine resortVabcd
