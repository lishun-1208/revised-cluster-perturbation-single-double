subroutine run(iCore, maxCore)
    use mod_orbit
    use mod_symm
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: irrepab, irrepcd, irrepci
    integer*8 :: numab, numcd, numci
    integer*8 :: first, nread, nleft
    integer*8 :: i, j, nSize
    integer*8, dimension(8) :: idx
    integer*8 :: i0, i1
    real*8 :: norm

    integer*8, external :: idsymsz, irpdso
    real*8, external :: dnrm2

    call readmodule

    i0 = 1

    !--------------------------------------
    ! ijka, i,j,k... for alpha, I,J,K... for beta
    ! need: <ij||ka>, <iJ||kA>, <Ij||Ka>, <IJ||KA>
    ! <iJ||Ka>=-<Ji||Ka>, <Ij||kA>=-<jI||kA>
    !--------------------------------------
    ! ABAB
    call clmoioso0(1, iCore(i0), popa, popb, 1, popa, vrtb, 1, f_ijka(1))
    nSize = idsymsz(1, 14, 18)
    call getall(iCore(i0), nSize, 1, 10)
    call putallso(iCore(i0), 1, f_ijka(1))
    write (6, *) 'norm of ijkaABAB:', dnrm2(nSize, iCore(i0), 1)
    if (UHF) then
        ! AAAA
        call clmoioso0(1, iCore(i0), popa, popa, 0, popa, vrta, 1, f_ijka(2))
        nSize = idsymsz(1, 3, 16)
        call getall(iCore(i0), nSize, 1, 7)
        call putallso(iCore(i0), 1, f_ijka(2))
        write (6, *) 'norm of ijkaAAAA:', dnrm2(nSize, iCore(i0), 1)
        ! BABA
        call clmoioso0(1, iCore(i0), popb, popa, 1, popb, vrta, 1, f_ijka(4))
        nSize = idsymsz(1, 14, 11)
        call getall(iCore(i0), nSize, 1, 9) ! <iJ||aK>
        i1 = i0 + nSize
        call isymtrso(iCore(i0), iCore(i1), popa, popb, 1, 1, vrta, popb, 1, 1, 1)
        call putallso(iCore(i1), 1, f_ijka(4))
        write (6, *) 'norm of ijkaBABA:', dnrm2(nSize, iCore(i1), 1)
        ! BBBB
        call clmoioso0(1, iCore(i0), popb, popb, 0, popb, vrtb, 1, f_ijka(3))
        nSize = idsymsz(1, 4, 17)
        call getall(iCore(i0), nSize, 1, 8)
        call putallso(iCore(i0), 1, f_ijka(3))
        write (6, *) 'norm of ijkaBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if

    !--------------------------------------
    ! iabj sort as biaj
    ! need: <ia||bj>, <iA||bJ>, <iA||Bj>, <Ia||bJ>, <IA||BJ>
    ! <Ia||Bj>=<Bj||Ia>=<jB||aI>
    !--------------------------------------
    ! AAAA->AAAA
    call clmoioso0(1, iCore(i0), vrta, popa, 1, vrta, popa, 1, f_iabj(1))
    nSize = idsymsz(1, 9, 9)
    call getall(iCore(i0), nSize, 1, 23) ! <ia||jb> sort as biaj
    call dscal(nSize, -1.0d0, iCore(i0), 1) ! <ia||bj>=-<ia||jb>
    call putallso(iCore(i0), 1, f_iabj(1))
    write (6, *) 'norm of iabjAAAA:', dnrm2(nSize, iCore(i0), 1)
    ! ABAB->AABB
    call clmoioso0(1, iCore(i0), vrta, popa, 1, vrtb, popb, 1, f_iabj(3))
    nSize = idsymsz(1, 9, 10)
    call getall(iCore(i0), nSize, 1, 18)
    call putallso(iCore(i0), 1, f_iabj(3))
    write (6, *) 'norm of iabjABAB:', dnrm2(nSize, iCore(i0), 1)
    ! BAAB->ABAB
    call clmoioso0(1, iCore(i0), vrta, popb, 1, vrta, popb, 1, f_iabj(2))
    nSize = idsymsz(1, 11, 11)
    call getall(iCore(i0), nSize, 1, 25)
    call dscal(nSize, -1.0d0, iCore(i0), 1)
    call putallso(iCore(i0), 1, f_iabj(2))
    write (6, *) 'norm of iabjBAAB:', dnrm2(nSize, iCore(i0), 1)
    if (UHF) then
        ! ABBA->BABA
        call clmoioso0(1, iCore(i0), vrtb, popa, 1, vrtb, popa, 1, f_iabj(5))
        nSize = idsymsz(1, 12, 12)
        call getall(iCore(i0), nSize, 1, 26)
        call dscal(nSize, -1.0d0, iCore(i0), 1)
        call putallso(iCore(i0), 1, f_iabj(5))
        write (6, *) 'norm of iabjABBA:', dnrm2(nSize, iCore(i0), 1)
        ! BBBB->BBBB
        call clmoioso0(1, iCore(i0), vrtb, popb, 1, vrtb, popb, 1, f_iabj(4))
        nSize = idsymsz(1, 10, 10)
        call getall(iCore(i0), nSize, 1, 24)
        call dscal(nSize, -1.0d0, iCore(i0), 1)
        call putallso(iCore(i0), 1, f_iabj(4))
        write (6, *) 'norm of iabjBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if

    !--------------------------------------
    ! abci
    ! need: <ab||ci>, <aB||cI>, <Ab||Ci>, <AB||CI>
    ! <aB||Ci>=-<Ba||Ci>, <Ab||cI>=-<bA||cI>
    !--------------------------------------
    call clmoioso0(1, iCore(i0), vrta, vrtb, 1, vrta, popb, 1, f_abci(1))
    if (UHF) then
        call clmoioso0(1, iCore(i0), vrta, vrta, 0, vrta, popa, 1, f_abci(2))
        call clmoioso0(1, iCore(i0), vrtb, vrta, 1, vrtb, popa, 1, f_abci(4))
        call clmoioso0(1, iCore(i0), vrtb, vrtb, 0, vrtb, popb, 1, f_abci(3))
    end if
    do irrepci = 1, nirrep
        irrepab = dirprd(irrepci, 1)
        ! ABAB
        numab = irpdso(irrepab, vrta, vrtb, 1)
        numci = irpdso(irrepci, vrta, popb, 1)
        do i = 1, numci
            call getlst(iCore(i0), i, 1, 1, irrepci, 30)
            call putlstso(iCore(i0), i, 1, irrepci, 1, f_abci(1))
        end do
        if (UHF) then
            ! AAAA
            numab = irpdso(irrepab, vrta, vrta, 0)
            numci = irpdso(irrepci, vrta, popa, 1)
            do i = 1, numci
                call getlst(iCore(i0), i, 1, 1, irrepci, 27)
                call putlstso(iCore(i0), i, 1, irrepci, 1, f_abci(2))
            end do
            ! BABA
            numab = irpdso(irrepab, vrtb, vrta, 1)
            numci = irpdso(irrepci, vrtb, popa, 1)
            do i = 1, numci
                call getlst(iCore(i0), i, 1, 1, irrepci, 29) ! <aB||iC>
                i1 = i0 + numab
                ! aB->Ba
                call isymtrso(iCore(i0), iCore(i1), vrta, vrtb, 1, 1, zorb, zorb, 1, 0, irrepab)
                ! iC->Ci
                call getidx(i, idx, popa, vrtb, 1, zorb, zorb, 1, irrepci)
                call putidx(j, idx(2), idx(1), 1, 1, idx(6), idx(5), 1, 1, &
                            vrtb, popa, 1, zorb, zorb, 1, irrepci)
                call putlstso(iCore(i1), j, 1, irrepci, 1, f_abci(4))
            end do
            ! BBBB
            numab = irpdso(irrepab, vrtb, vrtb, 0)
            numci = irpdso(irrepci, vrtb, popb, 1)
            do i = 1, numci
                call getlst(iCore(i0), i, 1, 1, irrepci, 28)
                call putlstso(iCore(i0), i, 1, irrepci, 1, f_abci(3))
            end do
        end if
    end do
    call calLargeNorm(norm, f_abci(1), vrta, vrtb, 1, vrta, popb, 1, 1, iCore, maxCore)
    write (6, *) 'norm of abciABAB:', norm
    if (UHF) then
        call calLargeNorm(norm, f_abci(2), vrta, vrta, 0, vrta, popa, 1, 1, iCore, maxCore)
        write (6, *) 'norm of abciAAAA:', norm
        call calLargeNorm(norm, f_abci(4), vrtb, vrta, 1, vrtb, popa, 1, 1, iCore, maxCore)
        write (6, *) 'norm of abciBABA:', norm
        call calLargeNorm(norm, f_abci(3), vrtb, vrtb, 0, vrtb, popb, 1, 1, iCore, maxCore)
        write (6, *) 'norm of abciBBBB:', norm
    end if

    !--------------------------------------
    ! abij
    ! need: <ab||ij>, <aB||iJ>, <AB||IJ>
    ! <aB||Ij>=-<aB||jI>, <Ab||Ij>=<bA||jI>, <Ab||iJ>=-<bA||iJ>
    !--------------------------------------
    ! ABAB
    call clmoioso0(1, iCore(i0), vrta, vrtb, 1, popa, popb, 1, f_abij(1))
    nSize = idsymsz(1, 13, 14)
    call getall(iCore(i0), nSize, 1, 16)
    call putallso(iCore(i0), 1, f_abij(1))
    write (6, *) 'norm of abijABAB', dnrm2(nSize, iCore(i0), 1)
    if (UHF) then
        ! AAAA
        call clmoioso0(1, iCore(i0), vrta, vrta, 0, popa, popa, 0, f_abij(2))
        nSize = idsymsz(1, 1, 3)
        call getall(iCore(i0), nSize, 1, 14)
        call putallso(iCore(i0), 1, f_abij(2))
        write (6, *) 'norm of abijAAAA:', dnrm2(nSize, iCore(i0), 1)
        ! BBBB
        call clmoioso0(1, iCore(i0), vrtb, vrtb, 0, popb, popb, 0, f_abij(3))
        nSize = idsymsz(1, 2, 4)
        call getall(iCore(i0), nSize, 1, 15)
        call putallso(iCore(i0), 1, f_abij(3))
        write (6, *) 'norm of abijBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if

    !--------------------------------------
    ! ijkl
    ! need: <ij||kl>, <iJ||kL>, <IJ||Kl>
    ! <iJ||Kl>=-<iJ||lK>, <Ij||Kl>=<jI||lK>, <Ij||kL>=-<jI||kL>
    !--------------------------------------
    ! ABAB
    call clmoioso0(1, iCore(i0), popa, popb, 1, popa, popb, 1, f_ijkl(1))
    nSize = idsymsz(1, 14, 14)
    call getall(iCore(i0), nSize, 1, 13)
    call putallso(iCore(i0), 1, f_ijkl(1))
    write (6, *) 'norm of ijklABAB:', dnrm2(nSize, iCore(i0), 1)
    if (UHF) then
        ! AAAA
        call clmoioso0(1, iCore(i0), popa, popa, 0, popa, popa, 0, f_ijkl(2))
        nSize = idsymsz(1, 3, 3)
        call getall(iCore(i0), nSize, 1, 11)
        call putallso(iCore(i0), 1, f_ijkl(2))
        write (6, *) 'norm of ijklAAAA:', dnrm2(nSize, iCore(i0), 1)
        ! BBBB
        call clmoioso0(1, iCore(i0), popb, popb, 0, popb, popb, 0, f_ijkl(3))
        nSize = idsymsz(1, 4, 4)
        call getall(iCore(i0), nSize, 1, 12)
        call putallso(iCore(i0), 1, f_ijkl(3))
        write (6, *) 'norm of ijklBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if

    !--------------------------------------
    ! abcd
    ! need: <ab||cd>, <aB||cD>, <AB||CD>
    ! <aB||Cd>=-<aB||dC>, <Ab||Cd>=<bA||dC>, <Ab||cD>=-<bA||cD>
    !--------------------------------------
    call clmoioso0(1, iCore(i0), vrta, vrtb, 1, vrta, vrtb, 1, f_abcd(1))
    if (UHF) then
        call clmoioso0(1, iCore(i0), vrta, vrta, 0, vrta, vrta, 0, f_abcd(2))
        call clmoioso0(1, iCore(i0), vrtb, vrtb, 0, vrtb, vrtb, 0, f_abcd(3))
    end if
    do irrepcd = 1, nirrep
        irrepab = dirprd(irrepcd, 1)
        ! ABAB
        numab = irpdso(irrepab, vrta, vrtb, 1)
        numcd = irpdso(irrepcd, vrta, vrtb, 1)
        first = 1
        nleft = numcd
        nread = min(maxCore/numab, numcd)
        do
            call getlst(iCore(i0), first, nread, 1, irrepcd, 233)
            call putlstso(iCore(i0), first, nread, irrepcd, 1, f_abcd(1))
            first = first + nread
            nleft = nleft - nread
            if (nleft < nread) then
                if (nleft == 0) exit
                nread = nleft
            end if
        end do
        if (UHF) then
            ! AAAA
            numab = irpdso(irrepab, vrta, vrta, 0)
            numcd = irpdso(irrepcd, vrta, vrta, 0)
            first = 1
            nleft = numcd
            nread = min(maxCore/numab, numcd)
            do
                call getlst(iCore(i0), first, nread, 1, irrepcd, 231)
                call putlstso(iCore(i0), first, nread, irrepcd, 1, f_abcd(2))
                first = first + nread
                nleft = nleft - nread
                if (nleft < nread) then
                    if (nleft == 0) exit
                    nread = nleft
                end if
            end do
            ! BBBB
            numab = irpdso(irrepab, vrtb, vrtb, 0)
            numcd = irpdso(irrepcd, vrtb, vrtb, 0)
            first = 1
            nleft = numcd
            nread = min(maxCore/numab, numcd)
            do
                call getlst(iCore(i0), first, nread, 1, irrepcd, 232)
                call putlstso(iCore(i0), first, nread, irrepcd, 1, f_abcd(3))
                first = first + nread
                nleft = nleft - nread
                if (nleft < nread) then
                    if (nleft == 0) exit
                    nread = nleft
                end if
            end do
        end if
    end do
    call calLargeNorm(norm, f_abcd(1), vrta, vrtb, 1, vrta, vrtb, 1, 1, iCore, maxCore)
    write (6, *) 'norm of abcdABAB:', norm
    if (UHF) then
        call calLargeNorm(norm, f_abcd(2), vrta, vrta, 0, vrta, vrta, 0, 1, iCore, maxCore)
        write (6, *) 'norm of abcdAAAA:', norm
        call calLargeNorm(norm, f_abcd(3), vrtb, vrtb, 0, vrtb, vrtb, 0, 1, iCore, maxCore)
        write (6, *) 'norm of abcdBBBB:', norm
    end if

    write (6, *) 'finish 2e-integral'

    !---------------------------------
    ! f(i,j), not contain orbital energy, 0 for canical HF
    ! need: f(i,j), f(I,J)
    !---------------------------------
    ! AA
    call clmoioso0(1, iCore(i0), popa, popa, 1, zorb, zorb, 1, f_fij(1))
    nSize = irpdso(1, popa, popa, 1)
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call putallso(iCore(i0), 1, f_fij(1))
    write (6, *) 'norm of fijAA:', 0.0d0
    if (UHF) then
        ! BB
        call clmoioso0(1, iCore(i0), popb, popb, 1, zorb, zorb, 1, f_fij(2))
        nSize = irpdso(1, popb, popb, 1)
        call dscal(nSize, 0.0d0, iCore(i0), 1)
        call putallso(iCore(i0), 1, f_fij(2))
        write (6, *) 'norm of fijBB:', 0.0d0
    end if

    !---------------------------------
    ! f(a,b), not contain orbital energy, 0 for canical HF
    ! need: f(a,b), f(A,B)
    !---------------------------------
    ! AA
    call clmoioso0(1, iCore(i0), vrta, vrta, 1, zorb, zorb, 1, f_fab(1))
    nSize = irpdso(1, vrta, vrta, 1)
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call putallso(iCore(i0), 1, f_fab(1))
    write (6, *) 'norm of fabAA:', 0.0d0
    if (UHF) then
        ! BB
        call clmoioso0(1, iCore(i0), vrtb, vrtb, 1, zorb, zorb, 1, f_fab(2))
        nSize = irpdso(1, vrtb, vrtb, 1)
        call dscal(nSize, 0.0d0, iCore(i0), 1)
        call putallso(iCore(i0), 1, f_fab(2))
        write (6, *) 'norm of fabBB:', 0.0d0
    end if

    !---------------------------------
    ! f(i,a) sort as (a,i), 0 for canical HF
    ! need f(i,a), f(I,A)
    !---------------------------------
    ! AA
    call clmoioso0(1, iCore(i0), vrta, popa, 1, zorb, zorb, 1, f_fia(1))
    nSize = irpdso(1, vrta, popa, 1)
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call putallso(iCore(i0), 1, f_fia(1))
    write (6, *) 'norm of fiaAA:', 0.0d0
    if (UHF) then
        ! BB
        call clmoioso0(1, iCore(i0), vrtb, popb, 1, zorb, zorb, 1, f_fia(2))
        nSize = irpdso(1, vrtb, popb, 1)
        call dscal(nSize, 0.0d0, iCore(i0), 1)
        call putallso(iCore(i0), 1, f_fia(2))
        write (6, *) 'norm of fiaBB:', 0.0d0
    end if

    write (6, *) 'finish fock matrix'

end subroutine run
