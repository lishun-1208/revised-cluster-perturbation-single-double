! generate pure alpha 2e-integral when use RHF as reference
! 2025.12.15
subroutine generate_pure_alpha(iCore)
    use mod_orbit
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore

    integer*8 :: i0, i1, i2, nSize
    integer*8 :: idsymsoc, irpdso

    i0 = 1

    !--- ijkl
    call clmoioso0(1, iCore(i0), popa, popa, iop, popa, popa, iop, file_ijklAAAA)
    nSize = idsymsoc(1, popa, popb, 1, popa, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijklABAB)
    nSize = idsymsoc(1, popa, popa, iop, popa, popa, iop)
    i2 = i1 + nSize
    call p_abab2aaaa(iCore(i0), iCore(i1), popa, popa, iop, popa, popa, iop, 1, iCore(i2), 'pq')
    call putallso(iCore(i1), 1, file_ijklAAAA)

    !--- ijka
    call clmoioso0(1, iCore(i0), popa, popa, iop, popa, vrta, 1, file_ijkaAAAA)
    nSize = idsymsoc(1, popa, popb, 1, popa, vrtb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijkaABAB)
    nSize = idsymsoc(1, popa, popa, iop, popa, vrta, 1)
    i2 = i1 + nSize
    call p_abab2aaaa(iCore(i0), iCore(i1), popa, popa, iop, popa, vrta, 1, 1, iCore(i2), 'pq')
    call putallso(iCore(i1), 1, file_ijkaAAAA)

    !--- abij
    call clmoioso0(1, iCore(i0), vrta, vrta, iop, popa, popa, iop, file_abijAAAA)
    nSize = idsymsoc(1, vrta, vrtb, 1, popa, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_abijABAB)
    nSize = idsymsoc(1, vrta, vrta, iop, popa, popa, iop)
    i2 = i1 + nSize
    call p_abab2aaaa(iCore(i0), iCore(i1), vrta, vrta, iop, popa, popa, iop, 1, iCore(i2), 'pq')
    call putallso(iCore(i1), 1, file_abijAAAA)

    !--- abci
    if (vvvo_incore) then
        call clmoioso0(1, iCore(i0), vrta, vrta, iop, vrta, popa, 1, file_abciAAAA)
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, popb, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abciABAB)
        nSize = idsymsoc(1, vrta, vrta, iop, vrta, popa, 1)
        i2 = i1 + nSize
        call p_abab2aaaa(iCore(i0), iCore(i1), vrta, vrta, iop, vrta, popa, 1, 1, iCore(i2), 'pq')
        call putallso(iCore(i1), 1, file_abciAAAA)
    else
        call p_abab2aaaa_h(file_abciABAB, file_abciAAAA, vrta, vrta, iop, vrta, popa, 1, 1, iCore(i0), 'pq')
    end if

    !--- abcd
    if (vvvv_incore) then
        call clmoioso0(1, iCore(i0), vrta, vrta, iop, vrta, vrta, iop, file_abcdAAAA)
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, vrtb, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abcdABAB)
        nSize = idsymsoc(1, vrta, vrta, iop, vrta, vrta, iop)
        i2 = i1 + nSize
        call p_abab2aaaa(iCore(i0), iCore(i1), vrta, vrta, iop, vrta, vrta, iop, 1, iCore(i2), 'pq')
        call putallso(iCore(i1), 1, file_abcdAAAA)
    else
        call p_abab2aaaa_h(file_abcdABAB, file_abcdAAAA, vrta, vrta, iop, vrta, vrta, iop, 1, iCore(i0), 'pq')
    end if

end subroutine generate_pure_alpha
