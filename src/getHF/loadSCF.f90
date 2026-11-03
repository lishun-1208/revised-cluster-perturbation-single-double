! load two-electron integral from HF calculation, and write in file
!------------------------------------
! modified at 2024.12.5
subroutine loadSCF(iCore, istart)
    use mod_orbit
    use mod_symm
    use mod_interface
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0, i1, i2
    integer*8 :: i, iEnd, nSize
    real*8 :: norm

    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2 ! provide by BLAS

    ! pointor for 2e-integral with different spin
    integer*8 :: iaaaa, ibbbb, iabab, ibaba
    integer*8 :: iaabb

    i0 = istart
    !------------- read in orbital energy
    i1 = i0 + noa + nva
    call read_node('orb_energy_alpha_sym', iCore(i0:i1 - 1), noa + nva)
    i2 = i1 + nob + nvb
    call read_node('orb_energy_beta_sym', iCore(i1:i2 - 1), nob + nvb)
    ! resort Eorb from (oa,va),(ob,vb) to (o,v)
    call resortEorb(iCore(iEorb), iCore(i0), iCore(i1))
    call printEorb_ns(iCore(iEorb))
    nSize = nVir + nOcc
    call add_node('orb_energy_no_spin', iCore(iEorb:iEorb + nSize - 1), nSize)

    !---------- read in Fock matrix
    !------ fij
    call clmoioso0(1, iCore(i0), popo, popo, 1, zorb, zorb, 1, file_fij)
    nSize = irpdso(1, popo, popo, 1)
    call dscal(nSize, 0.0d0, iCore(fij%idx), 1)
    ! AA
    nSize = irpdso(1, popa, popa, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_fijAA)
    write (6, *) 'norm of fijAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(fij%idx), iCore(i0), popa, popb, popa, popb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popo, popo, 1, zorb, zorb, 1, 1, 'AAAA')
    ! BB
    nSize = irpdso(1, popb, popb, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_fijBB)
        write (6, *) 'norm of fijBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(fij%idx), iCore(i0), popa, popb, popa, popb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popo, popo, 1, zorb, zorb, 1, 1, 'BBAA')
    nSize = irpdso(1, popo, popo, 1)
    norm = dnrm2(nSize, iCore(fij%idx), 1)
    write (6, *) "norm of fij =", norm
    call putallso(iCore(fij%idx), 1, file_fij)

    !------ fab
    call clmoioso0(1, iCore(i0), popv, popv, 1, zorb, zorb, 1, file_fab)
    nSize = irpdso(1, popv, popv, 1)
    call dscal(nSize, 0.0d0, iCore(fab%idx), 1)
    ! AA
    nSize = irpdso(1, vrta, vrta, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_fabAA)
    write (6, *) 'norm of fabAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(fab%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popv, popv, 1, zorb, zorb, 1, 1, 'AAAA')
    ! BB
    nSize = irpdso(1, vrtb, vrtb, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_fabBB)
        write (6, *) 'norm of fabBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(fab%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popv, popv, 1, zorb, zorb, 1, 1, 'BBAA')
    nSize = irpdso(1, popv, popv, 1)
    norm = dnrm2(nSize, iCore(fab%idx), 1)
    write (6, *) "norm of fab =", norm
    call putallso(iCore(fab%idx), 1, file_fab)

    !------ fia, sort as (a,i)
    call clmoioso0(1, iCore(i0), popv, popo, 1, zorb, zorb, 1, file_fia)
    nSize = irpdso(1, popv, popo, 1)
    call dscal(nSize, 0.0d0, iCore(fia%idx), 1)
    ! AA
    nSize = irpdso(1, vrta, popa, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_fiaAA)
    write (6, *) 'norm of fiaAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(fia%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popv, popo, 1, zorb, zorb, 1, 1, 'AAAA')
    ! BB
    nSize = irpdso(1, vrtb, popb, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_fiaBB)
        write (6, *) 'norm of fiaBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(fia%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      zorb, zorb, zorb, zorb, 1, 1, &
                      popv, popo, 1, zorb, zorb, 1, 1, 'BBAA')
    nSize = irpdso(1, popv, popo, 1)
    norm = dnrm2(nSize, iCore(fia%idx), 1)
    write (6, *) "norm of fia =", norm
    call putallso(iCore(fia%idx), 1, file_fia)

    !---------- read in 2e-integral
    !--------------------------------------------------------------------------------------------------------
    ! we store the 2e-integral as <pq||rs>=<pq|rs>-<pq|sr>, with six different spin condiations
    ! <aa||aa>, <bb||bb>, <ab||ab>=<ba||ba>=-<ab||ba>=-<ba||ab>=<ab|ab>
    !--------------------------------------------------------------------------------------------------------

    !--------------- ijkl, with aaaa,bbbb,abab
    call clmoioso0(1, iCore(i0), popo, popo, 0, popo, popo, 0, file_ijkl)
    nSize = idsymsoc(1, popo, popo, 0, popo, popo, 0)
    call dscal(nSize, 0.0d0, iCore(ijkl%idx), 1)
    ! get aaaa
    nSize = idsymsoc(1, popa, popa, 0, popa, popa, 0)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijklAAAA) ! i<j;k<l
    write (6, *) 'norm of ijklAAAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(ijkl%idx), iCore(i0), popa, popb, popa, popb, 0, &
                      popa, popb, popa, popb, 0, 1, &
                      popo, popo, 0, popo, popo, 0, 1, 'AAAA')
    ! get bbbb
    nSize = idsymsoc(1, popb, popb, 0, popb, popb, 0)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_ijklBBBB) ! i<j;k<l
        write (6, *) 'norm of ijklBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(ijkl%idx), iCore(i0), popa, popb, popa, popb, 0, &
                      popa, popb, popa, popb, 0, 1, &
                      popo, popo, 0, popo, popo, 0, 1, 'BBBB')
    ! get abab
    nSize = idsymsoc(1, popa, popb, 1, popa, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijklABAB) ! i,j;k,l
    write (6, *) 'norm of ijklABAB:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(ijkl%idx), iCore(i0), popa, popb, popa, popb, 1, &
                      popa, popb, popa, popb, 1, 1, &
                      popo, popo, 0, popo, popo, 0, 1, 'ABAB')

    nSize = idsymsoc(1, popo, popo, 0, popo, popo, 0)
    norm = dnrm2(nSize, iCore(ijkl%idx), 1)
    write (6, '(A,I2,A,F,I12)') "norm of ijkl with iop", 0, ' is:', norm, nSize
    call putallso(iCore(ijkl%idx), 1, file_ijkl)

    !---------------- ijka, with aaaa,bbbb,abab,baba
    call clmoioso0(1, iCore(i0), popo, popo, 0, popo, popv, 1, file_ijka)
    nSize = idsymsoc(1, popo, popo, 0, popo, popv, 1)
    call dscal(nSize, 0.0d0, iCore(ijka%idx), 1)
    ! get aaaa
    nSize = idsymsoc(1, popa, popa, 0, popa, vrta, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijkaAAAA) ! I<J;K,A
    write (6, *) 'norm of ijkaAAAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(ijka%idx), iCore(i0), popa, popb, popa, popb, 0, &
                      popa, popb, vrta, vrtb, 1, 1, &
                      popo, popo, 0, popo, popv, 1, 1, 'AAAA')
    ! get bbbb
    nSize = idsymsoc(1, popb, popb, 0, popb, vrtb, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_ijkaBBBB) ! i<j;k,a
        write (6, *) 'norm of ijkaBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(ijka%idx), iCore(i0), popa, popb, popa, popb, 0, &
                      popa, popb, vrta, vrtb, 1, 1, &
                      popo, popo, 0, popo, popv, 1, 1, 'BBBB')
    ! get abab
    nSize = idsymsoc(1, popa, popb, 1, popa, vrtb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_ijkaABAB) ! I,j;K,a
    write (6, *) 'norm of ijkaABAB:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(ijka%idx), iCore(i0), popa, popb, popa, popb, 1, &
                      popa, popb, vrta, vrtb, 1, 1, &
                      popo, popo, 0, popo, popv, 1, 1, 'ABAB')
    ! get baba
    nSize = idsymsoc(1, popb, popa, 1, popb, vrta, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_ijkaBABA) ! i,j;k,a
        write (6, *) 'norm of ijkaBABA:', dnrm2(nSize, iCore(i0), 1)
    end if
    ! <Ij||Ka>=-<jI||Ka>
    i2 = i1 + nSize
    call isymtrso(iCore(i0), iCore(i1), popb, popa, 1, 1, popb, vrta, 1, 0, 1)
    call dscal(nSize, -1.0d0, iCore(i1), 1)
    call combintegral(iCore(ijka%idx), iCore(i1), popa, popb, popa, popb, 1, &
                      popa, popb, vrta, vrtb, 1, 1, &
                      popo, popo, 0, popo, popv, 1, 1, 'ABBA')

    nSize = idsymsoc(1, popo, popo, 0, popo, popv, 1)
    norm = dnrm2(nSize, iCore(ijka%idx), 1)
    write (6, '(A,I2,A,F,I12)') "norm of ijka with iop", 0, ' is:', norm, nSize
    call putallso(iCore(ijka%idx), 1, file_ijka)

    !------------------ abij, with aaaa,bbbb,abab
    call clmoioso0(1, iCore(i0), popv, popv, 0, popo, popo, 0, file_abij)
    nSize = idsymsoc(1, popv, popv, 0, popo, popo, 0)
    call dscal(nSize, 0.0d0, iCore(abij%idx), 1)
    ! get aaaa
    nSize = idsymsoc(1, vrta, vrta, 0, popa, popa, 0)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_abijAAAA) ! A<B;I<J
    write (6, *) 'norm of abijAAAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(abij%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                      popa, popb, popa, popb, 0, 1, &
                      popv, popv, 0, popo, popo, 0, 1, 'AAAA')
    ! get bbbb
    nSize = idsymsoc(1, vrtb, vrtb, 0, popb, popb, 0)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_abijBBBB) ! A<B;I<J
        write (6, *) 'norm of abijBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(abij%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                      popa, popb, popa, popb, 0, 1, &
                      popv, popv, 0, popo, popo, 0, 1, 'BBBB')
    ! get abab
    nSize = idsymsoc(1, vrta, vrtb, 1, popa, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_abijABAB) ! A,b;I,j
    write (6, *) 'norm of abijABAB:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(abij%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 1, &
                      popa, popb, popa, popb, 1, 1, &
                      popv, popv, 0, popo, popo, 0, 1, 'ABAB')

    nSize = idsymsoc(1, popv, popv, 0, popo, popo, 0)
    norm = dnrm2(nSize, iCore(abij%idx), 1)
    write (6, '(A,I2,A,F,I12)') "norm of abij with iop", 0, ' is:', norm, nSize
    call putallso(iCore(abij%idx), 1, file_abij)

    !------------------- <ia||bj> would be sort as (bi,aj)
    !------------------- with aaaa,bbbb,abab,baba,aabb for (bi,aj)
    call clmoioso0(1, iCore(i0), popv, popo, 1, popv, popo, 1, file_iabj)
    nSize = idsymsoc(1, popv, popo, 1, popv, popo, 1)
    call dscal(nSize, 0.0d0, iCore(iabj%idx), 1)
    ! biaj
    nSize = idsymsoc(1, vrta, popa, 1, vrta, popa, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_iabjAAAA)
    write (6, *) 'norm of iabjAAAA:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(iabj%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'AAAA')
    ! BIAJ
    nSize = idsymsoc(1, vrtb, popb, 1, vrtb, popb, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_iabjBBBB)
        write (6, *) 'norm of iabjBBBB:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(iabj%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'BBBB')
    ! bIaJ=IabJ
    nSize = idsymsoc(1, vrta, popb, 1, vrta, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_iabjBAAB)
    write (6, *) 'norm of iabjBAAB:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(iabj%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'ABAB')
    ! BiAj=iABj
    nSize = idsymsoc(1, vrtb, popa, 1, vrtb, popa, 1)
    i1 = i0 + nSize
    if (.not. RHF) then
        call getallso(iCore(i0), 1, file_iabjABBA)
        write (6, *) 'norm of iabjABBA:', dnrm2(nSize, iCore(i0), 1)
    end if
    call combintegral(iCore(iabj%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'BABA')
    ! biAJ=iAbJ
    nSize = idsymsoc(1, vrta, popa, 1, vrtb, popb, 1)
    i1 = i0 + nSize
    call getallso(iCore(i0), 1, file_iabjABAB)
    write (6, *) 'norm of iabjABAB:', dnrm2(nSize, iCore(i0), 1)
    call combintegral(iCore(iabj%idx), iCore(i0), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'AABB')

    ! BIaj=IaBj->ajBI=jBaI
    i2 = i1 + nSize
    call itranspso(iCore(i0), iCore(i1), vrta, popa, 1, vrtb, popb, 1, 1)
    call combintegral(iCore(iabj%idx), iCore(i1), vrta, vrtb, popa, popb, 1, &
                      vrta, vrtb, popa, popb, 1, 1, &
                      popv, popo, 1, popv, popo, 1, 1, 'BBAA')

    nSize = idsymsoc(1, popv, popo, 1, popv, popo, 1)
    norm = dnrm2(nSize, iCore(iabj%idx), 1)
    write (6, *) "norm of iabj:", norm, nSize
    call putallso(iCore(iabj%idx), 1, file_iabj)

    !----------------- abci, with aaaa,bbbb,abab,baba
    if (abci%incore) then
        call clmoioso0(1, iCore(i0), popv, popv, 0, popv, popo, 1, file_abci)
        nSize = idsymsoc(1, popv, popv, 0, popv, popo, 1)
        call dscal(nSize, 0.0d0, iCore(abci%idx), 1)
        ! get aaaa
        nSize = idsymsoc(1, vrta, vrta, 0, vrta, popa, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abciAAAA) ! A<B;C,I
        write (6, *) 'norm of abciAAAA:', dnrm2(nSize, iCore(i0), 1)
        call combintegral(iCore(abci%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                          vrta, vrtb, popa, popb, 1, 1, &
                          popv, popv, 0, popv, popo, 1, 1, 'AAAA')
        ! get bbbb
        nSize = idsymsoc(1, vrtb, vrtb, 0, vrtb, popb, 1)
        i1 = i0 + nSize
        if (.not. RHF) then
            call getallso(iCore(i0), 1, file_abciBBBB) ! A<B;C,I
            write (6, *) 'norm of abciBBBB:', dnrm2(nSize, iCore(i0), 1)
        end if
        call combintegral(iCore(abci%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                          vrta, vrtb, popa, popb, 1, 1, &
                          popv, popv, 0, popv, popo, 1, 1, 'BBBB')
        ! get abab
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, popb, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abciABAB) ! A,b;C,i
        write (6, *) 'norm of abciABAB:', dnrm2(nSize, iCore(i0), 1)
        call combintegral(iCore(abci%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 1, &
                          vrta, vrtb, popa, popb, 1, 1, &
                          popv, popv, 0, popv, popo, 1, 1, 'ABAB')
        ! get baba
        nSize = idsymsoc(1, vrtb, vrta, 1, vrtb, popa, 1)
        i1 = i0 + nSize
        if (.not. RHF) then
            call getallso(iCore(i0), 1, file_abciBABA) ! a,b;c,i
            write (6, *) 'norm of abciBABA:', dnrm2(nSize, iCore(i0), 1)
        end if
        ! <Ab||Ci>=-<bA||Ci>
        i2 = i1 + nSize
        call isymtrso(iCore(i0), iCore(i1), vrtb, vrta, 1, 1, vrtb, popa, 1, 0, 1)
        call dscal(nSize, -1.0d0, iCore(i1), 1)
        call combintegral(iCore(abci%idx), iCore(i1), vrta, vrtb, vrta, vrtb, 1, &
                          vrta, vrtb, popa, popb, 1, 1, &
                          popv, popv, 0, popv, popo, 1, 1, 'ABBA')

        nSize = idsymsoc(1, popv, popv, 0, popv, popo, 1)
        norm = dnrm2(nSize, iCore(abci%idx), 1)
        write (6, '(A,I2,A,F,I12)') "norm of abci with iop", 0, ' is:', norm, nSize
        call putallso(iCore(abci%idx), 1, file_abci)
    else
        call clmoioso0(1, iCore(i0), popv, popv, 0, popv, popo, 1, file_abci)
        if (RHF) then
            call resortVabci(file_abci, file_abciAAAA, file_abciAAAA, file_abciABAB, file_abciABAB, &
                             vrta, vrta, vrta, popa, vrtb, vrtb, vrtb, popb, popv, popv, popv, popo, &
                             0, 1, 0, 1, 0, 1, iCore(i0), nCore - i0 + 1)
        else
            call resortVabci(file_abci, file_abciAAAA, file_abciBBBB, file_abciABAB, file_abciBABA, &
                             vrta, vrta, vrta, popa, vrtb, vrtb, vrtb, popb, popv, popv, popv, popo, &
                             0, 1, 0, 1, 0, 1, iCore(i0), nCore - i0 + 1)
        end if
        call calLargeNorm(norm, file_abciAAAA, vrta, vrta, 0, vrta, popa, 1, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, vrta, vrta, 0, vrta, popa, 1)
        write (6, *) "norm of abciAAAA:", norm, nSize
        if (.not. RHF) then
            call calLargeNorm(norm, file_abciBBBB, vrtb, vrtb, 0, vrtb, popb, 1, 1, iCore(i0), nCore - i0 + 1)
            nSize = idsymsoc(1, vrtb, vrtb, 0, vrtb, popb, 1)
            write (6, *) "norm of abciBBBB:", norm, nSize
        end if
        call calLargeNorm(norm, file_abciABAB, vrta, vrtb, 1, vrta, popb, 1, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, popb, 1)
        write (6, *) "norm of abciABAB:", norm, nSize
        if (.not. RHF) then
            call calLargeNorm(norm, file_abciBABA, vrtb, vrta, 1, vrtb, popa, 1, 1, iCore(i0), nCore - i0 + 1)
            nSize = idsymsoc(1, vrtb, vrta, 1, vrtb, popa, 1)
            write (6, *) "norm of abciBABA:", norm, nSize
        end if
        call calLargeNorm(norm, file_abci, popv, popv, 0, popv, popo, 1, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, popv, popv, 0, popv, popo, 1)
        write (6, '(A,I2,A,F,I12)') "norm of abci with iop", 0, ' is:', norm, nSize
    end if

    !--------------- abcd, with aaaa,bbbb,abab
    if (abcd%incore) then
        call clmoioso0(1, iCore(i0), popv, popv, 0, popv, popv, 0, file_abcd)
        nSize = idsymsoc(1, popv, popv, 0, popv, popv, 0)
        call dscal(nSize, 0.0d0, iCore(abcd%idx), 1)
        ! get aaaa
        nSize = idsymsoc(1, vrta, vrta, 0, vrta, vrta, 0)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abcdAAAA) ! A<B;C<D
        write (6, *) 'norm of abcdAAAA:', dnrm2(nSize, iCore(i0), 1)
        call combintegral(iCore(abcd%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                          vrta, vrtb, vrta, vrtb, 0, 1, &
                          popv, popv, 0, popv, popv, 0, 1, 'AAAA')
        ! get bbbb
        nSize = idsymsoc(1, vrtb, vrtb, 0, vrtb, vrtb, 0)
        i1 = i0 + nSize
        if (.not. RHF) then
            call getallso(iCore(i0), 1, file_abcdBBBB) ! a<b;c<d
            write (6, *) 'norm of abcdBBBB:', dnrm2(nSize, iCore(i0), 1)
        end if
        call combintegral(iCore(abcd%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 0, &
                          vrta, vrtb, vrta, vrtb, 0, 1, &
                          popv, popv, 0, popv, popv, 0, 1, 'BBBB')
        ! get abab
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, vrtb, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, file_abcdABAB) ! A,b;C,d
        write (6, *) 'norm of abcdABAB:', dnrm2(nSize, iCore(i0), 1)
        call combintegral(iCore(abcd%idx), iCore(i0), vrta, vrtb, vrta, vrtb, 1, &
                          vrta, vrtb, vrta, vrtb, 1, 1, &
                          popv, popv, 0, popv, popv, 0, 1, 'ABAB')

        nSize = idsymsoc(1, popv, popv, 0, popv, popv, 0)
        norm = dnrm2(nSize, iCore(abcd%idx), 1)
        write (6, '(A,I2,A,F,I12)') "norm of abcd with iop", 0, ' is:', norm, nSize
        call putallso(iCore(abcd%idx), 1, file_abcd)
    else
        call clmoioso0(1, iCore(i0), popv, popv, 0, popv, popv, 0, file_abcd)
        if (RHF) then
            call resortVabcd(file_abcd, file_abcdAAAA, file_abcdAAAA, file_abcdABAB, &
                             vrta, vrta, vrta, vrta, vrtb, vrtb, vrtb, vrtb, popv, popv, popv, popv, &
                             0, 0, 0, 0, 0, 0, iCore(i0), nCore - i0 + 1)
        else
            call resortVabcd(file_abcd, file_abcdAAAA, file_abcdBBBB, file_abcdABAB, &
                             vrta, vrta, vrta, vrta, vrtb, vrtb, vrtb, vrtb, popv, popv, popv, popv, &
                             0, 0, 0, 0, 0, 0, iCore(i0), nCore - i0 + 1)
        end if
        call calLargeNorm(norm, file_abcdAAAA, vrta, vrta, 0, vrta, vrta, 0, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, vrta, vrta, 0, vrta, vrta, 0)
        write (6, *) "norm of abcdAAAA:", norm, nSize
        if (.not. RHF) then
            call calLargeNorm(norm, file_abcdBBBB, vrtb, vrtb, 0, vrtb, vrtb, 0, 1, iCore(i0), nCore - i0 + 1)
            nSize = idsymsoc(1, vrtb, vrtb, 0, vrtb, vrtb, 0)
            write (6, *) "norm of abcdBBBB:", norm, nSize
        end if
        call calLargeNorm(norm, file_abcdABAB, vrta, vrtb, 1, vrta, vrtb, 1, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, vrta, vrtb, 1, vrta, vrtb, 1)
        write (6, *) "norm of abcdABAB:", norm, nSize
        call calLargeNorm(norm, file_abcd, popv, popv, 0, popv, popv, 0, 1, iCore(i0), nCore - i0 + 1)
        nSize = idsymsoc(1, popv, popv, 0, popv, popv, 0)
        write (6, '(A,I2,A,F,I12)') "norm of abcd with iop", 0, ' is:', norm, nSize
    end if

end subroutine loadSCF
