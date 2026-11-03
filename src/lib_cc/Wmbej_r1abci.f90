! 2026.7.6
!> Wmbej(em,bj) += fact*r(f,j)*V(fe,bm)
subroutine Wmbej_r1abci(Wmbej, r1, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wmbej
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abci%incore) then
        call vv2v_424(abci, r1, Wmbej, fact, 'febm', 'fj', 'embj', iCore, maxCore, 0)
    else
        write (6, *) 'ERROR: Wmbej_r1abci current not support outcore.'
        stop
        call p_Wmbej_r1abci_h(Wmbej%arr, r1%arr, abci%f, Wmbej%popq, Wmbej%popr, Wmbej%popp, Wmbej%pops, &
                              r1%popp, abci%ioppq, r1%sym, abci%sym, iCore, maxCore, fact, F)
    end if

end subroutine Wmbej_r1abci

! 2026.7.6
!> Wmbej(em,bj) += fact*r(f,j)*V(fe,bm)
subroutine Wmbej_r1abci_R(WmbejAAAA, WmbejABAB, WmbejBAAB, spin, r1AA, abciAAAA, abciABAB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB
    character(len=3) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abciAAAA, abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) += r(f,j)*V(fe,bm)
        call Wmbej_r1abci(WmbejAAAA, r1AA, abciAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) += r(F,J)*V(Fe,Bm) -> r(f,j)*V(fE,bM)
        call Wmbej_r1abci(WmbejABAB, r1AA, abciABAB, fact, iCore, maxCore)
    end if

    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) += r(F,J)*V(Fe,bM) -> -r(f,j)*V(eF,bM)
        if (abciABAB%incore) then
            call vv2v_424(abciABAB, r1AA, WmbejBAAB, -fact, 'efbm', 'fj', 'embj', iCore, maxCore, 0)
        else
            write (6, *) 'ERROR: Wmbej_r1abci_R current not support outcore.'
            stop
        end if
    end if

end subroutine Wmbej_r1abci_R

! 2026.7.6
!> Wmbej(em,bj) += fact*r(f,j)*V(fe,bm)
subroutine Wmbej_r1abci_U(WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA, spin, &
                          r1AA, r1BB, abciAAAA, abciABAB, abciBBBB, abciBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA
    character(len=6) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) += r(f,j)*V(fe,bm)
        call Wmbej_r1abci(WmbejAAAA, r1AA, abciAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) += r(F,J)*V(Fe,Bm)
        call Wmbej_r1abci(WmbejABAB, r1BB, abciBABA, fact, iCore, maxCore)
    end if

    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) += r(F,J)*V(Fe,bM) = -r(F,J)*V(eF,bM)
        if (abciABAB%incore) then
            call vv2v_424(abciABAB, r1BB, WmbejBAAB, -fact, 'eFbM', 'FJ', 'eMbJ', iCore, maxCore, 0)
        else
            write (6, *) 'ERROR: Wmbej_r1abci_U current not support outcore.'
            stop
        end if
    end if

    ! BBBB
    if (spin(4:4) == 'T') then
        ! WMBEJ(EM,BJ) += r(F,J)*V(FE,BM)
        call Wmbej_r1abci(WmbejBBBB, r1BB, abciBBBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(5:5) == 'T') then
        ! WMbEj(EM,bj) += r(f,j)*V(fE,bM)
        call Wmbej_r1abci(WmbejBABA, r1AA, abciABAB, fact, iCore, maxCore)
    end if

    ! ABBA
    if (spin(6:6) == 'T') then
        ! WmBEj(Em,Bj) += r(f,j)*V(fE,Bm) = -r(f,j)*V(Ef,Bm)
        if (abciBABA%incore) then
            call vv2v_424(abciBABA, r1AA, WmbejABBA, -fact, 'EfBm', 'fj', 'EmBj', iCore, maxCore, 0)
        else
            write (6, *) 'ERROR: Wmbej_r1abci_U current not support outcore.'
            stop
        end if
    end if

end subroutine Wmbej_r1abci_U

! modified at 2025.10.3
subroutine p_Wmbej_r1abci_h(Wmbej, r1, file_abci, pop_m, pop_b, pop_e, pop_j, pop_f, &
                            iop, symR, symV, iCore, maxCore, fact, exchange)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: Wmbej
    real*8, dimension(*), intent(in) :: r1
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_m, pop_b, pop_e, pop_j, pop_f
    integer*8, intent(in) :: iop ! iopef=iop
    integer*8, intent(in) :: symR, symV
    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact
    logical, intent(in) :: exchange !> whether exchange ef in <mB||Ef>

    target :: iCore, Wmbej
    real*8, dimension(:), pointer :: ptr
    real*8 :: fact2
    integer*8 :: symW
    integer*8 :: irrepm, numm, m
    integer*8 :: sym1, sym2
    integer*8 :: i0, i1, i2, i3, i4, i5, iEnd
    integer*8 :: nSize
    integer*8, external :: idsymsoc

    symW = dirprd(symR, symV)
    i0 = 1
    iEnd = i0

    !--- W(em,bj)->W(bj,em), i0
    nSize = idsymsoc(symW, pop_b, pop_j, 1, pop_e, pop_m, 1)
    i1 = i0 + nSize
    iEnd = max(iEnd, i1)
    call itranspso(Wmbej, iCore(i0), pop_e, pop_m, 1, pop_b, pop_j, 1, symW)
    ptr => iCore(i0:i0)

    do irrepm = 1, nirrep
        numm = pop_m(irrepm)
        sym1 = dirprd(irrepm, symW) ! symmetry of W(bj,e*)
        sym2 = dirprd(irrepm, symV) ! symmetry of V(fe,b*)
        do m = 1, numm
            !--- get W(bj,e*), i2
            nSize = idsymsoc(sym1, pop_b, pop_j, 1, pop_e, zorb, 1)
            i2 = i1 + nSize
            i3 = i2 + nSize
            iEnd = max(iEnd, i3)
            call getpqrso_incore(iCore(i2), m, irrepm, pop_b, pop_j, 1, &
                                 pop_e, pop_m, 1, ptr, symW)
            !--- W(bj,e*)->W(e*,bj), i1
            call itranspso(iCore(i2), iCore(i1), pop_b, pop_j, 1, pop_e, zorb, 1, sym1)

            !--- read in abci(fe,b*), i3
            nSize = idsymsoc(sym2, pop_f, pop_e, 1, pop_b, zorb, 1)
            i3 = i2 + nSize
            i4 = i3 + nSize
            iEnd = max(iEnd, i4)
            if (iop == 1) then
                call getpqrso(iCore(i3), m, irrepm, 1, file_abci)
            else if (iop == 0) then
                nSize = idsymsoc(sym2, pop_f, pop_e, 0, pop_b, zorb, 1)
                i5 = i4 + nSize
                iEnd = max(iEnd, i5)
                call getpqrso(iCore(i4), m, irrepm, 1, file_abci)
                call iexpso(iCore(i4), iCore(i3), pop_f, pop_e, 0, 1, pop_b, zorb, 1, 0, sym2)
            end if

            !--- abci(fe,b*)->abci(e*,bf), i2
            nSize = idsymsoc(sym2, pop_f, pop_e, 1, pop_b, zorb, 1)
            if (exchange) then
                fact2 = -1.0d0
                ! V(Ef,B*)->V(E*,Bf)
                call sstgenso(iCore(i3), iCore(i2), nSize, pop_e, pop_f, pop_b, zorb, &
                              iCore(iEnd), sym2, '1432')
            else
                fact2 = 1.0d0
                ! V(fe,b*)->V(e*,bf)
                call sstgenso(iCore(i3), iCore(i2), nSize, pop_f, pop_e, pop_b, zorb, &
                              iCore(i5), irrepm, '2431')
            end if

            !--- Wmbej(e*,bj) += abci(e*,bf)*t(f,j), i1
            call VpqrmVms_to_Vpqrs(iCore(i2), r1, sym2, symR, pop_f, &
                                   pop_e, zorb, 1, pop_b, pop_j, iCore(i1), sym1, fact*fact2, 1.0d0)

            !--- W(e*,bj)->W(bj,e*), i2
            call itranspso(iCore(i1), iCore(i2), pop_e, zorb, 1, pop_b, pop_j, 1, sym1)
            ! put W(bj,e*)
            call putpqrso_incore(iCore(i2), m, irrepm, pop_b, pop_j, 1, &
                                 pop_e, pop_m, 1, ptr, symR)
        end do
    end do

    ! W(bj,em)->W(em,bj)
    call itranspso(ptr, Wmbej, pop_b, pop_j, 1, pop_e, pop_m, 1, symW)

    if (iEnd >= maxCore) then
        write (6, *) 'iCore not enough for p_Wmbej_r1abci_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wmbej_r1abci_h

! 2025.12.17
subroutine Wmbej_r1abci_h_spin(WmbejAAAA, WmbejABAB, WmbejABBA, r1AA, file_abciAAAA, file_abciABAB, &
                               pop_mA, pop_mB, pop_bA, pop_bB, pop_eA, pop_eB, pop_jA, pop_jB, &
                               pop_fA, pop_fB, iop, symR, symV, &
                               iCore, maxCore, fact, AAAA, ABAB, ABBA)
    implicit none

    real*8, dimension(*), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejABBA
    real*8, dimension(*), intent(in) :: r1AA
    character(len=8), intent(in) :: file_abciAAAA, file_abciABAB
    integer*8, dimension(8), intent(in) :: pop_mA, pop_mB, pop_bA, pop_bB, pop_eA, pop_eB
    integer*8, dimension(8), intent(in) :: pop_jA, pop_jB, pop_fA, pop_fB
    integer*8, intent(in) :: iop
    integer*8, intent(in) :: symR, symV
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact
    logical, intent(in) ::  AAAA, ABAB, ABBA

    if (AAAA) then
        ! Wmbej += <mb||ef>*r(f,j)
        call p_Wmbej_r1abci_h(WmbejAAAA, r1AA, file_abciAAAA, pop_mA, pop_bA, pop_eA, pop_jA, &
                              pop_fA, iop, symR, symV, iCore, maxCore, fact, .false.)
    end if
    if (ABAB) then
        ! WmBeJ += <mB||eF>*r(F,J)
        call p_Wmbej_r1abci_h(WmbejABAB, r1AA, file_abciABAB, pop_mA, pop_bB, pop_eA, pop_jB, &
                              pop_fB, 1, symR, symV, iCore, maxCore, fact, .false.)
    end if
    if (ABBA) then
        ! WmBEj += <mB||Ef>*r(f,j)
        call p_Wmbej_r1abci_h(WmbejABBA, r1AA, file_abciABAB, pop_mA, pop_bB, pop_eB, pop_jA, &
                              pop_fA, 1, symR, symV, iCore, maxCore, fact, .true.)
    end if
end subroutine Wmbej_r1abci_h_spin
