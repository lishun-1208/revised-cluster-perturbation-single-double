! 2026.7.6
!> c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
subroutine c1_r2abci(c1, r2, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abci%incore) then
        if (r2%ioppq == 0) then
            call vv2v_442(abci, r2, c1, fact, 'feam', 'feim', 'ai', iCore, maxCore)
        else
            call vv2v_442(abci, r2, c1, 0.5d0*fact, 'feam', 'feim', 'ai', iCore, maxCore)
        end if
    else
        write (6, *) 'ERROR: c1_r2abci current not support outcore.'
        stop
        call p_c1_r2abci_h(c1%arr, r2%arr, abci%f, c1%popp, c1%popq, &
                           r2%popp, r2%popq, r2%popr, iCore, fact, r2%ioppq, r2%ioprs, r2%sym)
    end if

end subroutine c1_r2abci

! can use for abci(Ef,Am)*t(Ef,Im)
!-----------------------------
! modified at 2025.11.10
subroutine p_c1_r2abci_h(c1, r2, file_abci, pop_a, pop_i, pop_e, pop_f, pop_m, &
                         iCore, fact, iopef, iopim, symR)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c1
    real*8, dimension(*), intent(in) :: r2
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_i, pop_e, pop_f, pop_m
    real*8, dimension(*), intent(inout) :: iCore
    real*8, intent(in) :: fact
    integer*8, intent(in) :: iopef, iopim
    integer*8, intent(in) :: symR

    real*8 :: fact2
    integer*8 :: irrepm, numm, m
    integer*8 :: irrepT_m ! sym of t(ef,i*)
    integer*8 :: i0, i1, i2
    integer*8, external :: idsymsoc

    if (iopef == 1) then
        fact2 = 1.0d0
    else if (iopef == 0) then
        fact2 = 2.0d0
    end if

    i0 = 1
    ! cycle for irrepm
    do irrepm = 1, nirrep
        irrepT_m = dirprd(irrepm, symR)
        numm = pop_m(irrepm)

        i1 = i0 + idsymsoc(irrepm, pop_e, pop_f, iopef, pop_a, zorb, 1) ! i0 for v(ef,a*)
        ! cycle for m
        do m = 1, numm
            ! get abci(ef,a*), i0
            call getpqrso(iCore(i0), m, irrepm, 1, file_abci)
            ! get r2(ef,i*), i1
            call getpqrso_incore(iCore(i1), m, irrepm, pop_e, pop_f, iopef, &
                                 pop_i, pop_m, iopim, r2, symR)
            ! c1(*a,*i) += fact*0.5*abci(ef,*a)*r2(ef,*i)
            call VmnpqVmnrs_to_Vpqrs(iCore(i0), iCore(i1), irrepm, irrepT_m, pop_e, pop_f, iopef, &
                                     zorb, pop_a, 1, zorb, pop_i, 1, c1, symR, 0.5D0*fact*fact2, 1.0d0)
        end do
    end do

end subroutine p_c1_r2abci_h

! 2026.7.6
!> c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
subroutine c1_r2abci_R(c1AA, r2AAAA, r2ABAB, abciAAAA, abciABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: abciAAAA, abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
    call c1_r2abci(c1AA, r2AAAA, abciAAAA, fact, iCore, maxCore)
    ! c(a,i) += fact*r(fE,iM)*V(fE,aM)
    call c1_r2abci(c1AA, r2ABAB, abciABAB, 2.0d0*fact, iCore, maxCore)

end subroutine c1_r2abci_R

! 2026.7.6
!> c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
subroutine c1_r2abci_U(c1AA, c1BB, r2AAAA, r2ABAB, r2BBBB, abciAAAA, abciABAB, abciBBBB, abciBABA, &
                       fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1
    type(opt) :: tmp

    ! c(a,i) += fact*0.5*r(fe,im)*V(fe,am)
    call c1_r2abci(c1AA, r2AAAA, abciAAAA, fact, iCore, maxCore)
    ! c(a,i) += fact*r(fE,iM)*V(fE,aM)
    call c1_r2abci(c1AA, r2ABAB, abciABAB, 2.0d0*fact, iCore, maxCore)

    ! c(A,I) += fact*0.5*r(FE,IM)*V(FE,AM)
    call c1_r2abci(c1BB, r2BBBB, abciBBBB, fact, iCore, maxCore)

    !--- c(A,I) += fact*r(Fe,Im)*V(Fe,Am)
    ! r(eF,mI)->r(Fe,Im), i0
    i0 = 1
    call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                  1, 1, r2ABAB%sym, T, i0)
    i1 = i0 + tmp%len
    call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                  r2ABAB%popr, r2ABAB%pops, 1, 1, tmp%sym)
    call c1_r2abci(c1BB, tmp, abciBABA, 2.0d0*fact, iCore(i1), maxCore - tmp%len)

end subroutine c1_r2abci_U

! 2025.11.10
subroutine c1_r2abci_h_spin(c1AA, r2AAAA, r2ABAB, file_abciAAAA, file_abciABAB, pop_aA, pop_iA, &
                            pop_eA, pop_eB, pop_fA, pop_fB, pop_mA, pop_mB, &
                            iopfe, iopim, symR, iCore, maxCore, fact)
    use mod_orbit
    use mod_size
    implicit none

    real*8, dimension(*), intent(inout) :: c1AA
    real*8, dimension(*), intent(in) :: r2AAAA, r2ABAB
    character(len=8), intent(in) :: file_abciAAAA, file_abciABAB
    integer*8, dimension(8), intent(in) :: pop_aA, pop_iA
    integer*8, dimension(8), intent(in) :: pop_eA, pop_eB, pop_fA, pop_fB, pop_mA, pop_mB
    integer*8, intent(in) :: iopfe, iopim
    integer*8, intent(in) :: symR
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact

    ! c(a,i) += 0.5*<ja||bc>*r(bc,ji)
    call p_c1_r2abci_h(c1AA, r2AAAA, file_abciAAAA, pop_aA, pop_iA, &
                       pop_eA, pop_fA, pop_mA, iCore, fact, iopfe, iopim, symR)
    ! c(a,i) += <Ja||Bc>*r(Bc,Ji)
    call p_c1_r2abci_h(c1AA, r2ABAB, file_abciABAB, pop_aA, pop_iA, &
                       pop_eA, pop_fB, pop_mA, iCore, 2.0d0*fact, 1, 1, symR)

end subroutine c1_r2abci_h_spin
