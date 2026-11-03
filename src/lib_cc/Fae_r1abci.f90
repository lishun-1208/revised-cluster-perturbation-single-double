! 2026.7.4
!> F(a,e) += fact*r(f,m)*V(ef,am)
subroutine Fae_r1abci(Fae, r1, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fae
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abci%incore) then
        call vv2v_422(abci, r1, Fae, fact, 'efam', 'fm', 'ae', iCore, maxCore)
    else
        write (6, *) 'ERROR: Fae_r1abci current not support outcore.'
        stop
        !call p_Fae_r1abci_h(Fae%arr, r1%arr, abci%f, Fae%popp, Fae%popq, &
        !                    r1%popp, r1%popq, abci%ioppq, r1%sym, iCore, maxCore, fact)
    end if

end subroutine Fae_r1abci

! write by LiShun, 2024.5.28
! can use for VABAB*rBB
!---------------------------------
! modified 2024.10.21
subroutine p_Fae_r1abci_h(Fae, r1, file_abci, pop_a, pop_e, pop_f, pop_m, iopef, symR, &
                          iCore, maxCore, fact)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: Fae
    real*8, dimension(*), intent(in) :: r1
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_e, pop_f, pop_m
    integer*8, intent(in) :: iopef
    integer*8, intent(in) :: symR
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact

    integer*8 :: irrepm, numm, m
    integer*8 :: symR_m ! symmerty of t(f,*)
    integer*8 :: i0, i1, i2, i3
    integer*8 :: nSize
    integer*8, external :: idsymsoc, irpdso

    i0 = 1
    do irrepm = 1, nirrep
        numm = pop_m(irrepm)
        symR_m = dirprd(symR, irrepm)
        do m = 1, numm
            ! get abci(ef,a*), i0
            nSize = idsymsoc(irrepm, pop_e, pop_f, 1, pop_a, zorb, 1)
            i1 = i0 + nSize
            if (iopef == 1) then
                call getpqrso(iCore(i0), m, irrepm, 1, file_abci)
            else if (iopef == 0) then
                call getpqrso(iCore(i1), m, irrepm, 1, file_abci)
                call iexpso(iCore(i1), iCore(i0), pop_e, pop_f, 0, 1, &
                            pop_a, zorb, 1, 0, irrepm) ! expand a<b to a,b
            end if
            ! get t(**,f*), i1
            i2 = i1 + pop_f(symR_m)
            call getpqrso_incore(iCore(i1), m, irrepm, zorb, zorb, 1, &
                                 pop_f, pop_m, 1, r1, symR)
            ! abci(ef,a*)->abci(f*,ae), i2
            nSize = idsymsoc(irrepm, pop_e, pop_f, 1, pop_a, zorb, 1)
            i3 = i2 + nSize
            call sstgenso(iCore(i0), iCore(i2), nSize, pop_e, pop_f, pop_a, zorb, &
                          iCore(i3), irrepm, '2431')
            ! Fae += v(f*,ae)*t(f*,**)
            call VmnpqVmnrs_to_Vpqrs(iCore(i2), iCore(i1), irrepm, symR_m, pop_f, zorb, 1, &
                                     pop_a, pop_e, 1, zorb, zorb, 1, Fae, symR, fact, 1.0D0)

        end do
    end do

end subroutine p_Fae_r1abci_h

! 2026.7.6
!> F(a,e) += fact*r(f,m)*V(ef,am)
subroutine Fae_r1abci_R(FaeAA, r1AA, abciAAAA, abciABAB, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: FaeAA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abciAAAA, abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(a,e) += r(f,m)*V(ef,am)
    call Fae_r1abci(FaeAA, r1AA, abciAAAA, fact, iCore, maxCore)
    ! F(a,e) += r(F,M)*V(eF,aM) -> r(f,m)*V(eF,aM)
    call Fae_r1abci(FaeAA, r1AA, abciABAB, fact, iCore, maxCore)

end subroutine Fae_r1abci_R

! 2026.7.6
!> F(a,e) += fact*r(f,m)*V(ef,am)
subroutine Fae_r1abci_U(FaeAA, FaeBB, r1AA, r1BB, abciAAAA, abciABAB, abciBBBB, abciBABA, &
                        fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: FaeAA, FaeBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(a,e) += r(f,m)*V(ef,am)
    call Fae_r1abci(FaeAA, r1AA, abciAAAA, fact, iCore, maxCore)
    ! F(a,e) += r(F,M)*V(eF,aM)
    call Fae_r1abci(FaeAA, r1BB, abciABAB, fact, iCore, maxCore)

    ! F(A,E) += r(f,m)*V(Ef,Am)
    call Fae_r1abci(FaeBB, r1AA, abciBABA, fact, iCore, maxCore)
    ! F(A,E) += r(F,M)*V(EF,AM)
    call Fae_r1abci(FaeBB, r1BB, abciBBBB, fact, iCore, maxCore)

end subroutine Fae_r1abci_U
