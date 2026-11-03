! 2026.7.10
!> W(ij,kb) -= fact*U(kl,ij)*r(b,l)
subroutine Wkbij_r1Uklij(Wkbij, r1, Uklij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wkbij
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: Uklij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(Uklij, r1, Wkbij, -fact, 'klij', 'bl', 'ijkb', iCore, maxCore, 0)

end subroutine Wkbij_r1Uklij

! 2026.7.10
!> W(ij,kb) -= fact*U(kl,ij)*r(b,l)
subroutine Wkbij_r1Uklij_R(WkbijABAB, r1AA, UklijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: UklijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(iJ,kB) = -U(kL,iJ)*r(B,L) -> -U(kL,iJ)*r(b,l)
    call Wkbij_r1Uklij(WkbijABAB, r1AA, UklijABAB, fact, iCore, maxCore)

end subroutine Wkbij_r1Uklij_R

! 2026.7.10
!> W(ij,kb) -= fact*U(kl,ij)*r(b,l)
subroutine Wkbij_r1Uklij_U(WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA, spin, &
                           r1AA, r1BB, UklijAAAA, UklijABAB, UklijBBBB, &
                           fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA
    character(len=4) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: UklijAAAA, UklijABAB, UklijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ij,kb) = -U(kl,ij)*r(b,l)
        call Wkbij_r1Uklij(WkbijAAAA, r1AA, UklijAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(iJ,kB) = -U(kL,iJ)*r(B,L)
        call Wkbij_r1Uklij(WkbijABAB, r1BB, UklijABAB, fact, iCore, maxCore)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(IJ,KB) = -U(KL,IJ)*r(B,L)
        call Wkbij_r1Uklij(WkbijBBBB, r1BB, UklijBBBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(1:1) == 'T') then
        ! W(Ij,Kb) = -U(Kl,Ij)*r(b,l) = -U(lK,jI)*r(b,l)
        call vv2v_424(UklijABAB, r1AA, WkbijBABA, -fact, 'lKjI', 'bl', 'IjKb', iCore, maxCore, 0)
    end if

end subroutine Wkbij_r1Uklij_U
