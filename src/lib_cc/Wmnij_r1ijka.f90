! 2026.7.6
!> W(mn,ij) += fact*P(ij){r(e,j)*V(mn,ie)}
subroutine Wmnij_r1ijka(Wmnij, r1, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wmnij
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(ijka, r1, Wmnij, fact, 'mnie', 'ej', 'mnij', iCore, maxCore, 2)

end subroutine Wmnij_r1ijka

! 2026.7.6
!> W(mn,ij) += fact*P(ij){r(e,j)*V(mn,ie)}
subroutine Wmnij_r1ijka_R(WmnijABAB, r1AA, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmnijABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(mN,iJ) = r(E,J)*V(mN,iE)
    call vv2v_424(ijkaABAB, r1AA, WmnijABAB, fact, 'mNiE', 'EJ', 'mNiJ', iCore, maxCore, 0)
    ! W(mN,iJ) =-r(e,i)*V(mN,Je) = r(e,i)*V(Nm,Je) -> r(e,i)*V(nM,jE)
    call vv2v_424(ijkaABAB, r1AA, WmnijABAB, fact, 'nmje', 'ei', 'mnij', iCore, maxCore, 0)

end subroutine Wmnij_r1ijka_R

! 2026.7.6
!> W(mn,ij) += fact*P(ij){r(e,j)*V(mn,ie)}
subroutine Wmnij_r1ijka_U(WmnijAAAA, WmnijABAB, WmnijBBBB, spin, &
                          r1AA, r1BB, ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmnijAAAA, WmnijABAB, WmnijBBBB
    character(len=3), intent(in) :: spin
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(mn,ij) = P(ij){r(e,j)*V(mn,ie)}
        call Wmnij_r1ijka(WmnijAAAA, r1AA, ijkaAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(mN,iJ) = r(E,J)*V(mN,iE)
        call vv2v_424(ijkaABAB, r1BB, WmnijABAB, fact, 'mNiE', 'EJ', 'mNiJ', iCore, maxCore, 0)
        ! W(mN,iJ) = -r(e,i)*V(mN,Je) = r(e,i)*V(Nm,Je)
        call vv2v_424(ijkaBABA, r1AA, WmnijABAB, fact, 'NmJe', 'ei', 'mNiJ', iCore, maxCore, 0)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(MN,IJ) = P(IJ){r(E,J)*V(MN,IE)}
        call Wmnij_r1ijka(WmnijBBBB, r1BB, ijkaBBBB, fact, iCore, maxCore)
    end if

end subroutine Wmnij_r1ijka_U
