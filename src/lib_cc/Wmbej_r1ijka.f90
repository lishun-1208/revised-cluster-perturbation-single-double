! 2026.7.6
!> Wmbej(em,bj) -= fact*r(b,n)*V(nm,je)
subroutine Wmbej_r1ijka(Wmbej, r1, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wmbej
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(ijka, r1, Wmbej, -fact, 'nmje', 'bn', 'embj', iCore, maxCore, 0)

end subroutine Wmbej_r1ijka

! 2026.7.6
!> Wmbej(em,bj) -= fact*r(b,n)*V(nm,je)
subroutine Wmbej_r1ijka_R(WmbejAAAA, WmbejABAB, WmbejBAAB, spin, &
                          r1AA, ijkaAAAA, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB
    character(len=3) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) = -r(b,n)*V(nm,je)
        call Wmbej_r1ijka(WmbejAAAA, r1AA, ijkaAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) = -r(B,N)*V(Nm,Je) -> -r(b,n)*V(nM,jE)
        call Wmbej_r1ijka(WmbejABAB, r1AA, ijkaABAB, fact, iCore, maxCore)
    end if
    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) = -r(b,n)*V(nM,Je) = r(b,n)*V(Mn,Je) -> r(b,n)*V(mN,jE)
        call vv2v_424(ijkaABAB, r1AA, WmbejBAAB, fact, 'mnje', 'bn', 'embj', iCore, maxCore, 0)
    end if

end subroutine Wmbej_r1ijka_R

! 2026.7.6
!> Wmbej(em,bj) -= fact*r(b,n)*V(nm,je)
subroutine Wmbej_r1ijka_U(WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA, spin, &
                          r1AA, r1BB, ijkaAAAA, ijkaABAB, ijkaBABA, ijkaBBBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA
    character(len=6) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBABA, ijkaBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) = -r(b,n)*V(nm,je)
        call Wmbej_r1ijka(WmbejAAAA, r1AA, ijkaAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) = -r(B,N)*V(Nm,Je)
        call Wmbej_r1ijka(WmbejABAB, r1BB, ijkaBABA, fact, iCore, maxCore)
    end if
    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) = -r(b,n)*V(nM,Je) = r(b,n)*V(Mn,Je)
        call vv2v_424(ijkaBABA, r1AA, WmbejBAAB, fact, 'MnJe', 'bn', 'eMbJ', iCore, maxCore, 0)
    end if
    ! BBBB
    if (spin(4:4) == 'T') then
        ! WMBEJ(EM,BJ) = -r(B,N)*V(NM,JE)
        call Wmbej_r1ijka(WmbejBBBB, r1BB, ijkaBBBB, fact, iCore, maxCore)
    end if
    ! BABA
    if (spin(5:5) == 'T') then
        ! WMbEj(EM,bj) = -r(b,n)*V(nM,jE)
        call Wmbej_r1ijka(WmbejBABA, r1AA, ijkaABAB, fact, iCore, maxCore)
    end if
    ! ABBA
    if (spin(6:6) == 'T') then
        ! WmBEj(Em,Bj) = -r(B,N)*V(Nm,jE) = r(B,N)*V(mN,jE)
        call vv2v_424(ijkaABAB, r1BB, WmbejABBA, fact, 'mNjE', 'BN', 'EmBj', iCore, maxCore, 0)
    end if

end subroutine Wmbej_r1ijka_U
