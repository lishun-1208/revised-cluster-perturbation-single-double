! 2026.7.6
!> Wmbej(em,bj) += fact*0.5*r(bf,jn)*V(fe,nm)
subroutine Wmbej_r2abij(Wmbej, r2, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wmbej
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_444(abij, r2, Wmbej, 0.5d0*fact, 'fenm', 'bfjn', 'embj', iCore, maxCore, 0)

end subroutine Wmbej_r2abij

! 2026.7.6
!> Wmbej(em,bj) += fact*0.5*r(bf,jn)*V(fe,nm)
subroutine Wmbej_r2abij_R(WmbejAAAA, WmbejABAB, WmbejBAAB, spin, &
                          r2AAAA, r2ABAB, abijAAAA, abijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB
    character(len=3) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: abijAAAA, abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) = 0.5*r(bf,jn)*V(fe,nm)
        call Wmbej_r2abij(WmbejAAAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
        ! Wmbej(em,bj) = 0.5*r(bF,jN)*V(Fe,Nm) -> 0.5*r(bF,jN)*V(fE,nM)
        call Wmbej_r2abij(WmbejAAAA, r2ABAB, abijABAB, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) = 0.5*r(Bf,Jn)*V(fe,nm) -> 0.5*r(bF,jN)*V(fe,nm)
        call Wmbej_r2abij(WmbejABAB, r2ABAB, abijAAAA, fact, iCore, maxCore)
        ! WmBeJ(em,BJ) = 0.5*r(BF,JN)*V(Fe,Nm) -> 0.5*r(bf,jn)*V(fE,nM)
        call Wmbej_r2abij(WmbejABAB, r2AAAA, abijABAB, fact, iCore, maxCore)
    end if
    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) = 0.5*r(bF,Jn)*V(Fe,nM) = 0.5*r(bF,nJ)*V(eF,nM)
        call vv2v_444(abijABAB, r2ABAB, WmbejBAAB, 0.5d0*fact, 'efnm', 'bfnj', 'embj', &
                      iCore, maxCore, 0)
    end if

end subroutine Wmbej_r2abij_R

! 2026.7.6
!> Wmbej(em,bj) += fact*0.5*r(bf,jn)*V(fe,nm)
subroutine Wmbej_r2abij_U(WmbejAAAA, WmbejABAB, WmbejBAAB, &
                          WmbejBBBB, WmbejBABA, WmbejABBA, spin, &
                          r2AAAA, r2ABAB, r2BBBB, abijAAAA, abijABAB, abijBBBB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA
    character(len=6) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! Wmbej(em,bj) = 0.5*r(bf,jn)*V(fe,nm)
        call Wmbej_r2abij(WmbejAAAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
        ! Wmbej(em,bj) = 0.5*r(bF,jN)*V(Fe,Nm) = 0.5*r(bF,jN)*V(eF,mN)
        call vv2v_444(abijABAB, r2ABAB, WmbejAAAA, 0.5d0*fact, 'eFmN', 'bFjN', 'embj', &
                      iCore, maxCore, 0)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! WmBeJ(em,BJ) = 0.5*r(Bf,Jn)*V(fe,nm) = 0.5*r(fB,nJ)*V(fe,nm)
        call vv2v_444(abijAAAA, r2ABAB, WmbejABAB, 0.5d0*fact, 'fenm', 'fBnJ', 'emBJ', &
                      iCore, maxCore, 0)
        ! WmBeJ(em,BJ) = 0.5*r(BF,JN)*V(Fe,Nm) = 0.5*r(BF,JN)*V(eF,mN)
        call vv2v_444(abijABAB, r2BBBB, WmbejABAB, 0.5d0*fact, 'eFmN', 'BFJN', 'emBJ', &
                      iCore, maxCore, 0)
    end if
    ! BAAB
    if (spin(3:3) == 'T') then
        ! WMbeJ(eM,bJ) = 0.5*r(bF,Jn)*V(Fe,nM) = 0.5*r(bF,nJ)*V(eF,nM)
        call vv2v_444(abijABAB, r2ABAB, WmbejBAAB, 0.5d0*fact, 'eFnM', 'bFnJ', 'eMbJ', &
                      iCore, maxCore, 0)
    end if
    ! BBBB
    if (spin(4:4) == 'T') then
        ! WMBEJ(EM,BJ) = 0.5*r(BF,JN)*V(FE,NM)
        call Wmbej_r2abij(WmbejBBBB, r2BBBB, abijBBBB, fact, iCore, maxCore)
        ! WMBEJ(EM,BJ) = 0.5*r(Bf,Jn)*V(fE,nM) = 0.5*r(fB,nJ)*V(fE,nM)
        call vv2v_444(abijABAB, r2ABAB, WmbejBBBB, 0.5d0*fact, 'fEnM', 'fBnJ', 'EMBJ', &
                      iCore, maxCore, 0)
    end if
    ! BABA
    if (spin(5:5) == 'T') then
        ! WMbEj(EM,bj) = 0.5*r(bf,jn)*V(fE,nM)
        call Wmbej_r2abij(WmbejBABA, r2AAAA, abijABAB, fact, iCore, maxCore)
        ! WMbEj(EM,bj) = 0.5*r(bF,jN)*V(FE,NM)
        call Wmbej_r2abij(WmbejBABA, r2ABAB, abijBBBB, fact, iCore, maxCore)
    end if
    ! ABBA
    if (spin(6:6) == 'T') then
        ! WmBEj(Em,Bj) = 0.5*r(Bf,jN)*V(fE,Nm) = 0.5*r(fB,jN)*V(fE,mN)
        call vv2v_444(abijABAB, r2ABAB, WmbejABBA, 0.5d0*fact, 'fEmN', 'fBjN', 'EmBj', &
                      iCore, maxCore, 0)
    end if

end subroutine Wmbej_r2abij_U
