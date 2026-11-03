! 2026.7.6
!> W(lk,jd) += fact*r(c,j)*V(cd,lk)
subroutine Wlkjd_r1abij(Wlkjd, r1, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wlkjd
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(abij, r1, Wlkjd, fact, 'cdlk', 'cj', 'lkjd', iCore, maxCore, 0)

end subroutine Wlkjd_r1abij

! 2026.7.6
!> W(lk,jd) += fact*r(c,j)*V(cd,lk)
subroutine Wlkjd_r1abij_R(WlkjdAAAA, WlkjdABAB, spin, r1AA, abijAAAA, abijABAB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WlkjdAAAA, WlkjdABAB
    character(len=2) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abijAAAA, abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(lk,jd) += fact*r(c,j)*V(cd,lk)
        call Wlkjd_r1abij(WlkjdAAAA, r1AA, abijAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(lK,jD) += fact*r(c,j)*V(cD,lK)
        call Wlkjd_r1abij(WlkjdABAB, r1AA, abijABAB, fact, iCore, maxCore)
    end if

end subroutine Wlkjd_r1abij_R

! 2026.7.6
!> W(lk,jd) += fact*r(c,j)*V(cd,lk)
subroutine Wlkjd_r1abij_U(WlkjdAAAA, WlkjdABAB, WlkjdBBBB, WlkjdBABA, spin, &
                          r1AA, r1BB, abijAAAA, abijABAB, abijBBBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WlkjdAAAA, WlkjdABAB, WlkjdBBBB, WlkjdBABA
    character(len=2) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(lk,jd) += fact*r(c,j)*V(cd,lk)
        call Wlkjd_r1abij(WlkjdAAAA, r1AA, abijAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(lK,jD) += fact*r(c,j)*V(cD,lK)
        call Wlkjd_r1abij(WlkjdABAB, r1AA, abijABAB, fact, iCore, maxCore)
    end if
    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(LK,JD) += fact*r(C,J)*V(CD,LK)
        call Wlkjd_r1abij(WlkjdBBBB, r1BB, abijBBBB, fact, iCore, maxCore)
    end if
    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Lk,Jd) += fact*r(C,J)*V(Cd,Lk) = r(C,J)*V(dC,kL)
        call vv2v_424(abijABAB, r1BB, WlkjdBABA, fact, 'dCkL', 'CJ', 'LkJd', iCore, maxCore, 0)
    end if

end subroutine Wlkjd_r1abij_U
