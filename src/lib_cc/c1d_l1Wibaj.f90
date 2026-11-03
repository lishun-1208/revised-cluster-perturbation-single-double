! 2026.7.10
!> c(a,i) += fact*Wibaj(ai,bj)*l(b,j)
subroutine c1d_l1Wibaj(c1d, l1, Wibaj, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: l1
    type(opt), intent(in) :: Wibaj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(Wibaj, l1, c1d, fact, 'aibj', 'bj', 'ai', iCore, maxCore)

end subroutine c1d_l1Wibaj

! 2026.7.10
!> c(a,i) += fact*Wibaj(ai,bj)*l(b,j)
subroutine c1d_l1Wibaj_R(c1dAA, l1AA, WibajAAAA, WibajABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA
    type(opt), intent(in) :: l1AA
    type(opt), intent(in) :: WibajAAAA, WibajABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = Wibaj(ai,bj)*l(b,j)
    !         +WiBaJ(ai,BJ)*l(B,J) -> WiBaJ(ai,BJ)*l(b,j)
    call c1d_l1Wibaj(c1dAA, l1AA, WibajAAAA, fact, iCore, maxCore)
    call c1d_l1Wibaj(c1dAA, l1AA, WibajABAB, fact, iCore, maxCore)

end subroutine c1d_l1Wibaj_R

! 2026.7.10
!> c(a,i) += fact*Wibaj(ai,bj)*l(b,j)
subroutine c1d_l1Wibaj_U(c1dAA, c1dBB, l1AA, l1BB, &
                         WibajAAAA, WibajABAB, WibajBBBB, WibajBABA, &
                         fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA, c1dBB
    type(opt), intent(in) :: l1AA, l1BB
    type(opt), intent(in) :: WibajAAAA, WibajABAB, WibajBBBB, WibajBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = Wibaj(ai,bj)*l(b,j)
    !         +WiBaJ(ai,BJ)*l(B,J)
    call c1d_l1Wibaj(c1dAA, l1AA, WibajAAAA, fact, iCore, maxCore)
    call c1d_l1Wibaj(c1dAA, l1BB, WibajABAB, fact, iCore, maxCore)

    ! c(A,I) = WIbAj(AI,bj)*l(b,j)
    !         +WIBAJ(AI,BJ)*l(B,J)
    call c1d_l1Wibaj(c1dBB, l1AA, WibajBABA, fact, iCore, maxCore)
    call c1d_l1Wibaj(c1dBB, l1BB, WibajBBBB, fact, iCore, maxCore)

end subroutine c1d_l1Wibaj_U
