! 2026.7.6
!> c(a,i) += fact*r(e,i)*F(a,e)
subroutine c1_r1Fae(c1, r1, Fae, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: Fae
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(Fae, r1, c1, fact, 'ae', 'ei', 'ai', iCore, maxCore)

end subroutine c1_r1Fae

! 2026.7.6
!> c(a,i) += fact*r(e,i)*F(a,e)
subroutine c1_r1Fae_R(c1AA, r1AA, FaeAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: FaeAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += fact*r(e,i)*F(a,e)
    call c1_r1Fae(c1AA, r1AA, FaeAA, fact, iCore, maxCore)

end subroutine c1_r1Fae_R

! 2026.7.6
!> c(a,i) += fact*r(e,i)*F(a,e)
subroutine c1_r1Fae_U(c1AA, c1BB, r1AA, r1BB, FaeAA, FaeBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: FaeAA, FaeBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += fact*r(e,i)*F(a,e)
    call c1_r1Fae(c1AA, r1AA, FaeAA, fact, iCore, maxCore)
    ! c(A,I) += fact*r(E,I)*F(A,E)
    call c1_r1Fae(c1BB, r1BB, FaeBB, fact, iCore, maxCore)

end subroutine c1_r1Fae_U
