! 2026.7.4
!> F(a,e) -= fact*f(e,l)*r(a,l)
subroutine Fae_r1fia(Fae, r1, fia, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fae
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: fia
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(r1, fia, Fae, -fact, 'al', 'el', 'ae', iCore, maxCore)

end subroutine Fae_r1fia

! 2026.7.6
!> F(a,e) -= fact*f(e,l)*r(a,l)
subroutine Fae_r1fia_R(FaeAA, r1AA, fiaAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FaeAA
    type(opt), intent(in) :: r1AA, fiaAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    ! F(a,e) -= f(e,l)*r(a,l)
    call Fae_r1fia(FaeAA, r1AA, fiaAA, fact, iCore, maxCore)

end subroutine Fae_r1fia_R

! 2026.7.4
!> F(a,e) -= fact*f(e,l)*r(a,l)
subroutine Fae_r1fia_U(FaeAA, FaeBB, r1AA, r1BB, fiaAA, fiaBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FaeAA, FaeBB
    type(opt), intent(in) :: r1AA, r1BB, fiaAA, fiaBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    ! F(a,e) -= f(e,l)*r(a,l)
    call Fae_r1fia(FaeAA, r1AA, fiaAA, fact, iCore, maxCore)
    ! F(A,E) -= f(E,L)*r(A,L)
    call Fae_r1fia(FaeBB, r1BB, fiaBB, fact, iCore, maxCore)

end subroutine Fae_r1fia_U
