! 2026.7.6
!> F(m,i) += fact*r(d,i)*f(d,m)
subroutine Fmi_r1fia(Fmi, r1, fia, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fmi
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: fia
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(fia, r1, Fmi, fact, 'dm', 'di', 'mi', iCore, maxCore)

end subroutine Fmi_r1fia

! 2026.7.6
!> F(m,i) += fact*r(d,i)*f(d,m)
subroutine Fmi_r1fia_R(FmiAA, r1AA, fiaAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmiAA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: fiaAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*r(d,i)*f(d,m)
    call Fmi_r1fia(FmiAA, r1AA, fiaAA, fact, iCore, maxCore)

end subroutine Fmi_r1fia_R

! 2026.7.6
!> F(m,i) += fact*r(d,i)*f(d,m)
subroutine Fmi_r1fia_U(FmiAA, FmiBB, r1AA, r1BB, fiaAA, fiaBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmiAA, FmiBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: fiaAA, fiaBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*r(d,i)*f(d,m)
    call Fmi_r1fia(FmiAA, r1AA, fiaAA, fact, iCore, maxCore)
    ! F(M,I) += fact*r(D,I)*f(D,M)
    call Fmi_r1fia(FmiBB, r1BB, fiaBB, fact, iCore, maxCore)

end subroutine Fmi_r1fia_U
