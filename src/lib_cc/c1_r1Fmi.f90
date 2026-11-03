! 2026.7.6
!> c(a,i) -= fact*r(a,m)*F(m,i)
subroutine c1_r1Fmi(c1, r1, Fmi, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: Fmi
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(r1, Fmi, c1, -fact, 'am', 'mi', 'ai', iCore, maxCore)

end subroutine c1_r1Fmi

! 2026.7.6
!> c(a,i) -= fact*r(a,m)*F(m,i)
subroutine c1_r1Fmi_R(c1AA, r1AA, FmiAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: FmiAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) -= fact*r(a,m)*F(m,i)
    call c1_r1Fmi(c1AA, r1AA, FmiAA, fact, iCore, maxCore)

end subroutine c1_r1Fmi_R

! 2026.7.6
!> c(a,i) -= fact*r(a,m)*F(m,i)
subroutine c1_r1Fmi_U(c1AA, c1BB, r1AA, r1BB, FmiAA, FmiBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: FmiAA, FmiBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) -= fact*r(a,m)*F(m,i)
    call c1_r1Fmi(c1AA, r1AA, FmiAA, fact, iCore, maxCore)
    ! c(A,I) -= fact*r(A,M)*F(M,I)
    call c1_r1Fmi(c1BB, r1BB, FmiBB, fact, iCore, maxCore)

end subroutine c1_r1Fmi_U
