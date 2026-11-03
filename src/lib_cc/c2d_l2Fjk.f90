! 2026.7.10
!> c(ab,ij) -= fact*P(ij){F(j,k)*l(ab,ik)}
subroutine c2d_l2Fjk(c2d, l2, Fjk, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Fjk
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(l2, Fjk, c2d, -fact, 'abik', 'jk', 'abij', iCore, maxCore, 2)

end subroutine c2d_l2Fjk

! 2026.7.10
!> c(ab,ij) -= fact*P(ij){F(j,k)*l(ab,ik)}
subroutine c2d_l2Fjk_R(c2dABAB, l2ABAB, FjkAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dABAB
    type(opt), intent(in) :: l2ABAB
    type(opt), intent(in) :: FjkAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) = -F(J,K)*l(aB,iK) -> -F(j,k)*l(aB,iK)
    !            +F(i,k)*l(aB,Jk) = -F(i,k)*l(aB,kJ)
    call vv2v_424(l2ABAB, FjkAA, c2dABAB, -fact, 'abik', 'jk', 'abij', iCore, maxCore, 0)
    call vv2v_424(l2ABAB, FjkAA, c2dABAB, -fact, 'abkj', 'ik', 'abij', iCore, maxCore, 0)

end subroutine c2d_l2Fjk_R

! 2026.7.10
!> c(ab,ij) -= fact*P(ij){F(j,k)*l(ab,ik)}
subroutine c2d_l2Fjk_U(c2dAAAA, c2dABAB, c2dBBBB, l2AAAA, l2ABAB, l2BBBB, &
                       FjkAA, FjkBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dAAAA, c2dABAB, c2dBBBB
    type(opt), intent(in) :: l2AAAA, l2ABAB, l2BBBB
    type(opt), intent(in) :: FjkAA, FjkBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) = -P(ij){F(j,k)*l(ab,ik)}
    call c2d_l2Fjk(c2dAAAA, l2AAAA, FjkAA, fact, iCore, maxCore)

    ! c(aB,iJ) = -F(J,K)*l(aB,iK)
    !            +F(i,k)*l(aB,Jk) = -F(i,k)*l(aB,kJ)
    call vv2v_424(l2ABAB, FjkBB, c2dABAB, -fact, 'aBiK', 'JK', 'aBiJ', iCore, maxCore, 0)
    call vv2v_424(l2ABAB, FjkAA, c2dABAB, -fact, 'aBkJ', 'ik', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) = -P(IJ){F(J,K)*l(AB,IK)}
    call c2d_l2Fjk(c2dBBBB, l2BBBB, FjkBB, fact, iCore, maxCore)

end subroutine c2d_l2Fjk_U
