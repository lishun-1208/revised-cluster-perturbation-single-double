! 2026.7.10
!> c(ab,ij) += fact*P(ab){F(c,b)*l(ac,ij)}
subroutine c2d_l2Fcb(c2d, l2, Fcb, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Fcb
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(l2, Fcb, c2d, fact, 'acij', 'cb', 'abij', iCore, maxCore, 1)

end subroutine c2d_l2Fcb

! 2026.7.10
!> c(ab,ij) += fact*P(ab){F(c,b)*l(ac,ij)}
subroutine c2d_l2Fcb_R(c2dABAB, l2ABAB, FcbAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dABAB
    type(opt), intent(in) :: l2ABAB
    type(opt), intent(in) :: FcbAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) = F(C,B)*l(aC,iJ) -> F(c,b)*l(aC,iJ)
    !           -F(c,a)*l(Bc,iJ) = F(c,a)*l(cB,iJ)
    call vv2v_424(l2ABAB, FcbAA, c2dABAB, fact, 'acij', 'cb', 'abij', iCore, maxCore, 0)
    call vv2v_424(l2ABAB, FcbAA, c2dABAB, fact, 'cbij', 'ca', 'abij', iCore, maxCore, 0)

end subroutine c2d_l2Fcb_R

! 2026.7.10
!> c(ab,ij) += fact*P(ab){F(c,b)*l(ac,ij)}
subroutine c2d_l2Fcb_U(c2dAAAA, c2dABAB, c2dBBBB, l2AAAA, l2ABAB, l2BBBB, &
                       FcbAA, FcbBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dAAAA, c2dABAB, c2dBBBB
    type(opt), intent(in) :: l2AAAA, l2ABAB, l2BBBB
    type(opt), intent(in) :: FcbAA, FcbBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) = P(ab){F(c,b)*l(ac,ij)}
    call c2d_l2Fcb(c2dAAAA, l2AAAA, FcbAA, fact, iCore, maxCore)

    ! c(aB,iJ) = F(C,B)*l(aC,iJ)
    !           -F(c,a)*l(Bc,iJ) = F(c,a)*l(cB,iJ)
    call vv2v_424(l2ABAB, FcbBB, c2dABAB, fact, 'aCiJ', 'CB', 'aBiJ', iCore, maxCore, 0)
    call vv2v_424(l2ABAB, FcbAA, c2dABAB, fact, 'cBiJ', 'ca', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) = P(AB){F(C,B)*l(AC,IJ)}
    call c2d_l2Fcb(c2dBBBB, l2BBBB, FcbBB, fact, iCore, maxCore)

end subroutine c2d_l2Fcb_U
