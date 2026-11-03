! 2026.7.13
!> c(ab,ij) += fact*0.5*W(ij,kl)*l(ab,kl)
subroutine c2d_l2Wijkl(c2d, l2, Wijkl, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Wijkl
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (l2%ioprs == 1) then
        call vv2v_444(l2, Wijkl, c2d, 0.5d0*fact, 'abkl', 'ijkl', 'abij', iCore, maxCore, 0)
    else
        call vv2v_444(l2, Wijkl, c2d, fact, 'abkl', 'ijkl', 'abij', iCore, maxCore, 0)
    end if

end subroutine c2d_l2Wijkl

! 2026.7.13
!> c(ab,ij) += fact*0.5*W(ij,kl)*l(ab,kl)
subroutine c2d_l2Wijkl_R(c2dABAB, l2ABAB, WijklABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dABAB
    type(opt), intent(in) :: l2ABAB
    type(opt), intent(in) :: WijklABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) = W(iJ,kL)*l(aB,kL)
    call c2d_l2Wijkl(c2dABAB, l2ABAB, WijklABAB, 2.0d0*fact, iCore, maxCore)

end subroutine c2d_l2Wijkl_R

! 2026.7.13
!> c(ab,ij) += fact*0.5*W(ij,kl)*l(ab,kl)
subroutine c2d_l2Wijkl_U(c2dAAAA, c2dABAB, c2dBBBB, l2AAAA, l2ABAB, l2BBBB, &
                         WijklAAAA, WijklABAB, WijklBBBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dAAAA, c2dABAB, c2dBBBB
    type(opt), intent(in) :: l2AAAA, l2ABAB, l2BBBB
    type(opt), intent(in) :: WijklAAAA, WijklABAB, WijklBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) = 0.5*W(ij,kl)*l(ab,kl)
    call c2d_l2Wijkl(c2dAAAA, l2AAAA, WijklAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) = W(iJ,kL)*l(aB,kL)
    call c2d_l2Wijkl(c2dABAB, l2ABAB, WijklABAB, 2.0d0*fact, iCore, maxCore)

    ! c(AB,IJ) = 0.5*W(IJ,KL)*l(AB,KL)
    call c2d_l2Wijkl(c2dBBBB, l2BBBB, WijklBBBB, fact, iCore, maxCore)

end subroutine c2d_l2Wijkl_U
