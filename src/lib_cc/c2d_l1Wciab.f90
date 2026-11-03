! 2026.4.21
!> c(ab,ij) += fact*P(ij){W(ab,cj)*l(c,i)}
subroutine c2d_l1Wciab(c2d, l1, Wciab, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l1
    type(opt), intent(in) :: Wciab
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call c2_r1abci(c2d, l1, Wciab, fact, iCore, maxCore)

end subroutine c2d_l1Wciab
