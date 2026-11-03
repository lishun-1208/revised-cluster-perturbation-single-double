! 2026.4.21
!> c(a,i) -= 0.5*fact*W(kj,ib)*l(ab,kj)
subroutine c1d_l2Wibkj(c1d, l2, Wibkj, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Wibkj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call c1_r2ijka(c1d, l2, Wibkj, fact, iCore, maxCore)

end subroutine c1d_l2Wibkj
