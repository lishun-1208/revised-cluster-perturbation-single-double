! 2026.4.21
!> c(a,i) += fact*0.5*W(cb,aj)*l(cb,ij)
subroutine c1d_l2Wcbaj(c1d, l2, Wcbaj, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Wcbaj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call c1_r2abci(c1d, l2, Wcbaj, fact, iCore, maxCore)

end subroutine c1d_l2Wcbaj
