! 2026.4.21
!> c(ab,ij) += fact*0.5*V(cd,ab)*l(cd,ij)
!> need V(ab,cd)=V(cd,ab)
subroutine c2d_l2abcd(c2d, l2, abcd, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: abcd
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! V(cd,ab)=V(ab,cd)
    call c2_r2abcd(c2d, l2, abcd, fact, iCore, maxCore)

end subroutine c2d_l2abcd
