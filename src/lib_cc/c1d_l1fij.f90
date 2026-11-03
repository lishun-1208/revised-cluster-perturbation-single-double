! 2026.7.10
!> c(a,i) += fact*f(i,j)*l(a,j)
subroutine c1d_l1fij(c1d, l1, fij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: l1
    type(opt), intent(in) :: fij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(l1, fij, c1d, fact, 'aj', 'ij', 'ai', iCore, maxCore)

end subroutine c1d_l1fij

! 2026.7.10
!> c(a,i) += fact*f(i,j)*l(a,j)
subroutine c1d_l1fij_R(c1dAA, l1AA, fijAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA
    type(opt), intent(in) :: l1AA
    type(opt), intent(in) :: fijAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += f(i,j)*l(a,j)
    call c1d_l1fij(c1dAA, l1AA, fijAA, fact, iCore, maxCore)

end subroutine c1d_l1fij_R

! 2026.7.10
!> c(a,i) += fact*f(i,j)*l(a,j)
subroutine c1d_l1fij_U(c1dAA, c1dBB, l1AA, l1BB, fijAA, fijBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA, c1dBB
    type(opt), intent(in) :: l1AA, l1BB
    type(opt), intent(in) :: fijAA, fijBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += f(i,j)*l(a,j)
    call c1d_l1fij(c1dAA, l1AA, fijAA, fact, iCore, maxCore)

    ! c(A,I) += f(I,J)*l(A,J)
    call c1d_l1fij(c1dBB, l1BB, fijBB, fact, iCore, maxCore)

end subroutine c1d_l1fij_U
