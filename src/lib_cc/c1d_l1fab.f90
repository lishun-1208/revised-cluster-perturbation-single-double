! 2026.7.10
!> c(a,i) -= fact*f(b,a)*l(b,i)
subroutine c1d_l1fab(c1d, l1, fab, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: l1
    type(opt), intent(in) :: fab
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_222(fab, l1, c1d, -fact, 'ba', 'bi', 'ai', iCore, maxCore)

end subroutine c1d_l1fab

! 2026.7.10
!> c(a,i) -= fact*f(b,a)*l(b,i)
subroutine c1d_l1fab_R(c1dAA, l1AA, fabAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA
    type(opt), intent(in) :: l1AA
    type(opt), intent(in) :: fabAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -f(b,a)*l(b,i)
    call c1d_l1fab(c1dAA, l1AA, fabAA, fact, iCore, maxCore)

end subroutine c1d_l1fab_R

! 2026.7.10
!> c(a,i) -= fact*f(b,a)*l(b,i)
subroutine c1d_l1fab_U(c1dAA, c1dBB, l1AA, l1BB, fabAA, fabBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA, c1dBB
    type(opt), intent(in) :: l1AA, l1BB
    type(opt), intent(in) :: fabAA, fabBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -f(b,a)*l(b,i)
    call c1d_l1fab(c1dAA, l1AA, fabAA, fact, iCore, maxCore)

    ! c(A,I) = -f(B,A)*l(B,I)
    call c1d_l1fab(c1dBB, l1BB, fabBB, fact, iCore, maxCore)

end subroutine c1d_l1fab_U
