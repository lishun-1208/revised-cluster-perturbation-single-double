! write at 2025.7.4, lihsun
! calculate Cia(a,i) += fact*r(b,j)*l(ij,ab)
subroutine c1d_l2r1(c1d, r1, l2, pop_i, pop_a, pop_j, pop_b, iop, symL, symR, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: c1d
    real*8, dimension(*), intent(in) :: r1, l2
    integer*8, dimension(8), intent(in) :: pop_i, pop_a, pop_j, pop_b
    integer*8, intent(in) :: iop ! iop=iopab=iopij
    integer*8, intent(in) :: symL, symR
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore ! free work space
    integer*8, intent(in) :: maxCore

    ! C(a,i) += fact*r(b,j)*l(ab,ij)
    call c1_r2Fme(c1d, l2, r1, pop_a, pop_i, &
                  pop_j, pop_b, iop, symL, symR, iCore, maxCore, fact)

end subroutine c1d_l2r1
