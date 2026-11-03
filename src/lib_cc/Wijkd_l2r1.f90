! write at 2024.7.4, lishun
! calculate Wijkd += fact*r(c,k)*l(ij,cd)
! Lijcd sort as (cd,ij)
!> be deprecated, please use Wlkjd_r1abij
subroutine Wijkd_l2r1(Wijkd, l2, r1, pop_i, pop_j, pop_k, pop_d, &
                      pop_c, iop, symL, symR, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: Wijkd
    real*8, dimension(*), intent(in) :: l2, r1
    integer*8, dimension(8), intent(in) :: pop_i, pop_j, pop_k, pop_d, pop_c
    integer*8, intent(in) :: iop
    integer*8, intent(in) :: symL, symR
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1, i2, nSize
    integer*8, external :: irpdso, idsymsoc

    i0 = 1
    ! r(c,k)->r(k,c), i0
    nSize = irpdso(symR, pop_c, pop_k, 1)
    i1 = i0 + nSize
    call itranspso(r1, iCore(i0), pop_c, zorb, 1, pop_k, zorb, 1, symR)
    ! l(cd,ij)->l(ij,cd), i1
    nSize = idsymsoc(symL, pop_c, pop_d, iop, pop_i, pop_j, iop)
    i2 = i1 + nSize
    call itranspso(l2, iCore(i1), pop_c, pop_d, iop, pop_i, pop_j, iop, symL)
    ! Wijkd += fact*r(k,c)*l(ij,cd)
    call Wabdj_r2fld(Wijkd, iCore(i1), iCore(i0), pop_i, pop_j, pop_k, pop_d, &
                     pop_c, iop, symL, symR, -fact, iCore(i2), maxCore - i2 + 1)

end subroutine Wijkd_l2r1
