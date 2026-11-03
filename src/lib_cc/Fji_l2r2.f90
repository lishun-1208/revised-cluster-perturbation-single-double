! write at 2024.7.4, lishun
! calculate Fji += fact*0.5*t(dc,ik)*l(jk,dc)
! Ljkdc sort as (dc,jk)
!> ATTENTION: has been deprecated, plase use Fmi_r2abij
subroutine Fji_l2r2(Fji, l2, r2, pop_j, pop_i, &
                    pop_c, pop_d, pop_k, iop, symL, symR, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: Fji
    real*8, dimension(*), intent(in) :: l2, r2
    integer*8, dimension(8), intent(in) :: pop_j, pop_i, pop_c, pop_d, pop_k
    integer*8, intent(in) :: iop ! iop=iopdc=iopik=iopjk
    integer*8, intent(in) :: symL, symR
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symF
    integer*8 :: i0, i1, nSize
    integer*8, external :: irpdso, idsymsoc

    if (symR /= 1) then
        write (6, *) 'ERROR: Fji_l2r2 not support symR /= 1', symR
        stop
    end if

    i0 = 1
    symF = dirprd(symL, symR)
    ! tmp(i,j) = fact*0.5*t(dc,ik)*l(jk,dc)
    nSize = irpdso(symF, pop_i, pop_j, 1)
    i1 = i0 + nSize
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call Fmi_r2abij(iCore(i0), l2, r2, pop_i, pop_j, &
                    pop_d, pop_c, pop_k, fact, symL, iop, iop, iCore(i1), maxCore - nSize)
    ! tmp(i,j)->tmp(j,i), i1
    call itranspso(iCore(i0), iCore(i1), pop_i, zorb, 1, pop_j, zorb, 1, symF)
    ! Fji += tmp(j,i)
    call daxpy(nSize, 1.0d0, iCore(i1), 1, Fji, 1)

end subroutine Fji_l2r2
