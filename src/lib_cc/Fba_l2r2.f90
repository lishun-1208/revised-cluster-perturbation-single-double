! write at 2024.7.4, lishun
! calculate Fba -= fact*0.5*t(bc,lk)*l(lk,ac)
! Llkac sort as (ac,lk)
!> ATTENTION: has been deprecated, please use Fae_r2abij
subroutine Fba_l2r2(Fba, l2, r2, pop_b, pop_a, &
                    pop_c, pop_l, pop_k, iop, symL, symR, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: Fba
    real*8, dimension(*), intent(in) :: l2, r2
    integer*8, dimension(8), intent(in) :: pop_b, pop_a, pop_c, pop_l, pop_k
    integer*8, intent(in) :: iop ! iop=iopbc=ioplk=iopac
    integer*8, intent(in) :: symL, symR
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symF
    integer*8 :: i0, i1, nSize
    integer*8, external :: irpdso, idsymsoc

    if (symR /= 1) then
        write (6, *) 'ERROR: Fba_l2r2 not support symR /= 1', symR
        stop
    end if

    i0 = 1
    symF = dirprd(symL, symR)
    ! tmp(a,b) = -fact*0.5*t(bc,lk)*l(ac,lk)
    nSize = irpdso(symF, pop_a, pop_b, 1)
    i1 = i0 + nSize
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call Fae_r2abij(iCore(i0), l2, r2, pop_a, pop_b, &
                    pop_c, pop_l, pop_k, iCore(i1), maxCore - nSize, fact, iop, iop, symL)
    ! tmp(a,b)->tmp(b,a), i1
    call itranspso(iCore(i0), iCore(i1), pop_a, zorb, 1, pop_b, zorb, 1, symF)
    ! Fba += tmp(b,a)
    call daxpy(nSize, 1.0d0, iCore(i1), 1, Fba, 1)

end subroutine Fba_l2r2
