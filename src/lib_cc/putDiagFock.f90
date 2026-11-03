! write by lishun, 2024.9.30
!> put orbital energy as diagonal element of Fock matrix
!> f(p,p) = fact*e(p)
!---------------------------
! modified at 2024.11.13
subroutine putDiagFock(Fock, Eorb, pop_a, pop_i, type, fact)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: Fock
    real*8, dimension(*), intent(in) :: Eorb
    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    character(len=1), intent(in) :: type !> 'O' for fii, 'V' for faa
    real*8, intent(in) :: fact

    integer*8 :: irrep_p, nump, p
    integer*8, dimension(8) :: ioffo, ioffv, ioffp ! position of Eorb belong to an irrep
    integer*8, dimension(8) :: pop_p
    integer*8, dimension(3) :: idx ! f(r,s) belong to irreps,s,r
    integer*8 :: idxE
    integer*8, parameter :: sym = 1
    integer*8, external :: isymoffso

    ! give ioff
    ioffo(1) = 0
    do irrep_p = 1, nirrep - 1
        nump = pop_i(irrep_p)
        ioffo(irrep_p + 1) = ioffo(irrep_p) + nump
    end do
    ioffv(1) = ioffo(nirrep) + pop_i(nirrep)
    do irrep_p = 1, nirrep - 1
        nump = pop_a(irrep_p)
        ioffv(irrep_p + 1) = ioffv(irrep_p) + nump
    end do

    if (type == 'O') then
        ioffp = ioffo
        pop_p = pop_i
    else if (type == 'V') then
        ioffp = ioffv
        pop_p = pop_a
    end if

    do irrep_p = 1, nirrep
        nump = pop_p(irrep_p)
        idx(1) = isymoffso(irrep_p, sym, pop_p, pop_p, 1) - 1
        do p = 1, nump
            idxE = ioffp(irrep_p) + p
            idx(2) = idx(1) + nump*(p - 1) ! for f(r,p)
            idx(3) = idx(2) + p ! for f(p,p)
            Fock(idx(3)) = fact*Eorb(idxE)
        end do
    end do

end subroutine putDiagFock
