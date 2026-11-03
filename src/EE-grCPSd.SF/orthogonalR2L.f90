! write by lishun, 2024.11.25
! orthogonalize R to L
! |R> = |R>-(<L|R>/<L|L>)*|L>
! L1,L2 and R1,R2 need be continue
subroutine orthogonalR2L(R, L, symm)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: R
    real*8, dimension(*), intent(in) :: L
    integer*8, intent(in) :: symm

    real*8 :: dot, norm
    integer*8 :: nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: ddot, dnrm2

    nSize = irpdso(symm, popv, popo, 1) + idsymsoc(symm, popv, popv, iop, popo, popo, iop)
    ! <L|L>
    norm = dnrm2(nSize, L, 1)
    ! <L|R>
    dot = ddot(nSize, L, 1, R, 1)

    ! |R> = |R>-(<L|R>/<L|L>)*|L>
    call daxpy(nSize, -dot/norm, L, 1, R, 1)

    ! test <L|R>
    dot = ddot(nSize, L, 1, R, 1)
    if (abs(dot) > 1.0d-8) then
        write (6, *) 'WARN: <L^(0)|R^(n)> not be 0', dot
    end if

end subroutine orthogonalR2L
