! 2026.10.27
! iterative solve perturbation R equation
! [<1|[H+[V,MP2t], tau1]|0> <1|[V,tau2]|0>][R1] - w^(0)*R1 = [C1]
! [    <2|[V,tau1]|0>          e2*tau2    ][R2] - w^(0)*R2 = [C2]
subroutine solver(iCore, istart, w0, iL)
    use mod_variable
    use mod_common
    use mod_operator
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work spcae
    integer*8, intent(in) :: istart
    real*8, intent(in) :: w0
    integer*8, intent(in) :: iL !> L^(0)

    type(opt) :: pR1, pR2
    real*8 :: diisError
    integer*8 :: i0
    real*8 :: tt0, tt1
    integer*8 :: iter

    i0 = istart
    call copy_opt(R1, pR1, T, iCore, i0)
    i0 = i0 + pR1%len
    call copy_opt(R2, pR2, T, iCore, i0)
    i0 = i0 + pR2%len

    call cpu_time(tt0)
    iter = 0
    call dscal(C1%len, -1.0d0, C1%arr, 1)
    call dscal(C2%len, -1.0d0, C2%arr, 1)
    call dcopy(C1%len, C1%arr, 1, pR1%arr, 1)
    call dcopy(C2%len, C2%arr, 1, pR2%arr, 1)
    do
        iter = iter + 1
        write (6, *) '----------------- iteration:', iter

        !------------ claculate R1 and R2
        call denomT1(R1, pR1, iCore(iEorb))
        call denomT2(R2, pR2, iCore(iEorb))
        call orthogonalR2L(R1%arr, iCore(iL), R1%sym)

        !------------- DIIS
        call dosodiis2(iter, R1%arr, R1%len + R2%len, diis_dim, &
                       iCore(i0), nCore - i0 + 1, diisError)
        call orthogonalR2L(R1%arr, iCore(iL), R1%sym)
        write (6, '(A,F)') 'norm of R1:', R1%norm()
        write (6, '(A,F)') 'norm of R2:', R2%norm()

        !-------- convergence contral
        if (iter > 1 .and. abs(diisError) < diis_max_err) then
            write (6, *) 'DIIS convergenced'
            exit
        end if
        if (iter > diis_cyc) then
            write (6, *) 'WARN: DIIS not convergence'
            exit
        end if

        !-------- calculate pR1
        call dcopy(C1%len, C1%arr, 1, pR1%arr, 1)
        call daxpy(R1%len, -w0, R1%arr, 1, pR1%arr, 1)
        call calA0R1(iCore, i0, pR1, R1, R2)
        if (debug) then
            write (6, '(A,F)') 'norm of pR1:', pR1%norm()
        end if

        !------- calculate pR2
        call dcopy(C2%len, C2%arr, 1, pR2%arr, 1)
        call daxpy(R2%len, -w0, R2%arr, 1, pR2%arr, 1)
        call calA0R2(iCore, i0, pR2, R1)
        if (debug) then
            write (6, '(A,F)') 'norm of pR2:', pR2%norm()
        end if

    end do

    call cpu_time(tt1)
    write (6, '(A,F)') 'total cpu time for iteration:', tt1 - tt0

end subroutine solver
