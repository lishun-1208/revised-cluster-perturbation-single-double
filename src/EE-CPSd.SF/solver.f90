! 2026.10.27
! iterative solve perturbation R equation
! <1|[H,R1]|0> - w^(0)*R1 = C1
subroutine solver(iCore, istart, w0, iL)
    use mod_variable
    use mod_common
    use mod_operator
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work spcae
    integer*8, intent(in) :: istart
    real*8, intent(in) :: w0
    integer*8, intent(in) :: iL !> L^(0)

    type(opt) :: pR1
    real*8 :: diisError
    integer*8 :: i0
    real*8 :: tt0, tt1
    integer*8 :: iter

    i0 = istart
    call copy_opt(R1, pR1, T, iCore, i0)
    i0 = i0 + pR1%len

    call cpu_time(tt0)
    iter = 0
    call dscal(C1%len, -1.0d0, C1%arr, 1)
    call dcopy(C1%len, C1%arr, 1, pR1%arr, 1)
    do
        iter = iter + 1
        write (6, *) '----------------- iteration:', iter

        !------------ claculate R1 and R2
        call denomT1(R1, pR1, iCore(iEorb))
        call orthogonalR2L(R1%arr, iCore(iL), R1%sym)

        !------------- DIIS
        call dosodiis2(iter, R1%arr, R1%len, diis_dim, &
                       iCore(i0), nCore - i0 + 1, diisError)
        call orthogonalR2L(R1%arr, iCore(iL), R1%sym)
        write (6, '(A,F)') 'norm of R1:', R1%norm()

        !-------- convergence contral
        if (iter > 1 .and. abs(diisError) < diis_max_err) then
            write (6, *) 'DIIS convergenced'
            exit
        end if
        if (iter > diis_cyc) then
            write (6, *) 'WARN: DIIS not convergence'
            exit
        end if

        call dcopy(C1%len, C1%arr, 1, pR1%arr, 1)
        call daxpy(R1%len, -w0, R1%arr, 1, pR1%arr, 1)
        ! v(ja,bi)*r(b,j)
        call c1_r1iabj(pR1, R1, iabj, 1.0d0, iCore(i0), nCore - i0 + 1)
        if (debug) then
            write (6, '(A,F)') 'norm of pR1:', pR1%norm()
        end if

    end do

    call cpu_time(tt1)
    write (6, '(A,F)') 'total cpu time for iteration:', tt1 - tt0

end subroutine solver
