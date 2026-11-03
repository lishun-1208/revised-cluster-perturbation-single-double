! 2026.10.22
! iterative solve perturbation T equation
! [<1|[H+[V,MP2t], tau1]|0> <1|[V,tau2]|0>][T1] = [-C1]
! [    <2|[V,tau1]|0>          e2*tau2    ][T2] = [-C2]
subroutine solver(iCore, istart, order)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work spcae
    integer*8, intent(in) :: istart
    integer*8, intent(in) :: order !> k for T^(k)

    real*8 :: diisError
    integer*8 :: i0
    real*8 :: tt0, tt1
    integer*8 :: iter
    real*8 :: E, dE

    i0 = istart
    iter = 0
    E = 0.0d0
    dE = 0.0d0

    call cpu_time(tt0)
    call dcopy(C1%len, C1%arr, 1, R1%arr, 1)
    call dcopy(C2%len, C2%arr, 1, R2%arr, 1)
    do
        iter = iter + 1
        write (6, *) '----------------- iteration:', iter

        !------------ claculate T1 and T2
        call denomT1(T1, R1, iCore(iEorb))
        call denomT2(T2, R2, iCore(iEorb))

        !------------- DIIS
        call dosodiis2(iter, T1%arr, T1%len + T2%len, diis_dim, &
                       iCore(i0), nCore - i0 + 1, diisError)
        write (6, '(A,F)') 'norm of T1:', T1%norm()
        write (6, '(A,F)') 'norm of T2:', T2%norm()

        !-------- calculater 2 order correlation energy
        dE = E
        call calEcorr(E, iCore, i0, order)
        dE = E - dE
        write (6, '(A,2F)') 'correlation energy:', E, dE

        !-------- convergence contral
        if (iter > 1 .and. abs(diisError) < diis_max_err) then
            write (6, *) 'DIIS convergenced'
            Ecorr = Ecorr + E
            write (6, *) 'total correlation energy:', Ecorr
            exit
        end if
        if (iter > diis_cyc) then
            write (6, *) 'WARN: DIIS not convergence'
            Ecorr = Ecorr + E
            write (6, *) 'total correlation energy:', Ecorr, E
            exit
        end if

        !-------- calculate R1
        call dcopy(C1%len, C1%arr, 1, R1%arr, 1)
        call calA0R1(iCore, i0)
        write (6, *) 'norm of R1:', R1%norm()

        !------- calculate R2
        call dcopy(C2%len, C2%arr, 1, R2%arr, 1)
        call calA0R2(iCore, i0)
        write (6, '(A,F)') 'norm of R2:', R2%norm()

    end do

    call cpu_time(tt1)
    write (6, '(A,F)') 'total cpu time for iteration:', tt1 - tt0

end subroutine solver
