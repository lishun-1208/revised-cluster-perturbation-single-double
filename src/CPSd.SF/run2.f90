! write by lishun, 2024.12.5
! calculate second order perturbation T
! A^(0)*T^(2)+A^(1)*T^(1)=0
! seperate OED term:
! -e(a,i)*T1^(2) = <a,i|[V,T1^(2)]|0>+<a,i|[V,T2^(1)]|0>
! -e(ab,ij)*T2^(2) = <ab,ij|[V,T2^(1)]|0>
! E^(3) = 0.25*<ij||ab>*t(ab,ij)^(2), T1^(1)=0
subroutine run2(iCore, maxCore, max_cyc, maxError, maxDim)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), target, intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: max_cyc
    real*8, intent(in) :: maxError ! max DIIS error
    integer*8, intent(in) :: maxDim ! max DIIS dim

    integer*8 :: iter
    integer*8 :: i0, nSize
    real*8 :: Eorder ! E^(3)
    real*8 :: diisError
    real*8 :: tt0, tt1
    integer*8 :: idx

    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = iFree

    !------------------- calculate constant term
    ! c(a,i) = <a,i|[V,T2^(1)]|0>
    call calC1(iCore, maxCore, C1, i0, 2)
    write (6, *) 'norm of C1^(2):', C1%norm()
    ! c(ab,ij) = <ab,ij|[V,T^(1)]|0>
    call calC2(iCore, maxCore, C2, i0, 2)
    write (6, '(A,I2,A,F)') 'norm of C2^(2) with iop', iop, ':', C2%norm()

    !-------------- calculate T2
    ! -e(ab,ij)*T2^(2) = c(ab,ij)
    call denomT2(T2, C2, iCore(iEorb))
    write (6, *) 'norm of T2^(2):', T2%norm()

    !-------------------- solver T1 equation
    ! -e(a,i)*T1^(2) = <a,i|[V,T1^(2)]|0>+c(a,i)
    iter = 0
    ! initialize R1
    call dcopy(C1%len, C1%arr, 1, R1%arr, 1)
    do
        call cpu_time(tt0)
        iter = iter + 1
        write (6, *) '----------------- 2nd order iteration:', iter

        !------------ claculate T1
        call denomT1(T1, R1, iCore(iEorb))

        !------------- DIIS
        call dosodiis2(iter, T1%arr, T1%len, maxDim, iCore(i0), nCore - i0 + 1, diisError)
        write (6, *) 'norm of T1^(2):', T1%norm()

        !-------- convergence contral
        if (iter > 1 .and. abs(diisError) < maxError) then
            write (6, *) 'CPS(D) convergenced for 2nd order'
            exit
        end if
        if (iter > max_cyc) then
            write (6, *) 'CPS(D) not convergence for 2nd order', diisError
            exit
        end if

        !-------- calculate R1=<a,i|[V,T1^(2)]|0>+c(a,i)
        call dcopy(C1%len, C1%arr, 1, R1%arr, 1)
        call c1_r1iabj(R1, T1, iabj, 1.0d0, iCore(i0), nCore - i0 + 1)

        call cpu_time(tt1)
        write (6, '(A,I5,A,F)') 'cpu time for iteration', iter, ' :', tt1 - tt0
    end do

    !---------- calculate E^(3)
    call calEcorr(Eorder, iCore, maxCore, iFree, T2, 2)
    write (6, '(A,F)') '3rd order correlation energy', Eorder
    Ecorr = Ecorr + Eorder

    !-------- write in hard disk
    nSize = T1%len + T2%len
    idx = nSize + 1
    call IO_stream(0.0d0, T1%arr, idx, nSize, f_T, '23', 'old')

    ! print excitation information
    write (6, '(A)') '------ 2nd order T ------'
    write (6, '(A,F)') 'norm of T1:', T1%norm()
    call printEXinfo1_uhf(T1, 10)
    write (6, '(A,F)') 'norm of T2:', T2%norm()
    call printEXinfo2_uhf(T2, 10)

end subroutine run2
