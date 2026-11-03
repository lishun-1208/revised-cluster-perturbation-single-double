! write by lishun, 2024.10.21
! calculate second order perturbation T
! A^(0)*T^(3)+A^(1)*T^(2)+0.5*<mu|[[V,T^(1)],T^(1)]=0, seperate OED term:
! -eT(3) = (A^(0)-e)*T^(3)+A^(1)*T^(2)+0.5*<mu|[[V,T^(1)],T^(1)]
! -eT(3) = (A^(0)-e)*T^(3)+C
subroutine run3(iCore, istart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0
    real*8 :: tt0, tt1

    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart

    !------------------- calculate constant term
    !---- C = A^(1)*T^(2)+0.5*<mu|[[V,T*(1)],T^(1)]
    call cpu_time(tt0)
    call calC1(iCore, C1, i0, 3)
    write (6, '(A,F)') 'norm of C1^(3):', C1%norm()
    call calC2(iCore, C2, i0, 3)
    write (6, '(A,F)') 'norm of C2^(3):', C2%norm()
    call cpu_time(tt1)
    write (6, *) 'cpu time for constant term:', tt1 - tt0

    ! solve equation
    call solver(iCore, i0, 3)

    !-------- write in hard disk
    call IO_stream(0.0d0, T1%arr, (T1%len + T2%len)*2 + 1, &
                   T1%len + T2%len, f_T, '23', 'old')

    ! print excitation information
    write (6, '(A,F)') 'norm of T1^(3):', T1%norm()
    call printEXinfo1_uhf(T1, 10)
    write (6, '(A,F)') 'norm of T2^(3):', T2%norm()
    call printEXinfo2_uhf(T2, 10)

end subroutine run3
