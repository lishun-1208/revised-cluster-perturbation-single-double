! write by lishun, 2024.10.15
! calculate first order perturbation T
! A^(0)*T^(1)+<mu|H|0>=0, seperate OED term:
! -e*T^(1)=(A^(0)-e)*T^(1)+<mu|H|0>
subroutine run1(iCore, istart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0

    i0 = istart

    ! calculate C
    call C1%arr_clear
    call dcopy(C2%len, abij%arr, 1, C2%arr, 1)

    ! solve equation
    call solver(iCore, i0, 1)

    !-------- write in hard disk
    call IO_stream(0.0d0, T1%arr, 1, T1%len + T2%len, f_T, '23', 'replace')

    ! print excitation information
    write (6, '(A,F)') 'norm of T1^(1):', T1%norm()
    call printEXinfo1_uhf(T1, 10)
    write (6, '(A,F)') 'norm of T2^(1):', T2%norm()
    call printEXinfo2_uhf(T2, 10)

end subroutine run1
