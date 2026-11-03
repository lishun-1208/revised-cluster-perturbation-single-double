! write by lishun, 2024.10.18
! calculate second order perturbation T
! A^(0)*T^(2)+A^(1)*T^(1)=0, seperate OED term:
! -e*T^(2)=(A^(0)-e)*T^(2)+A^(1)*T^(1)
! -e*T^(2)=(A^(0)-e)*T^(2)+C
subroutine run2(iCore, istart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0
    real*8 :: tt0, tt1

    i0 = istart

    !------------------- calculate constant term
    ! C(a,i) = -[(Wkaci-V(ka,ci))*T(c,k)^(1)+Fli*T(a,l)^(1)-Fad*T(d,i)^(1)]
    ! C(ab,ij) = PijPab{V(kb,cj)*t(ac,ik)^(1)}+0.5*V(kl,ij)*t(ab,kl)^(1)+0.5*V(ab,cd)*t(cd,ij)^(1)
    call cpu_time(tt0)
    call calC1(iCore, C1, i0, 2)
    write (6, '(A,F)') 'norm of C1^(2):', C1%norm()
    call calC2(iCore, C2, i0, 2)
    write (6, '(A,F)') 'norm of C2^(2):', C2%norm()
    call cpu_time(tt1)
    write (6, *) 'cpu time for constant term:', tt1 - tt0

    ! solve equation
    call solver(iCore, i0, 2)

    !-------- write in hard disk
    call IO_stream(0.0d0, T1%arr, T1%len + T2%len + 1, &
                   T1%len + T2%len, f_T, '23', 'old')

    ! print excitation information
    write (6, '(A,F)') 'norm of T1^(2):', T1%norm()
    call printEXinfo1_uhf(T1, 10)
    write (6, '(A,F)') 'norm of T2^(2):', T2%norm()
    call printEXinfo2_uhf(T2, 10)

end subroutine run2
