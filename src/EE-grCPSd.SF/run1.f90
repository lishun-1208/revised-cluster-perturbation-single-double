! write by lishun, 2024.11.11
! calculate 1st order vector and 2nd order energy
subroutine run1(iCore, symm, iroot, nroot, istart)
    use mod_common
    use mod_variable
    use mod_orbit
    use mod_symm
    use mod_operator
    implicit none

    real*8, dimension(*), target, intent(inout) :: iCore
    integer*8, intent(in) :: symm
    integer*8, intent(in) :: iroot, nroot
    integer*8, intent(in) :: istart ! pointer of free iCore

    real*8 :: w0, w2
    integer*8 :: iL
    integer*8 :: i0, i1, nSize
    integer*8 :: idx
    real*8 :: tt0, tt1
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: ddot, dnrm2

    i0 = istart

    ! Bind state vectors after their workspace offsets are known
    call C1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, idx=i0)
    i0 = i0 + C1%len
    call C2%init(iCore, popv, popv, popo, popo, iop, iop, symm, idx=i0)
    i0 = i0 + C2%len
    call R1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, idx=i0)
    i0 = i0 + R1%len
    call R2%init(iCore, popv, popv, popo, popo, iop, iop, symm, idx=i0)
    i0 = i0 + R2%len
    iL = i0
    i0 = iL + R1%len + R2%len

    !----------- C = [w^(1)-J^(1)]*R^(0), J^(1) = <d|[V,\tau_d]|0>
    call calV_C(iCore, symm, iroot, i0, 1)
    write (6, '(A,F)') 'norm of vector C1^(1)', C1%norm()
    write (6, '(A,F)') 'norm of vector C2^(1):', C2%norm()

    !-------------- solve [J^(0)-w^(0)]*R = C
    !------- read in w^(0) and L^(0)
    idx = (1 + R1%len + R2%len)*(iroot - 1) + 1
    call IO_stream(w0, iCore(iL), idx, R1%len + R2%len, f_L0(symm), '11', 'old')

    call solver(iCore, i0, w0, iL)

    write (6, '(A,I2,A,I4)') '------ symmetry', symm, ' state', iroot
    write (6, '(A,F)') 'orthogonal to L^(0):', ddot(R1%len + R2%len, iCore(iL), 1, R1%arr, 1)
    write (6, '(A,F)') 'norm of R1^(1):', R1%norm()
    if (R1%norm() > 1.0d0) then
        write (6, '(A)') 'WARN'
    end if
    call printEXinfo1_uhf(R1, 10)
    write (6, '(A,F)') 'norm of R2^(1):', R2%norm()
    if (R2%norm() > 0.3d0) then
        write (6, '(A)') 'WARN'
    end if
    call printEXinfo2_uhf(R2, 10)
    ! write (*,R^(1))
    idx = (1 + R1%len + R2%len)*(nroot + iroot - 1) + 1
    call IO_stream(0.0d0, R1%arr, idx, R1%len + R2%len, f_R(symm), '21', 'old')

    !------------ calculate w^(2)
    ! Reuse the constant-term workspace to save C^(2) for the next R equation.
    call cpu_time(tt0)
    call calE_C(iCore, i0, symm, iroot, 2)
    write (6, '(A,F)') 'norm of energy C1^(2)', C1%norm()
    write (6, '(A,F)') 'norm of energy C2^(2):', C2%norm()
    w2 = ddot(C1%len + C2%len, iCore(iL), 1, C1%arr, 1)
    call cpu_time(tt1)
    write (6, '(A,F)') 'cpu time for calculate energy:', tt1 - tt0
    write (6, '(A,I2,A,I3,A,F)') '2nd order excitation energy for symmetry ', &
        symm, ' state', iroot, ':', w2

    !--- write (e^(2),*)
    idx = (1 + C1%len + C2%len)*(nroot + iroot - 1) + 1
    call IO_stream(w2, R1%arr, idx, 0, f_R(symm), '22', 'old')
    !--- write C^(2)
    idx = (C1%len + C2%len)*(iroot - 1) + 1
    call IO_stream(0.0d0, C1%arr, idx, C1%len + C2%len, f_C(symm), '23', 'old')
end subroutine run1
