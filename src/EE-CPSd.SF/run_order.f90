! write at 2025.1.13
! calculate n-th order vector and (n+1)-th order energy
subroutine run_order(iCore, symm, iroot, nroot, istart, order)
    use mod_variable
    use mod_common, only: iop, debug
    use mod_orbit, only: zorb
    implicit none

    real*8, dimension(*), target, intent(inout) :: iCore
    integer*8, intent(in) :: symm
    integer*8, intent(in) :: iroot, nroot
    integer*8, intent(in) :: istart ! pointer of free iCore
    integer*8, intent(in) :: order !> n for w^(n)

    real*8 :: diisError
    integer*8 :: iter
    real*8 :: w0, wn ! w^(0) and w^(n)
    integer*8 :: i0, i1, iL, nSize, idx, i
    real*8 :: tt0, tt1
    real*8 :: norm, dot1, dot2
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: ddot, dnrm2

    i0 = istart
    call C1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, idx=i0)
    i0 = i0 + C1%len
    call C2%init(iCore, popv, popv, popo, popo, iop, iop, symm, idx=i0)
    i0 = i0 + C2%len
    call R1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, idx=i0)
    i0 = i0 + R1%len
    call R2%init(iCore, popv, popv, popo, popo, iop, iop, symm, idx=i0)
    iL = i0 + R2%len
    i0 = iL + R1%len

    if (order /= 1) then
        !----------------------------------------
        ! calculate constant term C^(n)
        !----------------------------------------
        call calV_C(iCore, symm, iroot, i0, order - 1)
        write (6, '(A,I1,A,F)') 'norm of vector C1^(', order - 1, '):', C1%norm()
        write (6, '(A,I1,A,F)') 'norm of vector C2^(', order - 1, '):', C2%norm()

        !------ read in w^(0)=w^CIS and L^(0)=R^(0)=R^CIS
        idx = (1 + R1%len)*(iroot - 1) + 1
        call IO_stream(w0, iCore(iL), idx, R1%len, f_CIS(symm), '11', 'old')
        if (debug) then
            write (6, '(A,F)') 'CIS excitation energy', w0
        end if

        !-----------------------------------------------------------
        ! solve [J^(0)-w^(0)]*R^(n) = C^(n)
        !-----------------------------------------------------------
        !------ calculate r(ab,ij) = c(ab,ij)/[e(ab,ij)-w^(0)]
        call getallso(iCore(i0), 1, f_OED2(symm))
        !$OMP PARALLEL DO
        do i = 1, C2%len
            R2%arr(i) = C2%arr(i)/(iCore(i0 - 1 + i) - w0)
        end do
        !$OMP END PARALLEL DO

        if (order /= 2) then
            call solver(iCore, i0, w0, iL)
        else ! R1^(1) = 0
            call R1%arr_clear
        end if

        write (6, '(A,I2,A,I4)') '------ symmetry', symm, ' state', iroot
        write (6, '(A,F)') 'orthogonal to L^(0):', ddot(R1%len, iCore(iL), 1, R1%arr, 1)
        write (6, '(A,I1,A,F)') 'norm of R1^(', order - 1, '):', R1%norm()
        if (R1%norm() > 1.0d0) then
            write (6, '(A)') 'WARN'
        end if
        call printEXinfo1_uhf(R1, 10)
        write (6, '(A,I1,A,F)') 'norm of R2^(', order - 1, '):', R2%norm()
        if (R2%norm() > 0.3d0) then
            write (6, '(A)') 'WARN'
        end if
        call printEXinfo2_uhf(R2, 10)
        ! write (*,R^(n-1))
        idx = (1 + R1%len + R2%len)*((order - 1)*nroot + iroot - 1) + 1
        call IO_stream(0.0d0, R1%arr, idx, R1%len + R2%len, f_R(symm), '21', 'old')
    end if

    !------------ calculate w^(n)
    call cpu_time(tt0)
    call calE_C(iCore, i0, symm, iroot, order)
    write (6, '(A,I1,A,F)') 'norm of energy C1^(', order, ')', C1%norm()
    write (6, '(A,I1,A,F)') 'norm of energy C2^(', order, '):', C2%norm()
    wn = ddot(C1%len, iCore(iL), 1, C1%arr, 1)
    call cpu_time(tt1)
    write (6, '(A,F)') 'cpu time for calculate energy:', tt1 - tt0
    write (6, '(A,I2,A,I4,A,I2,A,F)') 'excitation energy for symmetry ', &
        symm, ' state', iroot, ' order', order, ':', wn

    !--- write (e^(n),*)
    idx = (1 + C1%len + C2%len)*((order - 1)*nroot + iroot - 1) + 1
    call IO_stream(wn, R1%arr, idx, 0, f_R(symm), '22', 'old')
    !--- write C^(n)
    idx = (C1%len + C2%len)*(iroot - 1) + 1
    if (order == 1 .and. iroot == 1) then
        call IO_stream(0.0d0, C1%arr, idx, C1%len + C2%len, f_C(symm), '23', 'replace')
    else
        call IO_stream(0.0d0, C1%arr, idx, C1%len + C2%len, f_C(symm), '23', 'old')
    end if
end subroutine run_order
