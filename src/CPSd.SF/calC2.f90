! write by lishun, 2024.12.5
! C2^(2) = <ab,ij|[V,T2^(1)]|0>
! C2^(3) = <ab,ij|[V,T^(2)]|0>+0.5*<ab,ij|[[V,T2^(1)],T2^(1)]|0>
subroutine calC2(iCore, maxCore, C2, istart, order)
    use mod_variable
    use mod_orbit, only: zorb
    implicit none

    real*8, dimension(maxCore), target, intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    type(opt), intent(inout) :: C2 ! C2^(n)
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order ! n

    integer*8 :: idx
    integer*8 :: i0, i1, i2, iEnd, nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart

    !----- initialize C2
    call dscal(C2%len, 0.0d0, C2%arr, 1)

    !-----------------------------------------
    ! C2^(n) = <ab,ij|[V,T^(n-1)]|0>
    !-----------------------------------------
    !---------- get T^(n-1)
    idx = (T1%len + T2%len)*(order - 2) + 1
    call IO_stream(0.0d0, T1%arr, idx, T1%len + T2%len, f_T, '13', 'old')
    write (6, '(A,I1,A,F)') 'norm of T1^(', order - 1, '):', T1%norm()
    write (6, '(A,I1,A,F)') 'norm of T2^(', order - 1, '):', T2%norm()

    call calC2_VR(iCore, maxCore, i0, C2, T1, T2, 1, 1.0d0, order)

    if (order == 3) then
        !----------------------------------------
        ! 0.5*<ab,ij|[[V,T^(1)],T^(1)]|0>
        !----------------------------------------
        !---------- get T^(1)
        ! T1
        call IO_stream(0.0d0, T1%arr, 1, T1%len+T2%len, f_T, '13', 'old')
        write (6, '(A,F)') 'norm of T1^(1):', T1%norm()
        write (6, '(A,F)') 'norm of T2^(1):', T2%norm()

        !------ T1^(1)=0
        call calC2_VTT(iCore, maxCore, i0, C2, T1, T2, 1, 1.0d0)
    end if

end subroutine calC2
