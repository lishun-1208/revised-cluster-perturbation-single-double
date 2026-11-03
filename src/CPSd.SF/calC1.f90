! write by lishun, 2024.10.18
! C1^(2) = <a,i|[V,T2^(1)]|0>
! C1^(3) = <a,i|[V,T2^(2)]|0>
subroutine calC1(iCore, maxCore, C1, istart, order)
    use mod_variable
    use mod_orbit
    implicit none

    real*8, dimension(maxCore), target, intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    type(opt), intent(inout) :: C1 !> C1^(n)
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> n, the order of deltaT need calculate

    integer*8 :: idx
    integer*8 :: i
    integer*8 :: i0, i1, i2, iEnd, nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart
    ! initialize C1
    call dscal(C1%len, 0.0d0, C1%arr, 1)
    !--------------------------------------
    ! C1^(n) = <a,i|[V,T2^(n-1)]|0>
    !--------------------------------------
    !---- get T2^(n-1), i0
    idx = (T1%len + T2%len)*(order - 2) + C1%len + 1
    call IO_stream(0.0d0, T2%arr, idx, T2%len, f_T, '13', 'old')
    write (6, '(A,I1,A,F)') 'norm of T2^(', order - 1, '):', T2%norm()

    call calC1_VR2(iCore, maxCore, i0, C1, T2, 1, 1.0d0)

end subroutine calC1
