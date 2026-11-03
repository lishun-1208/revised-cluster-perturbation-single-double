! Build the right-vector constant term from the constant vector already
! calculated with the excitation energy:
! Ctilde^(n) = sum_{k=1}^n w^(k) R^(n-k) - C^(n).
!> C1,C2 need be continue
subroutine calV_C(iCore, symm, iroot, istart, order)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: symm
    integer*8, intent(in) :: iroot
    integer*8, intent(in) :: istart
    integer*8, intent(in) :: order !> n

    integer*8 :: i0, idx, nroot
    real*8 :: wk
    integer*8 :: k

    nroot = estate_sym(symm)
    i0 = istart

    ! Ctilde^(n) = -C^(n)
    idx = (C1%len + C2%len)*(iroot - 1) + 1
    call IO_stream(0.0d0, C1%arr, idx, C1%len + C2%len, f_C(symm), '13', 'old')
    call dscal(C1%len, -1.0d0, C1%arr, 1)
    call dscal(C2%len, -1.0d0, C2%arr, 1)

    ! Ctilde^(n) += sum_{k=1}^n w^(k) R^(n-k), w^(1)=0
    do k = 2, order
        ! w^(k)
        idx = (1 + C1%len + C2%len)*(nroot*(k - 1) + iroot - 1) + 1
        call IO_stream(wk, iCore(i0), idx, 0, f_R(symm), '12', 'old')
        ! R^(n-k)
        idx = (1 + C1%len + C2%len)*(nroot*(order - k) + iroot - 1) + 2
        call IO_stream(0.0d0, iCore(i0), idx, C1%len + C2%len, f_R(symm), '13', 'old')

        call daxpy(C1%len + C2%len, wk, iCore(i0), 1, C1%arr, 1)
    end do

end subroutine calV_C
