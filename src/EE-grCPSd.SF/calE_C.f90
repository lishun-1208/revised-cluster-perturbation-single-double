! 2026.10.22
! calculate C^(n) = sum_(k=1,n){J^(k)*R(n-k)}
subroutine calE_C(iCore, istart, symR, iroot, order)
    use mod_variable
    use mod_operator
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    integer*8, intent(in) :: symR !> symmetry of state
    integer*8, intent(in) :: iroot
    integer*8, intent(in) :: order !> n for w^(n)

    type(opt) :: pR1, pR2
    type(opt) :: pT1, pT2
    integer*8 :: i0, idx
    integer*8 :: nroot

    i0 = istart
    call copy_opt(C1, pR1, T, iCore, i0)
    i0 = i0 + pR1%len
    call copy_opt(C2, pR2, T, iCore, i0)
    i0 = i0 + pR2%len
    call pT1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, idx=i0)
    i0 = i0 + pT1%len
    call pT2%init(iCore, popv, popv, popo, popo, iop, iop, 1, idx=i0)
    i0 = i0 + pT2%len
    nroot = estate_sym(symR)

    call C1%arr_clear
    call C2%arr_clear

    !------ C2^(n) = J^(1)*R^(n-1)
    idx = (1 + C1%len + C2%len)*(nroot*(order - 1) + iroot - 1) + (1 + C1%len) + 1
    call IO_stream(0.0d0, pR2%arr, idx, C2%len, f_R(symR), '13', 'old')
    call calC2_VR2(iCore, C2, pR2, symR, i0, 1.0d0)

    if (order == 1) then
        return
    end if

    !------- C2^(n) += J^(2)*R^(n-2)
    idx = (1 + C1%len + C2%len)*(nroot*(order - 2) + iroot - 1) + 1 + 1
    call IO_stream(0.0d0, pR1%arr, idx, C1%len + C2%len, f_R(symR), '13', 'old')
    call calC1_VMP2R1(iCore, C1, pR1, i0, 1.0d0)
    call calC_J2R(iCore, i0, C1, C2, pR1, pR2, 1.0d0)
    !--- <u|[[V,T^(1)],R]|0>
    call IO_stream(0.0d0, pT1%arr, 1, pT1%len + pT2%len, f_T, '13', 'old')
    call calC_VTR(iCore, i0, C1, C2, pT1, pT2, pR1, pR2, symR, 1.0d0)

    if (order == 2) then
        return
    end if

    !------- C2^(n) += J^(3)*R^(n-3)
    idx = (1 + C1%len + C2%len)*(nroot*(order - 3) + iroot - 1) + 1 + 1
    call IO_stream(0.0d0, pR1%arr, idx, C1%len + C2%len, f_R(symR), '13', 'old')
    call calC_J3R(iCore, i0, C1, C2, pR1, pR2, 1.0d0)
    !--- <u|[[V,T^(2)],R]|0>
    idx = pT1%len + pT2%len + 1
    call IO_stream(0.0d0, pT1%arr, idx, pT1%len + pT2%len, f_T, '13', 'old')
    call calC_VTR(iCore, i0, C1, C2, pT1, pT2, pR1, pR2, symR, 1.0d0)
    !--- <u|[[[V,T^(1)],T^(1)],R]|0>
    call IO_stream(0.0d0, pT1%arr, 1, pT1%len + pT2%len, f_T, '13', 'old')
    call calC_VTTR(iCore, i0, C2, pT1, pT2, pR1, pR2, symR, 1.0d0)

    if (order == 3) then
        return
    end if

end subroutine calE_C
