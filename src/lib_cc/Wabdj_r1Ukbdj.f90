! 2026.7.7
!> W(ab,dj) -= fact*P(ab){Ukbdj(dk,bj)*r(a,k)}
subroutine Wabdj_r1Ukbdj(Wabdj, r1, Ukbdj, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabdj
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: Ukbdj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wabdj%incore) then
        call vv2v_424(Ukbdj, r1, Wabdj, -fact, 'dkbj', 'ak', 'abdj', iCore, maxCore, 1)
    else
        call p_Wabdj_r1Ukbdj_h(Wabdj%f, r1%arr, Ukbdj%arr, Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                               r1%popq, r1%sym, Ukbdj%sym, Wabdj%ioppq, fact, iCore, maxCore)
    end if

end subroutine Wabdj_r1Ukbdj

! write at 2025.2.26
! Wabdj -= fact*Pab{Ukbdj*r(a,k)}
!       -= fact*Pab{U(dk,bj)*r(a,k)}
! not need U(dk,bj)=U(bj,dk)
subroutine p_Wabdj_r1Ukbdj_h(file_Wabdj, r1, Ukbdj, pop_a, pop_b, pop_d, pop_j, pop_k, &
                             symR, symU, iop, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    character(len=8), intent(in) :: file_Wabdj
    real*8, dimension(*), intent(in) :: r1, Ukbdj
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_k
    integer*8, intent(in) :: symR, symU ! symmetry of r1 and Ukbdj
    integer*8, intent(in) :: iop !> iop=iopab
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symW
    integer*8 :: irrep_j, num_j, j
    integer*8 :: irrep_abd
    integer*8 :: irrep_kbd
    integer*8 :: i0, i1, i2, i3, i4, i5, iEnd, nSize
    real*8 :: t0, t1
    integer*8, external :: idsymsoc, irpdso

    symW = dirprd(symR, symU)
    iEnd = 0
    i0 = 1

    ! U(dk,bj)->U(kb,dj), i0
    nSize = idsymsoc(symU, pop_d, pop_k, 1, pop_b, pop_j, 1)
    i1 = i0 + nSize
    call sstgenso(Ukbdj, iCore(i0), nSize, pop_d, pop_k, pop_b, pop_j, iCore(i1), symU, '2314')
    ! r(a,k)->r(k,a), i1
    nSize = irpdso(symR, pop_a, pop_k, 1)
    i2 = i1 + nSize
    call itranspso(r1, iCore(i1), zorb, pop_a, 1, zorb, pop_k, 1, symR)

    do irrep_j = 1, nirrep
        num_j = pop_j(irrep_j)
        irrep_abd = dirprd(irrep_j, symW)
        irrep_kbd = dirprd(irrep_j, symU)
        do j = 1, num_j
            ! get U(kb,d*), i2
            nSize = idsymsoc(irrep_kbd, pop_k, pop_b, 1, pop_d, zorb, 1)
            i3 = i2 + nSize
            call getpqrso_incore(iCore(i2), j, irrep_j, &
                                 pop_k, pop_b, 1, pop_d, pop_j, 1, iCore(i0), symU)
            ! tmp(ab,d*) = U(kb,d*)*r(k,a), i3
            nSize = idsymsoc(irrep_abd, pop_a, pop_b, 1, pop_d, zorb, 1)
            i4 = i3 + nSize
            call VmqrsVmp2Vpqrs(iCore(i2), iCore(i1), irrep_kbd, symR, pop_k, &
                                pop_a, pop_b, pop_d, zorb, 1, iCore(i3), irrep_abd, 1.0d0, 0.0d0)
            ! tmp(ab,d*)->tmp(ba,d*), i4
            i5 = i4 + nSize
            call isymtrso(iCore(i3), iCore(i4), &
                          pop_a, pop_b, 1, 1, pop_d, zorb, 1, 0, irrep_abd)
            ! tmp2(ab,d*) = tmp(ab,d*)-tmp(ba,d*), i3
            call daxpy(nSize, -1.0d0, iCore(i4), 1, iCore(i3), 1)
            if (iop == 0) then
                ! compress ab, i4
                call icompso(iCore(i3), iCore(i4), &
                             pop_a, pop_b, 1, 0, pop_d, zorb, 1, 1, irrep_abd)
                ! get W(ab,d*), i3
                call getpqrso(iCore(i3), j, irrep_j, 1, file_Wabdj)
                ! W(ab,d*) -= fact*tmp2(ab,d*)
                nSize = idsymsoc(irrep_abd, pop_a, pop_b, 0, pop_d, zorb, 1)
                call daxpy(nSize, -fact, iCore(i4), 1, iCore(i3), 1)
                ! write W(ab,d*)
                call putpqrso(iCore(i3), j, irrep_j, 1, file_Wabdj)
            else
                ! get W(ab,d*), i4
                call getpqrso(iCore(i4), j, irrep_j, 1, file_Wabdj)
                ! W(ab,d*) -= fact*tmp2(ab,d*)
                nSize = idsymsoc(irrep_abd, pop_a, pop_b, 1, pop_d, zorb, 1)
                call daxpy(nSize, -fact, iCore(i3), 1, iCore(i4), 1)
                ! write W(ab,d*)
                call putpqrso(iCore(i4), j, irrep_j, 1, file_Wabdj)
            end if
        end do
        iEnd = max(iEnd, i5)
    end do

    if ((iEnd - maxCore) >= 0) then
        write (6, *) 'ERROR: memory not enough for p_Wabdj_r1Ukbdj_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wabdj_r1Ukbdj_h

! 2026.7.7
!> W(ab,dj) -= fact*P(ab){Ukbdj(dk,bj)*r(a,k)}
subroutine Wabdj_r1Ukbdj_R(WabdjABAB, r1AA, UkbdjABAB, UkbdjBAAB, &
                           fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: UkbdjABAB, UkbdjBAAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,dJ) = -UkBdJ(dk,BJ)*r(a,k)
    !            +UKadJ(dK,aJ)*r(B,K) -> UKadJ(dK,aJ)*r(b,k)
    if (WabdjABAB%incore) then
        call vv2v_424(UkbdjABAB, r1AA, WabdjABAB, -fact, 'dkBJ', 'ak', 'aBdJ', iCore, maxCore, 0)
        call vv2v_424(UkbdjBAAB, r1AA, WabdjABAB, fact, 'dkaj', 'bk', 'abdj', iCore, maxCore, 0)
    else
        write (6, '(A)') 'ERROR: Wabdj_r1Ukbdj_R current not support outcore'
        stop
    end if

end subroutine Wabdj_r1Ukbdj_R

! 2026.7.7
!> W(ab,dj) -= fact*P(ab){Ukbdj(dk,bj)*r(a,k)}
subroutine Wabdj_r1Ukbdj_U(WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA, spin, &
                           r1AA, r1BB, UkbdjAAAA, UkbdjABAB, UkbdjBAAB, &
                           UkbdjBBBB, UkbdjBABA, UkbdjABBA, &
                           fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: UkbdjAAAA, UkbdjABAB, UkbdjBAAB
    type(opt), intent(in) :: UkbdjBBBB, UkbdjBABA, UkbdjABBA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ab,dj) = -P(ab){Ukbdj(dk,bj)*r(a,k)}
        call Wabdj_r1Ukbdj(WabdjAAAA, r1AA, UkbdjAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(aB,dJ) = -UkBdJ(dk,BJ)*r(a,k)
        !            +UKadJ(dK,aJ)*r(B,K)
        if (WabdjABAB%incore) then
            call vv2v_424(UkbdjABAB, r1AA, WabdjABAB, -fact, 'dkBJ', 'ak', 'aBdJ', iCore, maxCore, 0)
            call vv2v_424(UkbdjBAAB, r1BB, WabdjABAB, fact, 'dKaJ', 'BK', 'aBdJ', iCore, maxCore, 0)
        else
            write (6, '(A)') 'ERROR: Wabdj_r1Ukbdj_U current not support outcore'
            stop
        end if
    end if
    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(AB,DJ) = -P(AB){UKBDJ(DK,BJ)*r(A,K)}
        call Wabdj_r1Ukbdj(WabdjBBBB, r1BB, UkbdjBBBB, fact, iCore, maxCore)
    end if
    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ab,Dj) = -UKbDj(DK,bj)*r(A,K)
        !            +UkADj(Dk,Aj)*r(b,k)
        if (WabdjBABA%incore) then
            call vv2v_424(UkbdjBABA, r1BB, WabdjBABA, -fact, 'DKbj', 'AK', 'AbDj', iCore, maxCore, 0)
            call vv2v_424(UkbdjABBA, r1AA, WabdjBABA, fact, 'DkAj', 'bk', 'AbDj', iCore, maxCore, 0)
        else
            write (6, '(A)') 'ERROR: Wabdj_r1Ukbdj_U current not support outcore'
            stop
        end if
    end if

end subroutine Wabdj_r1Ukbdj_U
