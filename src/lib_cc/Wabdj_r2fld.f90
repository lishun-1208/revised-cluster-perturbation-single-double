! 2026.7.10
!> W(ab,dj) -= fact*F(d,l)*r(ab,lj)
subroutine Wabdj_r2fld(Wabdj, r2, Fld, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabdj
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Fld
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wabdj%incore) then
        call vv2v_424(r2, Fld, Wabdj, -fact, 'ablj', 'dl', 'abdj', iCore, maxCore, 0)
    else
        call p_Wabdj_r2fld_h(Wabdj%f, r2%arr, Fld%arr, Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                             Fld%popq, Wabdj%ioppq, r2%sym, Fld%sym, &
                             fact, iCore, maxCore)
    end if

end subroutine Wabdj_r2fld

! Wabdj -= fact*Fld*r(ab,lj)
! Fld store as (d,l)
subroutine p_Wabdj_r2fld_h(file_Wabdj, r2, Fld, pop_a, pop_b, pop_d, pop_j, pop_l, &
                           iop, symR, symF, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    character(len=8), intent(in) :: file_Wabdj
    real*8, dimension(*), intent(in) :: r2, Fld
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_l
    integer*8, intent(in) :: iop ! iop=iopab=ioplj
    integer*8, intent(in) :: symR, symF
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symW, symm, symm2
    integer*8 :: irrepj, numj, j
    integer*8 :: i0, i1, i2, i3, iEnd, nSize
    real*8 :: t0, t1
    integer*8, external :: idsymsoc, irpdso

    symW = dirprd(symR, symF)
    iEnd = 0
    i0 = 1
    do irrepj = 1, nirrep
        numj = pop_j(irrepj)
        symm = dirprd(symR, irrepj) ! symmetry of r(ab,l*)
        symm2 = dirprd(symW, irrepj) ! symmetry of W(ab,d*)
        do j = 1, numj
            ! get r(ab,l*), i0
            nSize = idsymsoc(symm, pop_a, pop_b, iop, pop_l, zorb, 1)
            i1 = i0 + nSize
            call getpqrso_incore(iCore(i0), j, irrepj, &
                                 pop_a, pop_b, iop, pop_l, pop_j, iop, r2, symR)
            ! tmp(ab,*d)=r(ab,*l)*F(d,l), i1
            nSize = idsymsoc(symm2, pop_a, pop_b, iop, pop_d, zorb, 1)
            i2 = i1 + nSize
            call VpqrmVsm_to_Vpqrs(iCore(i0), Fld, symm, symF, pop_l, &
                                   pop_a, pop_b, iop, zorb, pop_d, iCore(i1), symm2, 1.0d0, 0.0d0)
            ! get W(ab,d*), i2
            i3 = i2 + nSize
            iEnd = max(iEnd, i3)
            call getpqrso(iCore(i2), j, irrepj, 1, file_Wabdj)
            ! W(ab,d*) -= fact*tmp(ab,d*), i2
            call daxpy(nSize, -fact, iCore(i1), 1, iCore(i2), 1)
            ! put W(ab,d*)
            call putpqrso(iCore(i2), j, irrepj, 1, file_Wabdj)
        end do
    end do

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for p_Wabdj_r2fld_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wabdj_r2fld_h

! 2026.7.10
!> W(ab,dj) -= fact*F(d,l)*r(ab,lj)
subroutine Wabdj_r2fld_R(WabdjABAB, r2ABAB, FldAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: FldAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,dJ) = -F(d,l)*r(aB,lJ)
    call Wabdj_r2fld(WabdjABAB, r2ABAB, FldAA, fact, iCore, maxCore)

end subroutine Wabdj_r2fld_R

! 2026.7.10
!> W(ab,dj) -= fact*F(d,l)*r(ab,lj)
subroutine Wabdj_r2fld_U(WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA, spin, &
                         r2AAAA, r2ABAB, r2BBBB, FldAA, FldBB, &
                         fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: FldAA, FldBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ab,dj) = -F(d,l)*r(ab,lj)
        call Wabdj_r2fld(WabdjAAAA, r2AAAA, FldAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(aB,dJ) = -F(d,l)*r(aB,lJ)
        call Wabdj_r2fld(WabdjABAB, r2ABAB, FldAA, fact, iCore, maxCore)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(AB,DJ) = -F(D,L)*r(AB,LJ)
        call Wabdj_r2fld(WabdjBBBB, r2BBBB, FldBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ab,Dj) = -F(D,L)*r(Ab,Lj) = -F(D,L)*r(bA,jL)
        ! r(bA,jL)->r(Ab,Lk), i0
        i0 = 1
        call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                      1, 1, r2ABAB%sym, T, i0)
        call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                      r2ABAB%popr, r2ABAB%pops, 1, 1, r2ABAB%sym)
        call Wabdj_r2fld(WabdjBABA, tmp, FldBB, fact, iCore, maxCore)
    end if

end subroutine Wabdj_r2fld_U
