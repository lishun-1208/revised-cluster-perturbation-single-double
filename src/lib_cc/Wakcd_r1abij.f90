! 2026.7.10
!> W(cd,ak) -= fact*V(cd,lk)*r(a,l)
subroutine Wakcd_r1abij(Wakcd, r1, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wakcd
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wakcd%incore) then
        call vv2v_424(abij, r1, Wakcd, -fact, 'cdlk', 'al', 'cdak', iCore, maxCore, 0)
    else
        call p_Wakcd_r1abij_h(Wakcd%f, r1%arr, abij%arr, Wakcd%popr, Wakcd%pops, Wakcd%popp, Wakcd%popq, &
                              r1%popq, Wakcd%ioppq, r1%sym, fact, iCore, maxCore)
    end if

end subroutine Wakcd_r1abij

! write at 2025.7.1, lishun
! calculate Wakcd += -fact*V(lk,cd)*t(a,l)
! Wakcd sort as (cd,ak)
!-------------------------------
! modified at 2025.7.2
subroutine p_Wakcd_r1abij_h(file_Wakcd, r1, abij, pop_a, pop_k, pop_c, pop_d, &
                            pop_l, iop, symR, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    character(len=8), intent(in) :: file_Wakcd
    real*8, dimension(*), intent(in) :: r1, abij
    integer*8, dimension(8), intent(in) :: pop_a, pop_k, pop_c, pop_d, pop_l ! number of orb for each symmetry
    integer*8, intent(in) :: iop !> iop=0 for c<d, 1 for (c,d), iop=iopcd=iopab=iopij
    integer*8, intent(in) :: symR ! symmetry of r1
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1, i2, nSize, iEnd
    integer*8 :: irrepk, numk, k
    integer*8 :: sym, symW
    integer*8, external :: irpdso, idsymsoc

    symW = symR ! symmetry of Wakcd
    iEnd = 0
    i0 = 1
    ! -V(lk,cd)*t(a,l)=-V(cd,lk)*t(a,l)
    do irrepk = 1, nirrep
        sym = dirprd(irrepk, symW) ! symmetry of W(cd,a*)
        numk = pop_k(irrepk)
        do k = 1, numk
            ! get V(cd,l*), i0
            nSize = idsymsoc(irrepk, pop_c, pop_d, iop, pop_l, zorb, 1)
            i1 = i0 + nSize
            call getpqrso_incore(iCore(i0), k, irrepk, pop_c, pop_d, iop, pop_l, pop_k, iop, abij, 1)
            ! get W(cd,a*), i1
            nSize = idsymsoc(sym, pop_c, pop_d, iop, pop_a, zorb, 1)
            i2 = i1 + nSize
            call getpqrso(iCore(i1), k, irrepk, 1, file_Wakcd)
            ! tmp(cd,a*) = W(cd,*a)-fact*V(cd,*l)*t(a,l)
            call VpqrmVsm_to_Vpqrs(iCore(i0), r1, irrepk, symR, pop_l, &
                                   pop_c, pop_d, iop, zorb, pop_a, iCore(i1), sym, -fact, 1.0d0)
            ! put W(cd,a*)
            call putpqrso(iCore(i1), k, irrepk, 1, file_Wakcd)
        end do
        iEnd = max(i2, iEnd)
    end do

    if (iEnd > maxCore) then
        write (6, *) 'memory not enough for p_Wakcd_r1abij_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wakcd_r1abij_h

! 2026.7.10
!> W(cd,ak) -= fact*V(cd,lk)*r(a,l)
subroutine Wakcd_r1abij_R(WakcdABAB, r1AA, abijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WakcdABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(cD,aK) = -V(cD,lK)*r(a,l)
    call Wakcd_r1abij(WakcdABAB, r1AA, abijABAB, fact, iCore, maxCore)

end subroutine Wakcd_r1abij_R

! 2026.7.10
!> W(cd,ak) -= fact*V(cd,lk)*r(a,l)
subroutine Wakcd_r1abij_U(WakcdAAAA, WakcdABAB, WakcdBBBB, WakcdBABA, spin, &
                          r1AA, r1BB, abijAAAA, abijABAB, abijBBBB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WakcdAAAA, WakcdABAB, WakcdBBBB, WakcdBABA
    character(len=4), intent(in) :: spin !> AAAA, ABAB, BBBB, BABA
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(cd,ak) = -V(cd,lk)*r(a,l)
        call Wakcd_r1abij(WakcdAAAA, r1AA, abijAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(cD,aK) = -V(cD,lK)*r(a,l)
        call Wakcd_r1abij(WakcdABAB, r1AA, abijABAB, fact, iCore, maxCore)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(CD,AK) = -V(CD,LK)*r(A,L)
        call Wakcd_r1abij(WakcdBBBB, r1BB, abijBBBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Cd,Ak) = -V(Cd,Lk)*r(A,L) = -V(dC,kL)*r(A,L)
        if (WakcdBABA%incore) then
            call vv2v_424(abijABAB, r1BB, WakcdBABA, -fact, &
                          'dCkL', 'AL', 'CdAk', iCore, maxCore, 0)
        else
            write (6, '(A)') 'ERROR: Wakcd_r1abij_U current not support iostate=2'
            stop
        end if
    end if

end subroutine Wakcd_r1abij_U
