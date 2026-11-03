! 2026.7.10
!> W(ab,dj) += fact*0.5*V(lk,jd)*r(ab,kl)
subroutine Wabdj_r2ijka(Wabdj, r2, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabdj
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    i0 = 1
    if (Wabdj%incore) then
        ! V(lk,jd)->V(kl,dj), i0
        call tmp%init(iCore, ijka%popq, ijka%popp, ijka%pops, ijka%popr, &
                      ijka%ioppq, 1, ijka%sym, T, i0)
        i1 = i0 + tmp%len
        if (ijka%ioppq == 1) then
            call isymtrso(ijka%arr, tmp%arr, ijka%popp, ijka%popq, 1, 1, &
                          ijka%popr, ijka%pops, 1, 1, ijka%sym)
            ! W(ab,dj) += fact*0.5*V(kl,dj)*r(ab,kl)
            call vv2v_444(r2, tmp, Wabdj, 0.5d0*fact, 'abkl', 'kldj', 'abdj', &
                          iCore(i1), maxCore - tmp%len, 0)
        else if (ijka%ioppq == 0) then
            call isymtrso(ijka%arr, tmp%arr, ijka%popp, ijka%popq, 0, 0, &
                          ijka%popr, ijka%pops, 1, 1, ijka%sym)
            call dscal(tmp%len, -1.0d0, tmp%arr, 1)
            ! W(ab,dj) += fact*0.5*V(kl,dj)*r(ab,kl)
            call vv2v_444(r2, tmp, Wabdj, fact, 'abkl', 'kldj', 'abdj', &
                          iCore(i1), maxCore - tmp%len, 0)
        end if

    else
        call p_Wabdj_r2ijka_h(Wabdj%f, r2%arr, ijka%arr, Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                              r2%popr, r2%pops, Wabdj%ioppq, r2%sym, ijka%sym, &
                              fact, iCore, maxCore)
    end if

end subroutine Wabdj_r2ijka

! write at 2024.12.25
! Wabdj += fact*0.5*ijka(lk,jd)*r(ab,kl)
!       -= fact*0.5*ijka(kl,jd)*r(ab,kl)
! ab and kl is antisymmetric
subroutine p_Wabdj_r2ijka_h(file_Wabdj, r2, ijka, pop_a, pop_b, pop_d, pop_j, pop_k, pop_l, &
                            iop, symR, symV, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    character(len=8), intent(in) :: file_Wabdj
    real*8, dimension(*), intent(in) :: r2, ijka
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_k, pop_l
    integer*8, intent(in) :: iop ! iop=iopab=iopkl
    integer*8, intent(in) :: symR, symV
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    real*8 :: iopfact
    integer*8 :: symW, symm, symm2
    integer*8 :: irrepj, numj, j
    integer*8 :: i0, i1, i2, i3, i4, iEnd, nSize
    real*8 :: t0, t1
    integer*8, external :: idsymsoc

    if (iop == 0) then
        iopfact = 2.0d0
    else if (iop == 1) then
        iopfact = 1.0d0
    end if
    symW = dirprd(symR, symV)
    iEnd = 0
    i0 = 1
    ! ijka(kl,jd)->ijka(kl,dj), i0
    nSize = idsymsoc(symV, pop_k, pop_l, iop, pop_d, pop_j, 1)
    i1 = i0 + nSize
    call isymtrso(ijka, iCore(i0), pop_k, pop_l, iop, 0, pop_j, pop_d, 1, 1, symV)
    do irrepj = 1, nirrep
        numj = pop_j(irrepj)
        symm = dirprd(symV, irrepj) ! symmetry of ijka(kl,d*)
        symm2 = dirprd(symW, irrepj) ! symmetry of W(ab,d*)
        do j = 1, numj
            ! get ijka(kl,d*), i1
            nSize = idsymsoc(symm, pop_k, pop_l, iop, pop_d, zorb, 1)
            i2 = i1 + nSize
            call getpqrso_incore(iCore(i1), j, irrepj, &
                                 pop_k, pop_l, iop, pop_d, pop_j, 1, iCore(i0), symV)
            ! tmp(ab,d*)=r(ab,kl)*ijka(kl,d*), i2
            nSize = idsymsoc(symm2, pop_a, pop_b, iop, pop_d, zorb, 1)
            i3 = i2 + nSize
            call VpqmnVmnrs_to_Vpqrs(r2, iCore(i1), symR, symm, pop_k, pop_l, iop, &
                                     pop_a, pop_b, iop, pop_d, zorb, 1, iCore(i2), symm2, iopfact, 0.0d0)
            ! get W(ab,d*), i3
            i4 = i3 + nSize
            iEnd = max(iEnd, i4)
            call getpqrso(iCore(i3), j, irrepj, 1, file_Wabdj)
            ! W(ab,d*) -= fact*0.5*tmp(ab,d*), i3
            call daxpy(nSize, -fact*0.5d0, iCore(i2), 1, iCore(i3), 1)
            ! put W(ab,d*)
            call putpqrso(iCore(i3), j, irrepj, 1, file_Wabdj)
        end do
    end do

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for p_Wabdj_r2ijka_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wabdj_r2ijka_h

! 2026.7.10
!> W(ab,dj) += fact*0.5*V(lk,jd)*r(ab,kl)
subroutine Wabdj_r2ijka_R(WabdjABAB, r2ABAB, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,dJ) = V(Lk,Jd)*r(aB,kL) -> V(lK,jD)*r(aB,kL)
    call Wabdj_r2ijka(WabdjABAB, r2ABAB, ijkaABAB, 2.0d0*fact, iCore, maxCore)

end subroutine Wabdj_r2ijka_R

! 2026.7.10
!> W(ab,dj) += fact*0.5*V(lk,jd)*r(ab,kl)
subroutine Wabdj_r2ijka_U(WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA, spin, &
                          r2AAAA, r2ABAB, r2BBBB, ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ab,dj) = 0.5*V(lk,jd)*r(ab,kl)
        call Wabdj_r2ijka(WabdjAAAA, r2AAAA, ijkaAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(aB,dJ) = V(Lk,Jd)*r(aB,kL)
        call Wabdj_r2ijka(WabdjABAB, r2ABAB, ijkaBABA, 2.0d0*fact, iCore, maxCore)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(AB,DJ) = 0.5*V(LK,JD)*r(AB,KL)
        call Wabdj_r2ijka(WabdjBBBB, r2BBBB, ijkaBBBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ab,Dj) = V(lK,jD)*r(Ab,Kl) = V(lK,jD)*r(bA,lK)
        i0 = 1
        ! r(bA,lK)->r(Ab,Kl), i0
        call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                      1, 1, r2ABAB%sym, T, i0)
        i1 = i0 + tmp%len
        call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                      r2ABAB%popr, r2ABAB%pops, 1, 1, r2ABAB%sym)
        call Wabdj_r2ijka(WabdjBABA, tmp, ijkaABAB, 2.0d0*fact, iCore(i1), maxCore - tmp%len)
    end if

end subroutine Wabdj_r2ijka_U
