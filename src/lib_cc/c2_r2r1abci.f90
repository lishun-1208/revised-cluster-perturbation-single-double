! 2026.4.16
!> c(ab,ij) -= 0.5*fact*r(ef,ij)*Pab{r(b,m)*Vamef(ef,am)}
subroutine c2_r2r1abci(c2, r2, r1, abci, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2, r1
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abci%incore) then
        call p_c2_r2r1abci(c2%arr, r2%arr, r1%arr, abci%arr, &
                           c2%popp, c2%popq, c2%popr, c2%pops, r2%popp, r2%popq, r1%popq, &
                           r1%sym, r2%sym, c2%ioppq, fact, iCore, maxCore)
    else
        call p_c2_r2r1abci_h(c2%arr, r2%arr, r1%arr, abci%f, &
                             c2%popp, c2%popq, c2%popr, c2%pops, r2%popp, r2%popq, r1%popq, &
                             c2%ioppq, r1%sym, r2%sym, fact, iCore, maxCore, 1)
    end if

end subroutine c2_r2r1abci

! write by LiShun, 2024.5.27
! c2(ab,ij) -= 0.5*fact*r2(ef,ij)*P(ab){r1(b,m)*<am||ef>}
!           -= 0.5*fact*P(ab)*t(b,m)*t(ef,ij)*v(ef,am)
!           -= 0.5*fact*[t(b,m)*t(ef,ij)*v(ef,am)-t(a,m)t(ef,ij)*v(ef,bm)]
!           -= 0.5*fact*[t(b,m)*t(ef,ij)*v(ef,am)+t(a,m)t(ef,ji)*v(ef,bm)]
!           -= 0.5*fact*[t(b,m)*tv(ij,am)+t(a,m)tv(ji,bm)]
!           -= 0.5*fact*[ttv(ij,ab)+ttv(ji,ba)]
! tv(ij,am) = t(ef,ij)*v(ef,am)
! ttv(ij,ab) = tv(ij,am)*t(b,m), ab not antisymmtric
!----------------------------
! modified at 2024.10.26
subroutine p_c2_r2r1abci_h(c2, r2, r1, file_abci, pop_a, pop_b, pop_i, pop_j, &
                           pop_e, pop_f, pop_m, iop, symR1, symR2, fact, iCore, maxCore, type)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r2, r1
    character(8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_e, pop_f, pop_m
    integer*8, intent(in) :: iop !> need iopab=iopij=iopef=iop
    integer*8, intent(in) :: symR1, symR2
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(in) :: iCore
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: type
    ! 1 for r(ab,ij) -= 0.5*t(ef,ij)*P(ab){t(b,m)*<am||ef>}
    ! 2 for r(aB,iJ) -= 0.5*[t(B,M)*t(eF,iJ)*v(eF,aM)+t(a,m)*t(Ef,Ji)*v(Ef,Bm)]

    real*8 :: fact2
    integer*8 :: irrepm, numm, m
    integer*8 :: symTV, symC
    integer*8 :: i0, i1, i2, i3, iEnd, iT
    integer*8 :: nSize
    real*8 :: t0, t1
    integer*8, external :: idsymsoc

    call cpu_time(t0)

    symC = dirprd(symR1, symR2)
    if (iop == 0) then
        fact2 = 2.0d0 ! for e<f
    else if (iop == 1) then
        fact2 = 1.0d0
    end if
    i0 = 1
    if (type == 1) then ! iop can be 0 or 1
        nSize = idsymsoc(symR2, pop_i, pop_j, iop, pop_a, pop_m, 1)
        i1 = i0 + nSize ! i0 for tv(ij,am)
        do irrepm = 1, nirrep
            numm = pop_m(irrepm)
            symTV = dirprd(irrepm, symR2)
            do m = 1, numm
                ! get abci(ef,a*), i1
                nSize = idsymsoc(irrepm, pop_e, pop_f, iop, pop_a, zorb, 1)
                i2 = i1 + nSize
                call getpqrso(iCore(i1), m, irrepm, 1, file_abci)
                ! tv(ij,a*)=t(ef,ij)*abci(ef,a*), i2
                nSize = idsymsoc(symTV, pop_i, pop_j, iop, pop_a, zorb, 1)
                call VmnpqVmnrs_to_Vpqrs(r2, iCore(i1), symR2, irrepm, pop_e, pop_f, iop, &
                                         pop_i, pop_j, iop, pop_a, zorb, 1, iCore(i2), symTV, &
                                         fact2, 0.0d0)
                ! put tv(ij,a*), i0
                call putpqrso_incore(iCore(i2), m, irrepm, pop_i, pop_j, iop, &
                                     pop_a, pop_m, 1, iCore(i0), symR2)
            end do
        end do
        ! ttv(ij,ab)=tv(ij,am)*t(b,m), i1
        nSize = idsymsoc(symC, pop_i, pop_j, iop, pop_a, pop_b, 1)
        i2 = i1 + nSize
        call VpqrmVsm_to_Vpqrs(iCore(i0), r1, symR2, symR1, pop_m, &
                               pop_i, pop_j, iop, pop_a, pop_b, iCore(i1), symC, 1.0d0, 0.0d0)
        ! ttv(ij,ab)->ttv(ij,ba), i2
        iEnd = i2 + nSize
        call isymtrso(iCore(i1), iCore(i2), pop_i, pop_j, iop, 0, pop_a, pop_b, 1, 1, symC)
        ! tmp(ij,ab)=ttv(ij,ab)-ttv(ij,ba), i1
        call daxpy(nSize, -1.0d0, iCore(i2), 1, iCore(i1), 1)
        if (iop == 0) then
            ! compress ab, i2
            call icompso(iCore(i1), iCore(i2), pop_i, pop_j, 0, 1, pop_a, pop_b, 1, 0, symC)
            ! tmp(ij,ab)->tmp(ab,ij), i1
            call itranspso(iCore(i2), iCore(i1), pop_i, pop_j, 0, pop_a, pop_b, 0, symC)
            iT = i1
        else if (iop == 1) then
            ! tmp(ij,ab)->tmp(ab,ij), i2
            call itranspso(iCore(i1), iCore(i2), pop_i, pop_j, 1, pop_a, pop_b, 1, symC)
            iT = i2
        end if
        ! c2(ab,ij) -= fact*0.5*tmp(ab,ij)
        nSize = idsymsoc(symC, pop_a, pop_b, iop, pop_i, pop_j, iop)
        call daxpy(nSize, -0.5d0*fact, iCore(iT), 1, c2, 1)
    else if (type == 2) then ! iop must be 1
        ! i0 for tv(iJ,aM)
        nSize = idsymsoc(symR2, pop_i, pop_j, 1, pop_a, pop_m, 1)
        i1 = i0 + nSize
        do irrepm = 1, nirrep
            numm = pop_m(irrepm)
            symTV = dirprd(irrepm, symR2)
            do m = 1, numm
                ! get abci(eF,a*), i1
                nSize = idsymsoc(irrepm, pop_e, pop_f, 1, pop_a, zorb, 1)
                i2 = i1 + nSize
                call getpqrso(iCore(i1), m, irrepm, 1, file_abci)
                ! tv(iJ,a*)=t(eF,iJ)*abci(eF,a*), i2
                nSize = idsymsoc(symTV, pop_i, pop_j, 1, pop_a, zorb, 1)
                call VmnpqVmnrs_to_Vpqrs(r2, iCore(i1), symR2, irrepm, pop_e, pop_f, 1, &
                                         pop_i, pop_j, 1, pop_a, zorb, 1, iCore(i2), symTV, 1.0d0, 0.0d0)
                ! put tv(iJ,a*), i0
                call putpqrso_incore(iCore(i2), m, irrepm, pop_i, pop_j, 1, &
                                     pop_a, pop_m, 1, iCore(i0), symR2)
            end do
        end do
        ! ttv(iJ,aB)=tv(iJ,aM)*t(B,M), i1
        nSize = idsymsoc(symC, pop_i, pop_j, 1, pop_a, pop_b, 1)
        i2 = i1 + nSize
        call VpqrmVsm_to_Vpqrs(iCore(i0), r1, symR2, symR1, pop_m, &
                               pop_i, pop_j, 1, pop_a, pop_b, iCore(i1), symC, 1.0d0, 0.0d0)
        ! ttv(iJ,aB)->ttv(aB,iJ), i2
        iEnd = i2 + nSize
        call itranspso(iCore(i1), iCore(i2), pop_i, pop_j, 1, pop_a, pop_b, 1, symC)
        ! c2(aB,iJ) -= fact*0.5*ttv(aB,iJ)
        call daxpy(nSize, -0.5d0*fact, iCore(i2), 1, c2, 1)
        ! ttv(aB,iJ)->ttv(Ba,Ji), i1
        call isymtrso(iCore(i2), iCore(i1), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 1, symC)
        ! c2(aB,iJ) -= fact*0.5*ttv(Ba,Ji)
        call daxpy(nSize, -0.5d0*fact, iCore(i1), 1, c2, 1)
    end if

    call cpu_time(t1)
    write (6, '(A,F)') 'cpu_time for p_c2_r2r1abci_h(o3v3):', t1 - t0
    if (iEnd >= maxCore) then
        write (6, *) 'memory not enough for p_c2_r2r1abci_h', iEnd, maxCore
        stop
    end if

end subroutine p_c2_r2r1abci_h

! c(ab,ij) -= fact*0.5*r(ef,ij)*Pab{r(b,m)*<am||ef>}
!          -= fact*0.5*Pab{r(b,m)*r(ef,ij)*abci(ef,am)}
!          -= fact*0.5*Pab{r(b,m)*rV(ij,am)}
!------------------------------
! modified at 2024.11.22
subroutine p_c2_r2r1abci(c2, r2, r1, abci, pop_a, pop_b, pop_i, pop_j, pop_e, pop_f, pop_m, &
                         symR1, symR2, iop, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r2, r1, abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_e, pop_f, pop_m
    integer*8, intent(in) :: symR1, symR2 ! symmetry of R1 and R2
    integer*8, intent(in) :: iop !> iop=iopab=iopij=iopef
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1, i2, iT, iEnd, nSize
    integer*8 :: symR12 ! dirprd(symR1,symR2)
    real*8 :: iopfact
    real*8 :: t0, t1

    integer*8, external :: idsymsoc

    call cpu_time(t0)

    i0 = 1
    ! tmp1(ij,am)=r(ef,ij)*abci(ef,am), i0, o3v3
    nSize = idsymsoc(symR2, pop_i, pop_j, iop, pop_a, pop_m, 1)
    i1 = i0 + nSize
    if (iop == 0) then
        iopfact = 2.0d0 ! because e<f
    else if (iop == 1) then
        iopfact = 1.0d0
    end if
    call VmnpqVmnrs_to_Vpqrs(r2, abci, symR2, 1, pop_e, pop_f, iop, &
                             pop_i, pop_j, iop, pop_a, pop_m, 1, iCore(i0), symR2, iopfact, 0.0d0)
    ! tmp2(ij,ab)=tmp1(ij,am)*r1(b,m), i1, o3v2
    symR12 = dirprd(symR1, symR2)
    nSize = idsymsoc(symR12, pop_i, pop_j, iop, pop_a, pop_b, 1) ! ab not antisymm
    i2 = i1 + nSize
    call VpqrmVsm_to_Vpqrs(iCore(i0), r1, symR2, symR1, pop_m, &
                           pop_i, pop_j, iop, pop_a, pop_b, iCore(i1), symR12, 1.0d0, 0.0d0)
    ! tmp2(ij,ab)->tmp2(ij,ba), i2
    iEnd = i2 + nSize
    call isymtrso(iCore(i1), iCore(i2), pop_i, pop_j, iop, 0, pop_a, pop_b, 1, 1, symR12)
    ! tmp3(ij,ab)=tmp2(ij,ab)-tmp2(ij,ba), i1
    call daxpy(nSize, -1.0d0, iCore(i2), 1, iCore(i1), 1) ! ab is antisymm
    ! tmp3(ij,ab)->tmp3(ab,ij)
    if (iop == 0) then ! i1
        ! compress ab, i2
        call icompso(iCore(i1), iCore(i2), pop_i, pop_j, 0, 1, pop_a, pop_b, 1, 0, symR12)
        call itranspso(iCore(i2), iCore(i1), pop_i, pop_j, 0, pop_a, pop_b, 0, symR12)
        iT = i1
    else if (iop == 1) then ! i2
        call itranspso(iCore(i1), iCore(i2), pop_i, pop_j, 1, pop_a, pop_b, 1, symR12)
        iT = i2
    end if
    ! c(ab,ij) -= fact*0.5*tmp3(ab,ij)
    nSize = idsymsoc(symR12, pop_a, pop_b, iop, pop_i, pop_j, iop)
    call axpyso(nSize, -0.5d0*fact, iCore(iT), 1, c2, 1)

    call cpu_time(t1)
    write (6, '(A,F)') 'cpu_time for p_c2_r2r1abci(o3v3):', t1 - t0
    if (iEnd >= maxCore) then
        write (6, *) 'iCore not enough for p_c2_r2r1abci', iEnd, maxCore
        stop
    end if

end subroutine p_c2_r2r1abci
