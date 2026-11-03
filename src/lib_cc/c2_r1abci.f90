! 2026.7.6
!> c(ab,ij) += fact*P(ij){r(e,i)*V(ab,ej)}
subroutine c2_r1abci(c2, r1, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abci%incore) then
        call vv2v_424(abci, r1, c2, fact, 'abej', 'ei', 'abij', iCore, maxCore, 2)
    else
        write (6, *) 'ERROR: c2_r1abci current not support outcore.'
        stop
        call p_c2_r1abci_h(c2%arr, r1%arr, abci%f, c2%popp, c2%popq, c2%popr, c2%pops, &
                           r1%popp, c2%ioppq, c2%ioprs, abci%ioppq, r1%sym, abci%sym, &
                           iCore, maxCore, fact, 1, T)
    end if

end subroutine c2_r1abci

! 2026.7.6
!> c(ab,ij) += fact*P(ij){r(e,i)*V(ab,ej)}
subroutine c2_r1abci_R(c2ABAB, r1AA, abciABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) = r(e,i)*V(aB,eJ)
    !           +r(E,J)*V(Ba,Ei) -> r(e,j)*V(bA,eI)
    if (abciABAB%incore) then
        call vv2v_424(abciABAB, r1AA, c2ABAB, fact, 'aBeJ', 'ei', 'aBiJ', iCore, maxCore, 0)
        call vv2v_424(abciABAB, r1AA, c2ABAB, fact, 'baei', 'ej', 'abij', iCore, maxCore, 0)
    else
        write (6, *) 'ERROR: c2_r1abci_R current not support outcore.'
        stop
    end if

end subroutine c2_r1abci_R

! 2026.7.6
!> c(ab,ij) += fact*P(ij){r(e,i)*V(ab,ej)}
subroutine c2_r1abci_U(c2AAAA, c2ABAB, c2BBBB, r1AA, r1BB, &
                       abciAAAA, abciABAB, abciBBBB, abciBABA, &
                       fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) += fact*P(ij){r(e,i)*V(ab,ej)}
    call c2_r1abci(c2AAAA, r1AA, abciAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) += r(e,i)*V(aB,eJ)
    if (abciABAB%incore) then
        call vv2v_424(abciABAB, r1AA, c2ABAB, fact, 'aBeJ', 'ei', 'aBiJ', iCore, maxCore, 0)
    else
        write (6, *) 'ERROR: c2_r1abci_R current not support outcore.'
        stop
    end if
    ! c(aB,iJ) = -r(E,J)*V(aB,Ei) = r(E,J)*V(Ba,Ei)
    if (abciABAB%incore) then
        call vv2v_424(abciBABA, r1BB, c2ABAB, fact, 'BaEi', 'EJ', 'aBiJ', iCore, maxCore, 0)
    else
        write (6, *) 'ERROR: c2_r1abci_R current not support outcore.'
        stop
    end if

    ! c(AB,IJ) += fact*P(IJ)){r(E,I)*V(AB,EJ)}
    call c2_r1abci(c2BBBB, r1BB, abciBBBB, fact, iCore, maxCore)

end subroutine c2_r1abci_U

! write by LiShun, 2024.5.27
!-----------------------
! modified at 2025.11.25
subroutine p_c2_r1abci_h(c2, r1, file_abci, pop_a, pop_b, pop_i, pop_j, pop_e, &
                         iopCvv, iopCoo, iopVvv, symR, symV, &
                         iCore, maxCore, fact, type, antiSym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r1
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_e
    integer*8, intent(in) :: iopCvv, iopCoo, iopVvv
    integer*8, intent(in) :: symR, symV
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact
    ! 1 for abci*r1
    ! 2 for abciABAB*rAA
    integer*8, intent(in) :: type
    logical, intent(in) :: antiSym

    target :: iCore
    real*8, dimension(:), pointer :: ptr
    integer*8 :: symC, symm, symm2
    integer*8 :: irrepj, numj, j
    integer*8 :: nSize
    integer*8 :: i0, i1, i2, i3, iEnd
    integer*8, external :: idsymsoc, irpdso

    if (iopCvv /= iopVvv) then
        write (6, *) 'ERROR: p_c2_r1abci_h currently only support iopCvv=iopVvv.'
        stop
    end if
    i0 = 1
    iEnd = 1
    symC = dirprd(symR, symV)
    !--------------------------------------
    ! c(aB,iJ) = Pij{<aB||eJ>*r(e,i)}
    !          = abci(aB,eJ)*r(e,i)+abci(Ba,Ei)*r(E,J)
    !--------------------------------------
    nSize = idsymsoc(symC, pop_a, pop_b, iopCvv, pop_i, pop_j, 1)
    i1 = i0 + nSize ! i0 for tmp(aB,iJ)
    iEnd = max(iEnd, i1)
    do irrepj = 1, nirrep
        numj = pop_j(irrepj)
        symm = dirprd(irrepj, symV) ! symmetry of abci(aB,e*)
        symm2 = dirprd(irrepj, symC) ! symmetry of c(aB,i*)
        do j = 1, numj
            ! read abci(aB,e*), i1
            nSize = idsymsoc(symm, pop_a, pop_b, iopVvv, pop_e, zorb, 1)
            i2 = i1 + nSize
            iEnd = max(iEnd, i2)
            call getpqrso(iCore(i1), j, irrepj, 1, file_abci)
            ! tmp(aB,*i) = abci(aB,*e)*r(e,i), i2
            nSize = idsymsoc(symm2, pop_a, pop_b, iopVvv, pop_i, zorb, 1)
            i3 = i2 + nSize
            iEnd = max(iEnd, i3)
            call VpqrmVms_to_Vpqrs(iCore(i1), r1, symm, symR, pop_e, &
                                   pop_a, pop_b, iopVvv, zorb, pop_i, iCore(i2), symm2, 1.0d0, 0.0d0)
            ! put tmp(aB,i*), i0
            call putpqrso_incore(iCore(i2), j, irrepj, &
                                 pop_a, pop_b, iopCvv, pop_i, pop_j, 1, iCore(i0), symC)
        end do
    end do

    if (antiSym) then
        if (type == 1) then
            ! tmp(ab,ij)->tmp(ab,ji), i1
            nSize = idsymsoc(symC, pop_a, pop_b, iopCvv, pop_i, pop_j, 1)
            i2 = i1 + nSize
            iEnd = max(iEnd, i2)
            call isymtrso(iCore(i0), iCore(i1), pop_a, pop_b, iopCvv, 0, pop_i, pop_j, 1, 1, symC)
            ! tmp2(ab,ij)=tmp(ab,ij)-tmp(ab,ji), i0
            call daxpy(nSize, -1.0d0, iCore(i1), 1, iCore(i0), 1)
            ptr => iCore(i0:i0)
            if (iopCoo == 0) then
                ! compress ij, i1
                call icompso(iCore(i0), iCore(i1), pop_a, pop_b, iopCvv, 1, pop_i, pop_j, 1, 0, symC)
                ptr => iCore(i1:i1)
            end if
        else if (type == 2) then
            ! tmp(aB,iJ)->tmp(Ba,Ji), i1
            nSize = idsymsoc(symC, pop_a, pop_b, 1, pop_i, pop_j, 1)
            i2 = i1 + nSize
            iEnd = max(iEnd, i2)
            call isymtrso(iCore(i0), iCore(i1), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 1, symC)
            ! tmp2(aB,iJ)=tmp(aB,iJ)+tmp(Ba,Ji), i0
            call daxpy(nSize, 1.0d0, iCore(i1), 1, iCore(i0), 1)
            ptr => iCore(i0:i0)
        end if
    else ! for c(ab,*j) = v(ab,cj)*r(c,*)
        ptr => iCore(i0:i0)
    end if

    ! c2 += fact*tmp2(aB,iJ)
    nSize = idsymsoc(symC, pop_a, pop_b, iopCvv, pop_i, pop_j, iopCoo)
    call daxpy(nSize, fact, ptr, 1, c2, 1)

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for c1_r1abci_h', iEnd, maxCore
        stop
    end if

end subroutine p_c2_r1abci_h
