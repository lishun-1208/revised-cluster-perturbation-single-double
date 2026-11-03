! write at 2024.12.25
! c(ab,ij) -= fact*Pab{W(kb,ij)*r(a,k)}
! c(ab,ij) -= W(kb,ij)*r(a,k)-W(ka,ij)*r(b,k)
! consider spin:
! c(aB,iJ) -= W(kB,iJ)*r(a,k)+W(Ka,Ji)*r(B,K)
!--------------------------------
! Wkbij sort at (ij,kb)
!> be deprecated, please use c2_r1ijka
subroutine c2_r1Wkbij(c2, r1, Wkbij, pop_a, pop_b, pop_i, pop_j, pop_k, &
                      iop, symR, symW, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r1, Wkbij
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_k
    integer*8, intent(in) :: iop ! iop=iopab=iopij
    integer*8, intent(in) :: symR, symW
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symC
    integer*8 :: i0, i1, i2, i3, i4, iT, iEnd, nSize
    integer*8, external :: idsymsoc, irpdso

    symC = dirprd(symR, symW)
    i0 = 1
    ! W(ij,lb)->W(lb,ij), i0
    nSize = idsymsoc(symW, pop_k, pop_b, 1, pop_i, pop_j, iop)
    i1 = i0 + nSize
    call itranspso(Wkbij, iCore(i0), pop_i, pop_j, iop, pop_k, pop_b, 1, symW)
    ! r(a,l)->r(l,a), i1
    nSize = irpdso(symR, pop_k, pop_a, 1)
    i2 = i1 + nSize
    call itranspso(r1, iCore(i1), pop_a, zorb, 1, pop_k, zorb, 1, symR)
    ! tmp(ab,ij)=W(lb,ij)*r(l,a), i2
    nSize = idsymsoc(symC, pop_a, pop_b, 1, pop_i, pop_j, iop)
    i3 = i2 + nSize
    call VmqrsVmp2Vpqrs(iCore(i0), iCore(i1), symW, symR, pop_k, &
                        pop_a, pop_b, pop_i, pop_j, iop, iCore(i2), symC, 1.0d0, 0.0d0)
    ! tmp(ab,ij)->tmp(ba,ij), i3
    i4 = i3 + nSize
    iEnd = i4
    call isymtrso(iCore(i2), iCore(i3), pop_a, pop_b, 1, 1, pop_i, pop_j, iop, 0, symC)
    ! tmp2(ab,ij)=tmp(ab,ij)-tmp(ba,ij), i2
    call daxpy(nSize, -1.0d0, iCore(i3), 1, iCore(i2), 1)
    iT = i2
    if (iop == 0) then
        ! compress ab, i3
        call icompso(iCore(i2), iCore(i3), pop_a, pop_b, 1, 0, pop_i, pop_j, 0, 1, symC)
        iT = i3
    end if
    ! c(ab,ij) -= fact*tmp2(ab,ij)
    nSize = idsymsoc(symC, pop_a, pop_b, iop, pop_i, pop_j, iop)
    call daxpy(nSize, -fact, iCore(iT), 1, c2, 1)

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for c2_r1Wkbij', iEnd, maxCore
        stop
    end if
end subroutine c2_r1Wkbij

! c(ab,*j) -= Pab{W(kb,*j)*r(a,k)}
! consider spin:
! c(aB,*J) -= W(kB,*J)*r(a,k)-W(Ka,*J)*r(B,K)
!-------------------------------
! 2025.12.01
!> be depreciated, please use c2_r1ijka
subroutine EA2_c2_r1Wkbij(c2, r1, Wkbij, pop_a, pop_b, pop_j, pop_k, &
                          iopCvv, symR, symW, fact, iCore, maxCore, antiSym, exchange)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r1, Wkbij
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_j, pop_k
    integer*8, intent(in) :: iopCvv
    integer*8, intent(in) :: symR, symW
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    logical, intent(in) :: antiSym ! whether Pab
    logical, intent(in) :: exchange ! whether exchange ab for W(Ka,*J)*r(B,K)

    target :: iCore
    real*8, dimension(:), pointer :: ptr
    integer*8 :: symC
    integer*8 :: i0, i1, i2, i3, i4, iEnd, nSize
    real*8 :: fact2
    integer*8, external :: idsymsoc, irpdso

    symC = dirprd(symR, symW)
    i0 = 1
    iEnd = i0

    ! W(*j,kb)->W(kb,*j), i0
    nSize = idsymsoc(symW, zorb, pop_j, 1, pop_k, pop_b, 1)
    i1 = i0 + nSize
    iEnd = max(iEnd, i1)
    call itranspso(Wkbij, iCore(i0), zorb, pop_j, 1, pop_k, pop_b, 1, symW)

    ! r(a,k)->r(k,a), i1
    nSize = irpdso(symR, pop_a, pop_k, 1)
    i2 = i1 + nSize
    iEnd = max(iEnd, i2)
    call itranspso(r1, iCore(i1), pop_a, zorb, 1, pop_k, zorb, 1, symR)

    ! tmp(ab,*j)=W(kb,*j)*r(a,k), i2
    nSize = idsymsoc(symC, pop_a, pop_b, 1, zorb, pop_j, 1)
    i3 = i2 + nSize
    iEnd = max(iEnd, i3)
    call VmqrsVmp2Vpqrs(iCore(i0), iCore(i1), symW, symR, pop_k, &
                        pop_a, pop_b, zorb, pop_j, 1, iCore(i2), symC, 1.0d0, 0.0d0)

    if (antiSym) then
        fact2 = -1.0d0
        ! tmp(ab,*j)->tmp(ba,*j), i3
        i4 = i3 + nSize
        iEnd = max(iEnd, i4)
        call isymtrso(iCore(i2), iCore(i3), pop_a, pop_b, 1, 1, zorb, pop_j, 1, 0, symC)
        ! tmp2(ab,*j)=tmp(ab,*j)-tmp(ba,*j), i2
        call daxpy(nSize, -1.0d0, iCore(i3), 1, iCore(i2), 1)
        ptr => iCore(i2:i2)
        if (iopCvv == 0) then
            ! compress ab, i3
            call icompso(iCore(i2), iCore(i3), pop_a, pop_b, 1, 0, zorb, pop_j, 1, 1, symC)
            ptr => iCore(i3:i3)
        end if
    else
        if (exchange) then
            fact2 = 1.0d0
            ! tmp(Ba,*J)->tmp(aB,*J), i3
            i4 = i3 + nSize
            iEnd = max(iEnd, i4)
            call isymtrso(iCore(i2), iCore(i3), pop_b, pop_a, 1, 1, zorb, pop_j, 1, 0, symC)
            ptr => iCore(i3:i3)
        else
            fact2 = -1.0d0
            ptr => iCore(i2:i2)
        end if
    end if

    ! c(ab,*j) += tmp(ab,*j)
    nSize = idsymsoc(symC, pop_a, pop_b, iopCvv, zorb, pop_j, 1)
    call daxpy(nSize, fact*fact2, ptr, 1, c2, 1)

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for EA2_c2_r1Wkbij.', iEnd, maxCore
        stop
    end if

end subroutine EA2_c2_r1Wkbij

! 2025.12.01
subroutine EA2_c2_r1Wkbij_spin(c2ABB, r1AA, WkbijABB, WkbijBAB, pop_a, pop_b, pop_j, pop_k, &
                               symR, symW, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2ABB
    real*8, dimension(*), intent(in) :: r1AA, WkbijABB, WkbijBAB
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_j, pop_k
    integer*8, intent(in) :: symR, symW
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    ! W(kB,*J)*r(a,k)
    call EA2_c2_r1Wkbij(c2ABB, r1AA, WkbijABB, pop_a, pop_b, pop_j, pop_k, &
                        1, symR, symW, fact, iCore, maxCore, .false., .false.)
    ! W(Ka,*J)*r(B,K)
    call EA2_c2_r1Wkbij(c2ABB, r1AA, WkbijBAB, pop_a, pop_b, pop_j, pop_k, &
                        1, symR, symW, fact, iCore, maxCore, .false., .true.)

end subroutine EA2_c2_r1Wkbij_spin
