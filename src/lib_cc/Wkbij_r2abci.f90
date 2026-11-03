! 2024.7.6
!> Wkbij(ij,kb) += fact*0.5*V(dc,bk)*r(cd,ij), need V(dc,bk)=-V(cd,bk)
subroutine Wkbij_r2abci(Wkbij, r2, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wkbij
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1
    type(opt) :: tmp

    i0 = 1
    if (abci%incore) then
        ! Wkbij(ij,kb) = 0.5*V(dc,bk)*r(cd,ij) = -0.5*V(cd,bk)*r(cd,ij)
        ! tmp(ij,bk), i0
        call tmp%init(iCore, Wkbij%popp, Wkbij%popq, Wkbij%pops, Wkbij%popr, &
                      Wkbij%ioppq, 1, Wkbij%sym, T, i0)
        i1 = i0 + tmp%len
        call dscal(tmp%len, 0.0d0, tmp%arr, 1)
        if (r2%ioppq == 0) then
            call vv2v_444(r2, abci, tmp, -fact, 'cdij', 'cdbk', 'ijbk', &
                          iCore(i1), maxCore - tmp%len, 0)
        else
            call vv2v_444(r2, abci, tmp, -0.5d0*fact, 'cdij', 'cdbk', 'ijbk', &
                          iCore(i1), maxCore - tmp%len, 0)
        end if
        ! tmp(ij,bk)->tmp(ij,kb), i1
        call isymtrso(tmp%arr, iCore(i1), tmp%popp, tmp%popq, tmp%ioppq, 0, &
                      tmp%popr, tmp%pops, 1, 1, tmp%sym)
        call daxpy(Wkbij%len, 1.0d0, iCore(i1), 1, Wkbij%arr, 1)
    else
        write (6, *) 'ERROR: Wkbij_r2abci current not support outcore.'
        stop
        call p_Wkbij_r2abci_h(Wkbij%arr, r2%arr, abci%f, Wkbij%popr, Wkbij%pops, Wkbij%popp, Wkbij%popq, &
                              r2%popp, r2%popq, r2%sym, abci%sym, Wkbij%ioppq, r2%ioppq, r2%ioprs, abci%ioppq, &
                              fact, iCore, maxCore, F)
    end if

end subroutine Wkbij_r2abci

! 2024.7.6
!> Wkbij(ij,kb) += fact*0.5*V(dc,bk)*r(cd,ij), need V(dc,bk)=-V(cd,bk)
subroutine Wkbij_r2abci_R(WkbijABAB, r2ABAB, abciABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(iJ,kB) += V(Dc,Bk)*r(cD,iJ) -> V(dC,bK)*r(cD,iJ)
    if (abciABAB%incore) then
        call vv2v_444(r2ABAB, abciABAB, WkbijABAB, fact, 'cdij', 'dcbk', 'ijkb', iCore, maxCore, 0)
    else
        write (6, *) 'ERROR: Wkbij_r2abci current not support outcore.'
        stop
    end if

end subroutine Wkbij_r2abci_R

! 2024.7.6
!> W(ij,kb) += fact*0.5*V(dc,bk)*r(cd,ij), need V(dc,bk)=-V(cd,bk)
subroutine Wkbij_r2abci_U(WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA, spin, &
                          r2AAAA, r2ABAB, r2BBBB, abciAAAA, abciABAB, abciBBBB, abciBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA
    character(len=4) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8:: i0, i1
    type(opt) :: tmp

    i0 = 1

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ij,kb) += 0.5*V(dc,bk)*r(cd,ij)
        call Wkbij_r2abci(WkbijAAAA, r2AAAA, abciAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(iJ,kB) += V(Dc,Bk)*r(cD,iJ)
        if (abciBABA%incore) then
            call vv2v_444(r2ABAB, abciBABA, WkbijABAB, fact, 'cDiJ', 'DcBk', 'iJkB', iCore, maxCore, 0)
        else
            write (6, *) 'ERROR: Wkbij_r2abci current not support outcore.'
            stop
        end if
    end if
    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(IJ,KB) += 0.5*V(DC,BK)*r(CD,IJ)
        call Wkbij_r2abci(WkbijBBBB, r2BBBB, abciBBBB, fact, iCore, maxCore)
    end if
    ! BABA
    if (spin(4:4) == 'T') then
        !--- W(Ij,Kb) += V(dC,bK)*r(Cd,Ij)
        ! r(dC,jI)->r(Cd,Ij)
        call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                      1, 1, r2ABAB%sym, T, i0)
        i1 = i0 + tmp%len
        call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                      r2ABAB%popr, r2ABAB%pops, 1, 1, tmp%sym)

        if (abciABAB%incore) then
            call vv2v_444(tmp, abciABAB, WkbijBABA, fact, 'CdIj', 'dCbK', 'IjKb', &
                          iCore(i1), maxCore - tmp%len, 0)
        else
            write (6, *) 'ERROR: Wkbij_r2abci current not support outcore.'
            stop
        end if
    end if

end subroutine Wkbij_r2abci_U

! modified at 2025.10.4
subroutine p_Wkbij_r2abci_h(Wkbij, r2, file_abci, pop_k, pop_b, pop_i, pop_j, pop_c, pop_d, &
                            symR, symV, iopWoo, iopRvv, iopRoo, iopVvv, &
                            fact, iCore, maxCore, exchange)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: Wkbij
    real*8, dimension(*), intent(in) :: r2
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_k, pop_b, pop_i, pop_j, pop_c, pop_d
    integer*8, intent(in) :: symR, symV ! symmetry of r2 and abci
    integer*8, intent(in) :: iopWoo, iopRvv, iopRoo, iopVvv
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    logical, intent(in) :: exchange ! for W(kB,iJ), change Dc to cD

    target :: iCore
    real*8, dimension(:), pointer :: ptr
    integer*8 :: symW, symm, symm2
    integer*8 :: irrepk, numk, k
    integer*8 :: iW, i0, i1, i2, iEnd, nSize
    real*8 :: fact2
    real*8 :: t0, t1
    integer*8, external :: idsymsoc, irpdso

    if (iopRvv /= iopVvv) then
        write (6, *) 'ERROR: p_Wkbij_r2abci_h current not support iopRvv/=iopVvv.'
        stop
    end if
    if (iopRoo /= iopWoo) then
        write (6, *) 'ERROR: p_Wkbij_r2abci_h current not support iopRoo/=iopWoo.'
        stop
    end if

    if (iopVvv == 0) then
        fact2 = -2.0d0
    else if (iopVvv == 1) then
        fact2 = -1.0d0
    end if

    symW = dirprd(symV, symR)
    iW = 1 ! W(ij,bk)
    nSize = idsymsoc(symW, pop_i, pop_j, iopWoo, pop_b, pop_k, 1)
    i0 = iW + nSize
    iEnd = i0

    do irrepk = 1, nirrep
        numk = pop_k(irrepk)
        symm = dirprd(irrepk, symV) ! symmetry of V(cd,b*)
        symm2 = dirprd(irrepk, symW) ! symmetry of W(ij,b*)
        nSize = idsymsoc(symm, pop_c, pop_d, iopVvv, pop_b, zorb, 1)
        i1 = i0 + nSize
        iEnd = max(iEnd, i1)
        do k = 1, numk
            !--- get V(cd,b*), i0
            if (exchange) then
                fact2 = 1.0d0
                ! get V(Dc,B*), i1
                nSize = idsymsoc(symm, pop_d, pop_c, 1, pop_b, zorb, 1)
                i2 = i1 + nSize
                iEnd = max(iEnd, i2)
                call getpqrso(iCore(i1), k, irrepk, 1, file_abci)
                ! V(Dc,B*)->V(cD,B*), i0
                call isymtrso(iCore(i1), iCore(i0), pop_d, pop_c, 1, 1, pop_b, zorb, 1, 0, symm)
            else
                call getpqrso(iCore(i0), k, irrepk, 1, file_abci)
            end if

            ! tmp(ij,b*) = r(cd,ij)*V(cd,b*), i1
            nSize = idsymsoc(symm2, pop_i, pop_j, iopRoo, pop_b, zorb, 1)
            i2 = i1 + nSize
            iEnd = max(iEnd, i2)
            call VmnpqVmnrs_to_Vpqrs(r2, iCore(i0), symR, symm, pop_c, pop_d, iopVvv, &
                                     pop_i, pop_j, iopRoo, pop_b, zorb, 1, iCore(i1), symm2, 1.0d0, 0.0d0)

            ! put tmp(ij,bk)
            call putpqrso_incore(iCore(i1), k, irrepk, pop_i, pop_j, iopRoo, pop_b, pop_k, 1, &
                                 iCore(iW), symW)
        end do
    end do

    ! tmp(ij,bk)->tmp(ij,kb), i0
    nSize = idsymsoc(symW, pop_i, pop_j, iopRoo, pop_k, pop_b, 1)
    i1 = i0 + nSize
    iEnd = max(iEnd, i1)
    call isymtrso(iCore(iW), iCore(i0), pop_i, pop_j, iopRoo, 0, pop_b, pop_k, 1, 1, symW)

    ! W(ij,kb) = fact*0.5*tmp(ij,kb)
    call daxpy(nSize, 0.5d0*fact*fact2, iCore(i0), 1, Wkbij, 1)

    if (iEnd >= maxCore) then
        write (6, *) 'ERROR: memory not enough for p_Wkbij_r2abci_h', iEnd, maxCore
        stop
    end if

end subroutine p_Wkbij_r2abci_h
