! 2026.7.10
!> W(ab,dj) += fact*P(ab){V(dc,ak)*r(bc,jk)}
subroutine Wabdj_r2abci(Wabdj, r2, abci, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabdj
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abci
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wabdj%incore .and. abci%incore) then
        call vv2v_444(abci, r2, Wabdj, fact, 'dcak', 'bcjk', 'abdj', iCore, maxCore, 1)
    else if ((.not. Wabdj%incore) .and. (.not. abci%incore)) then
        call p_Wabdj_r2abci_h(Wabdj%f, r2%arr, abci%f, Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                              r2%popp, r2%popr, r2%sym, abci%sym, Wabdj%ioppq, &
                              fact, iCore, maxCore)
    else
        write (6, *) 'ERROR: Wabdj_r2abci IO state for Wabdj and abci need be the same.', &
            Wabdj%incore, abci%incore
        stop
    end if

end subroutine Wabdj_r2abci

! write at 2025.2.26
! Wabdj += fact*Pab{<ak||dc>*r(cb,kj)}
!       += fact*Pab{abci(dc,ak)*r(cb,kj)}
subroutine p_Wabdj_r2abci_h(file_Wabdj, r2, file_abci, pop_a, pop_b, pop_d, pop_j, pop_c, pop_k, &
                            symR, symV, iop, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    character(len=8), intent(in) :: file_Wabdj
    real*8, dimension(*), intent(in) :: r2
    character(len=8), intent(in) :: file_abci
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_c, pop_k
    integer*8, intent(in) :: symR, symV ! symmetry of r2 and iabj
    integer*8, intent(in) :: iop !> iop=iopab=iopdc=iopcb=iopkj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symW
    integer*8 :: irrep_k, num_k, k
    integer*8 :: irrep_j, num_j, j
    integer*8 :: irrep_dca, num_dca
    integer*8 :: irrep_cbk
    integer*8 :: irrep_cb, num_cb
    integer*8 :: irrep_abd
    integer*8 :: nBlock, iBlock, max_numk_block, numk_block, k_block
    integer*8 :: i0, i1, i2, i3, i4, i5, iT, itmp, iEnd, nSize
    real*8 :: t0, t1
    integer*8, external :: idsymsoc, irpdso

    integer*8 :: reserved_size = 1000000000 ! in mod_variable

    call cpu_time(t0)

    symW = dirprd(symR, symV)
    iEnd = 0
    i0 = 1

    !--------- cycle for k
    do irrep_k = 1, nirrep
        num_k = pop_k(irrep_k)
        irrep_dca = dirprd(irrep_k, symV)
        num_dca = idsymsoc(irrep_dca, pop_d, pop_c, 1, pop_a, zorb, 1)

        !-------------------------------------------
        ! I devide v(dc,ak) in some blocks,
        ! then cycle for eack block
        !-------------------------------------------
        max_numk_block = (maxCore - reserved_size)/num_dca ! max number of k for each block
        if (max_numk_block <= 0) then
            write (6, *) 'ERROR: memory not enough for p_Wabdj_r2abci_h 1', num_dca, maxCore
            stop
        end if
        nBlock = (num_k - 1 + max_numk_block)/max_numk_block ! number of block
        !if (debug) then
        !    write (6, *) 'max_numk_block =', max_numk_block
        !    write (6, *) 'nBlock =', nBlock
        !end if

        !--------- cycle for block
        do iBlock = 1, nBlock
            if (iBlock == nBlock) then
                numk_block = num_k - (nBlock - 1)*max_numk_block
            else
                numk_block = max_numk_block
            end if

            !------- read in v(dc,ak), i0
            i1 = i0 + numk_block*num_dca
            i2 = i1 + num_dca
            iEnd = max(iEnd, i2)
            do k_block = 1, numk_block
                k = k_block + (iBlock - 1)*max_numk_block
                iT = i0 + (k_block - 1)*num_dca
                ! get v(dc,a*), i1
                call getpqrso(iCore(i1), k, irrep_k, 1, file_abci)
                ! v(dc,a*)->v(ca,d*)
                if (iop == 0) then
                    ! expand dc, i2
                    i3 = i2 + num_dca
                    iEnd = max(iEnd, i3)
                    call iexpso(iCore(i1), iCore(i2), &
                                pop_d, pop_c, 0, 1, pop_a, zorb, 1, 0, irrep_dca)
                    call sstgenso(iCore(i2), iCore(iT), num_dca, &
                                  pop_d, pop_c, pop_a, zorb, iCore(i3), irrep_dca, '2314')
                else
                    call sstgenso(iCore(i1), iCore(iT), num_dca, &
                                  pop_d, pop_c, pop_a, zorb, iCore(i2), irrep_dca, '2314')
                end if
            end do

            !------- cycle for j
            do irrep_j = 1, nirrep
                num_j = pop_j(irrep_j)
                irrep_cbk = dirprd(irrep_j, symR)
                irrep_abd = dirprd(irrep_j, symW)
                do j = 1, num_j
                    ! get r(cb,k*), i1
                    nSize = idsymsoc(irrep_cbk, pop_c, pop_b, iop, pop_k, zorb, 1)
                    i2 = i1 + nSize
                    call getpqrso_incore(iCore(i1), j, irrep_j, &
                                         pop_c, pop_b, iop, pop_k, pop_j, iop, r2, symR)
                    !------ get r(cb,**), i2
                    irrep_cb = dirprd(irrep_cbk, irrep_k)
                    num_cb = irpdso(irrep_cb, pop_c, pop_b, 1)
                    i3 = i2 + numk_block*num_cb
                    do k_block = 1, numk_block
                        k = k_block + (iBlock - 1)*max_numk_block
                        iT = i2 + (k_block - 1)*num_cb
                        if (iop == 0) then
                            i4 = i3 + num_cb
                            call getpqrso_incore(iCore(i3), k, irrep_k, &
                                                 pop_c, pop_b, 0, zorb, pop_k, 1, iCore(i1), irrep_cbk)
                            ! expand cb
                            call iexpso(iCore(i3), iCore(iT), &
                                        pop_c, pop_b, 0, 1, zorb, zorb, 1, 0, irrep_cb)
                        else
                            call getpqrso_incore(iCore(iT), k, irrep_k, &
                                                 pop_c, pop_b, 1, zorb, pop_k, 1, iCore(i1), irrep_cbk)
                        end if
                    end do

                    !------ tmp(ba,d*) = v(ca,dk)*r(cb,k*), i3
                    nSize = idsymsoc(irrep_abd, pop_a, pop_b, 1, pop_d, zorb, 1)
                    i4 = i3 + nSize
                    call dscal(nSize, 0.0d0, iCore(i3), 1)
                    do k_block = 1, numk_block
                        iT = i0 + (k_block - 1)*num_dca ! pointer for v(ca,d*)
                        itmp = i2 + (k_block - 1)*num_cb ! pointer for r(cb,**)
                        call VmqrsVmp2Vpqrs(iCore(iT), iCore(itmp), irrep_dca, irrep_cb, pop_c, &
                                            pop_b, pop_a, pop_d, zorb, 1, iCore(i3), irrep_abd, 1.0d0, 1.0d0)
                    end do

                    ! tmp(ba,d*)->tmp(ab,d*), i4
                    i5 = i4 + nSize
                    iEnd = max(iEnd, i5)
                    call isymtrso(iCore(i3), iCore(i4), &
                                  pop_b, pop_a, 1, 1, pop_d, zorb, 1, 0, irrep_abd)
                    ! tmp2(ab,d*) = tmp(ab,d*)-tmp(ba,d*), i4
                    call daxpy(nSize, -1.0d0, iCore(i3), 1, iCore(i4), 1)
                    if (iop == 0) then
                        ! compree ab
                        call icompso(iCore(i4), iCore(i3), &
                                     pop_a, pop_b, 1, 0, pop_d, zorb, 1, 1, irrep_abd)
                        nSize = idsymsoc(irrep_abd, pop_a, pop_b, 0, pop_d, zorb, 1)
                        call dcopy(nSize, iCore(i3), 1, iCore(i4), 1)
                    end if

                    ! get W(ab,d*), i3
                    call getpqrso(iCore(i3), j, irrep_j, 1, file_Wabdj)
                    ! W(ab,d*) += fact*tmp(ab,d*)
                    nSize = idsymsoc(irrep_abd, pop_a, pop_b, iop, pop_d, zorb, 1)
                    call daxpy(nSize, fact, iCore(i4), 1, iCore(i3), 1)
                    ! put W(ab,d*)
                    call putpqrso(iCore(i3), j, irrep_j, 1, file_Wabdj)
                end do
            end do
        end do
    end do

    call cpu_time(t1)
    write (6, '(A,F)') 'cpu_time for p_Wabdj_r2abci_h(o2v4):', t1 - t0
    if ((iEnd - maxCore) >= 0) then
        write (6, *) 'ERROR: memory not enough for p_Wabdj_r2abci_h 2', iEnd, maxCore
        stop
    end if

end subroutine p_Wabdj_r2abci_h

! 2026.7.10
!> W(ab,dj) += fact*P(ab){V(dc,ak)*r(bc,jk)}
subroutine Wabdj_r2abci_R(WabdjABAB, r2AAAA, r2ABAB, abciAAAA, abciABAB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjABAB
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: abciAAAA, abciABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,dJ) = V(dc,ak)*r(Bc,Jk) = V(dc,ak)*r(cB,kJ)
    !           +V(dC,aK)*r(BC,JK) -> V(dC,aK)*r(bc,jk)
    !           -V(dC,Bk)*r(aC,Jk) = -V(Cd,Bk)*r(aC,kJ) -> -V(cD,bK)*r(aC,kJ)
    if (WabdjABAB%incore .and. abciAAAA%incore) then
        call vv2v_444(abciAAAA, r2ABAB, WabdjABAB, fact, 'dcak', 'cbkj', 'abdj', iCore, maxCore, 0)
    end if
    if (WabdjABAB%incore .and. abciABAB%incore) then
        call vv2v_444(abciABAB, r2AAAA, WabdjABAB, fact, 'dcak', 'bcjk', 'abdj', iCore, maxCore, 0)
        call vv2v_444(abciABAB, r2ABAB, WabdjABAB, -fact, 'cdbk', 'ackj', 'abdj', iCore, maxCore, 0)
    end if

end subroutine Wabdj_r2abci_R

! 2026.7.10
!> W(ab,dj) += fact*P(ab){V(dc,ak)*r(bc,jk)}
subroutine Wabdj_r2abci_U(WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA, spin, &
                          r2AAAA, r2ABAB, r2BBBB, &
                          abciAAAA, abciABAB, abciBBBB, abciBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abciAAAA, abciABAB, abciBBBB, abciBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    if ((.not. WabdjAAAA%incore) .or. &
        (.not. WabdjABAB%incore) .or. &
        (.not. WabdjBBBB%incore) .or. &
        (.not. WabdjBABA%incore) .or. &
        (.not. abciAAAA%incore) .or. &
        (.not. abciABAB%incore) .or. &
        (.not. abciBBBB%incore) .or. &
        (.not. abciBABA%incore)) then
        write (6, *) 'ERROR: Wabdj_r2abci_U current not support iostate=2.'
        stop
    end if

    i0 = 1

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ab,dj) = P(ab){V(dc,ak)*r(bc,jk)}
        !           +P(ab){V(dC,aK)*r(bC,jK)}
        call Wabdj_r2abci(WabdjAAAA, r2AAAA, abciAAAA, fact, iCore, maxCore)
        call Wabdj_r2abci(WabdjAAAA, r2ABAB, abciABAB, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(aB,dJ) = V(dc,ak)*r(Bc,Jk) = V(dc,ak)*r(cB,kJ)
        !           +V(dC,aK)*r(BC,JK)
        !           -V(dC,Bk)*r(aC,Jk) = -V(Cd,Bk)*r(aC,kJ)
        if (WabdjABAB%incore .and. abciAAAA%incore) then
            call vv2v_444(abciAAAA, r2ABAB, WabdjABAB, fact, 'dcak', 'cBkJ', 'aBdJ', iCore, maxCore, 0)
        end if
        if (WabdjABAB%incore .and. abciABAB%incore) then
            call vv2v_444(abciABAB, r2BBBB, WabdjABAB, fact, 'dCaK', 'BCJK', 'aBdJ', iCore, maxCore, 0)
        end if
        if (WabdjABAB%incore .and. abciBABA%incore) then
            call vv2v_444(abciBABA, r2ABAB, WabdjABAB, -fact, 'CdBk', 'aCkJ', 'aBdJ', iCore, maxCore, 0)
        end if
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(AB,DJ) = P(AB){V(DC,AK)*r(BC,JK)}
        !           +P(AB){V(Dc,Ak)*r(Bc,Jk)} = P(AB){V(Dc,Ak)*r(cB,kJ)}
        call Wabdj_r2abci(WabdjBBBB, r2BBBB, abciBBBB, fact, iCore, maxCore)
        ! r(cB,kJ)->r(Bc,Jk)
        call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                      1, 1, r2ABAB%sym, T, i0)
        i1 = i0 + tmp%len
        call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                      r2ABAB%popr, r2ABAB%pops, 1, 1, r2ABAB%sym)
        call Wabdj_r2abci(WabdjBBBB, tmp, abciBABA, fact, iCore(i1), maxCore - tmp%len)
    end if

    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ab,Dj) = V(Dc,Ak)*r(bc,jk)
        !           +V(DC,AK)*r(bC,jK)
        !           -V(Dc,bK)*r(Ac,jK) = -V(cD,bK)*r(cA,jK)
        if (WabdjBABA%incore .and. abciBABA%incore) then
            call vv2v_444(abciBABA, r2AAAA, WabdjBABA, fact, 'DcAk', 'bcjk', 'AbDj', iCore, maxCore, 0)
        end if
        if (WabdjBABA%incore .and. abciBBBB%incore) then
            call vv2v_444(abciBBBB, r2ABAB, WabdjBABA, fact, 'DCAK', 'bCjK', 'AbDj', iCore, maxCore, 0)
        end if
        if (WabdjBABA%incore .and. abciABAB%incore) then
            call vv2v_444(abciABAB, r2ABAB, WabdjBABA, -fact, 'cDbK', 'cAjK', 'AbDj', iCore, maxCore, 0)
        end if
    end if

end subroutine Wabdj_r2abci_U
