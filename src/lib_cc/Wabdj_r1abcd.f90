! 2026.7.7
!> W(ab,dj) += fact*V(ab,dc)*r(c,j)
subroutine Wabdj_r1abcd(Wabdj, r1, abcd, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabdj
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abcd
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wabdj%incore .and. abcd%incore) then
        call vv2v_424(abcd, r1, Wabdj, fact, 'abdc', 'cj', 'abdj', iCore, maxCore, 0)
    else if (Wabdj%incore .and. (.not. abcd%incore)) then
        call p_Wabdj_r1abcd_h_2(Wabdj%arr, r1%arr, abcd%f, &
                                Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                                r1%popp, r1%sym, abcd%sym, Wabdj%ioppq, fact, iCore, maxCore, 1)
    else if ((.not. Wabdj%incore) .and. (.not. abcd%incore)) then
        call p_Wabdj_r1abcd_h(Wabdj%f, r1%arr, abcd%f, Wabdj%popp, Wabdj%popq, Wabdj%popr, Wabdj%pops, &
                              r1%popp, r1%sym, abcd%sym, Wabdj%ioppq, fact, iCore, maxCore)
    else
        write (6, *) 'ERROR: Wabdj_r1abcd IO state error.', Wabdj%incore, abcd%incore
        stop
    end if

end subroutine Wabdj_r1abcd

! write at 2025.2.26
! can use for VABAB*rBB
subroutine p_Wabdj_r1abcd_h(file_Wabdj, r1, file_abcd, pop_a, pop_b, pop_d, pop_j, pop_c, &
                            symR, symV, iop, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    character(len=8), intent(in) :: file_Wabdj
    real*8, dimension(*), intent(in) :: r1
    character(len=8), intent(in) :: file_abcd
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_c
    integer*8, intent(in) :: symR, symV ! symmetry of r1 and iabj
    integer*8, intent(in) :: iop !> iop=iopab=iopdc
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symW
    integer*8 :: irrepc, numc, c
    integer*8 :: irrepj, numj, j
    integer*8 :: irrep_abd, num_abd
    integer*8 :: nBlock, iBlock, max_numc_block, numc_block, c_block
    integer*8 :: i0, i1, i2, i3, idx, iEnd
    integer*8, external :: idsymsoc, irpdso

    integer*8 :: reserved_size = 1000000000 ! 8GB

    symW = dirprd(symR, symV)
    iEnd = 0
    i0 = 1

    !--------- cycle for c
    do irrepc = 1, nirrep
        numc = pop_c(irrepc)
        irrepj = dirprd(irrepc, symR)
        numj = pop_j(irrepj)
        irrep_abd = dirprd(irrepj, symW)

        !---------------------------------------
        ! I devide v(ab,dc) in some blocks,
        ! each block have several v(ab,d*), and cycle for blocks
        !---------------------------------------
        num_abd = idsymsoc(irrep_abd, pop_a, pop_b, iop, pop_d, zorb, 1)
        max_numc_block = (maxCore - reserved_size)/num_abd ! max number of W(ab,d*) in a block
        if (max_numc_block <= 0) then
            write (6, *) 'ERROR: memory not enough for p_Wabdj_r1abcd_h 1', num_abd, maxCore
            stop
        end if
        nBlock = (numc - 1 + max_numc_block)/max_numc_block ! total number of block
        !if (debug) then
        !    write (6, *) 'max_numc_block =', max_numc_block
        !    write (6, *) 'nBlock =', nBlock
        !end if

        !--------- cycle for block
        do iBlock = 1, nBlock
            if (iBlock == nBlock) then
                ! for the last block, numc_block may less than max_numc_block
                numc_block = numc - max_numc_block*(nBlock - 1)
            else
                numc_block = max_numc_block
            end if

            ! get v(ab,d*), i0
            i1 = i0 + num_abd*numc_block
            do c_block = 1, numc_block
                c = (iBlock - 1)*max_numc_block + c_block
                idx = i0 + (c_block - 1)*num_abd
                call getpqrso(iCore(idx), c, irrepc, 1, file_abcd)
            end do

            !-------- cycle for j
            do j = 1, numj
                ! get r(c,j), i1
                i2 = i1 + numc
                call getpqrso_incore(iCore(i1), j, irrepj, &
                                     zorb, zorb, 1, pop_c, pop_j, 1, r1, symR)
                ! read in W(ab,d*), i2
                i3 = i2 + num_abd
                call getpqrso(iCore(i2), j, irrepj, 1, file_Wabdj)

                !------ W(ab,d*)=W(ab,d*)+fact*v(ab,dc)*r(c,*), i2
                !------ cycle for c
                do c_block = 1, numc_block
                    c = (iBlock - 1)*max_numc_block + c_block
                    idx = i0 + (c_block - 1)*num_abd
                    call daxpy(num_abd, fact*iCore(i1 - 1 + c), iCore(idx), 1, iCore(i2), 1)
                end do

                ! write W(ab,d*)
                call putpqrso(iCore(i2), j, irrepj, 1, file_Wabdj)
            end do
            iEnd = max(iEnd, i3)
        end do
    end do

    if ((iEnd - maxCore) >= 0) then
        write (6, *) 'ERROR: memory not enough for p_Wabdj_r1abcd_h 2', iEnd, maxCore
        stop
    end if

end subroutine p_Wabdj_r1abcd_h

! 2026.7.1
! modified at 2026.7.10
!> Wabdj += fact*<ab||dc>*r(c,j)
subroutine p_Wabdj_r1abcd_h_2(Wabdj, r1, f_abcd, pop_a, pop_b, pop_d, pop_j, pop_c, &
                              symR, symV, iop, fact, iCore, maxCore, type)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: Wabdj
    real*8, dimension(*), intent(in) :: r1
    character(len=8), intent(in) :: f_abcd
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_d, pop_j, pop_c
    integer*8, intent(in) :: symR, symV ! symmetry of r1 and abcd
    integer*8, intent(in) :: iop !> iop=iopab=iopdc
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: type !> 1 - (aB,dC), 2 - (bA,cD)

    integer*8 :: symW
    integer*8 :: irrepd, numd, d
    integer*8 :: sym1, sym2
    integer*8 :: i0, i1, i2, i3, nsize, idx, iEnd
    integer*8, external :: idsymsoc, irpdso

    symW = dirprd(symR, symV)
    i0 = 1
    ! cycle for symmetry of d
    do irrepd = 1, nirrep
        numd = pop_d(irrepd)
        sym1 = dirprd(irrepd, symV) ! symmetry of V(ab,*c)
        sym2 = dirprd(irrepd, symW) ! symmetry of W(ab,*j)
        ! cycle for d
        do d = 1, numd
            ! get W(ab,*j), i0
            nsize = idsymsoc(sym2, pop_a, pop_b, iop, zorb, pop_j, 1)
            i1 = i0 + nsize
            call getpqs_incore(iCore(i0), d, irrepd, pop_a, pop_b, iop, &
                               pop_d, pop_j, 1, Wabdj, symW)
            if (type == 1) then
                ! get V(ab,*c), i1
                nsize = idsymsoc(sym1, pop_a, pop_b, iop, zorb, pop_c, 1)
                i2 = i1 + nsize
                call getpqs(iCore(i1), d, irrepd, pop_a, pop_b, iop, &
                            pop_d, pop_c, iop, f_abcd, symV)
                ! W(ab,*j) = V(ab,*c)*r(c,j)
                call vv2v_dpd_424(iCore(i1), pop_a, pop_b, iop, zorb, pop_c, 1, sym1, 'ab*c', &
                                  r1, pop_c, pop_j, 1, symR, 'cj', &
                                  iCore(i0), pop_a, pop_b, iop, zorb, pop_j, 1, sym2, 'ab*j', &
                                  fact, 1.0d0)
            else
                ! get V(bA,c*), i1
                nsize = idsymsoc(sym1, pop_b, pop_a, iop, pop_c, zorb, 1)
                i2 = i1 + nsize
                call getpqrso(iCore(i1), d, irrepd, 1, f_abcd)
                ! W(Ab,*j) = V(bA,c*)*r(c,j)
                call vv2v_dpd_424(iCore(i1), pop_b, pop_a, iop, zorb, pop_c, 1, sym1, 'ba*c', &
                                  r1, pop_c, pop_j, 1, symR, 'cj', &
                                  iCore(i0), pop_a, pop_b, iop, zorb, pop_j, 1, sym2, 'ab*j', &
                                  fact, 1.0d0)
            end if
            ! put W(ab,*j), i0
            call putpqs_incore(iCore(i0), d, irrepd, pop_a, pop_b, iop, &
                               pop_d, pop_j, 1, Wabdj, symW)
        end do
    end do

end subroutine p_Wabdj_r1abcd_h_2

! 2026.7.7
!> W(ab,dj) += fact*V(ab,dc)*r(c,j)
subroutine Wabdj_r1abcd_R(WabdjABAB, r1AA, abcdABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abcdABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,dJ) = V(aB,dC)*r(C,J) -> V(aB,dC)*r(c,j)
    call Wabdj_r1abcd(WabdjABAB, r1AA, abcdABAB, fact, iCore, maxCore)

end subroutine Wabdj_r1abcd_R

! 2026.7.7
!> W(ab,dj) += fact*V(ab,dc)*r(c,j)
subroutine Wabdj_r1abcd_U(WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA, spin, &
                          r1AA, r1BB, abcdAAAA, abcdABAB, abcdBBBB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabdjAAAA, WabdjABAB, WabdjBBBB, WabdjBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abcdAAAA, abcdABAB, abcdBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ab,dj) = V(ab,dc)*r(c,j)
        call Wabdj_r1abcd(WabdjAAAA, r1AA, abcdAAAA, fact, iCore, maxCore)
    end if
    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(aB,dJ) = V(aB,dC)*r(C,J)
        call Wabdj_r1abcd(WabdjABAB, r1BB, abcdABAB, fact, iCore, maxCore)
    end if
    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(AB,DJ) = V(AB,DC)*r(C,J)
        call Wabdj_r1abcd(WabdjBBBB, r1BB, abcdBBBB, fact, iCore, maxCore)
    end if
    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ab,Dj) = V(Ab,Dc)*r(c,j) = V(bA,cD)*r(c,j)
        if (WabdjBABA%incore) then
            if (abcdABAB%incore) then
                call vv2v_424(abcdABAB, r1AA, WabdjBABA, fact, 'bAcD', 'cj', 'AbDj', &
                              iCore, maxCore, 0)
            else
                call p_Wabdj_r1abcd_h_2(WabdjBABA%arr, r1AA%arr, abcdABAB%f, &
                                        WabdjBABA%popp, WabdjBABA%popq, WabdjBABA%popr, WabdjBABA%pops, &
                                        r1AA%popp, r1AA%sym, abcdABAB%sym, 1, fact, iCore, maxCore, 2)
            end if
        else
            write (6, '(A)') 'ERROR: Wabdj_r1abcd_U current not support iostate=2.'
            stop
        end if
    end if

end subroutine Wabdj_r1abcd_U
