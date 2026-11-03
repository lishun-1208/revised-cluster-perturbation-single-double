! 2026.7.13
!> W(ab,cd) += fact*1/2*l(ab,ij)*r(cd,ij)
subroutine Wabcd_l2r2(Wabcd, l2, r2, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wabcd
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: r2
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wabcd%incore) then
        if (l2%ioprs == 1) then
            call vv2v_444(l2, r2, Wabcd, 0.5d0*fact, 'abij', 'cdij', 'abcd', iCore, maxCore, 0)
        else if (l2%ioprs == 0) then
            call vv2v_444(l2, r2, Wabcd, fact, 'abij', 'cdij', 'abcd', iCore, maxCore, 0)
        end if
    else
        call p_Wabcd_l2r2_h(Wabcd%f, l2%arr, r2%arr, &
                            Wabcd%popp, Wabcd%popq, Wabcd%popr, Wabcd%pops, &
                            r2%popr, r2%pops, Wabcd%ioppq, l2%sym, r2%sym, &
                            fact, iCore, maxCore)
    end if

end subroutine Wabcd_l2r2

! 2026.7.13
subroutine p_Wabcd_l2r2_h(f_Wabcd, l2, r2, pop_a, pop_b, pop_c, pop_d, &
                          pop_i, pop_j, iop, symL, symR, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    character(len=8), intent(in) :: f_Wabcd
    real*8, dimension(*), intent(in) :: l2, r2
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_c, pop_d, pop_i, pop_j
    integer*8, intent(in) :: iop !> iop=iopab=iopcd=iopij
    integer*8, intent(in) :: symL, symR
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    real*8 :: iopfact
    integer*8 :: symW, sym1, sym2
    integer*8 :: irrepcd, numcd, cd
    integer*8 :: i0, i1, i2, i3, nSize
    integer*8, external :: idsymsoc, irpdso

    if (iop == 0) then
        iopfact = 2.0d0
    else
        iopfact = 1.0d0
    end if
    symW = dirprd(symL, symR)
    i0 = 1
    ! r(cd,ij)->r(ij,cd), i0
    nSize = idsymsoc(symR, pop_c, pop_d, iop, pop_i, pop_j, iop)
    i1 = i0 + nSize
    call itranspso(r2, iCore(i0), pop_c, pop_d, iop, pop_i, pop_j, iop, symR)

    ! cycle for symmeyry of cd
    do irrepcd = 1, nirrep
        sym1 = dirprd(symW, irrepcd) ! symmetry of W(ab,**)
        sym2 = dirprd(symR, irrepcd) ! symmetry of r(ij,**)
        numcd = irpdso(irrepcd, pop_c, pop_d, iop)

        ! cycle for cd
        do cd = 1, numcd
            ! get W(ab,**), i1
            nSize = irpdso(sym1, pop_a, pop_b, iop)
            i2 = i1 + nSize
            call getlstso(iCore(i1), cd, 1, irrepcd, 1, f_Wabcd)
            ! get r(ij,**), i1
            nSize = irpdso(sym2, pop_i, pop_j, iop)
            i3 = i2 + nSize
            call getlstso_incore(iCore(i2), cd, 1, irrepcd, &
                                 pop_i, pop_j, iop, pop_c, pop_d, iop, iCore(i0), symR)
            ! W(ab,**) += fact*1/2*l(ab,ij)*r(ij,**)
            call vv2v_dpd_422(l2, pop_a, pop_b, iop, pop_i, pop_j, iop, symL, 'abij', &
                              iCore(i2), pop_i, pop_j, iop, sym2, 'ij', &
                              iCore(i1), pop_a, pop_b, iop, sym1, 'ab', &
                              0.5d0*fact*iopfact, 1.0d0)
            call putlstso(iCore(i1), cd, 1, irrepcd, 1, f_Wabcd)
        end do
    end do

end subroutine p_Wabcd_l2r2_h

! 2026.7.13
!> W(ab,cd) += fact*1/2*l(ab,ij)*r(cd,ij)
subroutine Wabcd_l2r2_R(WabcdABAB, l2ABAB, r2ABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabcdABAB
    type(opt), intent(in) :: l2ABAB
    type(opt), intent(in) :: r2ABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(aB,cD) = l(aB,iJ)*r(cD,iJ)
    call Wabcd_l2r2(WabcdABAB, l2ABAB, r2ABAB, 2.0d0*fact, iCore, maxCore)

end subroutine Wabcd_l2r2_R

! 2026.7.13
!> W(ab,cd) += fact*1/2*l(ab,ij)*r(cd,ij)
subroutine Wabcd_l2r2_U(WabcdAAAA, WabcdABAB, WabcdBBBB, &
                        l2AAAA, l2ABAB, l2BBBB, &
                        r2AAAA, r2ABAB, r2BBBB, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WabcdAAAA, WabcdABAB, WabcdBBBB
    type(opt), intent(in) :: l2AAAA, l2ABAB, l2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(ab,cd) = 1/2*l(ab,ij)*r(cd,ij)
    call Wabcd_l2r2(WabcdAAAA, l2AAAA, r2AAAA, fact, iCore, maxCore)

    ! W(aB,cD) = l(aB,iJ)*r(cD,iJ)
    call Wabcd_l2r2(WabcdABAB, l2ABAB, r2ABAB, 2.0d0*fact, iCore, maxCore)

    ! W(ab,cd) = 1/2*l(AB,IJ)*r(CD,IJ)
    call Wabcd_l2r2(WabcdBBBB, l2BBBB, r2BBBB, fact, iCore, maxCore)

end subroutine Wabcd_l2r2_U
