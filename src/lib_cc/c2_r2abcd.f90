! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ef,ij)*V(ab,ef)
subroutine c2_r2abcd(c2, r2, abcd, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abcd
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (abcd%incore) then
        if (r2%ioppq == 1) then
            call vv2v_444(abcd, r2, c2, 0.5d0*fact, 'abef', 'efij', 'abij', iCore, maxCore, 0)
        else
            call vv2v_444(abcd, r2, c2, fact, 'abef', 'efij', 'abij', iCore, maxCore, 0)
        end if
    else
        call p_c2_r2abcd_h(c2%arr, r2%arr, abcd%f, c2%popp, c2%popq, c2%popr, c2%pops, &
                           r2%popp, r2%popq, c2%ioppq, c2%ioprs, r2%ioppq, r2%ioprs, abcd%ioppq, abcd%ioprs, &
                           r2%sym, abcd%sym, fact, iCore, maxCore)
    end if

end subroutine c2_r2abcd

subroutine p_c2_r2abcd_h(c2, r2, f_abcd, pop_a, pop_b, pop_i, pop_j, pop_e, pop_f, &
                         iopCvv, iopCoo, iopRvv, iopRoo, iopV1, iopV2, &
                         symR, symV, fact, iCore, maxCore)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r2
    character(len=8), intent(in) :: f_abcd
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_e, pop_f
    integer*8, intent(in) :: iopCvv, iopCoo, iopRvv, iopRoo, iopV1, iopV2 ! iopV1=iopVab, iopV2=iopVef
    integer*8, intent(in) :: symR ! symmetry of r2
    integer*8, intent(in) :: symV ! symmetry of abcd
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free space
    integer*8, intent(in) :: maxCore

    real*8 :: fact2
    integer*8 :: symC
    integer*8 :: irrep_ab, num_ab, ab
    integer*8 :: sym1, sym2
    integer*8 :: i0, i1, i2, i3, iEnd
    integer*8 :: nSize
    integer*8, external :: idsymsoc, irpdso

    symC = dirprd(symR, symV)
    if (iopV1 /= iopCvv) then
        write (6, *) 'ERROR: p_c2_r2abcd_h current not support iopV1/=iopCvv.'
        stop
    end if
    if (iopRoo /= iopCoo) then
        write (6, *) 'ERROR: p_c2_r2abcd_h current not support iopRoo/=iopCoo.'
        stop
    end if
    if (iopV2 == 0) then
        fact2 = 1.0d0
    else if (iopV2 == 1) then
        fact2 = 0.5d0
    end if

    i0 = 1
    ! c(ab,ij)->c(ij,ab), i0
    nSize = idsymsoc(symC, pop_a, pop_b, iopCvv, pop_i, pop_j, iopCoo)
    i1 = i0 + nSize
    call itranspso(c2, iCore(i0), pop_a, pop_b, iopCvv, pop_i, pop_j, iopCoo, symC)

    ! cyclef for symmetry of ab
    do irrep_ab = 1, nirrep
        num_ab = irpdso(irrep_ab, pop_a, pop_b, iopCvv)
        sym1 = dirprd(symV, irrep_ab) ! symmerty of V(ef,**)
        sym2 = dirprd(symC, irrep_ab) ! symmetry of c(ij,**)
        ! cycle for ab
        do ab = 1, num_ab
            ! get V(ab,ef)=V(ef,ab) for given ab, i1
            nSize = irpdso(sym1, pop_e, pop_f, iopV2)
            i2 = i1 + nSize
            call getlstso(iCore(i1), ab, 1, irrep_ab, 1, f_abcd)
            ! get c(ij,ab) for given ab, i2
            nSize = irpdso(sym2, pop_i, pop_j, iopCoo)
            i3 = i2 + nSize
            call getlstso_incore(iCore(i2), ab, 1, irrep_ab, &
                                 pop_i, pop_j, iopCoo, pop_a, pop_b, iopCvv, iCore(i0), symC)
            ! c(ij,**) += fact*0.5*V(ef,**)*t(ef,ij)
            call VmnpqVmnrs_to_Vpqrs(r2, iCore(i1), symR, sym1, pop_e, pop_f, iopV1, &
                                     pop_i, pop_j, iopRoo, zorb, zorb, 1, iCore(i2), sym2, fact*fact2, 1.0d0)
            ! put c(ij,ab) for given ab
            call putlstso_incore(iCore(i2), ab, 1, irrep_ab, &
                                 pop_i, pop_j, iopCoo, pop_a, pop_b, iopCvv, iCore(i0), symC)
        end do
    end do

    ! c(ij,ab)->c(ab,ij)
    call itranspso(iCore(i0), c2, pop_i, pop_j, iopCoo, pop_a, pop_b, iopCvv, symC)

end subroutine p_c2_r2abcd_h

! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ef,ij)*V(ab,ef)
subroutine c2_r2abcd_R(c2ABAB, r2ABAB, abcdABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: abcdABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) += fact*r(eF,iJ)*V(aB,eF)
    call c2_r2abcd(c2ABAB, r2ABAB, abcdABAB, 2.0d0*fact, iCore, maxCore)

end subroutine c2_r2abcd_R

! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ef,ij)*V(ab,ef)
subroutine c2_r2abcd_U(c2AAAA, c2ABAB, c2BBBB, r2AAAA, r2ABAB, r2BBBB, &
                       abcdAAAA, abcdABAB, abcdBBBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abcdAAAA, abcdABAB, abcdBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) += 0.5*fact*r(ef,ij)*V(ab,ef)
    call c2_r2abcd(c2AAAA, r2AAAA, abcdAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) += fact*r(eF,iJ)*V(aB,eF)
    call c2_r2abcd(c2ABAB, r2ABAB, abcdABAB, 2.0d0*fact, iCore, maxCore)

    ! c(AB,IJ) += 0.5*fact*r(EF,IJ)*V(AB,EF)
    call c2_r2abcd(c2BBBB, r2BBBB, abcdBBBB, fact, iCore, maxCore)

end subroutine c2_r2abcd_U
