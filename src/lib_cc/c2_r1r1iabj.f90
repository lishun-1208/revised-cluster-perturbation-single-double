! 2026.4.16
!> c2(ab,ij) -= fact*PabPij{r(e,i)*r(a,m)*Vmbej(em,bj)}
subroutine c2_r1r1iabj(c2, r1, iabj, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: iabj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call p_c2_r1r1iabj(c2%arr, r1%arr, iabj%arr, c2%popp, c2%popq, c2%popr, c2%pops, &
                       r1%popp, r1%popq, c2%ioppq, r1%sym, iCore, maxCore, fact)

end subroutine c2_r1r1iabj

! c2(ab,ij) -= fact*PabPij{r1(e,i)*r1(a,m)*<mb||ej>}
! <ia||bj> store as (bi,aj)
!--------------------------
! modified at 2024.11.11
!
subroutine p_c2_r1r1iabj(c2, r1, iabj, pop_a, pop_b, pop_i, pop_j, pop_e, pop_m, &
                         iop, symR, iCore, maxCore, fact)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: c2
    real*8, dimension(*), intent(in) :: r1, iabj
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j, pop_e, pop_m
    integer*8, intent(in) :: iop !> iopab=iopij=iop
    integer*8, intent(in) :: symR
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    real*8, intent(in) :: fact

    integer*8 :: nSize
    integer*8 :: i0, i1, i2, i3, iEnd, itmp
    integer*8, external :: idsymsoc

    i0 = 1
    !------------ first sum e, then m
    !------------ tmp1=r1(e,i)*<mb||ej>=t(e,i)iabj(bj,em)
    ! iabj(bj,em)->iabj(bj,me), i0
    nSize = idsymsoc(1, pop_b, pop_j, 1, pop_e, pop_m, 1)
    i1 = i0 + nSize
    call isymtrso(iabj, iCore(i0), pop_b, pop_j, 1, 0, pop_e, pop_m, 1, 1, 1)
    ! tmp1(bj,mi)=(bj,me)*r1(e,i), i1, o3v2
    nSize = idsymsoc(1, pop_b, pop_j, 1, pop_m, pop_i, 1)
    i2 = i1 + nSize
    call VpqrmVms_to_Vpqrs(iCore(i0), r1, 1, symR, pop_e, &
                           pop_b, pop_j, 1, pop_m, pop_i, iCore(i1), symR, 1.0D0, 0.0D0)
    !------------ tmp2=r1(a,m)*tmp1(bj,mi)
    ! tmp1(bj,mi)->tmp1(bj,im), i2
    i3 = i2 + nSize
    call isymtrso(iCore(i1), iCore(i2), pop_b, pop_j, 1, 0, pop_m, pop_i, 1, 1, symR)
    ! tmp2(bj,ia)=tmp1(bj,im)*r1(a,m), i0, o3v2
    call VpqrmVsm_to_Vpqrs(iCore(i2), r1, symR, symR, pop_m, &
                           pop_b, pop_j, 1, pop_i, pop_a, iCore(i0), 1, 1.0D0, 0.0D0)
    !----------- c2
    ! tmp2(bj,ia)->tmp2(ab,ji), i1
    nSize = idsymsoc(1, pop_b, pop_j, 1, pop_i, pop_a, 1)
    i2 = i1 + nSize
    call sstgenso(iCore(i0), iCore(i1), nSize, pop_b, pop_j, pop_i, pop_a, iCore(i2), 1, "4123")
    ! tmp2(ab,ji)->tmp2(ab,ij), i0
    call isymtrso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_j, pop_i, 1, 1, 1)
    ! tmp3(ab,ij)=tmp2(ab,ij)-tmp2(ab,ji), i0
    call daxpy(nSize, -1.0D0, iCore(i1), 1, iCore(i0), 1)
    itmp = i0
    if (iop == 0) then
        ! compress ij, i1
        nSize = idsymsoc(1, pop_a, pop_b, 1, pop_i, pop_j, 0)
        i2 = i1 + nSize
        call icompso(iCore(i0), iCore(i1), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 0, 1)
        itmp = i1
    end if
    ! tmp3(ab,ij)->tmp3(ba,ij), i2
    call isymtrso(iCore(itmp), iCore(i2), pop_a, pop_b, 1, 1, pop_i, pop_j, iop, 0, 1)
    ! tmp4(ab,ij)=tmp3(ab,ij)-tmp3(ba,ij), itmp
    nSize = idsymsoc(1, pop_a, pop_b, 1, pop_i, pop_j, iop)
    call daxpy(nSize, -1.0D0, iCore(i2), 1, iCore(itmp), 1)
    if (iop == 0) then
        ! compress ab, i2
        call icompso(iCore(itmp), iCore(i2), pop_a, pop_b, 1, 0, pop_i, pop_j, 0, 1, 1)
        itmp = i2
    end if
    ! c2 -= fact*tmp4(ab,ij)
    nSize = idsymsoc(1, pop_a, pop_b, iop, pop_i, pop_j, iop)
    call daxpy(nSize, -fact, iCore(itmp), 1, c2, 1)

    iEnd = i3
    if (iEnd >= maxCore) then
        write (6, *) 'iCore not enough for p_c2_r1r1iabj'
        stop
    end if

end subroutine p_c2_r1r1iabj
