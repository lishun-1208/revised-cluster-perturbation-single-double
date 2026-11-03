! 2026.7.13
!> c(a,i) -= fact*W(ca,di)*F(c,d)
subroutine c1d_FcdWdica(c1d, Fcd, Wdica, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: Fcd
    type(opt), intent(in) :: Wdica
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Wdica%incore) then
        call vv2v_422(Wdica, Fcd, c1d, -fact, 'cadi', 'cd', 'ai', iCore, maxCore)
    else
        call p_c1d_FcdWdica_h(c1d%arr, Fcd%arr, Wdica%f, c1d%popq, c1d%popp, &
                              Fcd%popp, Fcd%popq, Wdica%ioppq, Fcd%sym, Wdica%sym, &
                              fact, iCore, maxCore)
    end if

end subroutine c1d_FcdWdica

! write at 2025.7.4, lishun
! calculate Cia -= fact*Wdica*Fcd
! Wdica and Cia sort as (ca,di) and (a,i)
subroutine p_c1d_FcdWdica_h(c1d, Fcd, file_Wdica, pop_i, pop_a, &
                            pop_c, pop_d, iop, symF, symW, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: c1d
    real*8, dimension(*), intent(in) :: Fcd
    character(len=8), intent(in) :: file_Wdica
    integer*8, dimension(8), intent(in) :: pop_i, pop_a, pop_c, pop_d
    integer*8, intent(in) :: iop ! iop=iopca
    integer*8, intent(in) :: symF, symW
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symC, sym1, sym2
    integer*8 :: irrepi, numi, i
    integer*8 :: i0, i1, i2, i3, i4, nSize, iEnd
    integer*8, external :: irpdso, idsymsoc

    i0 = 1
    iEnd = 0
    symC = dirprd(symF, symW)
    ! cycle for i
    do irrepi = 1, nirrep
        numi = pop_i(irrepi)
        sym1 = dirprd(symF, irrepi) ! symmetry if C(a,*)
        sym2 = dirprd(symW, irrepi) ! symmetry of W(ca,d*)
        do i = 1, numi
            ! get C(a,*), i0
            nSize = pop_a(sym1)
            i1 = i0 + nSize
            call getpqrso_incore(iCore(i0), i, irrepi, &
                                 zorb, zorb, 1, pop_a, pop_i, 1, c1d, symC)
            ! get W(ca,d*), i1
            nSize = idsymsoc(sym2, pop_c, pop_a, iop, pop_d, zorb, 1)
            i2 = i1 + nSize
            call getpqrso(iCore(i1), i, irrepi, 1, file_Wdica)
            ! W(ca,d*)->W(cd,a*), i2
            nSize = idsymsoc(sym2, pop_c, pop_a, 1, pop_d, zorb, 1)
            i3 = i2 + nSize
            if (iop == 0) then
                ! expand ca, i3
                i4 = i3 + nSize
                call iexpso(iCore(i1), iCore(i3), pop_c, pop_a, 0, 1, pop_d, zorb, 1, 0, sym2)
                call sstgenso(iCore(i3), iCore(i2), nSize, pop_c, pop_a, pop_d, zorb, &
                              iCore(i4), sym2, '1324')
                iEnd = max(iEnd, i4)
            else
                call sstgenso(iCore(i1), iCore(i2), nSize, pop_c, pop_a, pop_d, zorb, &
                              iCore(i3), sym2, '1324')
                iEnd = max(iEnd, i3)
            end if
            ! C(**,a*) -= fact*F(cd,**)*W(cd,a*)
            call VmnpqVmnrs_to_Vpqrs(Fcd, iCore(i2), symF, sym2, pop_c, pop_d, 1, &
                                     zorb, zorb, 1, pop_a, zorb, 1, iCore(i0), sym1, -fact, 1.0d0)
            ! put C(a,*)
            call putpqrso_incore(iCore(i0), i, irrepi, &
                                 zorb, zorb, 1, pop_a, pop_i, 1, c1d, symC)
        end do
    end do

    if (iEnd > maxCore) then
        write (6, *) 'ERROR: memory not enough for p_c1d_FcdWdica_h', iEnd, maxCore
        stop
    end if

end subroutine p_c1d_FcdWdica_h

! 2026.7.13
!> c(a,i) -= fact*W(ca,di)*F(c,d)
subroutine c1d_FcdWdica_R(c1dAA, FcdAA, WdicaAAAA, WdicaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA
    type(opt), intent(in) :: FcdAA
    type(opt), intent(in) :: WdicaAAAA, WdicaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -W(ca,di)*F(c,d)
    !          -W(Ca,Di)*F(C,D) -> -W(cA,dI)*F(C,D)
    call c1d_FcdWdica(c1dAA, FcdAA, WdicaAAAA, fact, iCore, maxCore)
    call c1d_FcdWdica(c1dAA, FcdAA, WdicaABAB, fact, iCore, maxCore)

end subroutine c1d_FcdWdica_R

! 2026.7.13
!> c(a,i) -= fact*W(ca,di)*F(c,d)
subroutine c1d_FcdWdica_U(c1dAA, c1dBB, FcdAA, FcdBB, &
                          WdicaAAAA, WdicaABAB, WdicaBBBB, WdicaBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA, c1dBB
    type(opt), intent(in) :: FcdAA, FcdBB
    type(opt), intent(in) :: WdicaAAAA, WdicaABAB, WdicaBBBB, WdicaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -W(ca,di)*F(c,d)
    !          -W(Ca,Di)*F(C,D)
    call c1d_FcdWdica(c1dAA, FcdAA, WdicaAAAA, fact, iCore, maxCore)
    call c1d_FcdWdica(c1dAA, FcdBB, WdicaBABA, fact, iCore, maxCore)

    ! c(A,I) = -W(cA,dI)*F(c,d)
    !          -W(CA,DI)*F(C,D)
    call c1d_FcdWdica(c1dBB, FcdAA, WdicaABAB, fact, iCore, maxCore)
    call c1d_FcdWdica(c1dBB, FcdBB, WdicaBBBB, fact, iCore, maxCore)

end subroutine c1d_FcdWdica_U
