! write by lishun, 2024.11.6
! calculate L1(a,i) += fact*<0|[L2,V]|a,i>
! <0|[L2,V]|a,i> = -0.5*v(bi,jk)*l(ba,jk)+0.5*v(bc,ja)*l(bc,ji)
subroutine calL1d_VL2(iCore, symm, pL1, pL2, iStart, fact)
    use mod_orbit
    use mod_symm
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: symm
    type(opt), intent(inout) :: pL1
    type(opt), intent(in) :: pL2
    integer*8, intent(in) :: iStart
    real*8, intent(in) :: fact

    ! -0.5*v(bi,jk)*l(ba,jk)
    call c1_r2ijka(pL1, pL2, ijka, fact, iCore(iStart), nCore - iStart + 1)

    ! 0.5*v(bc,ja)*l(bc,ji)
    call c1_r2abci(pL1, pL2, abci, fact, iCore(iStart), nCore - iStart + 1)

end subroutine calL1d_VL2
