! write by lishun, 2024.11.6
! calculate L2(ab,ij) += <0|[L1,V]|ab,ij>
! <0|[L1,V]|ab,ij> = -Pab{v(ij,kb)*l(k,a)}+Pij{v(cj,ab)*l(i,c)}
! L store as (v,o)
subroutine calL2d_VL1(iCore, symm, pL1, pL2, iStart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full space
    integer*8, intent(in) :: symm
    type(opt), intent(in) :: pL1
    type(opt), intent(inout) :: pL2
    integer*8, intent(in) :: iStart

    ! -Pab{v(ij,kb)*l(k,a)}
    call c2_r1ijka(pL2, pL1, ijka, 1.0d0, iCore(iStart), nCore - iStart + 1)

    ! Pij{v(cj,ab)*l(i,c)}
    call c2_r1abci(pL2, pL1, abci, 1.0d0, iCore(iStart), nCore - iStart + 1)

end subroutine calL2d_VL1
