! write by lishun, 2024.12.10
! calculate AR for Davidson iteraton
!----------------------------------------
! for CIS:
! AR(a,i) = <a,i|[H,R1]|0>=e(a,i)*r(a,i)+<a,i|[V,R1]|0>
subroutine calAR(iCore, maxCore, istart, iAR, iR, symR)
    use mod_variable
    use mod_orbit
    use mod_operator
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    integer*8, intent(in) :: iAR
    integer*8, intent(in) :: iR
    integer*8, intent(in) :: symR

    integer*8 :: i0, i1, i2, nSize

    i0 = istart
    call R1%init(iCore(iR), popv, popo, zorb, zorb, 1_8, 1_8, symR, idx=1)
    call AR1%init(iCore(iAR), popv, popo, zorb, zorb, 1_8, 1_8, symR, idx=1)

    ! e(a,i)*r(a,i)
    i1 = i0 + R1%len
    call getallso(iCore(i0), 1, f_OED1(symR))
    i2 = i1 + R1%len
    call vdmul(R1%len, iCore(i0), R1%arr, iCore(i1))
    call dcopy(R1%len, iCore(i1), 1, AR1%arr, 1)

    ! <a,i|[V,R1]|0>
    call c1_r1iabj(AR1, R1, iabj, 1.0d0, iCore(i0), nCore - i0 + 1)

end subroutine calAR
