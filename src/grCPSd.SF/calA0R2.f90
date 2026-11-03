! write by lishun, 2024.10.15
! calculate Au2v^(0)*tv += -Pab{V(bl,ji)*t(a,l)}+Pij{V(ba,jd)*t(d,i)}
!
subroutine calA0R2(iCore, istart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0, maxCore

    i0 = istart
    maxCore = nCore - i0 + 1

    !------ -Pab{V(bl,ji)*t(a,l)}
    call c2_r1ijka(R2, T1, ijka, 1.0d0, iCore(i0), maxCore)

    !----- Pij{V(ba,jd)*t(d,i)}
    call c2_r1abci(R2, T1, abci, 1.0d0, iCore(i0), maxCore)

end subroutine calA0R2
