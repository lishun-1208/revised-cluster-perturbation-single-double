! write by lishun, 2024.10.15
! calculate Au2v^(0)*tv += -Pab{V(bl,ji)*t(a,l)}+Pij{V(ba,jd)*t(d,i)}
!
subroutine calA0R2(iCore, istart, pR2, pT1)
    use mod_variable
    use mod_operator
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pR2
    type(opt), intent(in) :: pT1

    integer*8 :: i0, maxCore

    i0 = istart
    maxCore = nCore - i0 + 1

    !------ -Pab{V(bl,ji)*t(a,l)}
    call c2_r1ijka(pR2, pT1, ijka, 1.0d0, iCore(i0), maxCore)

    !----- Pij{V(ba,jd)*t(d,i)}
    call c2_r1abci(pR2, pT1, abci, 1.0d0, iCore(i0), maxCore)

end subroutine calA0R2
