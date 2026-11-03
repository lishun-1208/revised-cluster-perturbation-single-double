! write by lishun, 2024.12.6
! c(ab,ij) += fact*<ab,ij|[V,R]|0> = <ab,ij|[V,R1]|0>+<ab,ij|[V,R2]|0>
subroutine calC2_VR(iCore, maxCore, istart, C2, R1, R2, symR, fact, order)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    type(opt), intent(inout) :: C2
    type(opt), intent(in) :: R1, R2
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact
    integer*8, intent(in) :: order !> n, C^(n)

    if (order > 2) then
        call calC2_VR1(iCore, maxCore, istart, C2, R1, symR, fact)
    end if
    call calC2_VR2(iCore, maxCore, istart, C2, R2, symR, fact)

end subroutine calC2_VR
