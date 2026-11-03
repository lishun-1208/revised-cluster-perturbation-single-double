! write by lishun, 2024.10.22
! calculate correlation energy E^(k+1)
subroutine calEcorr(E, iCore, maxCore, istart, T2, order)
    use mod_variable
    implicit none

    real*8, intent(out) :: E !> E^(k+1)
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    type(opt), intent(in) :: T2 !> T2^(k)
    integer*8, intent(in) :: order !> k

    integer*8 :: i0
    real*8, external :: ddot

    i0 = istart
    if (order <= 3) then
        ! E^(k+1) = 0.25*V(kl,cd)*t(cd,kl)^(k)
        E = ddot(T2%len, T2%arr, 1, abij%arr, 1)
    end if

end subroutine calEcorr
