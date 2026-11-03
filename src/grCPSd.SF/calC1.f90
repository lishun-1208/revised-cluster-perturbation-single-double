! write by lishun, 2024.10.18
! calculate constant1 term
!
subroutine calC1(iCore, C1, istart, order)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    type(opt), intent(inout) :: C1 !> C1^(n)
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> n>=2, the order of deltaT need calculate

    integer*8 :: i
    integer*8 :: i0, i1, i2, iEnd, nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart

    ! initialize C1
    call dscal(C1%len, 0.0d0, C1%arr, 1)
    if (order == 3) then
        ! calculate Au1v^(2)*Tv^(n-2) term
        call calA2R1(iCore, C1, istart, order)

        call calL1_VTT(iCore, C1, istart, 1)
    end if

end subroutine calC1
