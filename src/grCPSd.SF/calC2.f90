! write by lishun, 2024.10.18
! calculate constant2 term
!
subroutine calC2(iCore, C2, istart, order)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    type(opt), intent(inout) :: C2 !> C^(n) constant doubles operator
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> n>=2, the order of deltaT need calculate

    integer*8 :: i0
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart

    !----- initialize C2
    call C2%arr_clear

    ! calculate Au2v^(1)*Tv(1) term
    call calA1R2(iCore, C2, i0, order)

    if (order == 3) then
        call calL2_VTT(iCore, C2, i0, 1)
    end if

end subroutine calC2
