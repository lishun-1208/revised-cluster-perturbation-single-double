! write at 2025.6.25, lishun
!> copy file generate by clmoioso0
subroutine fcopy(filein, fileout, popr, pops, ioprs, iCore)
    use mod_symm
    implicit none

    character(len=8), intent(in) :: filein, fileout
    integer*8, dimension(8), intent(in) :: popr, pops
    integer*8, intent(in) :: ioprs
    real*8, dimension(*), intent(inout) :: iCore !> free work space

    integer*8 :: irreprs, numrs, rs
    integer*8, external :: irpdso

    do irreprs = 1, nirrep
        numrs = irpdso(irreprs, popr, pops, ioprs)
        do rs = 1, numrs
            call getlstso(iCore, rs, 1, irreprs, 1, filein)
            call putlstso(iCore, rs, 1, irreprs, 1, fileout)
        end do
    end do
end subroutine fcopy
