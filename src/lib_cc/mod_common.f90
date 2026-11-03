! 2026.7.4
module mod_common
    use mod_operator
    use mod_orbit
    implicit none
    !--------------------------------------------
    !             input parameters
    !--------------------------------------------
    logical, private :: exist

    integer*8 :: iop = 0 !> iop for all tensors, current only support 0

    integer*8 :: memory = 10 !> memory with GiB

    integer*8 :: diis_cyc = 50 !> max DIIS cycle
    integer*8 :: diis_dim = 10 !> max DIIS dimension
    integer*8 :: diis_conv = 8 !> DIIS convergenced standard
    real*8 :: diis_max_err = 10.0d-8 !> max DIIS error

    logical :: incore = F !> be depreciated, whether all calculation run incore

    logical :: debug = F !> whether use debug module

    ! control the resevered size of iCore when read date from hard disk
    integer*8 :: reserved_size = 1000000000 !> must more than nv**3

    integer*8, dimension(8) :: estate_sym !> population of excited state
    ! file to sotore init guess R1
    character(len=10), dimension(8), parameter :: f_init_davidson = &
                                                  (/'initS_sym1', 'initS_sym2', 'initS_sym3', 'initS_sym4', &
                                                    'initS_sym5', 'initS_sym6', 'initS_sym7', 'initS_sym8'/)

    integer*8 :: davidson_cyc = 50 !> max cycle for Davidson method
    integer*8 :: davidson_dim = 20 !> max dimension of Davidson method
    integer*8 :: davidson_conv = 8 !> Davidson convergenced standard
    integer*8 :: lowest = F !> whether use the lowest eigenvalue in Davidson

    integer*8 :: iostate = 1 !> iostate for vvvo and vvvv
    logical :: vvvo_incore = T
    logical :: vvvv_incore = F

    logical :: CVS = F !> whether use CVS approach
    integer*8 :: CVS_orb = 0 !> which orbital be frozen, sort by spatial orbital energy
    integer*8, dimension(:), allocatable :: core_sym, core_idx ! the symmetry and idx of the core virtual orbital

contains
    subroutine get_variable
        use mod_interface
        implicit none

        character(len=20) :: string

        call grep_node('memory', exist)
        if (exist) then
            call read_node('memory', memory)
        end if

        call grep_node('diis_cyc', exist)
        if (exist) then
            call read_node('diis_cyc', diis_cyc)
        end if

        call grep_node('diis_dim', exist)
        if (exist) then
            call read_node('diis_dim', diis_dim)
        end if

        call grep_node('diis_conv', exist)
        if (exist) then
            call read_node('diis_conv', diis_conv)
            diis_max_err = 10.0d0**(-diis_conv)
        end if

        call grep_node('incore', exist)
        if (exist) then
            call read_node('incore', string)
            if (string(1:5) == 'false') then
                incore = F
            else if (string(1:4) == 'true') then
                incore = T
            end if
        end if

        call grep_node('debug', exist)
        if (exist) then
            call read_node('debug', string)
            if (string(1:5) == 'false') then
                debug = F
            else if (string(1:4) == 'true') then
                debug = T
            end if
        end if

        call grep_node('iop', exist)
        if (exist) then
            call read_node('iop', iop)
        end if

        call grep_node('reserved_size', exist)
        if (exist) then
            call read_node('reserved_size', reserved_size)
        end if

        call grep_node('estate_sym', exist)
        if (exist) then
            call read_node('estate_sym', estate_sym, 8)
        end if

        call grep_node('davidson_cyc', exist)
        if (exist) then
            call read_node('davidson_cyc', davidson_cyc)
        end if

        call grep_node('davidson_dim', exist)
        if (exist) then
            call read_node('davidson_dim', davidson_dim)
        end if

        call grep_node('davidson_conv', exist)
        if (exist) then
            call read_node('davidson_conv', davidson_conv)
        end if

        call grep_node('iostate', exist)
        if (exist) then
            call read_node('iostate', iostate)
            if (iostate == 0) then
                vvvo_incore = T
                vvvv_incore = T
            else if (iostate == 1) then
                vvvo_incore = T
                vvvv_incore = F
            else if (iostate == 2) then
                vvvo_incore = F
                vvvv_incore = F
            end if
        end if

        call grep_node('cvs', exist)
        if (exist) then
            call read_node('cvs', CVS_orb)
            CVS = T
        end if

        call grep_node('lowest', exist)
        if (exist) then
            call read_node('lowest', string)
            if (string(1:1) == 't') then
                lowest = T
            else
                lowest = F
            end if
        end if

    end subroutine get_variable

end module mod_common
