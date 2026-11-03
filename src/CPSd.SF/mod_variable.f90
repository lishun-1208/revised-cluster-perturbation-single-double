! 2026.7.4
module mod_variable
    use mod_operator
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
    integer*8 :: davidson_cyc = 50 !> max cycle for Davidson method
    integer*8 :: davidson_dim = 20 !> max dimension of Davidson method
    integer*8 :: davidson_conv = 8 !> Davidson convergenced standard
    integer*8 :: iostate = 1 !> iostate for vvvo and vvvv
    logical :: vvvo_incore = T
    logical :: vvvv_incore = F

    !--------------------------------------------
    !      common use variables
    !--------------------------------------------
    integer*8 :: nCore ! total size of work space iCore
    integer*8 :: iFree ! start position of free work space
    ! total number of (alpha/beta) virtual and occupied orbital
    integer*8 :: nVir, nOcc
    ! number of virtual and occupied orbital for irrep, alias of nv, no or vrta/vrtb, popa/popb in mod_orbit
    integer*8, dimension(8) :: popv, popo
    real*8 :: Escf ! SCF energy
    ! orbital energy
    integer*8 :: iEorb

    !------ Fock matrix, fia store as (a,i)
    type(opt) :: fij, fab, fia

    !------ 2e-integrals, iabj store as (bi,aj)
    ! for closed shell, I(pq,rs)=I(pQ,rS)+I(pQ,Rs)
    ! so, only ijklABAB, ijkaABAB, abijABAB, iabjABAB, iabjBAAB, abciABAB, and abcdABAB are necessary
    ! but for convient, I also use ijklAAAA, ijkaAAAA, abijAAAA, iabjAAAA, abciAAAA, and abcdAAAA
    !--------------------------------------------
    type(opt) :: ijkl, ijka, abij, iabj, abci, abcd

    !------ file to store integral
    ! fij
    character(len=8), parameter :: f_fij = 'fij_ns__'
    ! fab
    character(len=8), parameter :: f_fab = 'fab_ns__'
    ! fia, sort as (a,i)
    character(len=8), parameter :: f_fia = 'fia_ns__'
    ! ijkl
    character(len=8), parameter :: f_ijkl = 'ijkl_ns_'
    ! ijka
    character(len=8), parameter :: f_ijka = 'ijka_ns_'
    ! abij
    character(len=8), parameter :: f_abij = 'abij_ns_'
    ! iabj, sort as (bi,aj)
    character(len=8), parameter :: f_iabj = 'iabj_ns_'
    ! abci
    character(len=8), parameter :: f_abci = 'abci_ns_'
    ! abcd
    character(len=8), parameter :: f_abcd = 'abcd_ns_'

    !--------------------------------------------
    !       private
    !--------------------------------------------
    real*8 :: Ecorr
    type(opt) :: T1, T2, C1, C2, R1
    type(opt) :: Fli, Fad, Fld, Wkaci, Wklij
    ! file to store (T1,T2)
    character(len=8) :: f_T = 'CPSd_T__'

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

    end subroutine get_variable

end module mod_variable
