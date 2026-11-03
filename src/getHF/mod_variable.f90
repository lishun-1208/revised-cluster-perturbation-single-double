! 2026.7.1
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
    integer*8 :: iEorb, iEa, iEb

    !------ Fock matrix, fia store as (a,i)
    type(opt) :: fij, fab, fia
    type(opt) :: fijAA, fijBB
    type(opt) :: fabAA, fabBB
    type(opt) :: fiaAA, fiaBB

    !------ 2e-integrals, iabj store as (bi,aj)
    ! for closed shell, I(pq,rs)=I(pQ,rS)+I(pQ,Rs)
    ! so, only ijklABAB, ijkaABAB, abijABAB, iabjABAB, iabjBAAB, abciABAB, and abcdABAB are necessary
    ! but for convient, I also use ijklAAAA, ijkaAAAA, abijAAAA, iabjAAAA, abciAAAA, and abcdAAAA
    !--------------------------------------------
    type(opt) :: ijkl, ijka, abij, iabj, abci, abcd
    type(opt) :: ijklAAAA, ijklABAB, jklBBBB
    type(opt) :: ijkaAAAA, ijkaABAB, ijkaBABA, ijkaBBBB
    type(opt) :: abijAAAA, abijABAB, abijBBBB
    type(opt) :: iabjAAAA, iabjABAB, iabjBAAB, iabjABBA, iabjBBBB
    type(opt) :: abciAAAA, abciABAB, abciBABA, abciBBBB
    type(opt) :: abcdAAAA, abcdABAB, abcdBBBB

    !------ file to store integral
    ! fij
    character(len=8), parameter :: file_fij = 'fij_ns__'
    character(len=8), parameter :: file_fijAA = 'FIJAA___'
    character(len=8), parameter :: file_fijBB = 'FIJBB___'
    ! fab
    character(len=8), parameter :: file_fab = 'fab_ns__'
    character(len=8), parameter :: file_fabAA = 'FABAA___'
    character(len=8), parameter :: file_fabBB = 'FABBB___'
    ! fia, sort as (a,i)
    character(len=8), parameter :: file_fia = 'fia_ns__'
    character(len=8), parameter :: file_fiaAA = 'FAIAA___'
    character(len=8), parameter :: file_fiaBB = 'FAIBB___'
    ! ijkl
    character(len=8), parameter :: file_ijkl = 'ijkl_ns_'
    character(len=8), parameter :: file_ijklAAAA = 'IJKLAA__'
    character(len=8), parameter :: file_ijklABAB = 'IJKLAB__'
    character(len=8), parameter :: file_ijklBBBB = 'IJKLBB__'
    ! ijka
    character(len=8), parameter :: file_ijka = 'ijka_ns_'
    character(len=8), parameter :: file_ijkaAAAA = 'IJKAAA__'
    character(len=8), parameter :: file_ijkaABAB = 'IJKAAB__'
    character(len=8), parameter :: file_ijkaBABA = 'IJKABA__'
    character(len=8), parameter :: file_ijkaBBBB = 'IJKABB__'
    ! abij
    character(len=8), parameter :: file_abij = 'abij_ns_'
    character(len=8), parameter :: file_abijAAAA = 'ABIJAA__'
    character(len=8), parameter :: file_abijABAB = 'ABIJAB__'
    character(len=8), parameter :: file_abijBBBB = 'ABIJBB__'
    ! iabj, sort as (bi,aj)
    character(len=8), parameter :: file_iabj = 'iabj_ns_'
    character(len=8), parameter :: file_iabjAAAA = 'BIAJAAAA'
    character(len=8), parameter :: file_iabjABAB = 'BIAJAABB' ! iabjABAB
    character(len=8), parameter :: file_iabjBAAB = 'BIAJABAB' ! iabjBAAB
    character(len=8), parameter :: file_iabjABBA = 'BIAJBABA' ! iabjABBA
    character(len=8), parameter :: file_iabjBBBB = 'BIAJBBBB'
    ! abci
    character(len=8), parameter :: file_abci = 'abci_ns_'
    character(len=8), parameter :: file_abciAAAA = 'ABCIAA__'
    character(len=8), parameter :: file_abciABAB = 'ABCIAB__'
    character(len=8), parameter :: file_abciBABA = 'ABCIBA__'
    character(len=8), parameter :: file_abciBBBB = 'ABCIBB__'
    ! abcd
    character(len=8), parameter :: file_abcd = 'abcd_ns_'
    character(len=8), parameter :: file_abcdAAAA = 'ABCDAA__'
    character(len=8), parameter :: file_abcdABAB = 'ABCDAB__'
    character(len=8), parameter :: file_abcdBBBB = 'ABCDBB__'

    !--------------------------------------------
    !       private
    !--------------------------------------------
    logical :: spin_free = .true.
    logical :: RHF ! reference state whether be RHF state

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
