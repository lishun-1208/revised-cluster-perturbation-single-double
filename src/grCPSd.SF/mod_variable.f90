! write at 2025.2.28
module mod_variable
    use mod_operator
    implicit none

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
    character(len=8), parameter :: f_T = 'iCPSd_T_'
    real*8 :: Ecorr !> total correlation energy
    type(opt) :: T1, T2, R1, R2, C1, C2
    type(opt) :: MP2t, Wkaci, Fli, Fad
end module mod_variable
