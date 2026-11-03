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
    type(opt) :: R1, R2 ! eigenvector
    type(opt) :: C1, C2 ! energy and eigenvector consistent
    ! J^(2) intermediates, Wooov=0, Fov=0, Woooo not need for CPS(D-3)
    type(opt) :: J2_Foo, J2_Fvv, J2_Wovvo, J2_Wovoo
    ! J^(3) intermediates, Fov, Wovoo, Woooo, Wooov not need for CPS(D-3)
    type(opt) :: J3_Foo, J3_Fvv, J3_Wovvo

    !----------------- pointor of variables store in iCore
    ! file of perturbation T: (T^(1),T^(2),...)
    character(len=8) :: f_T = 'CPSd_T__'
    ! file of EE-CPS(D): (e^(k+1),R^(k))
    character(len=9), dimension(8) :: f_R = (/'CPSd_sym1', 'CPSd_sym2', 'CPSd_sym3', 'CPSd_sym4', &
                                              'CPSd_sym5', 'CPSd_sym6', 'CPSd_sym7', 'CPSd_sym8'/)
    ! file of CIS (e,R)
    character(len=8), dimension(8) :: f_CIS = (/'CIS_sym1', 'CIS_sym2', 'CIS_sym3', 'CIS_sym4', &
                                                'CIS_sym5', 'CIS_sym6', 'CIS_sym7', 'CIS_sym8'/)
    ! file of orbital energy difference(OED)
    character(len=8), dimension(8) :: f_OED1 = (/'OED1sym1', 'OED1sym2', 'OED1sym3', 'OED1sym4', &
                                                 'OED1sym5', 'OED1sym6', 'OED1sym7', 'OED1sym8'/)
    character(len=8), dimension(8) :: f_OED2 = (/'OED2sym1', 'OED2sym2', 'OED2sym3', 'OED2sym4', &
                                                 'OED2sym5', 'OED2sym6', 'OED2sym7', 'OED2sym8'/)
    ! energy constant C: (C_1^(k),C_2^(k),...)
    character(len=10), dimension(8) :: f_C = (/'CPSdC_sym1', 'CPSdC_sym2', 'CPSdC_sym3', 'CPSdC_sym4', &
                                               'CPSdC_sym5', 'CPSdC_sym6', 'CPSdC_sym7', 'CPSdC_sym8'/)

end module mod_variable
