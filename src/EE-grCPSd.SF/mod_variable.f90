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
    ! J^(2) intermediates
    type(opt) :: MP2t, Wkaci, Fli, Fad ! MP2 intermediates
    type(opt) :: J2_Foo, J2_Fvv, J2_Fov, J2_Wovoo, J2_Wovvo, J2_Woooo, J2_Wooov
    ! J^(3) intermediates
    type(opt) :: J3_Foo, J3_Fvv, J3_Fov, J3_Wovoo, J3_Wovvo, J3_Woooo, J3_Wooov

    !----------------- pointor of variables store in iCore
    ! file of perturbation T: (T^(1),T^(2),...)
    character(len=8) :: f_T = 'iCPSd_T_'
    ! file of EE-grCPS(D): (e^(k+1),R^(k))
    character(len=10), dimension(8) :: f_R = (/'iCPSd_sym1', 'iCPSd_sym2', 'iCPSd_sym3', 'iCPSd_sym4', &
                                               'iCPSd_sym5', 'iCPSd_sym6', 'iCPSd_sym7', 'iCPSd_sym8'/)
    ! file of CIS(Doo) (e,R)
    character(len=11), dimension(8) :: f_R0 = (/'CISdoo_sym1', 'CISdoo_sym2', 'CISdoo_sym3', 'CISdoo_sym4', &
                                                'CISdoo_sym5', 'CISdoo_sym6', 'CISdoo_sym7', 'CISdoo_sym8'/)
    ! file of CIS(Doo) (e,L)
    character(len=12), dimension(8) :: f_L0 = (/'CISdooL_sym1', 'CISdooL_sym2', 'CISdooL_sym3', 'CISdool_sym4', &
                                                'CISdooL_sym5', 'CISdooL_sym6', 'CISdooL_sym7', 'CISdool_sym8'/)
    ! file of orbital energy difference(OED)
    character(len=8), dimension(8) :: f_OED1 = (/'OED1sym1', 'OED1sym2', 'OED1sym3', 'OED1sym4', &
                                                 'OED1sym5', 'OED1sym6', 'OED1sym7', 'OED1sym8'/)
    character(len=8), dimension(8) :: f_OED2 = (/'OED2sym1', 'OED2sym2', 'OED2sym3', 'OED2sym4', &
                                                 'OED2sym5', 'OED2sym6', 'OED2sym7', 'OED2sym8'/)
    ! energy constant C: (C_1^(k),C_2^(k),...)
    character(len=11), dimension(8) :: f_C = (/'iCPSdC_sym1', 'iCPSdC_sym2', 'iCPSdC_sym3', 'iCPSdC_sym4', &
                                               'iCPSdC_sym5', 'iCPSdC_sym6', 'iCPSdC_sym7', 'iCPSdC_sym8'/)

end module mod_variable
