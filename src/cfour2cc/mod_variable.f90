! 2026.7.4
module mod_variable
    implicit none
    !--------------------------------------------
    !      common use variables
    !--------------------------------------------
    !------ 2e-integrals, iabj store as (bi,aj)
    ! for closed shell, I(pq,rs)=I(pQ,rS)+I(pQ,Rs)
    ! so, only ijklABAB, ijkaABAB, abijABAB, iabjABAB, iabjBAAB, abciABAB, and abcdABAB are necessary
    ! but for convient, I also use ijklAAAA, ijkaAAAA, abijAAAA, iabjAAAA, abciAAAA, and abcdAAAA
    !--------------------------------------------
    !> ij, ab, ia(a,i)
    character(len=2), dimension(2), parameter :: spin_1e = ['AA', & ! RHF
                                                            'BB'] ! UHF
    !> ijkl, abij, abcd
    character(len=4), dimension(3), parameter :: spin_2e_1 = ['ABAB', 'AAAA', & ! RHF
                                                              'BBBB'] ! UHF
    !> ijka, abci, kaij(ij,ka), ciab(ab,ci)
    character(len=4), dimension(4), parameter :: spin_2e_2 = ['ABAB', 'AAAA', & ! RHF
                                                              'BBBB', 'BABA'] ! UHF
    !> iabj(bi,aj), for iabj(bi,aj) /= jbai(aj,bi)
    character(len=4), dimension(6), parameter ::  spin_2e_3 = ['AAAA', 'ABAB', 'AABB', & ! RHF
                                                               'BBBB', 'BABA', 'BBAA'] ! UHF
    !------ file to store integral
    ! fij
    character(len=8), dimension(2), parameter :: f_fij = ['FIJAA___', 'FIJBB___']
    ! fab
    character(len=8), dimension(2), parameter :: f_fab = ['FABAA___', 'FABBB___']
    ! fia, sort as (a,i)
    character(len=8), dimension(2), parameter :: f_fia = ['FAIAA___', 'FAIBB___']
    ! ijkl
    character(len=8), dimension(3), parameter :: f_ijkl = ['IJKLAB__', 'IJKLAA__', &
                                                           'IJKLBB__']
    ! ijka
    character(len=8), dimension(4), parameter :: f_ijka = ['IJKAAB__', 'IJKAAA__', &
                                                           'IJKABB__', 'IJKABA__']
    ! abij
    character(len=8), dimension(3), parameter :: f_abij = ['ABIJAB__', 'ABIJAA__', &
                                                           'ABIJBB__']
    ! iabj, sort as (bi,aj)
    character(len=8), dimension(5), parameter :: f_iabj = ['BIAJAAAA', 'BIAJABAB', 'BIAJAABB', &
                                                           'BIAJBBBB', 'BIAJBABA']
    ! abci
    character(len=8), dimension(4), parameter :: f_abci = ['ABCIAB__', 'ABCIAA__', &
                                                           'ABCIBB__', 'ABCIBA__']
    ! abcd
    character(len=8), dimension(3), parameter :: f_abcd = ['ABCDAB__', 'ABCDAA__', &
                                                           'ABCDBB__']

    !--------------------------------------------
    !       private
    !--------------------------------------------
    logical :: UHF = .false. !> whether reference be UHF state
    integer*8 :: symX, symY, symZ

end module mod_variable
