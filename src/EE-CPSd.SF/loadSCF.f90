subroutine loadSCF(iCore)
    use mod_common
    use mod_variable
    use mod_interface
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work space

    integer*8 :: maxCore

    maxCore = nCore - iFree + 1

    !--- Eorb
    call read_node('orb_energy_no_spin', iCore(iEorb:iEorb + nVir + nOcc - 1), nVir + nOcc)
    call printEorb_ns(iCore(iEorb))
    !--- fock matrix
    ! fij
    call getallso(fij%arr, 1, fij%f)
    write (6, *) 'norm of fij:', fij%norm(iCore(iFree), maxCore)
    ! fab
    call getallso(fab%arr, 1, fab%f)
    write (6, *) 'norm of fab:', fab%norm(iCore(iFree), maxCore)
    ! fia
    call getallso(fia%arr, 1, fia%f)
    write (6, *) 'norm of fia:', fia%norm(iCore(iFree), maxCore)
    !--- 2e-integral
    write (6, *) 'iop =', iop
    ! ijkl
    call getallso(ijkl%arr, 1, ijkl%f)
    write (6, *) 'norm of ijkl:', ijkl%norm(iCore(iFree), maxCore)
    ! ijka
    call getallso(ijka%arr, 1, ijka%f)
    write (6, *) 'norm of ijka:', ijka%norm(iCore(iFree), maxCore)
    ! abij
    call getallso(abij%arr, 1, abij%f)
    write (6, *) 'norm of abij:', abij%norm(iCore(iFree), maxCore)
    ! iabj
    call getallso(iabj%arr, 1, iabj%f)
    write (6, *) 'norm of iabj:', iabj%norm(iCore(iFree), maxCore)
    ! abci
    if (abci%incore) then
        call getallso(abci%arr, 1, abci%f)
    end if
    write (6, *) 'norm of abci:', abci%norm(iCore(iFree), maxCore)
    ! abcd
    if (abcd%incore) then
        call getallso(abcd%arr, 1, abcd%f)
    end if
    write (6, *) 'norm of abcd:', abcd%norm(iCore(iFree), maxCore)

end subroutine loadSCF
