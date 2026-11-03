subroutine allocate_mem(iCore)
    use mod_variable
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work space

    integer*8 :: i0
    integer*8 :: maxCore

    !-------------- calculate the index of variables
    ! for orbital energy
    iEorb = 1
    i0 = iEorb + nVir + nOcc
    ! for Fock matrix
    call fij%init(iCore, popo, popo, zorb, zorb, 1, 1, 1, T, i0, f_fij)
    i0 = i0 + fij%len
    call fab%init(iCore, popv, popv, zorb, zorb, 1, 1, 1, T, i0, f_fab)
    i0 = i0 + fab%len
    call fia%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0, f_fia)
    i0 = i0 + fia%len
    ! for 2e-integral
    call ijkl%init(iCore, popo, popo, popo, popo, iop, iop, 1, T, i0, f_ijkl)
    i0 = i0 + ijkl%len
    call ijka%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0, f_ijka)
    i0 = i0 + ijka%len
    call abij%init(iCore, popv, popv, popo, popo, iop, iop, 1, T, i0, f_abij)
    i0 = i0 + abij%len
    call iabj%init(iCore, popv, popo, popv, popo, 1, 1, 1, T, i0, f_iabj) ! iabj(bi,aj)
    i0 = i0 + iabj%len
    ! for abci and abcd we store in hard disk
    if (vvvo_incore) then
        call abci%init(iCore, popv, popv, popv, popo, iop, 1, 1, T, i0, f_abci)
        i0 = i0 + abci%len
    else
        call abci%init(iCore, popv, popv, popv, popo, iop, 1, 1, F, file=f_abci)
    end if
    if (vvvv_incore) then
        call abcd%init(iCore, popv, popv, popv, popv, iop, iop, 1, T, i0, f_abcd)
        i0 = i0 + abcd%len
    else
        call abcd%init(iCore, popv, popv, popv, popv, iop, iop, 1, F, file=f_abcd)
    end if
    ! intermidates
    call T1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + T1%len
    call T2%init(iCore, popv, popv, popo, popo, iop, iop, 1, T, i0)
    i0 = i0 + T2%len
    call C1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + C1%len
    call C2%init(iCore, popv, popv, popo, popo, iop, iop, 1, T, i0)
    i0 = i0 + C2%len
    call R1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + R1%len
    
    iFree = i0
    write (6, '(A,I15)') 'iFree =', iFree
    if (iFree >= nCore) then
        write (6, '(A,2I15)') "ERROR: memory not enough for variables", iFree, nCore
        stop
    end if

    maxCore = nCore - iFree ! memory can be used for later calculate
    write (6, '(A,I12)') 'maxCore =', maxCore
    write (6, '(A,F15.6)') 'free memory(GiB):', maxCore*8/(1024.d0**3)

end subroutine allocate_mem
