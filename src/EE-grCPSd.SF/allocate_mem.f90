subroutine allocate_mem(iCore)
    use mod_orbit
    use mod_common
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work space

    integer*8 :: i0, maxCore

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
    call MP2t%init(iCore, popv, popv, popo, popo, iop, iop, 1, T, i0)
    i0 = i0 + MP2t%len
    call Wkaci%init(iCore, popv, popo, popv, popo, 1, 1, 1, T, i0)
    i0 = i0 + Wkaci%len
    call Fli%init(iCore, popo, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + Fli%len
    call Fad%init(iCore, popv, popv, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + Fad%len
    call J2_Foo%init(iCore, popo, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J2_Foo%len
    call J2_Fvv%init(iCore, popv, popv, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J2_Fvv%len
    call J2_Fov%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J2_Fov%len
    call J2_Woooo%init(iCore, popo, popo, popo, popo, iop, iop, 1, T, i0)
    i0 = i0 + J2_Woooo%len
    call J2_Wooov%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0)
    i0 = i0 + J2_Wooov%len
    call J2_Wovoo%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0)
    i0 = i0 + J2_Wovoo%len
    call J2_Wovvo%init(iCore, popv, popo, popv, popo, 1, 1, 1, T, i0)
    i0 = i0 + J2_Wovvo%len
    call J3_Foo%init(iCore, popo, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J3_Foo%len
    call J3_Fvv%init(iCore, popv, popv, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J3_Fvv%len
    call J3_Fov%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0)
    i0 = i0 + J3_Fov%len
    call J3_Woooo%init(iCore, popo, popo, popo, popo, iop, iop, 1, T, i0)
    i0 = i0 + J3_Woooo%len
    call J3_Wooov%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0)
    i0 = i0 + J3_Wooov%len
    call J3_Wovoo%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0)
    i0 = i0 + J3_Wovoo%len
    call J3_Wovvo%init(iCore, popv, popo, popv, popo, 1, 1, 1, T, i0)
    i0 = i0 + J3_Wovvo%len

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
