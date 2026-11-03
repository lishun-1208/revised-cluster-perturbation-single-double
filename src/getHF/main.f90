! write by lishun, 2024.11.15
! get fock matrix and 2e-integral in file
program getHF
    use mod_iop
    use mod_orbit
    use mod_symm
    use mod_interface
    use mod_variable
    implicit none

    real*8, dimension(:), allocatable :: iCore ! main work space
    integer*8 :: maxCore
    integer*8 :: status ! allocate stat
    ! cpu time
    real*8 :: tt0, tt1
    integer*8 :: i, irrep, i0
    logical :: exist
    character(len=5) :: string

    integer*8, external :: idsymsoc, irpdso

    call cpu_time(tt0)

    call prepare

    write (6, *) '-----------------------------------------'
    write (6, *) '           start allocate iCore          '
    write (6, *) '-----------------------------------------'
    allocate (iCore(nCore), stat=status)
    if (status /= 0) then
        write (6, *) "ERROR: allocate iCore error, memory not enough"
        stop
    end if

    !--- for RHF as reference
    call grep_node('ref_state', exist)
    if (exist) then
        call read_node('ref_state', string(1:3))
        if (string(1:3) == 'rhf') then
            RHF = .true.
        else if (string(1:3) == 'uhf') then
            RHF = .false.
        else
            write (6, *) 'ERROR: value of ref_state error, please use lower case.', string(1:3)
            stop
        end if
    else
        write (6, *) 'ERROR: please add ref_state=uhf/rhf in CCSOINP'
        stop
    end if

    if (RHF) then
        call generate_pure_alpha(iCore)
        write (6, *) 'finish generate pure alpha spin 2e-integral'
    end if

    inquire (file=trim(tmpdir)//'/'//file_fijAA, exist=exist)
    if (.not. exist) then
        call clmoioso0(1, iCore, popa, popa, 1, zorb, zorb, 1, file_fijAA)
    end if
    inquire (file=trim(tmpdir)//'/'//file_fabAA, exist=exist)
    if (.not. exist) then
        call clmoioso0(1, iCore, vrta, vrta, 1, zorb, zorb, 1, file_fabAA)
    end if
    inquire (file=trim(tmpdir)//'/'//file_fiaAA, exist=exist)
    if (.not. exist) then
        call clmoioso0(1, iCore, vrta, popa, 1, zorb, zorb, 1, file_fiaAA)
    end if
    if (.not. RHF) then
        inquire (file=trim(tmpdir)//'/'//file_fijBB, exist=exist)
        if (.not. exist) then
            call clmoioso0(1, iCore, popb, popv, 1, zorb, zorb, 1, file_fijBB)
        end if
        inquire (file=trim(tmpdir)//'/'//file_fabBB, exist=exist)
        if (.not. exist) then
            call clmoioso0(1, iCore, vrtb, vrtb, 1, zorb, zorb, 1, file_fabBB)
        end if
        inquire (file=trim(tmpdir)//'/'//file_fiaBB, exist=exist)
        if (.not. exist) then
            call clmoioso0(1, iCore, vrtb, popb, 1, zorb, zorb, 1, file_fiaBB)
        end if
    end if

    call grep_node('spin_free', exist)
    if (exist) then
        call read_node('spin_free', string)
        if (string(1:1) == 'f') then
            spin_free = .false.
        else if (string(1:1) == 't') then
            spin_free = .true.
        end if
    end if

    if (spin_free) then
        write (6, *) 'start generate no spin 2e_integral.'
        !-------------- calculate the index of variables
        ! for orbital energy
        iEorb = 1
        i0 = iEorb + nOcc + nVir
        ! for Fock matrix
        call fij%init(iCore, popo, popo, zorb, zorb, 1, 1, 1, T, i0, file_fij)
        i0 = i0 + fij%len
        call fab%init(iCore, popv, popv, zorb, zorb, 1, 1, 1, T, i0, file_fab)
        i0 = i0 + fab%len
        call fia%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, T, i0, file_fia) ! fia(a,i)
        i0 = i0 + fia%len
        ! for 2e-integral
        call ijkl%init(iCore, popo, popo, popo, popo, iop, iop, 1, T, i0, file_ijkl)
        i0 = i0 + ijkl%len
        call ijka%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0, file_ijka)
        i0 = i0 + ijka%len
        call abij%init(iCore, popv, popv, popo, popo, iop, iop, 1, T, i0, file_abij)
        i0 = i0 + abij%len
        call iabj%init(iCore, popv, popo, popv, popo, 1, 1, 1, T, i0, file_iabj) ! iabj(bi,aj)
        i0 = i0 + iabj%len
        ! for abci and abcd we store in hard disk
        if (vvvo_incore) then
            call abci%init(iCore, popv, popv, popv, popo, iop, 1, 1, T, i0, file_abci)
            i0 = i0 + abci%len
        else
            call abci%init(iCore, popv, popv, popv, popo, iop, 1, 1, F, file=file_abci)
        end if
        if (vvvv_incore) then
            call abcd%init(iCore, popv, popv, popv, popv, iop, iop, 1, T, i0, file_abcd)
            i0 = i0 + abcd%len
        else
            call abcd%init(iCore, popv, popv, popv, popv, iop, iop, 1, F, file=file_abcd)
        end if

        iFree = i0
        write (6, '(A,I15)') 'iFree =', iFree
        if (iFree >= nCore) then
            write (6, *) "ERROR: memory not enough for variables.", iFree, nCore
            stop
        end if

        maxCore = nCore - iFree ! memory can be used for later calculate
        write (6, '(A,I12)') 'maxCore =', maxCore
        write (6, '(A,F15.6)') 'free memory(GiB):', maxCore*8/(1024.d0**3)

        write (6, *) '-----------------------------------------'
        write (6, *) '               start loadSCF             '
        write (6, *) '-----------------------------------------'
        call loadSCF(iCore, iFree)
    end if

    call cpu_time(tt1)
    write (6, *) "total cpu time for getHF:", tt1 - tt0
end program getHF
