! lishun, 2015.10.11
program cfour2cc
    use mod_variable
    use mod_interface
    implicit none

    integer*8 :: icore
    common//icore(1)
    integer*8 :: iintln, ifltln, iintfp, ialone, ibitwd
    common/machsp/iintln, ifltln, iintfp, ialone, ibitwd
    integer*8 :: nocco, nvrto
    common/info/nocco(2), nvrto(2)
    integer*8 :: i0
    common/istart/i0
    integer*8 :: icrsiz, ichcsz, ioff
    common/iopos/icrsiz, ichcsz, ioff(2)
    integer :: symrx, symry, symrz ! symmerty of Rx, Ry and Rz
    common/symrxyz/symrx, symry, symrz
    integer*8 :: nstart, nirrep, irreps, dirprd
    common/syminf/nstart, nirrep, irreps(255, 2), dirprd(8, 8)
    integer*8 :: pop, vrt, nt, nfmi, nfea
    common/sym/pop(8, 2), vrt(8, 2), nt(2), nfmi(2), nfea(2)

    character(len=6) :: cmpgrp ! computational group
    integer*8 :: iuhf
    integer*8 :: i
    real*8, dimension(:), allocatable::Ea, Eb
    integer*8 :: noa, nva, nob, nvb
    real*8 :: Escf
    character(len=3) :: string
    logical :: exist

    call read_input_file
    call grep_node('ref_state', exist)
    if (exist) then
        call read_node('ref_state', string(1:3))
        if (string(1:3) == 'rhf') then
            UHF = .false.
        else if (string(1:3) == 'uhf') then
            UHF = .true.
        else
            write (6, *) 'ERROR: value of ref_state error, please use lower case.', string(1:3)
            stop
        end if
    else
        write (6, *) 'ERROR: please add ref_state=uhf/rhf in CCSOINP'
        stop
    end if

    call crapsi(icore, iuhf, 0) ! cfour main program

    write (6, '(A)') '------------------------------------'
    write (6, '(A)') '       data save in interface       '

    call add_node('iprcln', 32768)
    call add_node('iprcwd', 4096)
    call add_node('ichcsz', 4096)
    call add_node('frozen_orbit', 0)

    !--- determine group
    call getrec(20, 'JOBARC', 'COMPPTGP', 1, cmpgrp)
    write (6, '(2A)') 'computational group: ', cmpgrp
    write (6, '(A,I2)') 'number of irreducible representation:', nirrep
    write (6, '(A)') 'direct product table:'
    do i = 1, 8
        write (6, '(8I3)') dirprd(i, :)
    end do
    call add_node('group_n', cmpgrp, 6)
    call add_node('nirrep', nirrep)
    call add_node('dirprd', dirprd, 8, 8)

    !--- irrep name
    if (cmpgrp(1:3) == 'C1 ') then
        write (6, '(A)') 'irrep name: A'
        call add_node('irrep_name1', 'A', 1)
    else if (cmpgrp(1:3) == 'Cs ') then
        write (6, '(3A)') 'irrep name:', ' A''', ' A'''''
        call add_node('irrep_name1', 'A''', 2)
        call add_node('irrep_name2', 'A''''', 3)
    else if (cmpgrp(1:3) == 'Ci ') then
        write (6, '(A)') 'irrep name: Ag Au'
        call add_node('irrep_name1', 'Ag', 2)
        call add_node('irrep_name2', 'Au', 2)
    else if (cmpgrp(1:3) == 'C2 ') then
        write (6, '(A)') 'irrep name: A B'
        call add_node('irrep_name1', 'A', 1)
        call add_node('irrep_name2', 'B', 1)
    else if (cmpgrp(1:3) == 'C2v') then
        write (6, '(A)') 'irrep name: A1 B1 B2 A2'
        call add_node('irrep_name1', 'A1', 2)
        call add_node('irrep_name2', 'B1', 2)
        call add_node('irrep_name3', 'B2', 2)
        call add_node('irrep_name4', 'A2', 2)
    else if (cmpgrp(1:3) == 'C2h') then
        write (6, '(A)') 'irrep name: Ag Au Bu Bg'
        call add_node('irrep_name1', 'Ag', 2)
        call add_node('irrep_name2', 'Au', 2)
        call add_node('irrep_name3', 'Bu', 2)
        call add_node('irrep_name4', 'Bg', 2)
    else if (cmpgrp(1:3) == 'D2 ') then
        write (6, '(A)') 'irrep name: A B2 B1 B3'
        call add_node('irrep_name1', 'A', 1)
        call add_node('irrep_name2', 'B2', 2)
        call add_node('irrep_name3', 'B1', 2)
        call add_node('irrep_name4', 'B3', 2)
    else if (cmpgrp(1:3) == 'D2h') then
        write (6, '(A)') 'irrep name: Ag B2u B3u B1g B1u B3g B2g Au'
        call add_node('irrep_name1', 'Ag', 2)
        call add_node('irrep_name2', 'B2u', 3)
        call add_node('irrep_name3', 'B3u', 3)
        call add_node('irrep_name4', 'B1g', 3)
        call add_node('irrep_name5', 'B1u', 3)
        call add_node('irrep_name6', 'B3g', 3)
        call add_node('irrep_name7', 'B2g', 3)
        call add_node('irrep_name8', 'Au', 2)
    end if

    !--- determine symR
    if (cmpgrp(1:3) == 'C2 ' .or. cmpgrp(1:3) == 'Cs ') then
        symrx = 2
        symry = 2
        symrz = 1
    else if (cmpgrp(1:3) == 'C2v') then
        symrx = 3
        symry = 2
        symrz = 4
    end if
    if (cmpgrp(1:3) == 'C2h') then
        symrx = 4
        symry = 4
        symrz = 1
    else if (cmpgrp(1:3) == 'D2 ') then
        symrx = 4
        symry = 2
        symrz = 3
    else if (cmpgrp(1:3) == 'D2h') then
        symrx = 7
        symry = 6
        symrz = 4
    else if (cmpgrp(1:3) == 'C1 ' .or. cmpgrp(1:3) == 'Ci ') then
        symrx = 1
        symry = 1
        symrz = 1
    end if
    write (6, '(A,3I2)') 'symR(x,y,z):', symrx, symry, symrz
    call add_node('symrx', symrx)
    call add_node('symry', symry)
    call add_node('symrz', symrz)

    !--- determine symX
    if (cmpgrp(1:3) == 'C1 ') then
        symX = 1
        symY = 1
        symZ = 1
    else if (cmpgrp(1:3) == 'Cs ') then
        symX = 1
        symY = 1
        symZ = 2
    else if (cmpgrp(1:3) == 'Ci ') then
        symX = 2
        symY = 2
        symZ = 2
    else if (cmpgrp(1:3) == 'C2 ') then
        symX = 2
        symY = 2
        symZ = 1
    else if (cmpgrp(1:3) == 'C2v') then
        symX = 2
        symY = 3
        symZ = 1
    else if (cmpgrp(1:3) == 'C2h') then
        symX = 3
        symY = 3
        symZ = 2
    else if (cmpgrp(1:3) == 'D2 ') then
        symX = 4
        symY = 2
        symZ = 3
    else if (cmpgrp(1:3) == 'D2h') then
        symX = 3
        symY = 2
        symZ = 5
    end if
    write (6, '(A,3I2)') 'sym(x,y,z):', symrx, symry, symrz
    call add_node('symx', symX)
    call add_node('symy', symY)
    call add_node('symz', symZ)

    !--- population
    write (6, '(A)') 'pop:'
    write (6, '(8I5)') pop(:, 1)
    write (6, '(8I5)') pop(:, 2)
    write (6, '(A)') 'vrt:'
    write (6, '(8I5)') vrt(:, 1)
    write (6, '(8I5)') vrt(:, 2)
    call add_node('pop', pop, 8, 2)
    call add_node('vrt', vrt, 8, 2)

    !--- energy
    noa = sum(pop(1:nirrep, 1))
    nob = sum(pop(1:nirrep, 2))
    nva = sum(vrt(1:nirrep, 1))
    nvb = sum(vrt(1:nirrep, 2))
    allocate (Ea(noa + nva), Eb(nob + nvb))
    call getrec(20, 'JOBARC', 'SCFENEG ', iintfp, Escf)
    write (6, '(A,F18.10)') 'SCF energy:', Escf
    call getrec(20, 'JOBARC', 'SCFEVALA', iintfp*(noa + nva), Ea)
    write (6, '(A)') 'orbital energy alpha:'
    write (6, *) Ea
    call getrec(20, 'JOBARC', 'SCFEVALB', iintfp*(nob + nvb), Eb)
    write (6, '(A)') 'orbital energy beta:'
    write (6, *) Eb
    call add_node('Escf', Escf)
    call add_node('orb_energy_alpha_sym', Ea, noa + nva)
    call add_node('orb_energy_beta_sym', Eb, nob + nvb)

    write (6, '(A)') '------------------------------------'

    ! molecular integral
    write (6, '(A)') 'start gengeate ERI and Fock matrix file'
    call run(icore(i0), icrsiz)

    write (6, '(A)') 'cfour2cc finished'

end program cfour2cc
