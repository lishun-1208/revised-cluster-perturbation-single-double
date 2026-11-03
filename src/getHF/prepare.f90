! 2026.1.26
subroutine prepare
    use mod_orbit
    use mod_interface
    use mod_variable
    implicit none
    
    write (6, *) '-----------------------------------------'
    write (6, *) '          start read input file          '
    write (6, *) '-----------------------------------------'
    call read_input_file

    write (6, *) '-----------------------------------------'
    write (6, *) '              start readmodule           '
    write (6, *) '-----------------------------------------'
    call readmodule
    call get_variable

    !--- read mod_variable
    ! for UHF as reference and not consider spin
    nOcc = noa + nob
    nVir = nva + nvb
    popo = no
    popv = nv
    ! for RHF as reference and consider spin
!    nOcc = noa
!    nVir = nva
!    popo = popa
!    popv = vrta

    nCore = memory*1024**3/8
    call read_node('Escf', Escf)
    write (6, '(A,I5)') 'nOcc =', nOcc
    write (6, '(A,I5)') 'nVir =', nVir
    write (6, '(A,8I5)') 'popo =', popo
    write (6, '(A,8I5)') 'popv =', popv
    write (6, '(A,I)') 'nCore =', nCore
    write (6, '(A,F)') 'SCF energy =', Escf

    
end subroutine prepare