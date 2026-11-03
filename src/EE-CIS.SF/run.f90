! write by lishun, 2024.12.9
subroutine run(iCore, maxCore)
    use mod_symm
    use mod_orbit
    use mod_variable
    use mod_davidson
    use mod_common
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symm, iroot, nroot
    integer*8 :: i0, i1, nSize, idx
    real*8 :: exE ! excitation energy
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2
    external :: calAR

    i0 = iFree

    !----- calculate OED1
    do symm = 1, nirrep
        call clmoioso0(symm, iCore(i0), popv, popo, 1, zorb, zorb, 1, f_OED1(symm))
        call Eorb_diff1(iCore(i0), iCore(iEorb), popv, popo, symm)
        call putallso(iCore(i0), 1, f_OED1(symm))
    end do
    write (6, '(A)') 'finish calculate orbital energy difference'

    !----- prepare initial guess
    call preDavidson1_ns(iCore(iEorb), popv, popo, iCore(i0), nCore - i0 + 1)
    do symm = 1, nirrep
        if (estate_sym(symm) == 0) then
            cycle
        end if
        nSize = estate_sym(symm)*(1 + irpdso(symm, popv, popo, 1))
        call IO_stream(0.0d0, iCore(i0), 1, nsize, f_init_davidson(symm), '13', 'old')
        call IO_stream(0.0d0, iCore(i0), 1, nsize, f_R(symm), '23', 'replace')
    end do
    do symm = 1, nirrep
        write (6, '(A,I2,A)') '------ symmetry', symm, ' ------'
        call R1%init(iCore, popv, popo, zorb, zorb, 1_8, 1_8, symm, idx=i0)
        nroot = estate_sym(symm)
        do iroot = 1, nroot
            write (6, '(A,I3)') '------ state', iroot
            idx = (1 + R1%len)*(iroot - 1) + 1
            call IO_stream(exE, R1%arr, idx, R1%len, f_R(symm), '11', 'old')
            write (6, '(A,F12.8,A,F12.5,A)') 'excitation energy:', exE, 'a.u.', exE*27.21138d0, 'eV'
            call printEXinfo1_uhf(R1, 10)
        end do
    end do
    write (6, '(A)') 'finish prepare initial guess'

    write (6, *) '------------------------------------'
    write (6, *) '      start Davidson iteration      '
    write (6, *) '------------------------------------'
    iroot = 1
    do symm = 1, nirrep
        write (6, '(A,I2)') '------------ symmetry', symm
        ! get OED1 as diag elements, i0
        nSize = irpdso(symm, popv, popo, 1)
        i1 = i0 + nSize
        call getallso(iCore(i0), 1, f_OED1(symm))
!        call calDiag(iCore(i0), symm, iCore)

        nroot = estate_sym(symm)
        call Davidson(iCore, maxCore, i1, f_R(symm), &
                      iCore(i0), nSize, nroot, 1, symm, calAR)
        iroot = iroot + estate_sym(symm)
    end do

    write (6, *) '---------------------------------'
    write (6, *) '      excitation information     '
    write (6, *) '---------------------------------'
    do symm = 1, nirrep
        nroot = estate_sym(symm)
        if (nroot == 0) then
            cycle
        end if
        write (6, '(A,I2,A)') '------ symmetry', symm, ' ------'
        call R1%init(iCore, popv, popo, zorb, zorb, 1_8, 1_8, symm, idx=i0)
        do iroot = 1, nroot
            write (6, '(A,I3)') '------ state', iroot
            idx = (1 + R1%len)*(iroot - 1) + 1
            call IO_stream(exE, R1%arr, idx, R1%len, f_R(symm), '11', 'old')
            write (6, '(A,F12.8,A,F12.5,A)') 'excitation energy:', exE, 'a.u.', exE*27.21138d0, 'eV'
            call printEXinfo1_uhf(R1, 10)
        end do
    end do

end subroutine run
