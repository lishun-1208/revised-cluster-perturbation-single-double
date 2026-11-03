! write by lishun, 2024.11.5
subroutine run(iCore)
    use mod_symm
    use mod_orbit
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full

    real*8 :: tmp, tmp1
    integer*8 :: idx
    integer*8 :: symm
    integer*8 :: iroot
    integer*8 :: i0, nSize
    real*8 :: tt0, tt1, tt2, tt3
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    call cpu_time(tt0)
    i0 = iFree

    !----- calculate orbital energy difference, e(ab,ij)=e(a)+e(b)-e(i)-e(j)
    do symm = 1, nirrep
        ! OED1
        call clmoioso0(symm, iCore(i0), popv, popo, 1, zorb, zorb, 1, f_OED1(symm))
        call Eorb_diff1(iCore(i0), iCore(iEorb), popv, popo, symm)
        call putallso(iCore(i0), 1, f_OED1(symm))
        if (debug) then
            nSize = irpdso(symm, popv, popo, 1)
            write (6, '(A,I2,A,F)') 'norm of OED1 for symmetry', symm, ':', &
                dnrm2(nSize, iCore(i0), 1)
        end if
        ! OED2
        call clmoioso0(symm, iCore(i0), popv, popv, iop, popo, popo, iop, f_OED2(symm))
        call Eorb_diff2(iCore(i0), iCore(iEorb), popv, popv, popo, popo, symm, iop)
        call putallso(iCore(i0), 1, f_OED2(symm))
        if (debug) then
            nSize = idsymsoc(symm, popv, popv, iop, popo, popo, iop)
            write (6, '(A,I2,A,F)') 'norm of OED2 for symmetry', symm, ':', &
                dnrm2(nSize, iCore(i0), 1)
        end if
    end do

    !----- calculate MP2t, t(ab,ij) = -v(ab,ij)/e(ab,ij)
    call denomT2(MP2t, abij, iCore(iEorb))
    write (6, '(A,F)') 'norm of MP2t:', MP2t%norm()

    !-------- calculate Fli = -0.5*v(kl,cd)*MP2(cd,ki)
    call Fli%arr_clear
    call Fmi_r2abij(Fli, MP2t, abij, -1.0d0, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of Fli:', Fli%norm()
    end if

    !------- calculate Fad = 0.5*v(kl,cd)*MP2(ca,kl)
    call Fad%arr_clear
    call Fae_r2abij(Fad, MP2t, abij, -1.0d0, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of Fad:', Fad%norm()
    end if

    !------- calculate Wkaci = v(ka,ci)+v(kl,cd)*MP2(da,li), sort as (ck,ai)
    call dcopy(iabj%len, iabj%arr, 1, Wkaci%arr, 1)
    call Wmbej_r2abij(Wkaci, MP2t, abij, 2.0d0, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of Wkaci:', Wkaci%norm()
    end if

    call cpu_time(tt1)
    write (6, '(A,F)') 'cpu time for calculate intermidate:', tt1 - tt0

    ! calculate CIS(Doo) left eigenvector
    write (6, '(A)') '----------------------------------'
    write (6, '(A)') '             CIS(Doo)             '
    write (6, '(A)') '----------------------------------'
    call cpu_time(tt0)
    do symm = 1, nirrep
        do iroot = 1, estate_sym(symm)
            call cpu_time(tt2)
            write (6, '(A,I2,A,I4)') '----------- symmetry', symm, ' state', iroot
            call calEigenCISDoo(iCore, i0, symm, iroot)
            call cpu_time(tt3)
            write (6, '(A,F)') 'cpu time:', tt3 - tt2
        end do
    end do
    call cpu_time(tt1)
    write (6, '(A,F)') 'total cpu time for CIS(Doo):', tt1 - tt0

    ! EE-grCPS(D) calculate
    write (6, '(A)') '-----------------------------------------'
    write (6, '(A)') '       start EE-grCPS(D) calculate        '
    write (6, '(A)') '-----------------------------------------'
    !------ calculate J^(2) and J^(3)
    call cpu_time(tt0)
    call calJ2(iCore, i0)
    call calJ3(iCore, i0)
    call cpu_time(tt1)
    write (6, '(A,F)') 'cpu time for calculate J^(2) and J^(3)', tt1 - tt0

    call cpu_time(tt0)
    do symm = 1, nirrep
        do iroot = 1, estate_sym(symm)
            ! 2nd order energy
            call cpu_time(tt2)
            write (6, '(A,I2,A,I3,A)') '----------- symmetry', symm, ' state', iroot, ' order 2'
            call run1(iCore, symm, iroot, estate_sym(symm), i0)
            call cpu_time(tt3)
            write (6, '(A,F)') 'cpu time:', tt3 - tt2
            ! 3rd order energy
            call cpu_time(tt2)
            write (6, '(A,I2,A,I3,A)') '----------- symmetry', symm, ' state', iroot, ' order 3'
            call run2(iCore, symm, iroot, estate_sym(symm), i0)
            call cpu_time(tt3)
            write (6, '(A,F)') 'cpu time:', tt3 - tt2
        end do
    end do
    call cpu_time(tt1)
    write (6, '(A,F)') 'total cpu time for EE-grCPS(D) calculate', tt1 - tt0

    ! print result
    write (6, '(A)') '---------- EE-grCPS(D) excitation energy ----------'
    do symm = 1, nirrep
        do iroot = 1, estate_sym(symm)
            write (6, '(A,I2,A,I2,A)') '--- symmetry', symm, ' state', iroot, ':'
            nSize = 1 + irpdso(symm, popv, popo, 1) + idsymsoc(symm, popv, popv, iop, popo, popo, iop)

            idx = nSize*(iroot - 1) + 1
            call IO_stream(tmp, iCore(i0), idx, 0, f_L0(symm), '12', 'old')
            write (6, '(A,2F10.6)') 'CIS(Doo) (a.u./eV):', tmp, tmp*27.21138
            tmp1 = tmp

            idx = nSize*(iroot - 1) + 1
            call IO_stream(tmp, iCore(i0), idx, 0, f_R(symm), '12', 'old')
            write (6, '(A,F10.6,A)') '1st order:', tmp, 'a.u.'
            tmp1 = tmp1 + tmp
            write (6, '(A,2F10.6)') 'total excitation energy (a.u./eV):', tmp1, tmp1*27.21138

            idx = nSize*(estate_sym(symm) + iroot - 1) + 1
            call IO_stream(tmp, iCore(i0), idx, 0, f_R(symm), '12', 'old')
            write (6, '(A,F10.6,A)') '2nd order:', tmp, 'a.u.'
            tmp1 = tmp1 + tmp
            write (6, '(A,2F10.6)') 'total excitation energy (a.u./eV):', tmp1, tmp1*27.21138

            idx = nSize*(2*estate_sym(symm) + iroot - 1) + 1
            call IO_stream(tmp, iCore(i0), idx, 0, f_R(symm), '12', 'old')
            write (6, '(A,F10.6,A)') '3rd order:', tmp, 'a.u.'
            tmp1 = tmp1 + tmp
            write (6, '(A,2F10.6)') 'total excitation energy (a.u./eV):', tmp1, tmp1*27.21138

        end do
    end do
end subroutine run
