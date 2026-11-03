! write at 2025.2.28
! intermidate:
! Fac = fac-0.5*<kl||cd>*r(ad,kl)
! Fki = fki+0.5*<kl||cd>*r(cd,il)
! Wkaci = <ka||ci>+<kj||cb>*r(ba,ji)
subroutine run(iCore, maxCore)
    use mod_symm
    use mod_orbit
    use mod_variable
    use mod_davidson
    use mod_common
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: irrep
    integer*8 :: nSize, l1, l2, idx
    integer*8 :: i0, i1, iroot
    real*8 :: exE
    real*8 :: t0, t1, t2, t3
    real*8 :: norm
    integer*8, external :: idsymsoc, irpdso
    real*8, external :: dnrm2
    external :: calAR

    i0 = iFree

    call cpu_time(t0)
    !--------- calculate MP2 amplitudes, MP2t(ab,ij)=-<ab||ij>/e(ab,ij)
    call MP2t%arr_clear()
    call denomT2(MP2t, abij, iCore(iEorb))
    write (6, '(A,F)') 'norm of MP2 amplitudes:', MP2t%norm()

    !--------- Fac = fac-0.5*<kl||cd>*r(ad,kl)
    call Fac%arr_clear()
    call putDiagFock(Fac%arr, iCore(iEorb), popv, popo, 'V', 1.0d0)
    call Fae_r2abij(Fac, MP2t, abij, 1.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Fac:', Fac%norm()

    !--------- Fki = fki+0.5*<kl||cd>*r(cd,il)
    call Fki%arr_clear()
    call putDiagFock(Fki%arr, iCore(iEorb), popv, popo, 'O', 1.0d0)
    call Fmi_r2abij(Fki, MP2t, abij, 1.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Fki:', Fki%norm()

    !---------- Wkaci = <ka||ci>+<kj||cb>*r(ba,ji)
    call dcopy(iabj%len, iabj%arr, 1, Wkaci%arr, 1)
    call Wmbej_r2abij(Wkaci, MP2t, abij, 2.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Wkaci:', Wkaci%norm()
    call cpu_time(t1)
    write (6, *) 'cpu time for calculate intermidate MP2t, Fac, Fki and Wkaci:', t1 - t0

    call cpu_time(t0)
    !--------- prepare iniutial guess for Davidson
    call preDavidson1_ns(iCore(iEorb), popv, popo, iCore(i0), nCore - i0 + 1)
    do irrep = 1, nirrep
        if (estate_sym(irrep) == 0) then
            cycle
        end if
        l1 = irpdso(irrep, popv, popo, 1)
        l2 = idsymsoc(irrep, popv, popv, iop, popo, popo, iop)
        call dscal(1 + l1 + l2, 0.0d0, iCore(i0), 1)
        do iroot = 1, estate_sym(irrep)
            ! read
            idx = (1 + l1)*(iroot - 1) + 1
            call IO_stream(0.0d0, iCore(i0), idx, 1 + l1, f_init_davidson(irrep), '13', 'old')
            ! write
            idx = (1 + l1 + l2)*(iroot - 1) + 1
            if (iroot == 1) then
                call IO_stream(0.0d0, iCore(i0), idx, 1 + l1 + l2, f_R(irrep), '23', 'replace')
            else
                call IO_stream(0.0d0, iCore(i0), idx, 1 + l1 + l2, f_R(irrep), '23', 'old')
            end if
        end do
    end do

    !--------- calculate (OED1,OED2) as diagonal elements need by Davidson
    do irrep = 1, nirrep
        nSize = irpdso(irrep, popv, popo, 1)
        i1 = i0 + nSize
        ! OED1, i0
        call Eorb_diff1(iCore(i0), iCore(iEorb), popv, popo, irrep)
        ! OED2, i1
        call Eorb_diff2(iCore(i1), iCore(iEorb), popv, popv, popo, popo, irrep, iop)

        nSize = irpdso(irrep, popv, popo, 1) + idsymsoc(irrep, popv, popv, iop, popo, popo, iop)
        call IO_stream(0.0d0, iCore(i0), 1, nSize, f_Diag(irrep), '23', 'replace')
    end do
    call cpu_time(t1)
    write (6, '(A,F)') 'CPU time for prepare Davidson iteration:', t1 - t0

    call cpu_time(t0)
    write (6, '(A)') '--------------------------------------------------'
    write (6, '(A)') '                CIS(Doo) calculate                '
    write (6, '(A)') '--------------------------------------------------'
    do irrep = 1, nirrep
        write (6, '(A,I2)') '-------------- symmetry', irrep
        call cpu_time(t2)
        !--------- read in diagnoal elements
        nSize = irpdso(irrep, popv, popo, 1) + idsymsoc(irrep, popv, popv, iop, popo, popo, iop)
        i1 = i0 + nSize
        call IO_stream(0.0d0, iCore(i0), 1, nSize, f_Diag(irrep), '13', 'old')
        !-------- start Davidson iteration
        call Davidson(iCore, nCore, i1, f_R(irrep), iCore(i0), &
                      nSize, estate_sym(irrep), 1, irrep, calAR)

        call cpu_time(t3)
        write (6, '(A,I2,A,F)') 'CPU time for CIS(Doo) symmetry ', irrep, ':', t3 - t2
    end do
    call cpu_time(t1)
    write (6, '(A,F)') 'totaL CPU time for CIS(Doo) calculation:', t1 - t0

    !--------- print result
    write (6, '(A)') '----------------------------------------------'
    write (6, '(A)') '                CIS(Doo) result               '
    write (6, '(A)') '----------------------------------------------'
    do irrep = 1, nirrep
        if (estate_sym(irrep) == 0) then
            cycle
        end if
        write (6, '(A,I5)') 'excitation energy for symmetry:', irrep
        call R1%init(iCore, popv, popo, zorb, zorb, 1, 1, irrep, T, i0)
        i1 = i0 + R1%len
        call R2%init(iCore, popv, popv, popo, popo, iop, iop, irrep, T, i0)

        do iroot = 1, estate_sym(irrep)
            idx = (iroot - 1)*(1 + R1%len + R2%len) + 1
            call IO_stream(exE, R1%arr, idx, R1%len + R2%len, f_R(irrep), '11', 'old')
            write (6, '(A,I5,A,F12.6,A,F12.6,A)') 'state', iroot, ':', &
                exE, 'a.u.', exE*27.21138, 'eV'

            write (6, '(A,F)') 'percentage of sigle excitation:', R1%norm()**2/1.0d0
            call printEXinfo1_uhf(R1, 10)
            call printEXinfo2_uhf(R2, 10)
        end do
    end do

end subroutine run
