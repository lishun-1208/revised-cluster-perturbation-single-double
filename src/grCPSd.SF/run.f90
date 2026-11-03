! write by lishun, 2024.10.14
! 0 and 1 order Ecorr is 0
! n order Ecorr need (n-1) order T2 and (n-2) order T1
! E^(n)=0.25*V(cd,kl)*{t(cd,kl)+2*t(c,k)*t(l,d)}^(n-1)
!
subroutine run(iCore, istart)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart

    integer*8 :: i0
    real*8 :: tt0, tt1
    integer*8 :: nSize
    integer*8, external :: idsymsoc, irpdso
    real*8, external :: dnrm2

    write (6, *) 'max cycle for cc iteration:', diis_cyc
    write (6, *) 'cc convergence criteria:', diis_max_err
    write (6, *) 'max DIIS dimension:', diis_dim
    i0 = istart

    !------------- calculate intermidate
    call cpu_time(tt0)

    !-------- clculate MP2 amplitudes
    call denomT2(MP2t, abij, iCore(iEorb))
    write (6, '(A,F)') 'norm of MP2t:', MP2t%norm()
    !-------- calculate Fli = -0.5*v(kl,cd)*MP2(cd,ki)
    call Fli%arr_clear
    call Fmi_r2abij(Fli, MP2t, abij, -1.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Fli:', Fli%norm()
    !------- calculate Fad = 0.5*v(kl,cd)*MP2(ca,kl)
    call Fad%arr_clear
    call Fae_r2abij(Fad, MP2t, abij, -1.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Fad:', Fad%norm()
    !------- calculate Wkaci = V(ka,ci)+V(kl,cd)*MP2(da,li), sort as (ck,ai)
    call dcopy(Wkaci%len, iabj%arr, 1, Wkaci%arr, 1)
    call Wmbej_r2abij(Wkaci, MP2t, abij, 2.0d0, iCore(i0), nCore - i0 + 1)
    write (6, '(A,F)') 'norm of Wkaci:', Wkaci%norm()

    call cpu_time(tt1)
    write (6, *) 'finish calculate intermidate MP2 amplitudes, Fli, Fad and Wkaci'
    write (6, *) 'cpu time for calculate intermidate:', tt1 - tt0

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 1st order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run1(iCore, i0)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 1st order calculate:', tt1 - tt0
    write (6, *) '2nd total correlation energy:', Ecorr
    write (6, *) '2nd total energy:', Escf + Ecorr

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 2nd order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run2(iCore, i0)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 2nd order calculate:', tt1 - tt0
    write (6, *) '3rd total correlation energy:', Ecorr
    write (6, *) '3rd total energy:', Escf + Ecorr

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 3rd order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run3(iCore, i0)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 3rd order calculate:', tt1 - tt0
    write (6, *) '4th total correlation energy:', Ecorr
    write (6, *) '4th total energy:', Escf + Ecorr

end subroutine run
