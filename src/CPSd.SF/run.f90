! write by lishun, 2024.12.5
! 0 and 1 order correlation energy is 0
subroutine run(iCore, maxCore)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    ! iteration contral
    integer*8 :: max_cyc
    real*8 :: maxError
    integer*8 :: maxDim

    real*8 :: tt0, tt1
    integer*8 :: nSize
    integer*8, external :: idsymsoc, irpdso
    real*8, external :: dnrm2

    max_cyc = diis_cyc
    write (6, *) 'max cycle for cc iteration:', max_cyc
    maxError = diis_max_err
    write (6, *) 'cc mconvergence criteria:', maxError
    maxDim = diis_dim
    write (6, *) 'max DIIS dimension:', maxDim

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 1st order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run1(iCore, maxCore)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 1st order calculate:', tt1 - tt0
    write (6, *) 'total correlation energy:', Ecorr
    write (6, *) 'total energy:', Escf + Ecorr

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 2nd order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run2(iCore, maxCore, max_cyc, maxError, maxDim)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 2nd order calculate:', tt1 - tt0
    write (6, *) 'total correlation energy:', Ecorr
    write (6, *) 'total energy:', Escf + Ecorr

    write (6, *) '-----------------------------------------------------'
    write (6, *) '             start 3rd order calculation             '
    write (6, *) '-----------------------------------------------------'
    call cpu_time(tt0)
    call run3(iCore, maxCore, max_cyc, maxError, maxDim)
    call cpu_time(tt1)
    write (6, *) 'cpu time for 3rd order calculate:', tt1 - tt0
    write (6, *) 'total correlation energy:', Ecorr
    write (6, *) 'total energy:', Escf + Ecorr

end subroutine run
