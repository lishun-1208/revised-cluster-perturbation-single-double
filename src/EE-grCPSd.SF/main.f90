! write by lishun
! modified 2024.11.4
! spin free excited state iCPS(D) program
! need canonical Fock matrix
program EE_grCPSd_SF
    use mod_orbit
    use mod_symm
    use mod_variable
    implicit none

    real*8, dimension(:), allocatable, target :: iCore ! main work space
    integer*8 :: status ! allocate stat
    ! cpu time
    real*8 :: tt0, tt1
    real*8 :: t0, t1

    call cpu_time(tt0)
    call cpu_time(t0)
    call prepare

    write (6, '(A)') '-----------------------------------------'
    write (6, '(A)') '           start allocate iCore          '
    write (6, '(A)') '-----------------------------------------'
    allocate (iCore(nCore), stat=status)
    if (status /= 0) then
        write (6, '(A)') "ERROR: allocate iCore error, memory not enough"
        stop
    end if
    call allocate_mem(iCore)

    write (6, '(A)') '-----------------------------------------'
    write (6, '(A)') '               start loadSCF             '
    write (6, '(A)') '-----------------------------------------'
    call loadSCF(iCore)
    call cpu_time(t1)
    write (6, '(A,F)') 'cpu time for prepare:', t1 - t0

    call run(iCore)

    call cpu_time(tt1)
    write (6, '(A,F)') "total cpu time for program:", tt1 - tt0

end program EE_grCPSd_SF
