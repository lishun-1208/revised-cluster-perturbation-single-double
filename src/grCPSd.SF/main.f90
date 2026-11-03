! write by LiShun, 2024.10.14
! spin free grCPS(D) program
! need canonical Fock matrix
!
program grCPSd_SF
    use mod_variable
    implicit none

    real*8, dimension(:), allocatable :: iCore ! main work space
    integer*8 :: status ! allocate stat
    ! cpu time
    real*8 :: tt0, tt1, tt2, tt3

    call cpu_time(tt0)
    call cpu_time(tt2)
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
    call cpu_time(tt3)
    write (6, '(A,F)') 'cpu time for prepare:', tt3 - tt2

    call run(iCore, iFree)

    call cpu_time(tt1)
    write (6, '(A,F)') "total cpu time for program:", tt1 - tt0

end program grCPSd_SF
