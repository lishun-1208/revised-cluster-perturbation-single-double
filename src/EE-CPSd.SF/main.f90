! write at 2025.1.13
! excited state CPS(D) program without spin
! need canonical Fock matrix
program EE_CPSd_SF
    use mod_variable
    implicit none

    real*8, dimension(:), allocatable, target :: iCore ! main work space
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

    call run(iCore)

    call cpu_time(tt1)
    write (6, '(A,F)') "total cpu time for program:", tt1 - tt0

end program EE_CPSd_SF
