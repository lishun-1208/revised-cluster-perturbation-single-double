! write by lishun, 2024.11.29
! read in data from hard disk
! use i/o unit 6219
! sort as (val,vec)
!-----------------------------------
! modified at 2024.1.13
! suggest to use IO_stream
subroutine IO_direct(valinout, vecinout, recordLength, record, filerw, type, state)
    implicit none

    real*8, intent(inout) :: valinout
    real*8, dimension(*), intent(inout) :: vecinout
    integer*8, intent(in) :: recordLength ! length of every record
    integer*8, intent(in) :: record ! record number
    character(len=*), intent(in) :: filerw ! file name
    ! first: 1 for read, 2 for write
    ! second: 1 for only have vecinout, others for both vec and val
    !         2 for valinout and vecinout
    !         3 for only i/o valinout
    !         4 for only i/o vecinout
    character(len=2), intent(in) :: type
    character(len=*), intent(in) :: state ! 'old' or 'replace', only for write

    integer*8 :: i
    real*8 :: tmp

    if ((state /= 'old') .and. (state /= 'replace')) then
        write (6, *) 'ERROR: status in IO_direct error', state
        stop
    end if

    if (type == '12') then
        open (unit=6219, file=filerw, form='unformatted', action='read', &
              status='old', access='direct', recl=8*recordLength)
        read (6219, rec=record) valinout, (vecinout(i), i=1, recordLength - 1)
        close (6219)
    else if (type == '22') then
        open (unit=6219, file=filerw, form='unformatted', action='write', &
              status=state, access='direct', recl=8*recordLength)
        write (6219, rec=record) valinout, (vecinout(i), i=1, recordLength - 1)
        close (6219)
    else if (type == '11') then
        open (unit=6219, file=filerw, form='unformatted', action='read', &
              status='old', access='direct', recl=8*recordLength)
        read (6219, rec=record) (vecinout(i), i=1, recordLength)
        close (6219)
    else if (type == '21') then
        open (unit=6219, file=filerw, form='unformatted', action='write', &
              status=state, access='direct', recl=8*recordLength)
        write (6219, rec=record) (vecinout(i), i=1, recordLength)
        close (6219)
    else if (type == '13') then
        open (unit=6219, file=filerw, form='unformatted', action='read', &
              status='old', access='direct', recl=8*recordLength)
        read (6219, rec=record) valinout
        close (6219)
    else if (type == '23') then
        open (unit=6219, file=filerw, form='unformatted', action='write', &
              status=state, access='direct', recl=8*recordLength)
        write (6219, rec=record) valinout
        close (6219)
    else if (type == '14') then
        open (unit=6219, file=filerw, form='unformatted', action='read', &
              status=state, access='direct', recl=8*recordLength)
        read (6219, rec=record) tmp, (vecinout(i), i=1, recordLength - 1)
        close (6219)
    else if (type == '24') then
        ! read val temporary
        open (unit=6219, file=filerw, form='unformatted', action='read', &
              status=state, access='direct', recl=8*recordLength)
        read (6219, rec=record) tmp
        close (6219)
        ! write
        open (unit=6219, file=filerw, form='unformatted', action='write', &
              status=state, access='direct', recl=8*recordLength)
        write (6219, rec=record) tmp, (vecinout(i), i=1, recordLength - 1)
        close (6219)
    else
        write (6, *) 'ERROR: type in IO_direct error', type
        stop
    end if

end subroutine IO_direct

! write at 2025.3.11
! modified at 2025.11.7
subroutine IO_stream(val, vec, idx, nElements, filerw, type, state)
    use mod_iop
    implicit none

    real*8, intent(inout) :: val
    real*8, dimension(nelements), intent(inout) :: vec
    integer*8, intent(in) :: idx !> idx of first element in file, with 8 Byte(a real*8) as a unit
    integer*8, intent(in) :: nElements !> number of elements need be operate in vec
    character(len=*), intent(in) :: filerw ! file name
    ! first: 1 for read, 2 for write
    ! second: 1 for i/o both val and vec
    !         2 for only i/o val
    !         3 for only i/o vec
    character(len=2), intent(in) :: type
    character(len=*), intent(in) :: state !> 'old' or 'replace'

    character(len=:), allocatable :: file
    integer*8 :: i
    integer*8 :: position ! position of first element in file, with 1 Byte as a unit

    if ((state /= 'old') .and. (state /= 'replace')) then
        write (6, *) 'ERROR: status in IO_stream error', state
        stop
    end if

    if (nElements >= 2**31 - 1) then
        write (6, *) 'ERROR: nElements exceed max value of int32 in IO_stream', nElements
        stop
    end if

    file = trim(tmpdir)//filerw
    position = 8*(idx - 1) + 1
    if (type == '11') then
        open (unit=6219, file=file, form='unformatted', action='read', &
              status=state, access='stream')
        read (6219, pos=position) val, (vec(i), i=1, nElements)
        close (6219)
    else if (type == '12') then
        open (unit=6219, file=file, form='unformatted', action='read', &
              status=state, access='stream')
        read (6219, pos=position) val
        close (6219)
    else if (type == '13') then
        open (unit=6219, file=file, form='unformatted', action='read', &
              status=state, access='stream')
        read (6219, pos=position) (vec(i), i=1, nElements)
        close (6219)
    else if (type == '21') then
        open (unit=6219, file=file, form='unformatted', action='write', &
              status=state, access='stream')
        write (6219, pos=position) val, (vec(i), i=1, nElements)
        close (6219)
    else if (type == '22') then
        open (unit=6219, file=file, form='unformatted', action='write', &
              status=state, access='stream')
        write (6219, pos=position) val
        close (6219)
    else if (type == '23') then
        open (unit=6219, file=file, form='unformatted', action='write', &
              status=state, access='stream')
        write (6219, pos=position) (vec(i), i=1, nElements)
        close (6219)
    else
        write (6, *) 'ERROR: type in IO_stream error', type
        stop
    end if

end subroutine IO_stream
