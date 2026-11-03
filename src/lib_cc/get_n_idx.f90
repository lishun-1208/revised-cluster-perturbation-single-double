! WARN: no more use
! write by LiShun, 2024.9.27
! get the position of lowest elements in vecin
!----------------------------------
! rewrite at 2025.6.5
!> this function is deprecated, please use nmin
subroutine getLowest(vecin, nSize, nlowest, lowest)
    implicit none

    real*8, dimension(nSize), intent(inout) :: vecin ! input 1D array
    integer*8, intent(in) :: nSize ! size of vecin
    integer*8, intent(in) :: nlowest ! the number of lowest element
    integer*8, dimension(nlowest), intent(out) :: lowest ! position of lowest elements, length is nlowest

    integer*8 :: i, j, idx
    integer*8 :: start, end, mid

    if (nSize < nlowest) then
        write (6, '(A,2I)') 'ERROR: nSize < nlowest in getLowest', nSize, nlowest
        stop
    end if

    lowest = (/(i, i=1, nlowest)/)

    ! sort, s->l
    do i = 1, nlowest - 1
        do j = i + 1, nlowest
            if (vecin(lowest(i)) > vecin(lowest(j))) then
                idx = lowest(i)
                lowest(i) = lowest(j)
                lowest(j) = idx
            end if
        end do
    end do

    do i = nlowest + 1, nSize
        ! popup largest element
        if (vecin(i) < vecin(lowest(nlowest))) then
            ! find, I find the first element bigger than vecin(i)
            start = 1
            end = nlowest
            do
                ! for 2 element (s,e),m=s,vec(e)>val
                mid = start + (end - start)/2
                if (vecin(lowest(mid)) > vecin(i)) then
                    end = mid ! vec(s)>val, s
                else
                    start = mid + 1 ! vec(s)<val, e
                end if
                if (start == end) then
                    exit
                end if
            end do
            idx = end
            ! insert
            do j = 0, nlowest - idx - 1
                lowest(nlowest - j) = lowest(nlowest - j - 1)
            end do
            lowest(idx) = i
        end if
    end do

end subroutine getLowest

! write at 2025.6.5, lishun
!> get some max element in vecin
subroutine nmax(vecin, nSize, nhighest, highest)
    implicit none

    real*8, dimension(nSize), intent(inout) :: vecin !> input 1D array
    integer*8, intent(in) :: nSize !> size of vecin
    integer*8, intent(in) :: nhighest !> the number of max element
    integer*8, dimension(nhighest), intent(out) :: highest !> position of max elements, length is nhighest

    integer*8 :: i, j, idx
    integer*8 :: start, end, mid

    if (nSize < nhighest) then
        write (6, '(A,2I)') 'ERROR: nSize < nhighest in nmax', nSize, nhighest
        stop
    end if

    highest = (/(i, i=1, nhighest)/)

    ! sort, l->s
    do i = 1, nhighest - 1
        do j = i + 1, nhighest
            if (vecin(highest(i)) < vecin(highest(j))) then
                idx = highest(i)
                highest(i) = highest(j)
                highest(j) = idx
            end if
        end do
    end do

    do i = nhighest + 1, nSize
        ! popup smallest element
        if (vecin(i) > vecin(highest(nhighest))) then
            ! find, I find the first element smaller than vecin(i)
            start = 1
            end = nhighest
            do
                ! for 2 element (s,e),m=s,vec(e)<val
                mid = start + (end - start)/2
                if (vecin(highest(mid)) < vecin(i)) then
                    end = mid ! vec(s)<val, s
                else
                    start = mid + 1 ! vec(s)>val, e
                end if
                if (start == end) then
                    exit
                end if
            end do
            idx = end
            ! insert
            do j = 0, nhighest - idx - 1
                highest(nhighest - j) = highest(nhighest - j - 1)
            end do
            highest(idx) = i
        end if
    end do

end subroutine nmax

! write at 2025.6.5, lishun
!> get some absolute max element in vecin
subroutine namax(vecin, nSize, nhighest, highest)
    implicit none

    real*8, dimension(nSize), intent(inout) :: vecin !> input 1D array
    integer*8, intent(in) :: nSize !> size of vecin
    integer*8, intent(in) :: nhighest !> the number of absolute max element
    integer*8, dimension(nhighest), intent(out) :: highest !> position of absolute max elements, length is nhighest

    integer*8 :: i, j, idx
    integer*8 :: start, end, mid

    if (nSize < nhighest) then
        write (6, '(A,2I)') 'ERROR: nSize < nhighest in namax', nSize, nhighest
        stop
    end if

    highest = (/(i, i=1, nhighest)/)

    ! sort, l->s
    do i = 1, nhighest - 1
        do j = i + 1, nhighest
            if (abs(vecin(highest(i))) < abs(vecin(highest(j)))) then
                idx = highest(i)
                highest(i) = highest(j)
                highest(j) = idx
            end if
        end do
    end do

    do i = nhighest + 1, nSize
        ! popup smallest element
        if (abs(vecin(i)) > abs(vecin(highest(nhighest)))) then
            ! find, I find the first element smaller than vecin(i)
            start = 1
            end = nhighest
            do
                ! for 2 element (s,e),m=s,vec(e)<val
                mid = start + (end - start)/2
                if (abs(vecin(highest(mid))) < abs(vecin(i))) then
                    end = mid ! vec(s)<val, s
                else
                    start = mid + 1 ! vec(s)>val, e
                end if
                if (start == end) then
                    exit
                end if
            end do
            idx = end
            ! insert
            do j = 0, nhighest - idx - 1
                highest(nhighest - j) = highest(nhighest - j - 1)
            end do
            highest(idx) = i
        end if
    end do

end subroutine namax

! 2026.4.16
!> get some min element in vecin
subroutine nmin(vecin, nSize, nlowest, lowest)
    implicit none

    real*8, dimension(nSize), intent(inout) :: vecin !> input 1D array
    integer*8, intent(in) :: nSize !> size of vecin
    integer*8, intent(in) :: nlowest !> the number of lowest element
    integer*8, dimension(nlowest), intent(out) :: lowest !> position of lowest elements, length is nlowest

    call getLowest(vecin, nSize, nlowest, lowest)

end subroutine nmin

! 2026.5.7, lish
!> get some absolute min element in vecin
subroutine namin(vecin, nSize, nlowest, lowest)
    implicit none

    real*8, dimension(nSize), intent(inout) :: vecin !> input 1D array
    integer*8, intent(in) :: nSize !> size of vecin
    integer*8, intent(in) :: nlowest !> the number of absolute min element
    integer*8, dimension(nlowest), intent(out) :: lowest !> position of absolute min elements, length is nlowest

    integer*8 :: i, j, idx
    integer*8 :: start, end, mid

    if (nSize < nlowest) then
        write (6, '(A,2I)') 'ERROR: nSize < nlowest in namin', nSize, nlowest
        stop
    end if

    lowest = (/(i, i=1, nlowest)/)

    ! sort, s->l
    do i = 1, nlowest - 1
        do j = i + 1, nlowest
            if (abs(vecin(lowest(i))) > abs(vecin(lowest(j)))) then
                idx = lowest(i)
                lowest(i) = lowest(j)
                lowest(j) = idx
            end if
        end do
    end do

    do i = nlowest + 1, nSize
        ! popup lagest element
        if (abs(vecin(i)) < abs(vecin(lowest(nlowest)))) then
            ! find, I find the first element smaller than vecin(i)
            start = 1
            end = nlowest
            do
                ! for 2 element (s,e),m=s,vec(e)>val
                mid = start + (end - start)/2
                if (abs(vecin(lowest(mid))) > abs(vecin(i))) then
                    end = mid ! vec(s)>val, s
                else
                    start = mid + 1 ! vec(s)<val, e
                end if
                if (start == end) then
                    exit
                end if
            end do
            idx = end
            ! insert
            do j = 0, nlowest - idx - 1
                lowest(nlowest - j) = lowest(nlowest - j - 1)
            end do
            lowest(idx) = i
        end if
    end do

end subroutine namin

! 2026.5.13
!> get the number of elements which abs value lower than given threshold
subroutine namin_threshold(vecin, nSize, threshold, nlowest)
    implicit none

    real*8, dimension(nSize), intent(in) :: vecin !> input 1D array
    integer*8, intent(in) :: nSize !> size of vecin
    real*8, intent(in) :: threshold !> threshold
    integer*8, intent(out) :: nlowest

    integer*8 :: i, j

    nlowest = 0
    do i = 1, nSize
        if (abs(vecin(i)) < threshold) then
            nlowest = nlowest + 1
        end if
    end do

end subroutine namin_threshold