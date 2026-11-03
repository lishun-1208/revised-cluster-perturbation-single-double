! 2026.4.16
!> print excitation information for spin free t(a,i)
subroutine printEXinfo1_uhf(vecin, num)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num

    call p_printEXinfo1_uhf(vecin%arr, vecin%len, vecin%popp, vecin%popq, vecin%sym, num)

end subroutine printEXinfo1_uhf

! 2026.4.16
!> print excitation information for spin free t(ab,ij)
subroutine printEXinfo2_uhf(vecin, num)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num

    call p_printEXinfo2_uhf(vecin%arr, vecin%len, vecin%popp, vecin%popq, vecin%popr, vecin%pops, &
                            vecin%ioppq, vecin%ioprs, vecin%sym, num)

end subroutine printEXinfo2_uhf

! write by lishun, 2024.10.12
!-------------------------------
! modified at 2025.7.3
subroutine p_printEXinfo1_uhf(vecin, length, pop_a, pop_i, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num ! number of item need print

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item ! item is the number of items need print
    character(len=2), dimension(2) :: AB ! alpha,beta
    integer*8, dimension(2) :: idx ! i,a

    item = num
    AB(1) = ' a'
    AB(2) = ' b'

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 1 excitation information --------------'
    write (6, *) 'symi i ---> syma a'
    do i = 1, item
        idx = 1
        call getidx(maxIdx(i), info, pop_a, pop_i, 1, zorb, zorb, 1, sym)
        if (info(6) > popa(info(2))) then ! i
            idx(1) = 2
        end if
        if (info(5) > vrta(info(1))) then ! a
            idx(2) = 2
        end if

        write (6, '(I2,I4,A,A,I2,I4,A,F)') info(2), info(6), AB(idx(1)), &
            ' --->', info(1), info(5), AB(idx(2)), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printEXinfo1_uhf

! modified at 2025.7.3
subroutine p_printEXinfo2_uhf(vecin, length, pop_a, pop_b, pop_i, pop_j, iopab, iopij, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j
    integer*8, intent(in) :: iopab, iopij
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item
    character(len=2), dimension(2) :: AB ! alpha,beta
    integer*8, dimension(4) :: idx ! i,j,a,b

    item = num
    AB(1) = ' a'
    AB(2) = ' b'

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 2 excitation information --------------'
    write (6, *) 'symi i symj j ---> syma a symb b'
    do i = 1, item
        idx = 1
        call getidx(maxIdx(i), info, pop_a, pop_b, iopab, pop_i, pop_j, iopij, sym)
        if (info(7) > popa(info(3))) then ! i
            idx(1) = 2
        end if
        if (info(8) > popa(info(4))) then ! j
            idx(2) = 2
        end if
        if (info(5) > vrta(info(1))) then ! a
            idx(3) = 2
        end if
        if (info(6) > vrta(info(2))) then ! b
            idx(4) = 2
        end if

        write (6, '(I2,I4,2A,I2,I4,2A,I2,I4,2A,I2,I4,A,F)') info(3), info(7), AB(idx(1)), '|', &
            info(4), info(8), AB(idx(2)), ' --->', info(1), info(5), AB(idx(3)), '|', &
            info(2), info(6), AB(idx(4)), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printEXinfo2_uhf

! modified at 2025.7.3
!> print 1 and 2 excitation information for (t1,t2)
subroutine printEXinfo12_uhf(vecin, length, pop_a, pop_b, pop_i, pop_j, &
                             iopab, iopij, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j
    integer*8, intent(in) :: iopab, iopij
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item
    character(len=2), dimension(2) :: AB ! alpha,beta
    integer*8, dimension(4) :: idx ! i,j,a,b
    integer*8 :: nsingle
    integer*8, external :: irpdso, idsymsoc

    item = num
    AB(1) = ' a'
    AB(2) = ' b'

    !--- get maxIdx num elements
    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    nsingle = irpdso(sym, pop_a, pop_i, 1)
    write (6, *) '-------------- excitation information --------------'
    write (6, *) 'symi i symj j ---> syma a symb b'
    do i = 1, item
        if (maxIdx(i) <= nsingle) then
            !------ single
            idx = 1
            call getidx(maxIdx(i), info, pop_a, pop_i, 1, zorb, zorb, 1, sym)
            if (info(6) > popa(info(2))) then ! i
                idx(1) = 2
            end if
            if (info(5) > vrta(info(1))) then ! a
                idx(2) = 2
            end if
            write (6, '(I2,I4,A,A,I2,I4,A,F)') info(2), info(6), AB(idx(1)), &
                ' --->', info(1), info(5), AB(idx(2)), vecin(maxIdx(i))
        else
            !------ double
            idx = 1
            call getidx(maxIdx(i) - nsingle, info, pop_a, pop_b, iopab, pop_i, pop_j, iopij, sym)
            if (info(7) > popa(info(3))) then ! i
                idx(1) = 2
            end if
            if (info(8) > popa(info(4))) then ! j
                idx(2) = 2
            end if
            if (info(5) > vrta(info(1))) then ! a
                idx(3) = 2
            end if
            if (info(6) > vrta(info(2))) then ! b
                idx(4) = 2
            end if

            write (6, '(I2,I4,2A,I2,I4,2A,I2,I4,2A,I2,I4,A,F)') info(3), info(7), AB(idx(1)), '|', &
                info(4), info(8), AB(idx(2)), ' --->', info(1), info(5), AB(idx(3)), '|', &
                info(2), info(6), AB(idx(4)), vecin(maxIdx(i))
        end if
    end do
    write (6, *) '----------------------------------------------------'
    deallocate (maxIdx)

end subroutine printEXinfo12_uhf

! 2026.4.28
!> print excitation information for F(p,q)
subroutine printInfo1_uhf(vecin, num, type)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num
    character(len=2), intent(in) :: type !> 'O' or 'V'

    call p_printInfo1_uhf(vecin%arr, vecin%len, vecin%popp, vecin%popq, &
                          type(1:1), type(2:2), vecin%sym, num)

end subroutine printInfo1_uhf

! 2026.4.28
!> print excitation information for G(pq,rs)
subroutine printInfo2_uhf(vecin, num, type)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num
    character(len=4), intent(in) :: type !> 'O' or 'V'

    call p_printInfo2_uhf(vecin%arr, vecin%len, vecin%popp, vecin%popq, vecin%popr, vecin%pops, &
                          type(1:1), type(2:2), type(3:3), type(4:4), vecin%ioppq, vecin%ioprs, &
                          vecin%sym, num)

end subroutine printInfo2_uhf

! write at 2025.7.3
!> print excitation information for F(p,q)
subroutine p_printInfo1_uhf(vecin, length, popp, popq, typep, typeq, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: popp, popq
    character(len=1), intent(in) :: typep, typeq ! type of virtual of occupied, 'V','O'
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num ! number of item need print

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item ! item is the number of items need print
    character(len=2), dimension(2) :: AB ! alpha,beta
    integer*8, dimension(2) :: idx ! i,a

    item = num
    AB(1) = ' a'
    AB(2) = ' b'

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 1 excitation information --------------'
    write (6, *) 'symq q ---> symp p'
    do i = 1, item
        idx = 1
        call getidx(maxIdx(i), info, popp, popq, 1, zorb, zorb, 1, sym)
        ! p
        if (typep == 'O') then
            if (info(5) > popa(info(1))) then
                idx(2) = 2
            end if
        else if (typep == 'V') then
            if (info(5) > vrta(info(1))) then
                idx(2) = 2
            end if
        end if
        ! q
        if (typeq == 'O') then
            if (info(6) > popa(info(2))) then
                idx(1) = 2
            end if
        else if (typeq == 'V') then
            if (info(6) > vrta(info(2))) then
                idx(1) = 2
            end if
        end if

        write (6, '(I2,I4,A,A,I2,I4,A,F)') info(2), info(6), AB(idx(1)), &
            ' --->', info(1), info(5), AB(idx(2)), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printInfo1_uhf

! write at 2025.7.2
!> print excitation information for G(pq,rs)
subroutine p_printInfo2_uhf(vecin, length, popp, popq, popr, pops, &
                            typep, typeq, typer, types, ioppq, ioprs, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    character(len=1), intent(in) :: typep, typeq, typer, types
    integer*8, intent(in) :: ioppq, ioprs
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item
    character(len=2), dimension(2) :: AB ! alpha,beta
    integer*8, dimension(4) :: idx ! i,j,a,b

    item = num
    AB(1) = ' a'
    AB(2) = ' b'

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 2 excitation information --------------'
    write (6, *) 'symr r syms s ---> symp p symq q'
    do i = 1, item
        idx = 1
        call getidx(maxIdx(i), info, popp, popq, ioppq, popr, pops, ioprs, sym)
        ! p
        if (typep == 'O') then
            if (info(5) > popa(info(1))) then
                idx(3) = 2
            end if
        else if (typep == 'V') then
            if (info(5) > vrta(info(1))) then
                idx(3) = 2
            end if
        end if
        ! q
        if (typeq == 'O') then
            if (info(6) > popa(info(2))) then
                idx(4) = 2
            end if
        else if (typeq == 'V') then
            if (info(6) > vrta(info(2))) then
                idx(4) = 2
            end if
        end if
        ! r
        if (typer == 'O') then
            if (info(7) > popa(info(3))) then
                idx(1) = 2
            end if
        else if (typer == 'V') then
            if (info(7) > vrta(info(3))) then
                idx(1) = 2
            end if
        end if
        ! s
        if (types == 'O') then
            if (info(8) > popa(info(4))) then
                idx(2) = 2
            end if
        else if (types == 'V') then
            if (info(8) > vrta(info(4))) then
                idx(2) = 2
            end if
        end if

        write (6, '(I2,I4,2A,I2,I4,2A,I2,I4,2A,I2,I4,A,F)') info(3), info(7), AB(idx(1)), '|', &
            info(4), info(8), AB(idx(2)), ' --->', info(1), info(5), AB(idx(3)), '|', &
            info(2), info(6), AB(idx(4)), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printInfo2_uhf

! 2026.7.4
subroutine printEXinfo1(vecin, num)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num

    call p_printEXinfo1(vecin%arr, vecin%len, vecin%popp, vecin%popq, vecin%sym, num)

end subroutine printEXinfo1

! 2026.7.4
subroutine p_printEXinfo1(vecin, length, pop_a, pop_i, sym, num)
    use mod_orbit
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num ! number of item need print

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item ! item is the number of items need print

    item = num

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 1 excitation information --------------'
    write (6, *) 'symi i ---> syma a'
    do i = 1, item
        call getidx(maxIdx(i), info, pop_a, pop_i, 1, zorb, zorb, 1, sym)

        write (6, '(I2,I4,A,I2,I4,F)') info(2), info(6), &
            ' --->', info(1), info(5), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printEXinfo1

! 2026.7.4
subroutine printEXinfo2(vecin, num)
    use mod_operator
    implicit none

    type(opt), intent(in) :: vecin
    integer*8, intent(in) :: num

    call p_printEXinfo2(vecin%arr, vecin%len, vecin%popp, vecin%popq, vecin%popr, vecin%pops, &
                        vecin%ioppq, vecin%ioprs, vecin%sym, num)

end subroutine printEXinfo2

! 2026.7.4
subroutine p_printEXinfo2(vecin, length, pop_a, pop_b, pop_i, pop_j, iopab, iopij, sym, num)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(length), intent(in) :: vecin
    integer*8, intent(in) :: length ! length of vecin
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j
    integer*8, intent(in) :: iopab, iopij
    integer*8, intent(in) :: sym
    integer*8, intent(in) :: num

    integer*8, dimension(:), allocatable :: maxIdx
    integer*8, dimension(8) :: info ! irrepp,irrepq,irrepr,irreps,p,q,r,s
    integer*8 :: i, item

    item = num

    if (item > length) then
        item = length
    end if
    allocate (maxIdx(item))
    call namax(vecin, length, item, maxIdx)

    write (6, *) '-------------- 2 excitation information --------------'
    write (6, *) 'symi i symj j ---> syma a symb b'
    do i = 1, item
        call getidx(maxIdx(i), info, pop_a, pop_b, iopab, pop_i, pop_j, iopij, sym)

        write (6, '(I2,I4,A,I2,I4,A,I2,I4,A,I2,I4,F)') info(3), info(7), '|', &
            info(4), info(8), ' --->', info(1), info(5), '|', &
            info(2), info(6), vecin(maxIdx(i))
    end do
    write (6, *) '------------------------------------------------------'
    deallocate (maxIdx)
end subroutine p_printEXinfo2
