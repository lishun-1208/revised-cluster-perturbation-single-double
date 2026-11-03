! require libvv2v, libccsoint at /home/wangf1/CCSO_PROG
module mod_contract
    use mod_operator
    implicit none
    public

contains
    ! 2026.7.3
    !> Calculate v3 += fact*P1P2{v1*v2}
    subroutine vv2v_444(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore, anti)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=4), intent(in) :: type1, type2, type3 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore
        integer*8, intent(in) :: anti !> 0/1/2/3 for antisymetrize no/pair1/pair2/both in v3

        character(len=2), dimension(6) :: pair
        logical, dimension(6) :: flat !> whether idx need be flat
        integer*8, dimension(6) :: iop
        integer*8, dimension(6) :: iexpand !> whether idx need be expand, 1 - need, 0 - not
        integer*8 :: i
        target :: iCore
        real*8, dimension(:), pointer :: ptr1, ptr2, ptr3
        integer*8 :: i0, i1, i2, i3, nSize
        integer*8, external :: idsymsoc

        if ((.not. v1%incore) .or. (.not. v2%incore) .or. (.not. v3%incore)) then
            write (6, *) 'ERROR: vv2v_444 need all tensor be incore.', v1%incore, v2%incore, v3%incore
            stop
        end if
        if ((anti < 0) .or. (anti > 3)) then
            write (6, '(A,I2)') 'ERROR: vv2v_444 error in antisymmetrize.', anti
            stop
        end if

        i0 = 1
        ! check which pair of idx need be flat
        pair(1) = type1(1:2)
        pair(2) = type1(3:4)
        pair(3) = type2(1:2)
        pair(4) = type2(3:4)
        pair(5) = type3(1:2)
        pair(6) = type3(3:4)
        flat = .false.
        do i = 1, 6
            flat(i) = check_flat(pair(i), [pair(1:i - 1), pair(i + 1:6)], 5)
        end do
        ! get iop after expand
        iop(1) = v1%ioppq
        iop(2) = v1%ioprs
        iop(3) = v2%ioppq
        iop(4) = v2%ioprs
        iop(5) = v3%ioppq
        iop(6) = v3%ioprs
        iexpand = 0
        do i = 1, 6
            if (flat(i) .and. (iop(i) /= 1)) then
                iexpand(i) = 1
                iop(i) = 1
            end if
        end do

        !------ expand
        ! v3, i0
        ! if v3 need be expand, it means that also need antisymmetrize
        ! and need use a new space to store contract results
        if (anti /= 0) then
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, iop(5), v3%popr, v3%pops, iop(6))
            i1 = i0 + nSize
            call dscal(nSize, 0.0d0, iCore(i0), 1)
            ptr3(1:1) => iCore(i0:i0)
        else
            i1 = i0
            ptr3 => v3%arr
        end if
        ! v1, i1
        if ((iexpand(1) == 1) .or. (iexpand(2) == 1)) then
            nSize = idsymsoc(v1%sym, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2))
            i2 = i1 + nSize
            call iexpso(v1%arr, iCore(i1), v1%popp, v1%popq, v1%ioppq, iexpand(1), &
                        v1%popr, v1%pops, v1%ioprs, iexpand(2), v1%sym)
            ptr1(1:1) => iCore(i1:i1)
        else
            i2 = i1
            ptr1 => v1%arr
        end if
        ! v2, i2
        if ((iexpand(3) == 1) .or. (iexpand(4) == 1)) then
            nSize = idsymsoc(v2%sym, v2%popp, v2%popq, iop(3), v2%popr, v2%pops, iop(4))
            i3 = i2 + nSize
            call iexpso(v2%arr, iCore(i2), v2%popp, v2%popq, v2%ioppq, iexpand(3), &
                        v2%popr, v2%pops, v2%ioprs, iexpand(4), v2%sym)
            ptr2(1:1) => iCore(i2:i2)
        else
            i3 = i2
            ptr2 => v2%arr
        end if

        !------ calculate
        call vv2v_dpd_444(ptr1, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2), v1%sym, type1, &
                          ptr2, v2%popp, v2%popq, iop(3), v2%popr, v2%pops, iop(4), v2%sym, type2, &
                          ptr3, v3%popp, v3%popq, iop(5), v3%popr, v3%pops, iop(6), v3%sym, type3, &
                          fact, 1.0d0)

        !------ make permutation and compress
        !--- antisymetrize 1st pair
        if ((anti == 1) .or. (anti == 3)) then
            if ((iop(5) /= 1) .or. (.not. pop_eq(v3%popp, v3%popq))) then
                write (6, '(A,I2)') 'ERROR: vv2v_444 anti 1st pair error.', iop(5)
                write (6, '(8I5)') v3%popp
                write (6, '(8I5)') v3%popq
                stop
            end if
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, iop(5), v3%popr, v3%pops, iop(6))
            call isymtrso(ptr3, iCore(i1), v3%popp, v3%popq, iop(5), 1, &
                          v3%popr, v3%pops, iop(6), 0, v3%sym)
            call daxpy(nSize, -1.0d0, iCore(i1), 1, ptr3, 1)
        end if
        ! compress 1st pair
        if (iexpand(5) == 1) then
            i2 = i1 + nSize
            call icompso(ptr3, iCore(i1), v3%popp, v3%popq, iop(5), v3%ioppq, &
                         v3%popr, v3%pops, iop(6), 1, v3%sym)
            ptr3(1:1) => iCore(i1:i1)
        else
            i2 = i1
        end if
        !--- antisymetrize 2nd pair
        if ((anti == 2) .or. (anti == 3)) then
            if ((iop(6) /= 1) .or. (.not. pop_eq(v3%popr, v3%pops))) then
                write (6, '(A,I2)') 'ERROR: vv2v_444 anti 2nd pair error.', iop(6)
                write (6, '(8I5)') v3%popr
                write (6, '(8I5)') v3%pops
                stop
            end if
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, v3%ioppq, v3%popr, v3%pops, iop(6))
            call isymtrso(ptr3, iCore(i2), v3%popp, v3%popq, v3%ioppq, 0, &
                          v3%popr, v3%pops, iop(6), 1, v3%sym)
            call daxpy(nSize, -1.0d0, iCore(i2), 1, ptr3, 1)
        end if
        ! compress 2nd pair
        if (iexpand(6) == 1) then
            i3 = i2 + nSize
            call icompso(ptr3, iCore(i2), v3%popp, v3%popq, v3%ioppq, 1, &
                         v3%popr, v3%pops, iop(6), v3%ioprs, v3%sym)
            ptr3(1:1) => iCore(i2:i2)
        else
            i3 = i2
        end if

        !--- add
        if (anti /= 0) then
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, v3%ioppq, v3%popr, v3%pops, v3%ioprs)
            call daxpy(nSize, 1.0d0, ptr3, 1, v3%arr, 1)
        end if

    end subroutine vv2v_444

    ! 2026.7.3
    !> Calculate v3 += fact*v1*v2
    subroutine vv2v_222(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=2), intent(in) :: type1, type2, type3 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore

        call vv2v_444(v1, v2, v3, fact, type1//'=+', type2//'=-', type3//'+-', &
                      iCore, maxCore, 0)

    end subroutine vv2v_222

    ! 2026.7.3
    !> Calculate v3 += fact*v1*v2
    subroutine vv2v_422(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=4), intent(in) :: type1 ! contract type
        character(len=2), intent(in) :: type2, type3 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore

        call vv2v_444(v1, v2, v3, fact, type1, type2//'+-', type3//'+-', &
                      iCore, maxCore, 0)

    end subroutine vv2v_422

    ! 2026.7.4
    !> Calculate v3 += fact*P1P2{v1*v2}
    subroutine vv2v_224(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore, anti)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=2), intent(in) :: type1, type2 ! contract type
        character(len=4), intent(in) :: type3 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore
        integer*8, intent(in) :: anti !> 0/1/2/3 for antisymetrize no/pair1/pair2/both in v3

        call vv2v_444(v1, v2, v3, fact, type1//'+-', type2//'+-', type3, &
                      iCore, maxCore, anti)

    end subroutine vv2v_224

    ! 2026.7.4
    !> Calculate v3 = v3 + fact*P{v1*v2}
    subroutine vv2v_424(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore, anti)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=4), intent(in) :: type1, type3 ! contract type
        character(len=2), intent(in) :: type2 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore
        integer*8, intent(in) :: anti !> 0/1/2 for antisymetrize no/pair1/pair2 in v3

        character(len=2), dimension(4) :: pair
        logical, dimension(4) :: flat !> whether idx need be flat
        integer*8, dimension(4) :: iop
        integer*8, dimension(4) :: iexpand !> whether idx need be expand, 1 - need, 0 - not
        integer*8 :: i
        target :: iCore
        real*8, dimension(:), pointer :: ptr1, ptr2, ptr3
        integer*8 :: i0, i1, i2, i3, nSize
        integer*8, external :: idsymsoc

        if ((.not. v1%incore) .or. (.not. v2%incore) .or. (.not. v3%incore)) then
            write (6, *) 'ERROR: vv2v_424 need all tensor be incore.', v1%incore, v2%incore, v3%incore
            stop
        end if
        if ((anti < 0) .or. (anti > 2)) then
            write (6, '(A,I2)') 'ERROR: vv2v_424 error in antisymmetrize.', anti
            stop
        end if

        i0 = 1
        ! check which pair of idx need be flat
        pair(1) = type1(1:2)
        pair(2) = type1(3:4)
        pair(3) = type3(1:2)
        pair(4) = type3(3:4)
        flat = .false.
        do i = 1, 4
            flat(i) = check_flat(pair(i), [pair(1:i - 1), pair(i + 1:4)], 3)
        end do
        ! get iop after expand
        iop(1) = v1%ioppq
        iop(2) = v1%ioprs
        iop(3) = v3%ioppq
        iop(4) = v3%ioprs
        iexpand = 0
        do i = 1, 4
            if (flat(i) .and. (iop(i) /= 1)) then
                iexpand(i) = 1
                iop(i) = 1
            end if
        end do

        !------ expand
        ! v3, i0
        ! if v3 need be expand, it means that also need antisymmetrize
        ! and need use a new space to store contract results
        if (anti /= 0) then
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, iop(3), v3%popr, v3%pops, iop(4))
            i1 = i0 + nSize
            call dscal(nSize, 0.0d0, iCore(i0), 1)
            ptr3(1:1) => iCore(i0:i0)
        else
            i1 = i0
            ptr3 => v3%arr
        end if
        ! v1, i1
        if ((iexpand(1) == 1) .or. (iexpand(2) == 1)) then
            nSize = idsymsoc(v1%sym, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2))
            i2 = i1 + nSize
            call iexpso(v1%arr, iCore(i1), v1%popp, v1%popq, v1%ioppq, iexpand(1), &
                        v1%popr, v1%pops, v1%ioprs, iexpand(2), v1%sym)
            ptr1(1:1) => iCore(i1:i1)
        else
            i2 = i1
            ptr1 => v1%arr
        end if
        ! v2
        ptr2 => v2%arr

        !------ calculate
        call vv2v_dpd_424(ptr1, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2), v1%sym, type1, &
                          ptr2, v2%popp, v2%popq, v2%ioppq, v2%sym, type2, &
                          ptr3, v3%popp, v3%popq, iop(3), v3%popr, v3%pops, iop(4), v3%sym, type3, &
                          fact, 1.0d0)

        !------ make permutation and compress
        !--- antisymetrize 1st pair
        if (anti == 1) then
            if ((iop(3) /= 1) .or. (.not. pop_eq(v3%popp, v3%popq))) then
                write (6, '(A,I2)') 'ERROR: vv2v_424 anti 1st pair error.', iop(3)
                write (6, '(8I5)') v3%popp
                write (6, '(8I5)') v3%popq
                stop
            end if
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, iop(3), v3%popr, v3%pops, iop(4))
            call isymtrso(ptr3, iCore(i1), v3%popp, v3%popq, iop(3), 1, &
                          v3%popr, v3%pops, iop(4), 0, v3%sym)
            call daxpy(nSize, -1.0d0, iCore(i1), 1, ptr3, 1)
            ! compress 1st pair
            if (iexpand(3) == 1) then
                i2 = i1 + nSize
                call icompso(ptr3, iCore(i1), v3%popp, v3%popq, iop(3), v3%ioppq, &
                             v3%popr, v3%pops, iop(4), 1, v3%sym)
                ptr3(1:1) => iCore(i1:i1)
            else
                i2 = i1
            end if
        else if (iexpand(3) == 1) then
            write (6, '(A)') 'ERROR: vv2v_424 error for the 1st pair not antisymmertize but be expand.'
            stop
        end if
        !--- antisymetrize 2nd pair
        if (anti == 2) then
            if ((iop(4) /= 1) .or. (.not. pop_eq(v3%popr, v3%pops))) then
                write (6, '(A,I2)') 'ERROR: vv2v_424 anti 2nd pair error.', iop(4)
                write (6, '(8I5)') v3%popr
                write (6, '(8I5)') v3%pops
                stop
            end if
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, v3%ioppq, v3%popr, v3%pops, iop(4))
            call isymtrso(ptr3, iCore(i2), v3%popp, v3%popq, v3%ioppq, 0, &
                          v3%popr, v3%pops, iop(4), 1, v3%sym)
            call daxpy(nSize, -1.0d0, iCore(i2), 1, ptr3, 1)
            ! compress 2nd pair
            if (iexpand(4) == 1) then
                i3 = i2 + nSize
                call icompso(ptr3, iCore(i2), v3%popp, v3%popq, v3%ioppq, 1, &
                             v3%popr, v3%pops, iop(4), v3%ioprs, v3%sym)
                ptr3(1:1) => iCore(i2:i2)
            else
                i3 = i2
            end if
        else if (iexpand(4) == 1) then
            write (6, '(A)') 'ERROR: vv2v_424 error for the 2nd pair not antisymmertize but be expand.'
            stop
        end if
        !--- add
        if (anti /= 0) then
            nSize = idsymsoc(v3%sym, v3%popp, v3%popq, v3%ioppq, v3%popr, v3%pops, v3%ioprs)
            call daxpy(nSize, 1.0d0, ptr3, 1, v3%arr, 1)
        end if

    end subroutine vv2v_424

    ! 2026.7.4
    !> Calculate v3 = v3 + fact*v1*v2
    subroutine vv2v_442(v1, v2, v3, fact, type1, type2, type3, &
                        iCore, maxCore)
        implicit none

        type(opt), intent(in) :: v1, v2
        type(opt), intent(inout) :: v3
        real*8, intent(in) :: fact
        character(len=4), intent(in) :: type1, type2 ! contract type
        character(len=2), intent(in) :: type3 ! contract type
        real*8, dimension(maxCore), intent(inout) :: iCore !> free space
        integer*8, intent(in) :: maxCore

        character(len=2), dimension(4) :: pair
        logical, dimension(4) :: flat !> whether idx need be flat
        integer*8, dimension(4) :: iop
        integer*8, dimension(4) :: iexpand !> whether idx need be expand, 1 - need, 0 - not
        integer*8 :: i
        target :: iCore
        real*8, dimension(:), pointer :: ptr1, ptr2, ptr3
        integer*8 :: i0, i1, i2, i3, nSize
        integer*8, external :: idsymsoc

        if ((.not. v1%incore) .or. (.not. v2%incore) .or. (.not. v3%incore)) then
            write (6, *) 'ERROR: vv2v_442 need all tensor be incore.', v1%incore, v2%incore, v3%incore
            stop
        end if

        i0 = 1
        ! check which pair of idx need be flat
        pair(1) = type1(1:2)
        pair(2) = type1(3:4)
        pair(3) = type2(1:2)
        pair(4) = type2(3:4)
        flat = .false.
        do i = 1, 4
            flat(i) = check_flat(pair(i), [pair(1:i - 1), pair(i + 1:4)], 3)
        end do
        ! get iop after expand
        iop(1) = v1%ioppq
        iop(2) = v1%ioprs
        iop(3) = v2%ioppq
        iop(4) = v2%ioprs
        iexpand = 0
        do i = 1, 4
            if (flat(i) .and. (iop(i) /= 1)) then
                iexpand(i) = 1
                iop(i) = 1
            end if
        end do

        !------ expand
        ! v3, i0
        i1 = i0
        ptr3 => v3%arr
        ! v1, i1
        if ((iexpand(1) == 1) .or. (iexpand(2) == 1)) then
            nSize = idsymsoc(v1%sym, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2))
            i2 = i1 + nSize
            call iexpso(v1%arr, iCore(i1), v1%popp, v1%popq, v1%ioppq, iexpand(1), &
                        v1%popr, v1%pops, v1%ioprs, iexpand(2), v1%sym)
            ptr1(1:1) => iCore(i1:i1)
        else
            i2 = i1
            ptr1 => v1%arr
        end if
        ! v2, i2
        if ((iexpand(3) == 1) .or. (iexpand(4) == 1)) then
            nSize = idsymsoc(v2%sym, v2%popp, v2%popq, iop(3), v2%popr, v2%pops, iop(4))
            i3 = i2 + nSize
            call iexpso(v2%arr, iCore(i2), v2%popp, v2%popq, v2%ioppq, iexpand(3), &
                        v2%popr, v2%pops, v2%ioprs, iexpand(4), v2%sym)
            ptr2(1:1) => iCore(i2:i2)
        else
            i3 = i2
            ptr2 => v2%arr
        end if

        !------ calculate
        call vv2v_dpd_442(ptr1, v1%popp, v1%popq, iop(1), v1%popr, v1%pops, iop(2), v1%sym, type1, &
                          ptr2, v2%popp, v2%popq, iop(3), v2%popr, v2%pops, iop(4), v2%sym, type2, &
                          ptr3, v3%popp, v3%popq, v3%ioppq, v3%sym, type3, &
                          fact, 1.0d0)

    end subroutine vv2v_442

    ! 2026.7.3
    !> Given type used in sstgenso for given permutation
    subroutine permutation_map(type, permutation, result)
        implicit none

        character(len=4), intent(in) :: type
        character(len=2), intent(in) :: permutation
        character(len=4), intent(out) :: result

        character :: char
        integer*8 :: i
        integer*8 :: idx1, idx2

        result = '1234'
        idx1 = 0
        do i = 1, 4
            if (permutation(1:1) == type(i:i)) then
                idx1 = i
                exit
            end if
        end do
        idx2 = 0
        do i = 1, 4
            if (permutation(2:2) == type(i:i)) then
                idx2 = i
                exit
            end if
        end do
        char = result(idx1:idx1)
        result(idx1:idx1) = result(idx2:idx2)
        result(idx2:idx2) = char

    end subroutine permutation_map

    ! 2026.7.3
    !> Check whether a pair of index neeb be flat
    function check_flat(current, target, n) result(flat)
        implicit none

        character(len=2), intent(in) :: current
        character(len=2), dimension(n), intent(in) :: target
        integer*8, intent(in) :: n
        logical :: flat

        integer*8 :: i

        flat = .true.
        do i = 1, n
            if (current == target(i)) then
                flat = .false.
            end if
        end do

    end function check_flat

end module mod_contract
