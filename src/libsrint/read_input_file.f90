! read information from input file CCSOINP
! store in file interface as string, int64, and int64 array
! all keyword will change to lower case, but string value not change
!---------------------------------
! 2025.12.14, lish
subroutine read_input_file
    implicit none

    character(len=7), parameter :: input_f = 'CCSOINP'
    character(len=200), dimension(1000) :: input
    character(200) :: keywords, values
    integer*8 :: istat
    integer*8 :: i, imax
    logical :: status

    !------ read input file
    !--- test exist
    inquire (file=input_f, exist=status)
    if (.not. status) then
        write (6, *) 'WARN: file ', input_f, ' not exist'
        return
    end if
    !--- read
    write (6, '(A)') '--------------------------------------------------'
    write (6, '(3A)') '----------------- file ', input_f, ' -------------------'
    open (unit=56, file=input_f, form='formatted', action='read', status='old')
    rewind (56)
    imax = 0
    do
        ! read line
        imax = imax + 1
        read (56, '(A)', iostat=istat) input(imax)

        if (istat == 0) then
            write (6, *) input(imax)
        else if (istat == -1) then
            imax = imax - 1
            write (6, '(A)') '--------------- input file end -------------------'
            exit
        else if (istat > 0) then
            write (6, *) '---------- error in read input file --------------'
            stop
        end if
    end do
    close (56)
    write (6, '(A)') '--------------------------------------------------'

    !---------- get value and write in interface
    do i = 1, imax
        call read_line(input, i, keywords, values)
        if (len_trim(keywords) == 0 .or. len_trim(values) == 0) then
            cycle
        end if
        call write_in_interface(keywords, values)
    end do

end subroutine read_input_file

subroutine read_line(input, line, key, value)
    implicit none

    character(len=200), dimension(*), intent(in) :: input
    integer*8, intent(in) :: line
    character(len=*), intent(out) :: key, value

    integer*8 :: idx_equal ! index of '='

    idx_equal = index(input(line), '=')
    ! remove ' ' at left of string
    key = trim(adjustl(input(line) (1:idx_equal - 1)))
    value = trim(adjustl(input(line) (idx_equal + 1:)))

end subroutine read_line

subroutine write_in_interface(key, value)
    use mod_interface
    implicit none

    character(len=*), intent(in) :: key, value

    character(len=:), allocatable :: strValue
    integer*8 :: intValue
    real*8 :: floatValue
    integer*8, dimension(:), allocatable :: intVec
    integer*8 :: len, i, j, n
    integer*8, dimension(:), allocatable :: idxSep ! position of separator: '-', '/', ','
    logical :: int, float, str, sep
    logical, external :: check_integer, check_float

    len = len_trim(value)

    if (check_integer(value)) then
        read (value(1:len), '(I)') intValue
        write (6, '(2A,I)') trim(key), ' = ', intValue
        call add_node(trim(key), intValue)

    else if (check_float(value)) then
        read (value(1:len), '(F)') floatValue
        write (6, '(2A,F)') trim(key), ' = ', floatValue
        call add_node(trim(key), floatValue)

    else
        !--- determine whether have seperator
        n = 0 ! number of seperator
        do i = 1, len
            if (value(i:i) == '-' .or. value(i:i) == '/' .or. value(i:i) == ',') then
                n = n + 1
            end if
        end do

        if (n /= 0) then
            sep = .true.
            ! position of seperator
            allocate (idxSep(n + 2)) ! include (1,idx of seperator,len)
            idxSep(1) = 0
            j = 2
            do i = 1, len
                if (value(i:i) == '-' .or. value(i:i) == '/' .or. value(i:i) == ',') then
                    idxSep(j) = i
                    j = j + 1
                end if
            end do
            idxSep(n + 2) = len + 1
            allocate (intVec(n + 1))
        else
            sep = .false.
        end if

        if (.not. sep) then
            !--- pure string
            strValue = value(1:len)
            write (6, '(3A)') trim(key), ' = ', strValue
            call add_node(trim(key), strValue, len)
        else
            !--- integer vector
            do i = 1, n + 1
                read (value(idxSep(i) + 1:idxSep(i + 1) - 1), '(I)') intVec(i)
            end do
            write (6, *) trim(key), ' = ', intVec
            call add_node(trim(key), intVec, n + 1)
        end if
    end if

end subroutine write_in_interface

function up2low(str) result(newStr)
    implicit none

    character(len=*), intent(inout) :: str
    character(len=:), allocatable :: newStr

    integer*8 :: i, length
    character(len=1) :: c

    newStr = str
    length = len(str)
    do i = 1, length
        if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
            c = char(iachar(str(i:i)) + 32)
            newStr(i:i) = c
        end if
    end do
end function up2low

function check_integer(str) result(res)
    implicit none

    character(*), intent(in) :: str
    logical :: res
    integer :: i, n

    res = .true.
    n = len_trim(str)
    if (n == 0) then
        res = .false.
        return
    end if

    do i = 1, n
        ! 第一位可以是 + 或 -
        if (i == 1 .and. (str(i:i) == '+' .or. str(i:i) == '-')) then
            cycle
        end if

        ! 其他必须是 0-9
        if (str(i:i) < '0' .or. str(i:i) > '9') then
            res = .false.
            exit
        end if
    end do
end function check_integer

function check_float(str) result(res)
    implicit none

    character(*), intent(in) :: str
    logical :: res
    integer :: i, n, dot, exp
    logical :: valid

    n = len_trim(str)
    res = .true.
    dot = 0    ! 小数点计数
    exp = 0    ! 科学计数法 e/E 计数

    if (n == 0) then
        res = .false.
        return
    end if

    do i = 1, n
        valid = .false.

        ! 1. 数字：永远合法
        if (str(i:i) >= '0' .and. str(i:i) <= '9') valid = .true.

        ! 2. 正负号：只能在第1位 或 e/E 后面
        if (str(i:i) == '+' .or. str(i:i) == '-') then
            if (i == 1 .or. str(i - 1:i - 1) == 'e' .or. str(i - 1:i - 1) == 'E') then
                valid = .true.
            end if
        end if

        ! 3. 小数点：只能有1个，不能在 e/E 之后
        if (str(i:i) == '.') then
            dot = dot + 1
            if (dot <= 1 .and. exp == 0) valid = .true.
        end if

        ! 4. e/E：只能有1个
        if (str(i:i) == 'e' .or. str(i:i) == 'E') then
            exp = exp + 1
            if (exp <= 1) valid = .true.
        end if

        if (.not. valid) then
            res = .false.
            exit
        end if
    end do

    ! 浮点数必须包含 . 或 e/E
    if (dot == 0 .and. exp == 0) res = .false.
end function check_float
