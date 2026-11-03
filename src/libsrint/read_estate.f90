subroutine read_estate(estate_sym, input_string)
    implicit none
    character(len=100) :: input_string
    integer, dimension(8) :: estate_sym
    integer :: n, i, start, end, len_str
    character(len=100) :: temp_str, num_str
    integer :: num

    estate_sym = 0
    temp_str = trim(adjustl(input_string))
    len_str = len(temp_str)
    n = 1
    start = 1

    ! 循环遍历字符串
    do while (start <= len_str .and. n <= 8)
        ! 查找下一个分隔符的位置
        end = index(temp_str(start:len_str), "-")
        if (end == 0) then
            ! 没有找到分隔符，说明是最后一个数字
            num_str = temp_str(start:len_str)
            read(num_str, *) estate_sym(n)
            exit
        else
            ! 找到分隔符，提取数字
            num_str = temp_str(start:start+end-2)
            read(num_str, *) estate_sym(n)
            start = start + end  ! 更新起始位置
            n = n + 1
        endif
    end do

!   print *, '读取的数字为: ', estate_sym
return
end 
