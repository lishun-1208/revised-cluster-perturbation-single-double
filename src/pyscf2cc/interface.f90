! write at 2025.10.30, lishun
! 1. f2py required specified order of parameter, but I don't know why
!    so the length must intent before data
! 2. f2py will compile intent(out) in subroutine as the return value
!    so the read_node must use as a fucntion
!-----------------------------------------
! for different data type, need call correspond subroutine
! use of add_node and read_node:
!   _iv, integer*8 value
!   _dv, real*8 value
!   _sv, real*4 value
!   _i, integer*8 array
!   _d, real*8 array
!   _s, real*4 array
!   _i2, integer*8 2D-array
!   _i3, integer*8 3D-array
!   _str, string
! ATTENTION: read_node_str will return a bytes, not string, need be decode
! ie: read_node_str.decode()
!-----------------------------------------
subroutine add_node_i(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    integer*8, dimension(length), intent(in) :: data

    call add_node(name, data, length)

end subroutine add_node_i

subroutine add_node_d(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    real*8, dimension(length), intent(in) :: data

    call add_node(name, data, length)

end subroutine add_node_d

subroutine add_node_s(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    real*4, dimension(length), intent(in) :: data

    call add_node(name, data, length)

end subroutine add_node_s

subroutine add_node_iv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: data

    call add_node(name, data)

end subroutine add_node_iv

subroutine add_node_dv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    real*8, intent(in) :: data

    call add_node(name, data)

end subroutine add_node_dv

subroutine add_node_sv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    real*4, intent(in) :: data

    call add_node(name, data)

end subroutine add_node_sv

subroutine add_node_i2(name, data, dim1, dim2)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: dim1, dim2
    integer*8, dimension(dim1, dim2), intent(in) :: data

    call add_node(name, data, dim1, dim2)

end subroutine add_node_i2

subroutine add_node_i3(name, data, dim1, dim2, dim3)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: dim1, dim2, dim3
    integer*8, dimension(dim1, dim2, dim3), intent(in) :: data

    call add_node(name, data, dim1, dim2, dim3)

end subroutine add_node_i3

subroutine add_node_str(name, data, size)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: size ! character length of data
    character(len=size), intent(in) :: data

    call add_node(name, data, size)

end subroutine add_node_str

subroutine read_node_i(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    integer*8, dimension(length), intent(out) :: data

    call read_node(name, data, length)

end subroutine read_node_i

subroutine read_node_d(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    real*8, dimension(length), intent(out) :: data

    call read_node(name, data, length)

end subroutine read_node_d

subroutine read_node_s(name, data, length)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: length
    real*4, dimension(length), intent(out) :: data

    call read_node(name, data, length)

end subroutine read_node_s

subroutine read_node_iv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(out) :: data

    call read_node(name, data)

end subroutine read_node_iv

subroutine read_node_dv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    real*8, intent(out) :: data

    call read_node(name, data)

end subroutine read_node_dv

subroutine read_node_sv(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    real*4, intent(out) :: data

    call read_node(name, data)

end subroutine read_node_sv

subroutine read_node_i2(name, data, dim1, dim2)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: dim1, dim2
    integer*8, dimension(dim1, dim2), intent(out) :: data

    call read_node(name, data, dim1, dim2)

end subroutine read_node_i2

subroutine read_node_i3(name, data, dim1, dim2, dim3)
    use mod_interface

    character(len=*), intent(in) :: name
    integer*8, intent(in) :: dim1, dim2, dim3
    integer*8, dimension(dim1, dim2, dim3), intent(out) :: data

    call read_node(name, data, dim1, dim2, dim3)

end subroutine read_node_i3

subroutine read_node_str(name, data)
    use mod_interface

    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: data

    call read_node(name, data)

end subroutine read_node_str
