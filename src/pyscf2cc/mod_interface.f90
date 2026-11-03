! write at 2025.10.28, lishun
! this module is used to write and read interface file, not use for libsrint
! all keyword will be change to lower case
!-------------------------------------------
! structure of IOfile like this:
! metaPos,nNode,(data1,data2,...),(metadata1,metadata2,...)
!-------------------------------------------
! 2025.12.12
module mod_interface
    implicit none
    private

    character(len=9), parameter :: IOfile = 'interface'
    integer*8 :: metaPos ! position of meta data

    integer*8, parameter :: nameLen = 20
    type :: node
        character(len=nameLen) :: name ! data name
        integer*8 :: pos ! position of data in IOfile, unit=1byte
        integer*8 :: size ! length of datas in IOfile, unit=1byte
    end type node

    type(node), dimension(1000) :: dataset
    integer*8, parameter :: nBit = 36 ! size of per node, unit=1byte
    integer*8 :: nNode ! number of nodes

    interface add_node
        ! parameter be: add_node(name,data,size), for array
        module procedure add_node_i
        module procedure add_node_d
        module procedure add_node_s
        ! parameter be：add_node(name,data), for number
        module procedure add_node_iv
        module procedure add_node_dv
        module procedure add_node_sv
        ! parameter be: add_node(name,data,dim1,dim2,...), for nD-array
        module procedure add_node_i2
        module procedure add_node_i3
        ! parameter be: read_node(name,str,length), for string
        module procedure add_node_str
    end interface add_node

    interface read_node
        ! parameter be: read_node(name,data,size), for array
        module procedure read_node_i
        module procedure read_node_d
        module procedure read_node_s
        ! parameter be: read_node(name,data), for number
        module procedure read_node_iv
        module procedure read_node_dv
        module procedure read_node_sv
        ! parameter be: read_node(name,data,dim1,dim2,...), for nD-array
        module procedure read_node_i2
        module procedure read_node_i3
        ! parameter be: read_node(name,str,length), for string
        module procedure read_node_str
    end interface read_node

    public :: add_node, read_node, init_node, grep_node
contains
    !> change up to lower case
    subroutine up2low(str)
        implicit none

        character(len=*), intent(inout) :: str

        integer*8 :: i, length
        character(len=1) :: c

        length = len(str)
        do i = 1, length
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
                c = char(iachar(str(i:i)) + 32)
                str(i:i) = c
            end if
        end do
    end subroutine up2low

    !> initial variable in IOfile
    subroutine init_node
        implicit none

        metaPos = 17 ! 2*8+1
        nNode = 0

        ! write in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='replace', access='stream')
        write (6219, pos=1) metaPos, nNode
        close (6219)

    end subroutine init_node

    !> test whether file exist, if not create it
    subroutine file_exist
        implicit none

        logical :: status

        inquire (file=IOfile, exist=status)
        if (.not. status) then
            call init_node
        end if

    end subroutine file_exist

    !> read all meta data in IOfile
    subroutine read_dataset
        implicit none

        integer*8 :: i

        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=1) metaPos, nNode
        if (nNode > 0) then
            read (6219, pos=metaPos) (dataset(i), i=1, nNode)
        end if
        close (6219)

    end subroutine read_dataset

    !> write all meta data in IOfile;
    !> only use in add_node after find_node
    subroutine write_dataset
        implicit none

        integer*8 :: i

        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=1) metaPos, nNode
        write (6219, pos=metaPos) (dataset(i), i=1, nNode)
        close (6219)

    end subroutine write_dataset

    !> read_dataset and find index of node we want
    subroutine find_node(name, exist, idx)
        implicit none

        character(len=nameLen), intent(in) :: name !> node name
        logical, intent(out) :: exist
        integer*8, intent(out) :: idx !> index of node, if not exist, set as nNode+1

        integer*8 :: position
        integer*8 :: i

        exist = .false.
        do i = 1, nNode
            if (trim(dataset(i)%name) == trim(name)) then
                idx = i
                exist = .true.
                return
            end if
        end do
        if (.not. exist) then
            idx = nNode + 1
        end if

    end subroutine find_node

    !> find node whether exit
    subroutine grep_node(name, exist)
        implicit none

        character(len=*), intent(in) :: name
        logical, intent(out) :: exist

        character(len=nameLen) :: fullName
        integer*8 :: idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)

    end subroutine grep_node

    !> calculate values for node
    subroutine get_node(name, size, idx, exist, unit)
        implicit none

        character(len=nameLen), intent(in) :: name !> node name
        integer*8, intent(in) :: size !> number of data for node(idx)
        integer*8, intent(in) :: idx !> index of node in dataset
        logical, intent(in) :: exist
        integer*8, intent(in) :: unit !> length of per data, unit=byte

        if (.not. exist) then
            nNode = nNode + 1
            dataset(idx)%name = name
            dataset(idx)%pos = metaPos
            dataset(idx)%size = unit*size
            metaPos = dataset(idx)%pos + dataset(idx)%size
        else
            write (6, '(2A,I5)') 'old node: ', trim(name), idx
            if (dataset(idx)%size /= unit*size) then
                write (6, '(2A,I12,I12)') 'WARN: size of new data different from before. ', &
                    trim(name), dataset(idx)%size, unit*size

                if (dataset(idx)%size < unit*size) then
                    write (6, *) 'write new data in new space'
                    dataset(idx)%size = unit*size
                    dataset(idx)%pos = metaPos
                    metaPos = dataset(idx)%pos + dataset(idx)%size
                else
                    write (6, *) 'write new data in old space'
                    dataset(idx)%size = unit*size
                end if
            end if
        end if

    end subroutine get_node

    !> make string with len=nameLen, and change to lower case
    function full_name(name) result(newName)
        character(len=*), intent(in) :: name
        character(len=nameLen) :: newName

        integer*8 :: length
        integer*8 :: i

        length = len(name)
        if (length > nameLen) then
            write (6, *) 'ERROR: length of node name need <= ', nameLen, name
            stop
        else
            newName(1:length) = name
            do i = length + 1, nameLen
                newName(i:i) = ' '
            end do
        end if
        call up2low(newName)

    end function full_name

    !> add integer*8 array
    subroutine add_node_i(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(size), intent(in) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) (data(i), i=1, size)
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_i

    !> add real*8 array
    subroutine add_node_d(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        real*8, dimension(size), intent(in) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) (data(i), i=1, size)
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_d

    !> add real*4 array
    subroutine add_node_s(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        real*4, dimension(size), intent(in) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 4
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) (data(i), i=1, size)
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_s

    !> add integer*8 value
    subroutine add_node_iv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, intent(in) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) data
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_iv

    !> add real*8 value
    subroutine add_node_dv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        real*8, intent(in) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) data
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_dv

    !> add real*4 value
    subroutine add_node_sv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        real*4, intent(in) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 4
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) data
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_sv

    !> add integer*8 2d-array
    subroutine add_node_i2(name, data, dim1, dim2)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(dim1, dim2), intent(in) :: data
        integer*8, intent(in) :: dim1, dim2

        call add_node_i(name, reshape(data, [dim1*dim2]), dim1*dim2)

    end subroutine add_node_i2

    !> add integer*8 3d-array
    subroutine add_node_i3(name, data, dim1, dim2, dim3)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(dim1, dim2, dim3), intent(in) :: data
        integer*8, intent(in) :: dim1, dim2, dim3

        call add_node_i(name, reshape(data, [dim1*dim2*dim3]), dim1*dim2*dim3)

    end subroutine add_node_i3

    !> add string
    subroutine add_node_str(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, intent(in) :: size !> string length
        character(len=size), intent(in) :: data

        character(len=nameLen) :: fullName
        integer*8 :: unit = 1
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        ! decide value of node
        call get_node(fullName, size, idx, exist, unit)

        position = dataset(idx)%pos
        ! write data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='write', &
              status='old', access='stream')
        write (6219, pos=position) data
        close (6219)

        ! write meta data in IOfile
        call write_dataset

    end subroutine add_node_str

    !> read integer*8 array
    subroutine read_node_i(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(size), intent(out) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) (data(i), i=1, size)
        close (6219)

    end subroutine read_node_i

    !> read real*8 array
    subroutine read_node_d(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        real*8, dimension(size), intent(out) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) (data(i), i=1, size)
        close (6219)

    end subroutine read_node_d

    !> read real*4 array
    subroutine read_node_s(name, data, size)
        implicit none

        character(len=*), intent(in) :: name
        real*4, dimension(size), intent(out) :: data
        integer*8, intent(in) :: size

        character(len=nameLen) :: fullName
        integer*8 :: unit = 4
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) (data(i), i=1, size)
        close (6219)

    end subroutine read_node_s

    !> read integer*8 value
    subroutine read_node_iv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, intent(out) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) data
        close (6219)

    end subroutine read_node_iv

    !> read real*8 value
    subroutine read_node_dv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        real*8, intent(out) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 8
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) data
        close (6219)

    end subroutine read_node_dv

    !> read real*4 value
    subroutine read_node_sv(name, data)
        implicit none

        character(len=*), intent(in) :: name
        real*4, intent(out) :: data

        integer*8 :: size = 1
        character(len=nameLen) :: fullName
        integer*8 :: unit = 4
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        if (dataset(idx)%size /= unit*size) then
            write (6, '(2A,I12,I12)') 'ERROR: array len not correct. ', &
                name, dataset(idx)%size, unit*size
            stop
        end if
        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) data
        close (6219)

    end subroutine read_node_sv

    !> read integer*8 2d-array
    subroutine read_node_i2(name, data, dim1, dim2)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(dim1, dim2), intent(out) :: data
        integer*8, intent(in) :: dim1, dim2

        integer*8, dimension(dim1*dim2) :: tmp

        call read_node_i(name, tmp, dim1*dim2)
        data = reshape(tmp, [dim1, dim2])

    end subroutine read_node_i2

    !> read integer*8 3d-array
    subroutine read_node_i3(name, data, dim1, dim2, dim3)
        implicit none

        character(len=*), intent(in) :: name
        integer*8, dimension(dim1, dim2, dim3), intent(out) :: data
        integer*8, intent(in) :: dim1, dim2, dim3

        integer*8, dimension(dim1*dim2*dim3) :: tmp

        call read_node_i(name, tmp, dim1*dim2*dim3)
        data = reshape(tmp, [dim1, dim2, dim3])

    end subroutine read_node_i3

    !> read string
    subroutine read_node_str(name, data)
        implicit none

        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: data

        character(len=nameLen) :: fullName
        integer*8 :: unit = 1
        logical :: exist
        integer*8 :: position
        integer*8 :: i, idx
        character(len=:), allocatable :: tmp

        fullName = full_name(name)
        ! test file
        call file_exist
        ! read dataset
        call read_dataset
        ! test node be exist or not
        call find_node(fullName, exist, idx)
        if (.not. exist) then
            write (6, *) 'ERROR: node not exist ', fullName
            stop
        end if

        position = dataset(idx)%pos
        ! read data in IOfile
        open (unit=6219, file=IOfile, form='unformatted', action='read', &
              status='old', access='stream')
        read (6219, pos=position) (data(i:i), i=1, dataset(idx)%size)
        close (6219)

    end subroutine read_node_str

end module mod_interface

