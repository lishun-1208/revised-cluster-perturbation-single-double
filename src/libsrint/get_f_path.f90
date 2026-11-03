! modified at 2025.10.13, lishun
subroutine get_f_path(file_path)
    use mod_iop
    implicit none

    character(len=200), intent(out) :: file_path

    logical :: dir_exist

    file_path = trim(tmpdir//'/')

    inquire (file=trim(tmpdir), exist=dir_exist)

    if (.not. dir_exist) then
        call system('mkdir -p '//trim(tmpdir))
    end if
end subroutine get_f_path
