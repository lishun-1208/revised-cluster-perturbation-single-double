! write by lishun, 2024.10.22
! 2026.10.22
! calculate correlation energy
subroutine calEcorr(E, iCore, istart, order)
    use mod_variable
    use mod_orbit
    use mod_operator
    implicit none

    real*8, intent(out) :: E !> E^(k+1)
    real*8, dimension(*), intent(inout) :: iCore !> full work space
    integer*8, intent(in) :: istart
    integer*8, intent(in) :: order !> k

    type(opt) :: tmpT2, tmpT1_1, tmpT1_2
    integer*8 :: i, i0, i1, i2, nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: ddot

    i0 = istart

    call copy_opt(T2, tmpT2, T, iCore, i0)
    call dcopy(T2%len, T2%arr, 1, tmpT2%arr, 1)
    i0 = i0 + tmpT2%len

    ! E^(k+1) = 0.25*V(kl,cd)*t(cd,kl)^(k)
    if (order == 2) then
        ! E^(3) += 0.5*V(kl,cd)*t(c,k)^(1)*t(d,l)^(1)
        !---- read T(a,i)^(1)
        call copy_opt(T1, tmpT1_1, T, iCore, i0)
        i0 = i0 + T1%len
        call IO_stream(0.0d0, tmpT1_1%arr, 1, T1%len, f_T, '13', 'old')

        call effectT2(tmpT2, tmpT1_1, tmpT1_1, 1.0d0, iCore(i0), nCore - i0 + 1, 1_8)
    else if (order == 3) then
        ! E^(4) += V(kl,cd)*t(c,k)^(1)*t(d,l)^(2)
        !---- read T(a,i)^(2)
        call copy_opt(T1, tmpT1_1, T, iCore, i0)
        i0 = i0 + T1%len
        call IO_stream(0.0d0, tmpT1_1%arr, T1%len + T2%len + 1, T1%len, f_T, '13', 'old')
        !---- read in T(a,i)^(1)
        call copy_opt(T1, tmpT1_2, T, iCore, i0)
        i0 = i0 + T1%len
        call IO_stream(0.0d0, tmpT1_2%arr, 1, T1%len, f_T, '13', 'old')

        call effectT2(tmpT2, tmpT1_1, tmpT1_2, 1.0d0, iCore(i0), nCore - i0 + 1, 2_8)
    end if
    E = ddot(tmpT2%len, tmpT2%arr, 1, abij%arr, 1)

end subroutine calEcorr
