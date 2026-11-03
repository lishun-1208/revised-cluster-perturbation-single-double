! split integral and calculate norm
!-----------------------------------
! 2025.12.16
subroutine cal_spin_norm(Vpqrs, type, ioppq_ns, ioprs_ns, symV, spin, name, iCore, maxCore)
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(in) :: Vpqrs ! no spin integral
    character(len=4), intent(in) :: type ! 'o'/'v'/'n'
    integer*8, intent(in) :: ioppq_ns, ioprs_ns
    integer*8, intent(in) :: symV
    character(len=4), intent(in) :: spin ! such as 'A'/'B'
    character(len=*), intent(in) :: name
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: i
    integer*8, dimension(8, 4) :: pop_ns, pop_a, pop_b, pop_s
    integer*8 :: ioppq, ioprs
    integer*8 :: i0, i1, nSize
    real*8 :: norm
    integer*8, external :: idsymsoc
    real*8, external :: dnrm2

    !--- determine pop
    do i = 1, 4
        if (type(i:i) == 'o') then
            pop_ns(:, i) = no(:)
            pop_a(:, i) = popa(:)
            pop_b(:, i) = popb(:)
        else if (type(i:i) == 'v') then
            pop_ns(:, i) = nv(:)
            pop_a(:, i) = vrta(:)
            pop_b(:, i) = vrtb(:)
        else if (type(i:i) == 'n') then
            pop_ns(:, i) = zorb(:)
            pop_a(:, i) = zorb(:)
            pop_b(:, i) = zorb(:)
        end if
    end do
    do i = 1, 4
        if (spin(i:i) == 'A') then
            pop_s(:, i) = pop_a(:, i)
        else if (spin(i:i) == 'B') then
            pop_s(:, i) = pop_b(:, i)
        end if
    end do
    !--- determine iop
    if (ioppq_ns == 0 .and. spin(1:1) == spin(2:2)) then
        ioppq = 0
    else
        ioppq = 1
    end if
    if (ioprs_ns == 0 .and. spin(3:3) == spin(4:4)) then
        ioprs = 0
    else
        ioprs = 1
    end if

    i0 = 1
    nSize = idsymsoc(symV, pop_s(:, 1), pop_s(:, 2), ioppq, pop_s(:, 3), pop_s(:, 4), ioprs)
    i1 = i0 + nSize
    call dscal(nSize, 0.0d0, iCore(i0), 1)
    call splitintegral(Vpqrs, iCore(i0), &
                       pop_a(:, 1), pop_b(:, 1), pop_a(:, 2), pop_b(:, 2), ioppq, &
                       pop_a(:, 3), pop_b(:, 3), pop_a(:, 4), pop_b(:, 4), ioprs, &
                       symV, &
                       pop_ns(:, 1), pop_ns(:, 2), ioppq_ns, &
                       pop_ns(:, 3), pop_ns(:, 4), ioprs_ns, &
                       symV, spin, iCore(i1))
    norm = dnrm2(nSize, iCore(i0), 1)
    write (6, '(5A,F)') 'norm of ', name, ' with spin ', spin, ':', norm

end subroutine cal_spin_norm
