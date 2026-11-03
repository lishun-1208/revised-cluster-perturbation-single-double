! write by LiShun, 2024.9.19
! modified at 2026.5.15
!-------------------------------------------
!> prepare initial guess for spin free calculation
! sort as (e1,Rs1),(e2,Rs2),...
subroutine preDavidson1_ns(Eorb, pop_a, pop_i, iCore, maxCore)
    use mod_orbit
    use mod_symm
    use mod_interface
    use mod_operator
    use mod_common
    implicit none

    real*8, dimension(*), intent(in) :: Eorb !> spin free orbital energy
    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: iEa, iEb ! alpha/beta orbital energy
    integer*8 :: i0
    character(len=5) :: string
    logical :: exist, specified_ex
    integer*8, dimension(5) :: EX_specified1 = 0 !> (sym,irrepi,i,irrepa,a)
    character(len=2) :: spin_flip = 'no'
    integer*8 :: ex_mult = 1
    logical :: RHF

    !------ determine reference state
    call grep_node('shell', exist)
    if (exist) then
        call read_node('shell', string)
        if (string(1:5) == 'close') then
            RHF = .true.
        else if (string(1:4) == 'open') then
            RHF = .false.
        else
            write (6, *) 'ERROR: value of shell error, please use lower case.', string(1:3)
            stop
        end if
    else
        write (6, *) 'ERROR: please add shell=open/close in CCSOINP'
        write (6, *) 'this keyword use different form ref_state'
        stop
    end if
    if (.not. RHF .and. CVS) then
        write (6, *) 'ERROR: currently only support close shell CVS'
        stop
    end if

    !------ determine whether specified excitation
    call grep_node('ex_specified1', specified_ex)
    if (specified_ex) then
        call read_node('ex_specified1', EX_specified1, 5)
    end if
    if (specified_ex .and. CVS) then
        write (6, *) 'ERROR: currently not support CVS with specified excitation'
        stop
    end if

    !------ determine multiple of RHF excitated state
    call grep_node('ex_mult', exist)
    if (exist) then
        call read_node('ex_mult', ex_mult)
    end if

    !------ determine spin_flip
    call grep_node('spin_flip', exist)
    if (exist) then
        call read_node('spin_flip', spin_flip)
    end if

    !------ seperate orbital energy
    iEa = 1
    iEb = iEa + nva + noa
    i0 = iEb + nvb + nob
    call seperateEorb(Eorb, iCore(iEa), iCore(iEb))

    !------ for spin flip calculation
    if (spin_flip == 'ab') then
        call p_pre_spin_flip_ns(pop_a, pop_i, iCore(iEa), iCore(iEb), &
                                iCore(i0), maxCore - i0 + 1, .true.)
        return
    else if (spin_flip == 'ba') then
        if (CVS) then
            write (*, *) 'ERROR: current only support spin flip a->b for CVS'
            stop
        end if
        call p_pre_spin_flip_ns(pop_a, pop_i, iCore(iEa), iCore(iEb), &
                                iCore(i0), maxCore - i0 + 1, .false.)
        return
    else if (spin_flip /= 'no') then
        write (6, *) 'ERROR: value of keyword spin_flip error'
        stop
    end if

    if (RHF) then
        if (ex_specified1(1) /= 0) then
            !------- for specified single excitation
            call p_pre_close_specified_ns(pop_a, pop_i, iCore(iEa), iCore(iEb), &
                                          iCore(i0), maxCore - i0 + 1, EX_specified1, ex_mult)
        else
            write (6, *) 'Set the 1e lowest OED as initial guess'
            call p_pre_close_ns(pop_a, pop_i, iCore(iEa), iCore(iEb), &
                                iCore(i0), maxCore - i0 + 1, EX_specified1, ex_mult)
        end if
    else if (.not. RHF) then
        if (ex_specified1(1) /= 0) then
            write (6, '(A)') 'ERROR: current not support open shell specified excited state calculation'
            stop
        else
            call p_pre_open_ns(pop_a, pop_i, iCore(iEa), iCore(iEb), &
                               iCore(i0), maxCore - i0 + 1)
        end if
    end if

end subroutine preDavidson1_ns

! 2026.5.15
!> prepare initial guess for spin flip
subroutine p_pre_spin_flip_ns(pop_a, pop_i, Ea, Eb, iCore, maxCore, ab)
    use mod_symm
    use mod_orbit
    use mod_operator
    use mod_common
    implicit none

    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    real*8, dimension(*), intent(inout) :: Ea !> orbital energy alpha
    real*8, dimension(*), intent(inout) :: Eb !> orbital energy beta
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    logical, intent(in) :: ab !> whether be alphe->beta

    integer*8 :: irrepex !> symmerty of excited state
    integer*8 :: spatial_sym !> spatial symmetry of excited state
    integer*8 :: nroot, iroot
    integer*8, dimension(:), allocatable :: lowest_e !> index of lowest OED1
    real*8 :: exE !> initial excitation energy
    integer*8 :: i0, i1, i2, nSize, idx
    integer*8 :: v_len !> length of the eigenvector (r0,R1,R2,...)
    integer*8, dimension(8) :: exInfor !> (irrepp,iirepq,irrepr,irreps,p,q,r,s)
    type(opt) :: tmp
    integer*8, external :: irpdso, idsymsoc

    i0 = 1
    do irrepex = 1, nirrep
        spatial_sym = dirprd(irrepex, symRy)
        write (6, '(A,2I2)') '-------- symmetry of excitation state :', irrepex, spatial_sym
        nroot = estate_sym(irrepex) ! number of root
        if (nroot == 0) then
            cycle
        end if
        allocate (lowest_e(nroot))

        v_len = irpdso(irrepex, pop_a, pop_i, 1)

        ! calculate OED(a,i)=Eb(a)-Ea(i), i0
        if (ab) then ! a -> b
            nSize = irpdso(spatial_sym, vrtb, popa, 1)
            i1 = i0 + nSize
            call Eorb_diff_ab(iCore(i0), Eb, vrtb, popb, &
                              Ea, vrta, popa, spatial_sym)
            if (CVS) then
                call getpqrso_incore(iCore(i1), core_idx(CVS_orb), core_sym(CVS_orb), &
                                     zorb, zorb, 1, vrtb, popa, 1, iCore(i0), irrepex)
                iCore(i0:i1 - 1) = 1.0d10
                call putpqrso_incore(iCore(i1), core_idx(CVS_orb), core_sym(CVS_orb), &
                                     zorb, zorb, 1, vrtb, popa, 1, iCore(i0), irrepex)
            end if
        else ! b -> a
            nSize = irpdso(spatial_sym, vrta, popb, 1)
            i1 = i0 + nSize
            call Eorb_diff_ab(iCore(i0), Ea, vrta, popa, &
                              Eb, vrtb, popb, spatial_sym)
        end if

        ! get position of lowest OED
        call nmin(iCore(i0), nSize, nroot, lowest_e)

        ! creat file
        call IO_stream(0.0d0, iCore(i0), 1, &
                       0, f_init_davidson(irrepex), '22', 'replace')

        ! cycle for root
        do iroot = 1, nroot
            write (6, '(A,I2,I3,A)') '------ root', irrepex, iroot, ' ------'
            exE = iCore(i0 - 1 + lowest_e(iroot))
            write (6, '(A,F)') 'initial excitation energy:', exE
            !-------- get information of initiall guess energy
            if (ab) then
                call getidx(lowest_e(iroot), exInfor, &
                            vrtb, popa, 1, zorb, zorb, 1, spatial_sym)
            else
                call getidx(lowest_e(iroot), exInfor, &
                            vrta, popb, 1, zorb, zorb, 1, spatial_sym)
            end if

            !------- initialize eigenvector, i1
            i2 = i1 + v_len
            call dscal(v_len, 0.0d0, iCore(i1), 1)

            if (ab) then
                exInfor(1) = dirprd(exInfor(1), symRy) ! a be beta
                exInfor(5) = vrta(exInfor(1)) + exInfor(5)
            else
                exInfor(2) = dirprd(exInfor(2), symRy) ! i be beta
                exInfor(6) = popa(exInfor(2)) + exInfor(6)
            end if
            call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                        exInfor(5), exInfor(6), 1, 1, &
                        pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
            iCore(i1 + idx - 1) = 1.0d0

            call IO_stream(exE, iCore(i1), (v_len + 1)*(iroot - 1) + 1, &
                           v_len, f_init_davidson(irrepex), '21', 'old')

            call tmp%init(iCore, pop_a, pop_i, zorb, zorb, 1, 1, irrepex, T, i1)
            call printEXinfo1_uhf(tmp, 4)
        end do
        deallocate (lowest_e)
    end do

end subroutine p_pre_spin_flip_ns

! 2026.5.15
!> prepare initial guess for open shell excitation
subroutine p_pre_open_ns(pop_a, pop_i, Ea, Eb, iCore, maxCore)
    use mod_symm
    use mod_orbit
    use mod_operator
    use mod_common
    implicit none

    !> file to store eigen value and vector, like (E,r0,R1,R2)
    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    real*8, dimension(*), intent(inout) :: Ea !> orbital energy alpha
    real*8, dimension(*), intent(inout) :: Eb !> orbital energy beta
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: irrepex !> symmerty of excited state
    integer*8 :: spatial_sym !> spatial symmetry of excited state
    integer*8 :: nroot, iroot
    integer*8, dimension(:), allocatable :: lowest_e !> index of lowest OED1
    real*8 :: exE !> initial excitation energy
    integer*8 :: i0, i1, i2, i3, nSize, idx
    integer*8 :: v_len !> length of the eigenvector (r0,R1,R2,...)
    integer*8, dimension(8) :: exInfor !> (irrepp,iirepq,irrepr,irreps,p,q,r,s)
    type(opt) :: tmp
    integer*8, external :: irpdso, idsymsoc

    i0 = 1
    do irrepex = 1, nirrep
        write (6, '(A,I5)') '-------- symmetry of excitation state :', irrepex
        nroot = estate_sym(irrepex)
        if (nroot == 0) then
            cycle
        end if
        allocate (lowest_e(nroot))

        v_len = irpdso(irrepex, pop_a, pop_i, 1)

        !------ calculate OED(a,i)=E(a)-E(i), i0
        ! alpha, i0
        nSize = irpdso(irrepex, vrta, popa, 1)
        i1 = i0 + nSize
        call Eorb_diff1(iCore(i0), Ea, vrta, popa, irrepex)
        ! beta, i1
        nSize = irpdso(irrepex, vrtb, popb, 1)
        i2 = i1 + nSize
        call Eorb_diff1(iCore(i1), Eb, vrtb, popb, irrepex)

        ! get position of lowest OED
        nSize = irpdso(irrepex, vrta, popa, 1) + irpdso(irrepex, vrtb, popb, 1)
        call nmin(iCore(i0), nSize, nroot, lowest_e)

        !------ create file
        call IO_stream(0.0d0, iCore(i2), 1, &
                       0, f_init_davidson(irrepex), '22', 'replace')

        ! cycle for root
        do iroot = 1, nroot
            write (6, '(A,I2,I3,A)') '------ root', irrepex, iroot, ' ------'
            exE = iCore(i0 - 1 + lowest_e(iroot))
            write (6, '(A,F)') 'initial excitation energy:', exE

            !-------- get information of initiall guess energy
            nSize = irpdso(irrepex, vrta, popa, 1)
            if (lowest_e(iroot) <= nSize) then
                ! alpha
                call getidx(lowest_e(iroot), exInfor, &
                            vrta, popa, 1, zorb, zorb, 1, irrepex)
            else
                ! beta
                idx = lowest_e(iroot) - nSize
                call getidx(idx, exInfor, &
                            vrtb, popb, 1, zorb, zorb, 1, irrepex)
                exInfor(1) = dirprd(exInfor(1), symry) ! irrepa
                exInfor(2) = dirprd(exInfor(2), symry) ! irrepi
                exInfor(5) = exInfor(5) + vrta(exInfor(1)) ! a
                exInfor(6) = exInfor(6) + popa(exInfor(2)) ! i
            end if

            !------- initialize eigenvector, i2
            i3 = i2 + v_len
            call dscal(v_len, 0.0d0, iCore(i2), 1)

            call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                        exInfor(5), exInfor(6), 1, 1, &
                        pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
            iCore(i2 + idx - 1) = 1.0d0

            call IO_stream(exE, iCore(i2), (v_len + 1)*(iroot - 1) + 1, &
                           v_len, f_init_davidson(irrepex), '21', 'old')

            call tmp%init(iCore, pop_a, pop_i, zorb, zorb, 1, 1, irrepex, T, i2)
            call printEXinfo1_uhf(tmp, 4)
        end do
        deallocate (lowest_e)
    end do

end subroutine p_pre_open_ns

! 2026.5.15
!> prepare initial guess for closed shell specified excitation
subroutine p_pre_close_specified_ns(pop_a, pop_i, Ea, Eb, &
                                    iCore, maxCore, ex_specified1, ex_mult)
    use mod_symm
    use mod_orbit
    use mod_operator
    use mod_common
    implicit none

    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    real*8, dimension(*), intent(inout) :: Ea !> orbital energy alpha
    real*8, dimension(*), intent(inout) :: Eb !> orbital energy beta
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    integer*8, dimension(5), intent(in) :: EX_specified1 !> (sym,irrepi,i,irrepa,a)
    integer*8, intent(in) :: ex_mult !> multiple of excited state

    integer*8 :: irrepex !> symmerty of excited state
    integer*8 :: spatial_sym !> spatial symmetry of excited state
    integer*8, dimension(:), allocatable :: lowest_e !> index of lowest OED1
    real*8 :: exE !> initial excitation energy
    integer*8 :: i0, i1, i2, nSize, idx
    integer*8 :: v_len !> length of the eigenvector (r0,R1,R2,...)
    integer*8, dimension(8) :: exInfor !> (irrepp,iirepq,irrepr,irreps,p,q,r,s)
    type(opt) :: tmp
    integer*8, external :: irpdso, idsymsoc

    i0 = 1

    v_len = irpdso(irrepex, pop_a, pop_i, 1)

    irrepex = EX_specified1(1) ! symmetry of state, only consider spatial symmetry
    exInfor(1) = EX_specified1(4) ! irrepa
    exInfor(2) = EX_specified1(2) ! irrepi
    exInfor(5) = EX_specified1(5) ! a
    exInfor(6) = EX_specified1(3) ! i

    if (irrepex /= dirprd(exInfor(1), exInfor(2))) then
        write (6, *) 'ex_specified1 error', exInfor(2), exInfor(1), irrepex
        stop
    end if

    lowest = .false.
    write (6, *) 'ATTENTATION: Specified excited state, during Davidson iterations not use the lowest eigenvalue'

    ! calculate alpha OED for irrepex, i0
    nSize = irpdso(irrepex, vrta, popa, 1)
    i1 = i0 + nSize
    call Eorb_diff1(iCore(i0), Ea, vrta, popa, irrepex)

    ! get position(idx) of alpha OED we need at OED1
    call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                exInfor(5), exInfor(6), 1, 1, &
                vrta, popa, 1, zorb, zorb, 1, irrepex)
    exE = iCore(i0 - 1 + idx)
    write (6, *) 'initial excitation energy:', exE

    !------- initialize eigenvector, i1
    i2 = i1 + v_len
    call dscal(v_len, 0.0d0, iCore(i1), 1)
    ! alpha->alpha
    call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                exInfor(5), exInfor(6), 1, 1, &
                pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
    iCore(i1 + idx - 1) = 1.0d0/sqrt(2.0d0)
    ! beta->beta
    exInfor(1) = dirprd(exInfor(1), symry) ! irrepa
    exInfor(2) = dirprd(exInfor(2), symry) ! irrepi
    exInfor(5) = vrta(exInfor(1)) + exInfor(5) ! a
    exInfor(6) = popa(exInfor(2)) + exInfor(6) ! i
    call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                exInfor(5), exInfor(6), 1, 1, &
                pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
    if (ex_mult == 1) then
        iCore(i1 + idx - 1) = 1.0d0/sqrt(2.0d0)
    else if (ex_mult == 3) then
        iCore(i1 + idx - 1) = -1.0d0/sqrt(2.0d0)
    else
        write (6, *) 'ERROR: ex_mult error', ex_mult
        stop
    end if

    !-------- write in initial eigenvalue and eigenvector
    call IO_stream(exE, iCore(i1), 1, v_len, f_init_davidson(irrepex), '21', 'replace')

    call tmp%init(iCore, pop_a, pop_i, zorb, zorb, 1, 1, irrepex, T, i1)
    call printEXinfo1_uhf(tmp, 4)

end subroutine p_pre_close_specified_ns

! 2026.5.15
!> prepare initial guess for closed shell excitation
subroutine p_pre_close_ns(pop_a, pop_i, Ea, Eb, &
                          iCore, maxCore, ex_specified1, ex_mult)
    use mod_symm
    use mod_orbit
    use mod_operator
    use mod_common
    implicit none

    integer*8, dimension(8), intent(in) :: pop_a, pop_i
    real*8, dimension(*), intent(inout) :: Ea !> orbital energy alpha
    real*8, dimension(*), intent(inout) :: Eb !> orbital energy beta
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    integer*8, dimension(5), intent(in) :: EX_specified1 !> (sym,irrepi,i,irrepa,a)
    integer*8, intent(in) :: ex_mult !> multiple of excited state

    integer*8 :: irrepex !> symmerty of excited state
    integer*8 :: spatial_sym !> spatial symmetry of excited state
    integer*8, dimension(:), allocatable :: lowest_e !> index of lowest OED1
    real*8 :: exE !> initial excitation energy
    integer*8 :: i0, i1, i2, nSize, idx
    integer*8 :: nroot, iroot
    integer*8 :: v_len !> length of the eigenvector (r0,R1,R2,...)
    integer*8, dimension(8) :: exInfor !> (irrepp,iirepq,irrepr,irreps,p,q,r,s)
    type(opt) :: tmp
    integer*8, external :: irpdso, idsymsoc

    i0 = 1

    do irrepex = 1, nirrep
        ! create file
        call IO_stream(0.0d0, iCore(i0), 1, 0, f_init_davidson(irrepex), '22', 'replace')
        write (6, '(A,I5)') '-------- symmetry of excitation state :', irrepex
        nroot = estate_sym(irrepex)
        if (nroot == 0) then
            cycle
        end if
        allocate (lowest_e(nroot))

        v_len = irpdso(irrepex, pop_a, pop_i, 1)

        !------- calculate alpha OED1 for irrepex, i0
        nSize = irpdso(irrepex, vrta, popa, 1)
        i1 = i0 + nSize
        call Eorb_diff1(iCore(i0), Ea, vrta, popa, irrepex)
        if (CVS) then
            call getpqrso_incore(iCore(i1), core_idx(CVS_orb), core_sym(CVS_orb), &
                                 zorb, zorb, 1, vrta, popa, 1, iCore(i0), irrepex)
            iCore(i0:i1 - 1) = 1.0d10
            call putpqrso_incore(iCore(i1), core_idx(CVS_orb), core_sym(CVS_orb), &
                                 zorb, zorb, 1, vrta, popa, 1, iCore(i0), irrepex)
        end if
        ! get position of lowest OED for irrepex
        call nmin(iCore(i0), nSize, nroot, lowest_e)

        do iroot = 1, nroot
            write (6, '(A,I2,I3,A)') '------ root', irrepex, iroot, ' ------'
            exE = iCore(i0 - 1 + lowest_e(iroot))
            write (6, '(A,F)') 'initial excitation energy:', exE
            !-------- get information of inituall guess energy
            call getidx(lowest_e(iroot), exInfor, vrta, popa, 1, zorb, zorb, 1, irrepex)

            !------- initialize eigenvector, i1
            i2 = i1 + v_len
            call dscal(v_len, 0.0d0, iCore(i1), 1)
            ! alpha->alpha
            call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                        exInfor(5), exInfor(6), 1, 1, &
                        pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
            iCore(i1 + idx - 1) = 1.0d0/sqrt(2.0d0)
            ! beta->beta
            exInfor(1) = dirprd(exInfor(1), symry) ! irrepa
            exInfor(2) = dirprd(exInfor(2), symry) ! irrepi
            exInfor(5) = vrta(exInfor(1)) + exInfor(5) ! a
            exInfor(6) = popa(exInfor(2)) + exInfor(6) ! i
            call putidx(idx, exInfor(1), exInfor(2), 1, 1, &
                        exInfor(5), exInfor(6), 1, 1, &
                        pop_a, pop_i, 1, zorb, zorb, 1, irrepex)
            if (ex_mult == 1) then
                iCore(i1 + idx - 1) = 1.0d0/sqrt(2.0d0)
            else if (ex_mult == 3) then
                iCore(i1 + idx - 1) = -1.0d0/sqrt(2.0d0)
            else
                write (6, *) 'ERROR: ex_mult error'
                stop
            end if

            !-------- write in initial eigenvalue and eigenvector
            idx = (1 + v_len)*(iroot - 1) + 1
            call IO_stream(exE, iCore(i1), idx, v_len, &
                           f_init_davidson(irrepex), '21', 'old')

            call tmp%init(iCore, pop_a, pop_i, zorb, zorb, 1, 1, irrepex, T, i1)
            call printEXinfo1_uhf(tmp, 4)
        end do
        deallocate (lowest_e)
    end do

end subroutine p_pre_close_ns
