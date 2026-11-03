! modified by LiShun, 2024.9.14
!> read variables in interface and calculate some variables need in module
!---------------------------
! 2025.12.15
subroutine readmodule
    use mod_symm
    use mod_orbit
    use mod_size
    use mod_iop
    use mod_interface
    implicit none

    integer*8 :: irrep
    logical :: exist

    integer*8 :: i, j
    integer*8, external :: idsymsoc, irpdso

    !--- mod_sym
    call read_node('symrx', symrx)
    call read_node('symry', symry)
    call read_node('symrz', symrz)
    call read_node('dirprd', dirprd, 8, 8)
    call read_node('nirrep', nirrep)

    write (6, '(A,3I2)') 'symR:', symRx, symRy, symRz
    write (6, '(A,I2)') 'nirrep:', nirrep
    write (6, '(A)') 'dirprd:'
    do i = 1, 8
        write (6, '(8I2)') (dirprd(i, j), j=1, 8)
    end do

    !--- mod_orbit
    call read_node('pop', pop, 8, 2)
    call read_node('vrt', vrt, 8, 2)
    popa(:) = pop(:, 1)
    popb(:) = pop(:, 2)
    vrta(:) = vrt(:, 1)
    vrtb(:) = vrt(:, 2)
    do irrep = 1, nirrep
        nv(irrep) = vrta(irrep) + vrtb(dirprd(irrep, symry))
        no(irrep) = popa(irrep) + popb(dirprd(irrep, symry))
    end do
    noa = sum(popa(1:nirrep))
    nob = sum(popb(1:nirrep))
    nva = sum(vrta(1:nirrep))
    nvb = sum(vrtb(1:nirrep))

    write (6, 100) "vrta(8) =", vrta(1:8)
    write (6, 100) "vrtb(8) =", vrtb(1:8)
    write (6, 100) "popa(8) =", popa(1:8)
    write (6, 100) "popb(8) =", popb(1:8)

    !--- mod_size
    call read_node('iprcln', iprcln)
    call read_node('iprcwd', iprcwd)
    call read_node('ichcsz', ichcsz)

    write (6, '(A,I)') 'iprcln', iprcln
    write (6, '(A,I)') 'iprcwd', iprcwd
    write (6, '(A,I)') 'ichcsz', ichcsz

    !--- mod_iop, only double precision
    prec = 8
    isd = 2
    call grep_node('tmpdir', exist)
    if (exist) then
        call read_node('tmpdir', tmpdir)
    end if

    ! the number of integral with different occupy, in libsrint
    nvvvv = max(idsymsoc(1, nv, nv, 1, nv, nv, 1), idsymsoc(symrz, nv, nv, 1, nv, nv, 1)) ! why symrz?
    write (6, 102) "nvvvv(ab,cd) =", nvvvv
    write (6, 103) "memory for nvvvv(MiB) =", nvvvv*8/(1024.d0**2)
    write (6, 103) "memory for nvvvv(GiB) =", nvvvv*8/(1024.d0**3)
    nvvvo = max(idsymsoc(1, nv, nv, 1, nv, no, 1), idsymsoc(symrz, nv, nv, 1, nv, no, 1))
    write (6, 102) "nvvvo(ab,ci) =", nvvvo
    write (6, 103) "memory for nvvvo(MiB) =", nvvvo*8/(1024.d0**2)
    write (6, 103) "memory for nvvvo(GiB) =", nvvvo*8/(1024.d0**3)
    nvvoo = max(idsymsoc(1, nv, nv, 1, no, no, 1), idsymsoc(symrz, nv, nv, 1, no, no, 1))
    write (6, 102) "nvvoo(ab,ij) =", nvvoo
    write (6, 103) "memory for nvvoo(MiB) =", nvvoo*8/(1024.d0**2)
    write (6, 103) "memory for nvvoo(GiB) =", nvvoo*8/(1024.d0**3)
    nvooo = max(idsymsoc(1, nv, no, 1, no, no, 1), idsymsoc(symrz, nv, no, 1, no, no, 1))
    write (6, 102) "nvooo(ai,jk) =", nvooo
    write (6, 103) "memory for nvooo(MiB) =", nvooo*8/(1024.d0**2)
    noooo = max(idsymsoc(1, no, no, 1, no, no, 1), idsymsoc(symrz, no, no, 1, no, no, 1))
    write (6, 102) "noooo(ij,kl) =", noooo
    write (6, 103) "memory for noooo(MiB) =", noooo*8/(1024.d0**2)

100 format(A, 8I5)
102 format(A, I12)
103 format(A, F15.6)
end subroutine readmodule
