! 2025.12.10
!> print orbital energy for not consider spin case
subroutine printEorb_ns(Eorb)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(in) :: Eorb

    integer*8 :: irrep
    integer*8 :: i, j, idx

    write (6, *) 'number of orbital:', noa + nob + nva + nvb
    write (6, *) '------------------ orbital energy -----------------'

    idx = 0
    write (6, *) 'occupied orbital'
    do irrep = 1, nirrep
        write (6, '(A,I5)') 'symmetry:', irrep
        do i = 1, popa(irrep)
            write (6, '(I5,I5,A,F)') idx + i, i, '   alpha', Eorb(idx + i)
        end do
        do i = 1, popb(dirprd(irrep, symry))
            j = popa(irrep) + i
            write (6, '(I5,I5,A,F)') idx + j, j, '   beta ', Eorb(idx + j)
        end do
        idx = idx + no(irrep)
    end do

    write (6, *) 'virtual orbital'
    do irrep = 1, nirrep
        write (6, '(A,I5)') 'symmetry:', irrep
        do i = 1, vrta(irrep)
            write (6, '(I5,I5,A,F)') idx + i, i, '   alpha', Eorb(idx + i)
        end do
        do i = 1, vrtb(dirprd(irrep, symry))
            j = vrta(irrep) + i
            write (6, '(I5,I5,A,F)') idx + j, j, '   beta ', Eorb(idx + j)
        end do
        idx = idx + nv(irrep)
    end do

    write (6, *) '---------------------------------------------------'

end subroutine printEorb_ns

! 2025.12.16
!> print orbital for use the given pop
subroutine printEorb(Eorb, nOcc, nVir, popo, popv)
    use mod_symm
    implicit none

    real*8, dimension(*), intent(in) :: Eorb
    integer*8, intent(in) :: nOcc, nVir
    integer*8, dimension(8), intent(in) :: popo, popv

    integer*8 :: irrep
    integer*8 :: i, idx

    write (6, '(A,I)') 'number of orbital:', nOcc + nVir

    write (6, *) '----------------- orbital energy -----------------'

    idx = 0
    write (6, *) 'occupied orbital'
    do irrep = 1, nirrep
        write (6, '(A,I5)') 'symmetry:', irrep
        do i = 1, popo(irrep)
            write (6, '(I5,I5,F)') idx + i, i, Eorb(idx + i)
        end do
        idx = idx + popo(irrep)
    end do

    write (6, *) 'virtual orbital'
    do irrep = 1, nirrep
        write (6, '(A,I5)') 'symmetry:', irrep
        do i = 1, popv(irrep)
            write (6, '(I5,I5,F)') idx + i, i, Eorb(idx + i)
        end do
        idx = idx + popv(irrep)
    end do

    write (6, *) '--------------------------------------------------'

end subroutine printEorb
