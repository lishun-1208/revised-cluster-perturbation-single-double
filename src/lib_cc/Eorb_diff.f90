! write by LiShun, 2024.5.14
! calculate orbital energy difference
! e(a,i)=e(a)-e(i)
!-------------------------------
! modified at 2024.12.10
subroutine Eorb_diff1(OED1, Eorb, pop_a, pop_i, sym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: OED1 ! store OED
    real*8, dimension(*), intent(in) :: Eorb
    integer*8, dimension(8), intent(in) :: pop_a, pop_i ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of OED1(a,i)

    integer*8 :: i, a
    integer*8 :: irrepi, irrepa ! represnetation index
    integer*8 :: numa, numi ! number of orbitals beloning to one representation
    integer*8, dimension(8) :: ioff_a, ioff_i ! start index in Eorb of every representation
    integer*8 :: idxi, idxa ! start index in Eorb
    integer*8, dimension(3) :: idxT ! index of T belong to irrepi,i,a
    real*8 :: dtmp ! Dia=fii-faa

    integer*8, external :: isymoffso

    !------------- give ioff_i
    ioff_i(1) = 0
    do irrepi = 1, nirrep - 1
        ioff_i(irrepi + 1) = ioff_i(irrepi) + pop_i(irrepi)
    end do
    !------------- give ioff_a
    ioff_a(1) = ioff_i(nirrep) + pop_i(nirrep)
    do irrepa = 1, nirrep - 1 ! nirrep in mod_symm, is the number of representation
        ioff_a(irrepa + 1) = ioff_a(irrepa) + pop_a(irrepa)
    end do

    ! cycle for all representation of i to get representation of a, then calculate t
    do irrepi = 1, nirrep
        irrepa = dirprd(irrepi, sym)
        numa = pop_a(irrepa)
        numi = pop_i(irrepi)
        ! position of T belong to irrepi
        idxT(1) = isymoffso(irrepi, sym, pop_a, pop_i, 1) - 1
        ! for position in Eorb
        idxa = ioff_a(irrepa)
        idxi = ioff_i(irrepi)
        ! cycle for all i,a belong to irrepi and irrepa
        do i = 1, numi
            idxT(2) = idxT(1) + numa*(i - 1)
            !$OMP PARALLEL DO
            do a = 1, numa
                idxT(3) = idxT(2) + a
                ! D(a,i)=fii-faa
                dtmp = Eorb(idxa + a) - Eorb(idxi + i)
                ! t(a,i)=r(a,i)/D(a,i)
                OED1(idxT(3)) = dtmp
            end do
            !$OMP END PARALLEL DO
        end do
    end do

    return
end subroutine Eorb_diff1

! e(ab,ij)=e(a)+e(b)-e(i)-e(j)
!---------------------------------------
! modified at 2025.2.28
subroutine Eorb_diff2(OED2, Eorb, pop_a, pop_b, pop_i, pop_j, sym, iop)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: OED2 ! store orbital energy difference
    real*8, dimension(*), intent(in) :: Eorb
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of OED2(ab,ij)
    integer*8, intent(in) :: iop ! 0 for (a<b,i<j), 1 for (ab,ij), 2 for (a<=b,i<=j)

    integer*8 :: i, j, a, b
    integer*8 :: irrepi, irrepj, irrepa, irrepb
    integer*8 :: numi, numj, numa, numb
    integer*8 :: irrepab, irrepij
    integer*8 :: numab, numij
    integer*8, dimension(8) :: ioffv, ioffo ! start index in Eorb of every representation
    integer*8, dimension(8) :: iofft2_ij ! start index of T2 belong to specfic reprensentation of ij
    integer*8 :: idxi, idxj, idxa, idxb ! start index of Eorb
    integer*8, dimension(7) :: idxT ! pointer of OED2(ab,ij) belong to irrepij,irrepj,j,i,irrepb,b,a
    real*8 :: dtmp ! Dijab=fii+fjj-faa-fbb

    integer*8, external :: irpdso, isymoffso

    !------------- give ioffo
    ioffo(1) = 0
    do irrepi = 1, nirrep - 1
        ioffo(irrepi + 1) = ioffo(irrepi) + pop_i(irrepi)
    end do
    !------------- give ioffv
    ioffv(1) = ioffo(nirrep) + pop_i(nirrep)
    do irrepa = 1, nirrep - 1
        ioffv(irrepa + 1) = ioffv(irrepa) + pop_a(irrepa)
    end do
    !------------- give iofft2_ij
    iofft2_ij(1) = 1
    do irrepij = 1, nirrep - 1 ! cycle for all representation of ij to get representation of ab
        irrepab = dirprd(irrepij, sym)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        iofft2_ij(irrepij + 1) = iofft2_ij(irrepij) + numij*numab
    end do

    ! cycle for all representation of ij to get representation of ab
    do irrepij = 1, nirrep
        irrepab = dirprd(irrepij, sym)
        numab = irpdso(irrepab, pop_a, pop_a, iop)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        idxT(1) = iofft2_ij(irrepij) - 1
        if (iop == 1) then
            ! cycle for all representation of j to get representation of i
            do irrepj = 1, nirrep
                irrepi = dirprd(irrepij, irrepj)
                numi = pop_i(irrepi)
                numj = pop_j(irrepj)
                idxi = ioffo(irrepi)
                idxj = ioffo(irrepj)
                idxT(2) = idxT(1) + numab*(isymoffso(irrepj, irrepij, pop_i, pop_j, 1) - 1)
                ! cycle for j
                do j = 1, numj
                    idxT(3) = idxT(2) + numab*numi*(j - 1)
                    ! cycle for i
                    do i = 1, numi
                        idxT(4) = idxT(3) + numab*(i - 1)
                        ! cycle for irrepb
                        do irrepb = 1, nirrep
                            irrepa = dirprd(irrepab, irrepb)
                            numa = pop_a(irrepa)
                            numb = pop_b(irrepb)
                            idxa = ioffv(irrepa)
                            idxb = ioffv(irrepb)
                            idxT(5) = idxT(4) + isymoffso(irrepb, irrepab, pop_a, pop_b, 1) - 1
                            ! cycle for b
                            do b = 1, numb
                                idxT(6) = idxT(5) + numa*(b - 1)
                                ! cycle for a
                                !$OMP PARALLEL DO
                                do a = 1, numa
                                    idxT(7) = idxT(6) + a
                                    dtmp = Eorb(idxa + a) + Eorb(idxb + b) &
                                           - Eorb(idxi + i) - Eorb(idxj + j)
                                    OED2(idxT(7)) = dtmp
                                end do
                                !$OMP END PARALLEL DO
                            end do
                        end do
                    end do
                end do
            end do
        else if (iop == 0) then
            ! cycle for irrepj
            do irrepj = 1, nirrep
                irrepi = dirprd(irrepij, irrepj)
                if (irrepi > irrepj) then ! iop=0 need irrepi<=irrepj, the same to ab
                    cycle
                end if
                numj = pop_j(irrepj)
                idxi = ioffo(irrepi)
                idxj = ioffo(irrepj)
                idxT(2) = idxT(1) + numab*(isymoffso(irrepj, irrepij, pop_i, pop_j, 0) - 1)
                ! cycle for j
                do j = 1, numj
                    if (irrepi == irrepj) then ! for irrepi=irrepj, need i<j
                        numi = j - 1
                        idxT(3) = idxT(2) + numab*((j - 1)*(j - 2)/2)
                    else if (irrepi < irrepj) then
                        numi = pop_i(irrepi)
                        idxT(3) = idxT(2) + numab*numi*(j - 1)
                    end if
                    ! cycle for i
                    do i = 1, numi
                        idxT(4) = idxT(3) + numab*(i - 1)
                        ! cycle for irrepb
                        do irrepb = 1, nirrep
                            irrepa = dirprd(irrepab, irrepb)
                            if (irrepa > irrepb) then
                                cycle
                            end if
                            numb = pop_b(irrepb)
                            idxa = ioffv(irrepa)
                            idxb = ioffv(irrepb)
                            idxT(5) = idxT(4) + isymoffso(irrepb, irrepab, pop_a, pop_b, 0) - 1
                            ! cycle for b
                            do b = 1, numb
                                if (irrepa == irrepb) then
                                    numa = b - 1
                                    idxT(6) = idxT(5) + (b - 1)*(b - 2)/2
                                else if (irrepa < irrepb) then
                                    numa = pop_a(irrepa)
                                    idxT(6) = idxT(5) + numa*(b - 1)
                                end if
                                ! cycle for a
                                !$OMP PARALLEL DO
                                do a = 1, numa
                                    idxT(7) = idxT(6) + a
                                    dtmp = Eorb(idxa + a) + Eorb(idxb + b) &
                                           - Eorb(idxi + i) - Eorb(idxj + j)
                                    OED2(idxT(7)) = dtmp
                                end do
                                !$OMP END PARALLEL DO
                            end do
                        end do
                    end do
                end do
            end do
        else if (iop == 2) then
            write (6, *) 'Eorb_diff2 current not support iop=2'
            stop
        end if
    end do

    return
end subroutine Eorb_diff2

! write at 2025.5.26
!-------------------------------------
!> e(a,i)=e1(a)-e2(i)
subroutine Eorb_diff_ab(OED1, E1, pop_a1, pop_i1, E2, pop_a2, pop_i2, sym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: OED1
    real*8, dimension(*), intent(in) :: E1
    integer*8, dimension(8), intent(in) :: pop_a1, pop_i1
    real*8, dimension(*), intent(in) :: E2
    integer*8, dimension(8), intent(in) :: pop_a2, pop_i2
    integer*8, intent(in) :: sym ! symmerty of e(a,i)

    integer*8 :: i, a
    integer*8 :: irrepi, irrepa ! represnetation of occupy and virtual orbital
    integer*8 :: numa, numi ! number of orbitals beloning to a representation
    integer*8, dimension(8) :: ioff_a1, ioff_i2 ! start pointer of E1/E2 for a representation
    integer*8 :: idxi, idxa ! start pointer in Eorb
    integer*8, dimension(3) :: idx ! index of e(a,i) belong to irrepi,i,a
    real*8 :: tmp

    integer*8, external :: isymoffso

    !------------- give ioff_i2
    ioff_i2(1) = 0
    do irrepi = 1, nirrep - 1
        ioff_i2(irrepi + 1) = ioff_i2(irrepi) + pop_i2(irrepi)
    end do
    !------------- give ioff_a1
    i = 0
    do irrepi = 1, nirrep
        i = i + pop_i1(irrepi)
    end do
    ioff_a1(1) = i
    do irrepa = 1, nirrep - 1 ! nirrep in mod_symm, is the number of representation
        ioff_a1(irrepa + 1) = ioff_a1(irrepa) + pop_a1(irrepa)
    end do

    ! cycle for all representation of i
    do irrepi = 1, nirrep
        irrepa = dirprd(irrepi, sym)
        numa = pop_a1(irrepa)
        numi = pop_i2(irrepi)
        ! position of e(a,i) belong to irrepi
        idx(1) = isymoffso(irrepi, sym, pop_a1, pop_i2, 1) - 1
        ! position of E1/E2
        idxa = ioff_a1(irrepa)
        idxi = ioff_i2(irrepi)
        ! cycle for all i,a
        do i = 1, numi
            ! position of e(a,i) belong to an i
            idx(2) = idx(1) + numa*(i - 1)
            !$OMP PARALLEL DO
            do a = 1, numa
                ! position of e(a,i) belong to an a
                idx(3) = idx(2) + a
                ! e(a,i)=e1(a)-e2(i)
                tmp = E1(idxa + a) - E2(idxi + i)
                OED1(idx(3)) = tmp
            end do
            !$OMP END PARALLEL DO
        end do
    end do

end subroutine Eorb_diff_ab

! 2026.7.14
!> e(aB,iJ)=e(a)+e(B)-e(i)-e(J)
subroutine Eorb_diff2_ab(OED2, Ea, Eb, pop_a, pop_b, pop_i, pop_j, sym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: OED2 ! store orbital energy difference
    real*8, dimension(*), intent(in) :: Ea, Eb
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of OED2(ab,ij)

    integer*8 :: iop = 1
    integer*8 :: i, j, a, b
    integer*8 :: irrepi, irrepj, irrepa, irrepb
    integer*8 :: numi, numj, numa, numb
    integer*8 :: irrepab, irrepij
    integer*8 :: numab, numij
    integer*8, dimension(8) :: iofft2_ij ! start index of T2 belong to specfic reprensentation of ij
    integer*8 :: idxi, idxj, idxa, idxb ! start index of Eorb
    integer*8, dimension(7) :: idxT ! pointer of OED2(ab,ij) belong to irrepij,irrepj,j,i,irrepb,b,a
    real*8 :: dtmp

    integer*8, external :: irpdso, isymoffso

    !------------- give iofft2_ij
    iofft2_ij(1) = 0
    do irrepij = 1, nirrep - 1 
        irrepab = dirprd(irrepij, sym)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        iofft2_ij(irrepij + 1) = iofft2_ij(irrepij) + numij*numab
    end do

    do irrepij = 1, nirrep
        irrepab = dirprd(irrepij, sym)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        idxT(1) = iofft2_ij(irrepij)

        do irrepj = 1, nirrep
            irrepi = dirprd(irrepij, irrepj)
            numi = pop_i(irrepi)
            numj = pop_j(irrepj)
            idxi = sum(pop_i(1:irrepi - 1))
            idxj = sum(pop_j(1:irrepj - 1))
            idxT(2) = idxT(1) + numab*(isymoffso(irrepj, irrepij, pop_i, pop_j, iop) - 1)
            ! cycle for j
            do j = 1, numj
                idxT(3) = idxT(2) + numab*numi*(j - 1)
                ! cycle for i
                do i = 1, numi
                    idxT(4) = idxT(3) + numab*(i - 1)
                    ! cycle for irrepb
                    do irrepb = 1, nirrep
                        irrepa = dirprd(irrepab, irrepb)
                        numa = pop_a(irrepa)
                        numb = pop_b(irrepb)
                        idxa = noa + sum(pop_a(1:irrepa - 1))
                        idxb = nob + sum(pop_b(1:irrepb - 1))
                        idxT(5) = idxT(4) + isymoffso(irrepb, irrepab, pop_a, pop_b, iop) - 1
                        ! cycle for b
                        do b = 1, numb
                            idxT(6) = idxT(5) + numa*(b - 1)
                            ! cycle for a
                            !$OMP PARALLEL DO
                            do a = 1, numa
                                idxT(7) = idxT(6) + a
                                dtmp = Ea(idxa + a) + Eb(idxb + b) - &
                                       Ea(idxi + i) - Eb(idxj + j)
                                OED2(idxT(7)) = dtmp
                            end do
                            !$OMP END PARALLEL DO
                        end do
                    end do
                end do
            end do
        end do
    end do

    return
end subroutine Eorb_diff2_ab
