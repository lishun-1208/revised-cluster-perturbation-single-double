! 2026.4.16
!> t(a,i)=r(a,i)*D(a,i), D(a,i)=e(i)-e(a)=-e(a,i)
subroutine epsilonT1(t1, r1, Eorb)
    use mod_operator
    type(opt), intent(out) :: t1
    type(opt), intent(in) :: r1
    real*8, dimension(*), intent(out) :: Eorb

    call p_epsilonT1(t1%arr, r1%arr, Eorb, t1%popp, t1%popq, t1%sym)

end subroutine epsilonT1

! 2026.4.16
!> t(ab,ij)=r(ab,ij)*D(ab,ij), D(ab,ij)=e(i)+e(j)-e(a)-e(b)=-e(ab,ij)
subroutine epsilonT2(t2, r2, Eorb)
    use mod_operator
    type(opt), intent(out) :: t2
    type(opt), intent(in) :: r2
    real*8, dimension(*), intent(out) :: Eorb

    call p_epsilonT2(t2%arr, r2%arr, Eorb, t2%popp, t2%popq, t2%popr, t2%pops, t2%sym, t2%ioppq)

end subroutine epsilonT2

! write by LiShun, 2024.5.14
! denomT is used to calculate t=r*D=-r*e
! t(a,i)=r(a,i)*D(a,i), D(a,i)=e(i)-e(a)=-e(a,i)
! t(ab,ij)=r(ab,ij)*D(ab,ij), D(ab,ij)=e(i)+e(j)-e(a)-e(b)=-e(ab,ij)
! current only support sym=1
!--------------------------------------
! modified at 2024.11.12
!
subroutine p_epsilonT1(t1, r1, Eorb, pop_a, pop_i, sym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: t1
    real*8, dimension(*), intent(in) :: r1
    real*8, dimension(*), intent(in) :: Eorb
    integer*8, dimension(8), intent(in) :: pop_a, pop_i ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of t1(a,i)

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
                dtmp = Eorb(idxi + i) - Eorb(idxa + a)
                ! t(a,i)=r(a,i)*D(a,i)
                t1(idxT(3)) = r1(idxT(3))*dtmp
            end do
            !$OMP END PARALLEL DO
        end do
    end do

    return
end subroutine p_epsilonT1

! t2 and r2 need with the same iop
subroutine p_epsilonT2(t2, r2, Eorb, pop_a, pop_b, pop_i, pop_j, sym, iop)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: t2
    real*8, dimension(*), intent(in) :: r2
    real*8, dimension(*), intent(in) :: Eorb
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of t2(ab,ij)
    integer*8, intent(in) :: iop ! iop=iopab=iopij

    integer*8 :: i, j, a, b
    integer*8 :: irrepi, irrepj, irrepa, irrepb
    integer*8 :: numi, numj, numa, numb
    integer*8 :: irrepab, irrepij
    integer*8 :: numab, numij
    integer*8, dimension(8) :: ioffv, ioffo ! start index in Eorb of every representation
    integer*8, dimension(8) :: iofft2_ij ! start index of T2 belong to specfic reprensentation of ij
    integer*8 :: idxi, idxj, idxa, idxb ! start index of Eorb
    integer*8, dimension(7) :: idxT ! pointer of t2(ab,ij) belong to irrepij,irrepj,j,i,irrepb,b,a
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
    ! cycle for all representation of ij to get representation of ab
    do irrepij = 1, nirrep - 1
        irrepab = dirprd(irrepij, sym)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        iofft2_ij(irrepij + 1) = iofft2_ij(irrepij) + numij*numab
    end do

    if (iop == 1) then
        ! cycle for all representation of ij to get representation of ab
        do irrepij = 1, nirrep
            irrepab = dirprd(irrepij, sym)
            numab = irpdso(irrepab, pop_a, pop_a, iop)
            numij = irpdso(irrepij, pop_i, pop_j, iop)
            idxT(1) = iofft2_ij(irrepij) - 1
            ! cycle for all representation of j to get representation of i
            do irrepj = 1, nirrep
                irrepi = dirprd(irrepij, irrepj)
                numi = pop_i(irrepi)
                numj = pop_j(irrepj)
                idxi = ioffo(irrepi)
                idxj = ioffo(irrepj)
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
                            idxa = ioffv(irrepa)
                            idxb = ioffv(irrepb)
                            idxT(5) = idxT(4) + isymoffso(irrepb, irrepab, pop_a, pop_b, iop) - 1
                            ! cycle for b
                            do b = 1, numb
                                idxT(6) = idxT(5) + numa*(b - 1)
                                ! cycle for a
                                !$OMP PARALLEL DO
                                do a = 1, numa
                                    idxT(7) = idxT(6) + a
                                    dtmp = Eorb(idxi + i) + Eorb(idxj + j) - Eorb(idxa + a) - Eorb(idxb + b)
                                    t2(idxT(7)) = r2(idxT(7))*dtmp
                                end do
                                !$OMP END PARALLEL DO
                            end do
                        end do
                    end do
                end do
            end do
        end do
    else if (iop == 0) then
        ! cycle for irrepij
        do irrepij = 1, nirrep
            irrepab = dirprd(irrepij, sym)
            numab = irpdso(irrepab, pop_a, pop_a, iop)
            numij = irpdso(irrepij, pop_i, pop_j, iop)
            idxT(1) = iofft2_ij(irrepij) - 1
            ! cycle for irrepj
            do irrepj = 1, nirrep
                irrepi = dirprd(irrepij, irrepj)
                if (irrepi > irrepj) then ! iop=0 need irrepi<=irrepj, the same to ab
                    cycle
                end if
                numj = pop_j(irrepj)
                idxi = ioffo(irrepi)
                idxj = ioffo(irrepj)
                idxT(2) = idxT(1) + numab*(isymoffso(irrepj, irrepij, pop_i, pop_j, iop) - 1)
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
                            idxT(5) = idxT(4) + isymoffso(irrepb, irrepab, pop_a, pop_b, iop) - 1
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
                                    dtmp = Eorb(idxi + i) + Eorb(idxj + j) - Eorb(idxa + a) - Eorb(idxb + b)
                                    t2(idxT(7)) = r2(idxT(7))*dtmp
                                end do
                                !$OMP END PARALLEL DO
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if

end subroutine p_epsilonT2

! 2026.7.10
!> t(aB,iJ)=r(aB,iJ)*D(aB,iJ), D(aB,iJ)=e(i)+e(J)-e(a)-e(B)=-e(aB,iJ)
subroutine epsilonT2_ab(t2, r2, Ea, Eb)
    use mod_operator
    type(opt), intent(out) :: t2
    type(opt), intent(in) :: r2
    real*8, dimension(*), intent(out) :: Ea, Eb

    if ((t2%ioppq == 0) .or. (t2%ioprs == 0) .or. (r2%ioppq == 0) .or. (r2%ioprs == 0)) then
        write (6, '(A,4I2)') 'ERROR: iop error in epsilonT2_ab.', &
            t2%ioppq, t2%ioprs, r2%ioppq, r2%ioprs
        stop
    end if

    call p_epsilonT2_ab(t2%arr, r2%arr, Ea, Eb, t2%popp, t2%popq, t2%popr, t2%pops, t2%sym)

end subroutine epsilonT2_ab

subroutine p_epsilonT2_ab(t2, r2, Ea, Eb, pop_a, pop_b, pop_i, pop_j, sym)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: t2
    real*8, dimension(*), intent(in) :: r2
    real*8, dimension(*), intent(in) :: Ea, Eb
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: sym ! irrep of t2(ab,ij)

    integer*8 :: iop = 1
    integer*8 :: i, j, a, b
    integer*8 :: irrepi, irrepj, irrepa, irrepb
    integer*8 :: numi, numj, numa, numb
    integer*8 :: irrepab, irrepij
    integer*8 :: numab, numij
    integer*8, dimension(8) :: iofft2_ij ! start index of T2 belong to specfic reprensentation of ij
    integer*8 :: idxi, idxj, idxa, idxb ! start index of Eorb
    integer*8, dimension(7) :: idxT ! pointer of t2(ab,ij) belong to irrepij,irrepj,j,i,irrepb,b,a
    real*8 :: dtmp ! Dijab=fii+fjj-faa-fbb

    integer*8, external :: irpdso, isymoffso

    !------------- give iofft2_ij
    iofft2_ij(1) = 1
    ! cycle for all representation of ij to get representation of ab
    do irrepij = 1, nirrep - 1
        irrepab = dirprd(irrepij, sym)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        iofft2_ij(irrepij + 1) = iofft2_ij(irrepij) + numij*numab
    end do

    ! cycle for all representation of ij to get representation of ab
    do irrepij = 1, nirrep
        irrepab = dirprd(irrepij, sym)
        numab = irpdso(irrepab, pop_a, pop_b, iop)
        numij = irpdso(irrepij, pop_i, pop_j, iop)
        idxT(1) = iofft2_ij(irrepij) - 1
        ! cycle for all representation of j to get representation of i
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
                                dtmp = Ea(idxi + i) + Eb(idxj + j) - Ea(idxa + a) - Eb(idxb + b)
                                t2(idxT(7)) = r2(idxT(7))*dtmp
                            end do
                            !$OMP END PARALLEL DO
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine p_epsilonT2_ab
