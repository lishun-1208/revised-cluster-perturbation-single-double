! write by lishun
! resort Eorb from (oa,va),(ob,vb) to (o(oa,ob),v(va,vb))
! modified at 2024.11.13
subroutine resortEorb(Eorb, Ea, Eb)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: Eorb ! without spin
    real*8, dimension(*), intent(in) :: Ea, Eb ! with spin

    integer*8 :: irrep ! index of representation
    integer*8, dimension(8) :: ioffoa, ioffva, ioffob, ioffvb ! start index in Ea&Eb for different representation
    integer*8, dimension(8) :: ioffo, ioffv ! start index of Eorb belong to different representation
    integer*8 :: irrepa, irrepb ! representation for a&b
    integer*8 :: numa, numb ! number of orbital for representation of a&b
    integer*8 :: idxa, idxb, idx ! index used in Ea,Eb,Eorb

    !--------- give ioffoa, ioffob, ioffo
    ioffoa(1) = 1
    ioffob(1) = 1
    ioffo(1) = 1
    do irrep = 1, nirrep - 1
        ioffoa(irrep + 1) = ioffoa(irrep) + popa(irrep)
        ioffob(irrep + 1) = ioffob(irrep) + popb(irrep)
        ioffo(irrep + 1) = ioffo(irrep) + no(irrep)
    end do
    !-------- give ioffva, ioffvb, ioffv
    ioffva(1) = ioffoa(nirrep) + popa(nirrep)
    ioffvb(1) = ioffob(nirrep) + popb(nirrep)
    ioffv(1) = ioffo(nirrep) + no(nirrep)
    do irrep = 1, nirrep - 1
        ioffva(irrep + 1) = ioffva(irrep) + vrta(irrep)
        ioffvb(irrep + 1) = ioffvb(irrep) + vrtb(irrep)
        ioffv(irrep + 1) = ioffv(irrep) + nv(irrep)
    end do

    !----------- cycle for all representation
    do irrep = 1, nirrep
        irrepa = irrep
        ! why sympy
        irrepb = dirprd(irrepa, symry)

        !---------- for occ
        numa = popa(irrepa)
        numb = popb(irrepb)
        ! for a
        idxa = ioffoa(irrepa)
        idx = ioffo(irrep)
        call dcopy(numa, Ea(idxa), 1, Eorb(idx), 1)
        ! for b
        idxb = ioffob(irrepb)
        idx = ioffo(irrep) + numa
        call dcopy(numb, Eb(idxb), 1, Eorb(idx), 1)

        !------------- for vrt
        numa = vrta(irrepa)
        numb = vrtb(irrepb)
        ! for a
        idxa = ioffva(irrepa)
        idx = ioffv(irrep)
        call dcopy(numa, Ea(idxa), 1, Eorb(idx), 1)
        ! for b
        idxb = ioffvb(irrepb)
        idx = ioffv(irrep) + numa
        call dcopy(numb, Eb(idxb), 1, Eorb(idx), 1)
    end do

end subroutine resortEorb
