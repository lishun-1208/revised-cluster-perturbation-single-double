! write by lishun, 2024.11.25
! test CIS(Doo) left eigenvector
subroutine testEigenL(iCore, pL1, pL2, exE, symm, istart, maxError)
    use mod_common, only: iop
    use mod_variable
    use mod_orbit
    use mod_symm
    use mod_operator
    implicit none

    real*8, dimension(*), target, intent(inout) :: iCore ! full work space
    type(opt), intent(in) :: pL1, pL2
    real*8, intent(in) :: exE ! CIS(Doo) excitation energy
    integer*8, intent(in) :: symm
    integer*8, intent(in) :: istart
    real*8, intent(in) :: maxError

    type(opt) :: pC1, pC2
    real*8 :: norm, error
    integer*8 ::i0, i1, nSize, i
    integer*8, external :: irpdso, idsymsoc, idamax
    real*8, external :: dnrm2

    i0 = istart
    call pC1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, T, i0)
    i0 = i0 + pC1%len
    call pC2%init(iCore, popv, popv, popo, popo, iop, iop, symm, T, i0)
    i0 = i0 + pC2%len

    !--------------- Ls*Ass+Ld*Ads-w*Ls=0
    !-------- -w*Ls
    call dcopy(pL1%len, pL1%arr, 1, pC1%arr, 1)
    call dscal(pL1%len, -exE, pC1%arr, 1)

    !-------- Ls*Ass = <0|[L1,H]|a,i>+<0|[L1,[V,MP2t]]|a,i>
    !                = Wibaj*l(j,b)-Gba*l(i,b)+Gij*l(j,a)+e(a,i)*Ls
    ! Gba = Fad(b,a)
    ! Gij = Fli(i,j)
    ! Wibaj = Wkaci(ai,bj), sort as (ai,bj)
    !---- Wibaj*l(b,j)
    call c1d_l1Wibaj(pC1, pL1, Wkaci, 1.0d0, iCore(i0), nCore - i0 + 1)
    !---- -G(b,a)*l(i,b)
    call c1d_l1fab(pC1, pL1, Fad, 1.0d0, iCore(i0), nCore - i0 + 1)
    !---- G(i,j)*l(j,a)
    call c1d_l1fij(pC1, pL1, Fli, 1.0d0, iCore(i0), nCore - i0 + 1)
    !---------- e(a,i)*l(a,i)
    call getallso(iCore(i0), 1, f_OED1(symm))
    !$OMP PARALLEL DO
    do i = 1, pL1%len
        pC1%arr(i) = pC1%arr(i) + iCore(i0 - 1 + i)*pL1%arr(i)
    end do
    !$OMP END PARALLEL DO

    !--------- Ld*Ads
    call calL1d_VL2(iCore, symm, pC1, pL2, i0, 1.0d0)

    error = pC1%arr(idamax(pC1%len, pC1%arr, 1))
    write (6, *) 'max L1 error:', error
    if (abs(error) > maxError) then
        write (6, *) 'WARN: L1 may not be eigenvector'
    else
        write (6, *) "L1 pass test"
    end if

    !--------------- Ls*Asd+Ld*Add-w*Ld=0
    !------ -w*Ld
    call dcopy(pC2%len, pL2%arr, 1, pC2%arr, 1)
    call dscal(pC2%len, -exE, pC2%arr, 1)

    !------ Ls*Asd
    call calL2d_VL1(iCore, symm, pL1, pC2, i0)

    !------ Ld*Add
    call getallso(iCore(i0), 1, f_OED2(symm))
    !$OMP PARALLEL DO
    do i = 1, pL2%len
        pC2%arr(i) = pC2%arr(i) + iCore(i0 - 1 + i)*pL2%arr(i)
    end do
    !$OMP END PARALLEL DO

    error = pC2%arr(idamax(pC2%len, pC2%arr, 1))
    write (6, *) 'max L2 error:', error
    if (abs(error) > maxError) then
        write (6, *) 'WARN: L2 may not be eigenvector'
    else
        write (6, *) "L2 pass test"
    end if

end subroutine testEigenL
