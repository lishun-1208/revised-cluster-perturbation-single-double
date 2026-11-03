! write by lishun, 2024.11.5
! calculete eigenvector of CIS(Doo), and 1st order pertubation excitation energy
subroutine calEigenCISDoo(iCore, iStart, symm, iroot)
    use mod_symm
    use mod_orbit
    use mod_variable
    use mod_operator
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full work space
    integer*8, intent(in) :: iStart ! start pointer of free iCore
    integer*8, intent(in) :: symm ! the symmetry of root
    integer*8, intent(in) :: iroot

    real*8 :: exE, w1 ! excitation energy
    type(opt) :: L1, L2
    integer*8 :: idEmw ! -(Add^(0)-w)^(-1)
    integer*8 :: iter
    real*8 :: errdiis
    integer*8 :: i0, i1, iEnd, nSize
    integer*8 :: idx
    integer*8 :: i
    real*8 :: norm
    character(len=7) :: fileStatus ! for file, 'replace' or 'old'
    real*8 :: tt0, tt1
    integer*8, external :: idsymsoc
    real*8, external :: ddot

    if (iroot == 1) then
        fileStatus = 'replace'
    else
        fileStatus = 'old'
    end if

    i0 = iStart
    call R1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, T, i0)
    i0 = i0 + R1%len
    call R2%init(iCore, popv, popv, popo, popo, iop, iop, symm, T, i0)
    i0 = i0 + R2%len
    call L1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, T, i0)
    i0 = i0 + L1%len
    call L2%init(iCore, popv, popv, popo, popo, iop, iop, symm, T, i0)
    i0 = i0 + L2%len
    call C1%init(iCore, popv, popo, zorb, zorb, 1, 1, symm, T, i0)
    i0 = i0 + C1%len
    call C2%init(iCore, popv, popv, popo, popo, iop, iop, symm, T, i0)
    i0 = i0 + C2%len
    idEmw = i0
    i0 = idEmw + idsymsoc(symm, popv, popv, iop, popo, popo, iop)

    !-----------------------------------------
    !           CIS(Doo) R
    !-----------------------------------------
    idx = (1 + C1%len + C2%len)*(iroot - 1) + 1
    call IO_stream(exE, R1%arr, idx, C1%len + C2%len, f_R0(symm), '11', 'old')

    write (6, '(A,I2,A,I3,A,F)') 'CIS(Doo) excitation energy for symmerty', &
        symm, ' root', iroot, ":", exE
    norm = sqrt(R1%norm()**2 + R2%norm()**2)
    if (abs(norm - 1.0d0) > diis_max_err) then
        write (6, '(A,F)') 'WARN: CIS(Doo) R may not be normalized', abs(norm - 1.0d0)
        call dscal(nSize - 1, 1.0d0/norm, R1%arr, 1)
    end if

    !--- write EE-grCPS(D) (*,R^(0))
    call IO_stream(exE, R1%arr, idx, R1%len + R2%len, f_R(symm), '21', fileStatus)

    !-------------------------------------------------
    !         solve L1 equation
    !-------------------------------------------------
    !------------ calculate dEmw(ab,ij) = -(Add^(0)-w)^(-1)
    call getallso(iCore(idEmw), 1, f_OED2(symm)) ! get OED(ab,ij)=Add^(0)
    nSize = idsymsoc(symm, popv, popv, iop, popo, popo, iop)
    !$OMP PARALLEL DO
    do i = 1, nSize
        iCore(idEmw - 1 + i) = 1.0d0/(exE - iCore(idEmw - 1 + i))
    end do
    !$OMP END PARALLEL DO

    !------------------------------------------------
    ! solve L1 equation:
    !   Ls*[Ass-Asd*(Add-w)^(-1)*Ads-w] = 0
    ! use:
    !   D(a,i)*l(i,a) = Ls*[(Ass-e(a,i))-Asd*(Add-w)^(-1)*Ads-w]
    !------------------------------------------------
    call dcopy(R1%len, R1%arr, 1, L1%arr, 1)
    iter = 0
    do
        iter = iter + 1

        !------ DIIS
        ! normalize l(a,i) to r(a,i)
        norm = ddot(L1%len, L1%arr, 1, R1%arr, 1)
        call dscal(L1%len, 1.0d0/norm, L1%arr, 1)
        call dosodiis2(iter, L1%arr, L1%len, diis_dim, iCore(i0), nCore - i0 + 1, errdiis)
        ! normalize l(a,i) to r(a,i)
        norm = ddot(L1%len, L1%arr, 1, R1%arr, 1)
        call dscal(L1%len, 1.0d0/norm, L1%arr, 1)

        if ((abs(errdiis) < diis_max_err) .and. (iter > 1)) then
            write (6, '(A)') 'DIIS convergenced for CIS(Doo) L1'
            exit
        else if (iter > diis_cyc) then
            write (6, '(A,F)') 'WARN: DIIS not convergenc for CIS(Doo) L1', errdiis
            exit
        end if

        !------- -w*Ls
        call dcopy(L1%len, L1%arr, 1, C1%arr, 1)
        call dscal(L1%len, -exE, C1%arr, 1)

        !-------- Ls*Ass = (<0|[L1,H]|a,i>-e(a,i)*Ls)+<0|[L1,[V,MP2t]]|a,i>
        !                = Wibaj*l(j,b)-Gba*l(i,b)+Gij*l(j,a)
        ! Gba = Fad(b,a)
        ! Gij = Fli(i,j)
        ! Wibaj = Wkaci(ai,bj), sort as (ai,bj)
        !---- Wibaj*l(b,j)
        call c1d_l1Wibaj(C1, L1, Wkaci, 1.0d0, iCore(i0), nCore - i0 + 1)
        !---- -G(b,a)*l(i,b)
        call c1d_l1fab(C1, L1, Fad, 1.0d0, iCore(i0), nCore - i0 + 1)
        !---- G(i,j)*l(j,a)
        call c1d_l1fij(C1, L1, Fli, 1.0d0, iCore(i0), nCore - i0 + 1)

        !------ -Ls*Asd*(Add-w)^(-1)*Ads = Ls*Asd*dEmw*Ads
        call dscal(L2%len, 0.0d0, L2%arr, 1)
        !--- C(ab,ij) = Ls*Asd, i0
        call calL2d_VL1(iCore, symm, L1, L2, i0)

        !--- C = C*dEmw, i0
        call vdmul(L2%len, iCore(idEmw), L2%arr, iCore(i0))
        call dcopy(L2%len, iCore(i0), 1, L2%arr, 1)

        !--- D*Ads
        call calL1d_VL2(iCore, symm, C1, L2, i0, 1.0d0)

        !--------- l(a,i)=-c(a,i)/e(a,i)
        call denomT1(L1, C1, iCore(iEorb))
    end do

    !---------------------------------------
    ! calculate L2
    ! Ld = Ls*Asd*(w-Add)^(-1)=Ls*Asd*dEmw
    !---------------------------------------
    call dscal(L2%len, 0.0d0, L2%arr, 1)
    ! Ls*Asd
    call calL2d_VL1(iCore, symm, L1, L2, i0)
    ! [Ls*Asd]*dEmw
    call vdmul(L2%len, iCore(idEmw), L2%arr, iCore(i0))
    call dcopy(L2%len, iCore(i0), 1, L2%arr, 1)

    !------- normalize L to R
    norm = 1.0d0 + ddot(L2%len, L2%arr, 1, R2%arr, 1)
    call dscal(L1%len + L2%len, 1.0d0/norm, L1%arr, 1)

    !----------- test L
    call testEigenL(iCore, L1, L2, exE, symm, i0, diis_max_err)

    write (6, '(A,2I4)') 'CIS(Doo) left vector for state:', symm, iroot
    call printEXinfo1_uhf(L1, 10)
    call printEXinfo2_uhf(L2, 10)
    write (6, *) 'normalize:', ddot(L1%len + L2%len, L1%arr, 1, R1%arr, 1)
    write (6, *) 'R single fact:', R1%norm()**2
    write (6, *) 'L single fact:', L1%norm()**2/(L1%norm()**2 + L2%norm()**2)

    !--- write CIS(Doo) (e,L)
    call IO_stream(exE, L1%arr, idx, L1%len + L2%len, f_L0(symm), '21', fileStatus)

    !--------------------------------------------
    ! calculate 1st order perturbation excitation energy
    ! w^(1) = L*J^(1)*R, J^(1) = <d|[V,\tau_d]|0>
    !--------------------------------------------
    !------- J^(1)*R = <ab,ij|[V,R2]|0>
    call cpu_time(tt0)
    call calE_C(iCore, i0, symm, iroot, 1)
    write (6, '(A,F)') 'norm of energy C1^(1)', C1%norm()
    write (6, '(A,F)') 'norm of energy C2^(1):', C2%norm()
    !------ L*[J^(1)*R]
    w1 = ddot(C1%len + C2%len, L1%arr, 1, C1%arr, 1)
    call cpu_time(tt1)
    write (6, '(A,F)') 'cpu time for calculate energy:', tt1 - tt0
    write (6, '(A,I2,A,I3,A,F)') '1st order excitation energy for symmetry', &
        symm, ' root', iroot, ':', w1

    !--- write e^(1)
    idx = (1 + C1%len + C2%len)*(iroot - 1) + 1
    call IO_stream(w1, iCore(i0), idx, 0, f_R(symm), '22', 'old')
    !--- write C^(1)
    idx = (C1%len + C2%len)*(iroot - 1) + 1
    call IO_stream(0.0d0, C1%arr, idx, C1%len + C2%len, f_C(symm), '23', fileStatus)

end subroutine calEigenCISDoo
