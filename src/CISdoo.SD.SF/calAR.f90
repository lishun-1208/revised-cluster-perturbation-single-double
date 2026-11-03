! write by lishun, 2025.3.11
! calculate AR for Davidson iteraton
!--------------------------------------
! for CIS(Doo), let C=AR, then:
! c(a,i) = Fac*r(c,i)-Fki*r(a,k)+Wkaci*r(c,k)
!         +0.5*<ak||dc>*r(dc,ik)-0.5*<lk||ic>*r(ac,lk)
! c(ab,ij) = Pij{<ab||cj>*r(c,i)}-Pab{<kb||ij>*r(a,k)}+e(ab,ij)*r(ab,ij)
! e(ab,ij) = e(a)+e(b)-e(i)-e(j)
subroutine calAR(iCore, maxCore, istart, iAR, iR, symR)
    use mod_variable
    use mod_orbit
    use mod_common
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    integer*8, intent(in) :: iAR, iR
    integer*8, intent(in) :: symR

    integer*8 :: iR1, iR2 ! position of R1 and R2 in R(also AR)
    integer*8 :: i0, i1, nSize, idx
    real*8 :: t0, t1
    integer*8, external :: irpdso, idsymsoc

    i0 = istart
    call C1%init(iCore, popv, popo, zorb, zorb, 1_8, 1_8, symR, idx=iAR)
    call C2%init(iCore, popv, popv, popo, popo, iop, iop, symR, idx=iAR + C1%len)
    call R1%init(iCore, popv, popo, zorb, zorb, 1_8, 1_8, symR, idx=iR)
    call R2%init(iCore, popv, popv, popo, popo, iop, iop, symR, idx=iR + R1%len)

    !--------------------------------------------
    ! AR(a,i) = Fac*r(c,i)
    !          -Fki*r(a,k)
    !          +Wkaci*r(c,k)
    !          +0.5*<ak||dc>*r(dc,ik)
    !          -0.5*<lk||ic>*r(ac,lk)
    !--------------------------------------------
    call C1%arr_clear

    !--------- Fac*r(c,i)
    call c1_r1Fae(C1, R1, Fac, 1.0d0, iCore(i0), nCore - i0 + 1)

    !--------- -Fki*r(a,k)
    call c1_r1Fmi(C1, R1, Fki, 1.0d0, iCore(i0), nCore - i0 + 1)

    !--------- Wkaci*r(c,k)
    call c1_r1iabj(C1, R1, Wkaci, 1.0d0, iCore(i0), nCore - i0 + 1)

    !--------- 0.5*<ak||dc>*r(dc,ik)
!    call cpu_time(t0)
    call c1_r2abci(C1, R2, abci, 1.0d0, iCore(i0), nCore - i0 + 1)
!    call cpu_time(t1)
!    write (6, '(A,F)') 'CPU time for c1_r2abci(o2v3):', t1 - t0

    !--------- -0.5*<lk||ic>*r(ac,lk)
    call c1_r2ijka(C1, R2, ijka, 1.0d0, iCore(i0), nCore - i0 + 1)

    !--------------------------------------------
    ! AR(ab,ij) = Pij{<ab||cj>*r(c,i)}
    !            -Pab{<kb||ij>*r(a,k)}
    !            +e(ab,ij)*r(ab,ij)
    !--------------------------------------------
    !--------- get e(ab,ij) from file_Diag, i0
    i1 = i0 + R1%len
    idx = 1 + R1%len
    call IO_stream(0.0d0, iCore(i0), idx, R2%len, f_Diag(symR), '13', 'old')
    !--------- e(ab,ij)*r(ab,ij)
    call vdmul(R2%len, iCore(i0), R2%arr, C2%arr)

    !--------- Pij{<ab||cj>*r(c,i)}
!    call cpu_time(t0)
    call c2_r1abci(C2, R1, abci, 1.0d0, iCore(i0), nCore - i0 + 1)
!    call cpu_time(t1)
!    write (6, '(A,F)') 'CPU time for c2_r1abci(o2v3):', t1 - t0

    !--------- -Pab{<kb||ij>*r(a,k)}
    call c2_r1ijka(C2, R1, ijka, 1.0d0, iCore(i0), nCore - i0 + 1)

end subroutine calAR
