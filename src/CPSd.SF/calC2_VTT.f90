! write by lishun, 2024.10.21
! calculate c(ab,ij) += fact*0.5*<ab,ij|[[V,T],T]|0>
! c(ab,ij) += -PijPab{<ak||ic>*t(c,j)*t(b,k)}
!             +0.5*Pij{<ab||cd>*t(c,i)*t(d,j)}
!             +0.5*Pab{<kl||ij>*t(a,k)*t(b,l)}
!             +Pab{Fad*t(db,ij)}
!             -Pij{Fli*t(ab,lj)}
!             +PijPab{Wkaci*t(cb,kj)}
!             +Wklij*t(ab,kl)
!             -0.5*Pab{<kb||cd>*t(cd,ij)*t(a,k)}
! Fad = <ka||cd>*t(c,k)-0.5*<kl||cd>*t(ca,kl)
! Fli = <kl||ci>*t(c,k)+0.5*<kl||cd>*t(cd,ki)
! Wkaci = -<kl||ci>*t(a,l)+<ka||cd>*t(d,i)+0.5*<lk||dc>*t(ad,il)
! Wklij = 0.5*Pij{<kl||cj>*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
!----------------------------------------
! modified at 2024.12.9
subroutine calC2_VTT(iCore, maxCore, istart, C2, T1, T2, symT, fact)
    use mod_variable
    use mod_orbit
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    type(opt), intent(inout) :: C2
    type(opt), intent(in) :: T1, T2
    integer*8, intent(in) :: symT
    real*8, intent(in) :: fact

    integer*8 :: i, i0, i1, i2, nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart
    if (debug) then
        write (6, *) '********** start C2_VTT **********'
        write (6, *) 'norm of T2:', T2%norm()
        write (6, *) 'norm of C2 0:', C2%norm()
    end if
    !---- -PijPab{<ak||ic>*t(c,j)*t(b,k)}
!    calL c2_r1r1iabj(C2, T1, iCore(iiabj), popv, popv, popo, popo, popv, popo, &
!                     iop, symT, iCore(i0), nCore - i0 + 1, fact)
!    if (debug) then
!        nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
!        write (6, *) 'norm of C2 1:', iop, dnrm2(nSize, C2%arr, 1)
!    end if

    !---- 0.5*Pij{<ab||cd>*t(c,i)*t(d,j)}
!    ! tau(cd,ij) = 0.5*Pij{t(c,i)*t(d,j)}, i0
!    nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
!    i1 = i0 + nSize
!    call dscal(nSize, 0.0d0, iCore(i0), 1)
!    call effectT2(iCore(i0), T1, T1, iCore(i1), 0.5d0, nCore - i1 + 1, &
!                  popv, popv, popo, popo, iop, symT, symT, 1, 1)
!    ! <ab||cd>*tau(cd,ij)
!    if (incore) then
!        call c2_r2abcd(C2, iCore(i0), iCore(iabcd), popv, popv, popo, popo, popv, popv, &
!                       iop, 1, fact*2.0d0)
!    else
!        call c2_r2abcd_h(C2, iCore(i0), file_abcd, popv, popv, popo, popo, popv, popv, &
!                         iop, 1, 2.0d0*fact, iCore(i1), nCore - i1 + 1)
!    end if
!    if (debug) then
!        nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
!        write (6, *) 'norm of C2 2:', iop, dnrm2(nSize, C2%arr, 1)
!    end if

    !---- 0.5*Pab{<kl||ij>*t(a,k)*t(b,l)}
!    call c2_r1r1ijkl(C2, T1, iCore(iijkl), popv, popv, popo, popo, popo, popo, &
!                     iop, symT, iCore(i0), nCore - i0 + 1, fact)
!    if (debug) then
!        nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
!        write (6, *) 'norm of C2 3:', iop, dnrm2(nSize, C2%arr, 1)
!    end if

    !------------------------------------------
    ! -Pij{Fli*t(ab,lj)}
    ! Fli = <kl||ci>*t(c,k)+0.5*<kl||cd>*t(cd,ki)
    !------------------------------------------
    call Fli%init(iCore, popo, popo, zorb, zorb, 1_8, 1_8, symT, T, idx=i0)
    call Fli%arr_clear
    i1 = i0 + Fli%len
    ! Fli = <kl||ci>*t(c,k)+0.5*<kl||cd>*t(cd,ki)
!    call Fmi_r1ijka(iCore(i0), T1, iCore(iijka), popo, popo, popv, popo, &
!                    symT, iop, 1.0d0, iCore(i1), nCore - i1 + 1)
!    if (debug) then
!        write (6, *) 'norm of Fli:', dnrm2(nSize, iCore(i0), 1)
!    end if

    call Fmi_r2abij(Fli, T2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of Fli:', Fli%norm()
    end if

    ! -Pij{Fli*t(ab,lj)}
    call c2_r2Fmi(C2, T2, Fli, fact, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of C2 4:', C2%norm()
    end if

    !--------------------------------------------
    ! Pab{Fad*t(db,ij)}
    ! Fad = <ka||cd>*t(c,k)-0.5*<kl||cd>*t(ca,kl)
    !--------------------------------------------
    call Fad%init(iCore, popv, popv, zorb, zorb, 1_8, 1_8, symT, T, idx=i0)
    call Fad%arr_clear
    i1 = i0 + Fad%len
    ! Fad = <ka||cd>*t(c,k)-0.5*<kl||cd>*t(ca,kl)
!    if (incore) then
!        call Fae_r1abci(iCore(i0), T1, iCore(iabci), popv, popv, popv, popo, &
!                        iop, symT, 1.0d0, iCore(i1), nCore - i1 + 1)
!    else
!        call Fae_r1abci_h(iCore(i0), T1, file_abci, popv, popv, popv, popo, &
!                          iop, symT, iCore(i1), nCore - i1 + 1, 1.0d0)
!    end if
!    if (debug) then
!        write (6, *) 'norm of Fad:', dnrm2(nSize, iCore(i0), 1)
!    end if

    call Fae_r2abij(Fad, T2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of Fad:', Fad%norm()
    end if

    ! Pab{Fad*t(db,ij)}
    call c2_r2Fae(C2, T2, Fad, fact, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of C2 5:', C2%norm()
    end if

    !----------------------------------------------
    ! PijPab{Wkaci*t(cb,kj)}
    ! Wkaci = -<kl||ci>*t(a,l)+<ka||cd>*t(d,i)+0.5*<lk||dc>*t(ad,il)
    !----------------------------------------------
    call Wkaci%init(iCore, popv, popo, popv, popo, 1_8, 1_8, symT, idx=i0)
    call Wkaci%arr_clear
    i1 = i0 + Wkaci%len
    ! Wkaci = -<kl||ci>*t(a,l)+<ka||cd>*t(d,i)+0.5*<lk||dc>*t(ad,il)
!    call Wmbej_r1ijka(iCore(i0), T1, iCore(iijka), popo, popv, popv, popo, popo, &
!                      iop, symT, iCore(i1), nCore - i1 + 1, 1, 1.0d0)
!    if (debug) then
!        write (6, *) 'norm of Wkaci:', dnrm2(nSize, iCore(i0), 1)
!    end if
!
!    if (incore) then
!        call Wmbej_r1abci(iCore(i0), T1, iCore(iabci), popo, popv, popv, popo, popv, &
!                          symT, iop, 1.0d0, iCore(i1), nCore - i1 + 1)
!    else
!        call Wmbej_r1abci_h(iCore(i0), T1, file_abci, popo, popv, popv, popo, popv, &
!                            iop, symT, iCore(i1), nCore - i1 + 1, 1, 1.0d0)
!    end if
!    if (debug) then
!        write (6, *) 'norm of Wkaci:', dnrm2(nSize, iCore(i0), 1)
!    end if

    call Wmbej_r2abij(Wkaci, T2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of Wkaci:', dnrm2(nSize, iCore(i0), 1)
    end if

    ! PijPab{Wkaci*t(cb,kj)}
    call c2_r2Wmbej(C2, T2, Wkaci, fact, iCore(i1), nCore - i1 + 1)
    if (debug) then
        nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
        write (6, *) 'norm of C2 6:', iop, dnrm2(nSize, C2%arr, 1)
        call cal_spin_norm(C2%arr, 'vvoo', 0, 0, 1, 'ABAB', 'C2 6', iCore(istart), nCore - istart + 1)
    end if

    !---- -0.5*Pab{<kb||cd>*t(cd,ij)*t(a,k)}
!    if (incore) then
!        call c2_r2r1abci(C2, T2, T1, iCore(iabci), popv, popv, popo, popo, popv, popv, popo, &
!                         symT, symT, iop, fact, iCore(i0), nCore - i0 + 1)
!    else
!        call c2_r2r1abci_h(C2, T2, T1, file_abci, popv, popv, popo, popo, popv, popv, popo, &
!                           iop, symT, symT, fact, iCore(i0), nCore - i0 + 1, 1)
!    end if
!    if (debug) then
!        nSize = idsymsoc(1, popv, popv, iop, popo, popo, iop)
!        write (6, *) 'norm of C2 7:', iop, dnrm2(nSize, C2%arr, 1)
!    end if

    !-----------------------------------------------
    ! Wklij*t(ab,kl)
    ! Wklij = 0.5*Pij{<kl||cj>*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
    !-----------------------------------------------
    call Wklij%init(iCore, popo, popo, popo, popo, iop, iop, symT, T, idx=i0)
    call Wklij%arr_clear
    i1 = i0 + Wklij%len
    ! Wklij = 0.5*Pij{<kl||cj>*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
!    call Wmnij_r1ijka(iCore(i0), T1, iCore(iijka), popo, popo, popo, popo, popv, &
!                      iop, symT, iCore(i1), nCore - i1 + 1, 0.5d0, 1)
!    if (debug) then
!        write (6, *) 'norm of Wklij:', dnrm2(nSize, iCore(i0), 1)
!    end if

    call Wmnij_r2abij(Wklij, T2, abij, 0.5d0, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of Wklij:', Wklij%norm()
    end if

    ! Wklij*t(ab,kl)
    call c2_r2Wmnij(C2, T2, Wklij, 2.0d0*fact, iCore(i1), nCore - i1 + 1)
    if (debug) then
        write (6, *) 'norm of C2 8:', C2%norm()
        write (6, *) '*********** end C2_VTT ***********'
    end if

end subroutine calC2_VTT
