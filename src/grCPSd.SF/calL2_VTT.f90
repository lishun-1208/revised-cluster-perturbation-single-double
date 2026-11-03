! write by lishun, 2024.10.21
! calculate L(ab,ij) += 0.5*<ab,ij|[[V,T^(k)],T^(k)]|0>, total order is 2k+1
! L(ab,ij) += -PijPab{V(kb,cj)*t(c,i)*t(a,k)}
!             +0.5*Pij{V(ab,cd)*t(c,i)*t(d,j)}
!             +0.5*Pab{V(kl,ij)*t(a,k)*t(b,l)}
!             +Pab{Gad*t(db,ij)}
!             -Pij{Gli*t(ab,lj)}
!             +PijPab{Wkaci*t(bc,jk)}
!             +Wklij*t(ab,kl)
!             -0.5*Pab{V(kb,cd)*t(cd,ij)*t(a,k)}
! Gli = V(kl,ci)*t(c,k)+0.5*<kl||cd>*t(cd,ki)
! Gad = V(ka,cd)*t(c,k)-0.5*<kl||cd>*t(ca,kl)
! Wkaci = -V(kl,ci)*t(a,l)+V(ka,cd)*t(d,i)+0.5*<lk||dc>*t(ad,il)
! Wklij = 0.5*Pij{V(kl,cj)*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
subroutine calL2_VTT(iCore, L2, istart, order)
    use mod_variable
    use mod_orbit
    use mod_operator
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    type(opt), intent(inout) :: L2 !> L2^(n)
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> k for T^(k)

    integer*8 :: i0, i1, i2, nSize, idx
    integer*8 :: maxCore
    type(opt) :: tmpT1, tmpT2, tmp
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart
    call copy_opt(T1, tmpT1, T, iCore, i0)
    i0 = i0 + T1%len
    call copy_opt(T2, tmpT2, T, iCore, i0)
    i0 = i0 + T2%len
    maxCore = nCore - i0 + 1

    !-------- read in T1^(k) and T2^(k)
    idx = (T1%len + T2%len)*(order - 1) + 1
    call IO_stream(0.0d0, tmpT1%arr, idx, T1%len + T2%len, f_T, '13', 'old')

    !-------------------------------------------------
    ! -PijPab{V(kb,cj)*t(c,i)*t(a,k)}
    !-------------------------------------------------
    call tmp%init(iCore, popo, popo, popo, popv, iop, 1, 1, T, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    call Wkbij_r1iabj(tmp, tmpT1, iabj, 1.0d0, iCore(i1), nCore - i1 + 1)
    call c2_r1ijka(L2, tmpT1, tmp, -1.0d0, iCore(i1), nCore - i1 + 1)

    !--------------------------------------------
    ! 0.5*Pij{V(ab,cd)*t(c,i)*t(d,j)}
    !--------------------------------------------
    call copy_opt(abci, tmp, vvvo_incore, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    call Wabdj_r1abcd(tmp, tmpT1, abcd, 1.0d0, iCore(i1), nCore - i1 + 1)
    call c2_r1abci(L2, tmpT1, tmp, 0.5d0, iCore(i1), nCore - i1 + 1)

    !------------------------------------
    ! 0.5*Pab{V(kl,ij)*t(a,k)*t(b,l)}
    !------------------------------------
    call copy_opt(ijka, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    call Wkbij_r1Uklij(tmp, tmpT1, ijkl, 1.0d0, iCore(i1), nCore - i1 + 1)
    call c2_r1ijka(L2, tmpT1, tmp, 0.5d0, iCore(i1), nCore - i1 + 1)

    !-----------------------------------
    ! -Pij{Gli*t(ab,lj)}
    ! Gli = V(kl,ci)*t(c,k)+0.5*<kl||cd>*t(cd,ki)
    !-----------------------------------
    call copy_opt(fij, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    ! Gli = V(kl,ci)*t(c,k)+0.5*<kl||cd>*t(cd,ki)
    call Fmi_r1ijka(tmp, tmpT1, ijka, 1.0d0, iCore(i1), nCore - i1 + 1)
    call Fmi_r2abij(tmp, tmpT2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    ! -Pij{Gli*t(ab,lj)}
    call c2_r2Fmi(L2, tmpT2, tmp, 1.0d0, iCore(i1), nCore - i1 + 1)

    !-----------------------------
    ! Pab{Gad*t(db,ij)}
    ! Gad = V(ka,cd)*t(c,k)-0.5*<kl||cd>*t(ca,kl)
    !-----------------------------
    call copy_opt(fab, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    ! Gad = V(ka,cd)*t(c,k)-0.5*<kl||cd>*t(ca,kl)
    call Fae_r1abci(tmp, tmpT1, abci, 1.0d0, iCore(i1), nCore - i1 + 1)
    call Fae_r2abij(tmp, tmpT2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    ! Pab{Gad*t(db,ij)}
    call c2_r2Fae(L2, tmpT2, tmp, 1.0d0, iCore(i1), nCore - i1 + 1)

    !-------------------------------
    ! PijPab{Wkbcj*t(ac,ik)}
    ! Wkaci = -V(kl,ci)*t(a,l)+V(ka,cd)*t(d,i)+0.5*<lk||dc>*t(ad,il)
    !-------------------------------
    call copy_opt(iabj, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    ! Wkaci = -V(kl,ci)*t(a,l)+V(ka,cd)*t(d,i)+0.5*<lk||dc>*t(ad,il)
    call Wmbej_r1ijka(tmp, tmpT1, ijka, 1.0d0, iCore(i1), nCore - i1 + 1)
    call Wmbej_r1abci(tmp, tmpT1, abci, 1.0d0, iCore(i1), nCore - i1 + 1)
    call Wmbej_r2abij(tmp, tmpT2, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    ! PijPab{Wkbcj*t(ac,ik)}
    call c2_r2Wmbej(L2, tmpT2, tmp, 1.0d0, iCore(i1), nCore - i1 + 1)

    !-----------------------------------------
    ! -0.5*Pab{V(kb,cd)*t(cd,ij)*t(a,k)}
    !-----------------------------------------
    call copy_opt(ijka, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    call Wkbij_r2abci(tmp, tmpT2, abci, 1.0d0, iCore(i1), nCore - i1 + 1)
    call c2_r1ijka(L2, tmpT1, tmp, 1.0d0, iCore(i1), nCore - i1 + 1)

    !-----------------------------
    ! Wklij*t(ab,kl)
    ! Wklij = 0.5*Pij{V(kl,cj)*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
    !-----------------------------
    call copy_opt(ijkl, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call tmp%arr_clear
    ! Wklij = 0.5*Pij{V(kl,cj)*t(c,i)}+0.25*<kl||cd>*t(cd,ij)
    call Wmnij_r1ijka(tmp, tmpT1, ijka, 0.5d0, iCore(i1), nCore - i1 + 1)
    call Wmnij_r2abij(tmp, tmpT2, abij, 0.5d0, iCore(i1), nCore - i1 + 1)
    ! Wklij*t(ab,kl)
    call c2_r2Wmnij(L2, tmpT2, tmp, 2.0d0, iCore(i1), nCore - i1 + 1)

end subroutine calL2_VTT
