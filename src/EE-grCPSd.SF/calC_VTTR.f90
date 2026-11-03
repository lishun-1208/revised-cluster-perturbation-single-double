! 2026.10.27
!> calculate no intermediate term of C += fact*<u|[[[V,T],T],R]|0>
subroutine calC_VTTR(iCore, istart, pC2, pT1, pT2, pR1, pR2, symR, fact)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC2
    type(opt), intent(in) :: pT1, pT2
    type(opt), intent(in) :: pR1, pR2
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact

    integer*8 :: i0, i1, i2, i3, maxCore
    type(opt) :: tmp1, tmp2

    i0 = istart

    !--- Gld = V(kl,cd)*r(c,k)
    call tmp1%init(iCore, popv, popo, zorb, zorb, 1, 1, symR, idx=i0)
    call tmp1%arr_clear
    i1 = i0 + tmp1%len
    call Fme_r1abij(tmp1, pR1, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    !--- F(a,d) -= Gld*t(a,l)
    call tmp2%init(iCore, popv, popv, zorb, zorb, 1, 1, symR, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Fae_r1fia(tmp2, pT1, tmp1, 1.0d0, iCore(i2), nCore - i2 + 1)
    !--- c(ab,ij) += F(a,d)*t(db,ij)
    call c2_r2Fae(pC2, pT2, tmp2, fact, iCore(i2), nCore - i2 + 1)

    !--- F(l,i) = Gld*t(d,i)
    call tmp2%init(iCore, popo, popo, zorb, zorb, 1, 1, symR, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Fmi_r1fia(tmp2, pT1, tmp1, 1.0d0, iCore(i2), nCore - i2 + 1)
    !--- F(l,d) = V(kl,cd)*t(c,k)
    call tmp1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, idx=i2)
    call tmp1%arr_clear
    i3 = i2 + tmp1%len
    call Fme_r1abij(tmp1, pT1, abij, 1.0d0, iCore(i3), nCore - i3 + 1)
    !--- F(l,i) += Fld*r(d,i)
    call Fmi_r1fia(tmp2, pR1, tmp1, 1.0d0, iCore(i3), nCore - i3 + 1)
    !--- c(ab,ij) -= F(l,i)*t(ab,lj)
    call c2_r2Fmi(pC2, pT2, tmp2, fact, iCore(i2), nCore - i2 + 1)

    !--- c(ab,ij) += 0.5*V(kl,cd)*r(c,i)*t(d,j)*t(ab,kl)
    ! W(kl,id) = V(kl,cd)*r(c,i)
    call tmp1%init(iCore, popo, popo, popo, popv, iop, 1, symR, idx=i0)
    call tmp1%arr_clear
    i1 = i0 + tmp1%len
    call Wlkjd_r1abij(tmp1, pR1, abij, 1.0d0, iCore(i1), nCore - i1 + 1)
    ! W(kl,ij) = W(kl,id)*t(d,j)
    call tmp2%init(iCore, popo, popo, popo, popo, iop, iop, symR, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Wmnij_r1ijka(tmp2, pT1, tmp1, 1.0d0, iCore(i2), nCore - i2 + 1)
    ! c(ab,ij) += 0.5*W(kl,cd)*t(cd,ij)
    call c2_r2Wmnij(pC2, pT2, tmp2, fact, iCore(i2), nCore - i2 + 1)

    !--- W(kb,ij) = -0.5*W(kl,ij)*t(b,l)+W(kb,cj)*r(c,i)+W(kl,id)*t(db,lj)
    call tmp1%init(iCore, popo, popo, popo, popv, iop, 1, symR, idx=i0)
    call tmp1%arr_clear
    i1 = i0 + tmp1%len
    ! W(kl,ij) = 0.5*V(kl,cd)*r(cd,ij)+V(lk,jc)*r(c,i)
    call tmp2%init(iCore, popo, popo, popo, popo, iop, iop, symR, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Wmnij_r2abij(tmp2, pR2, abij, 1.0d0, iCore(i2), nCore - i2 + 1)
    call Wmnij_r1ijka(tmp2, pR1, ijka, 1.0d0, iCore(i2), nCore - i2 + 1)
    call Wkbij_r1Uklij(tmp1, pT1, tmp2, 0.5d0, iCore(i2), nCore - i2 + 1)
    ! W(kb,cj) = V(bk,dc)*t(d,j)
    call tmp2%init(iCore, popv, popo, popv, popo, 1, 1, 1, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Wmbej_r1abci(tmp2, pT1, abci, 1.0d0, iCore(i2), nCore - i2 + 1)
    call Wkbij_r1iabj(tmp1, pR1, tmp2, -1.0d0, iCore(i2), nCore - i2 + 1)
    ! W(kl,id) = V(kl,cd)*r(c,i)
    call tmp2%init(iCore, popo, popo, popo, popv, iop, 1, symR, idx=i1)
    call tmp2%arr_clear
    i2 = i1 + tmp2%len
    call Wlkjd_r1abij(tmp2, pR1, abij, 1.0d0, iCore(i2), nCore - i2 + 1)
    call Wkbij_r2ijka(tmp1, pT2, tmp2, 1.0d0, iCore(i2), nCore - i2 + 1)
    !-- c(ab,ij) -= W(kb,ij)*t(a,k)
    call c2_r1ijka(pC2, pT1, tmp1, fact, iCore(i1), nCore - i1 + 1)

end subroutine calC_VTTR
