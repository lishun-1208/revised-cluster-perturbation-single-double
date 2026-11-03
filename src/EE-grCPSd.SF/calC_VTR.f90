! 2026.10.27
!> calculate no intermediate term of C += fact*<u|[[V,T],R]|0>
subroutine calC_VTR(iCore, istart, pC1, pC2, pT1, pT2, pR1, pR2, symR, fact)
    use mod_variable
    use mod_common
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC1, pC2
    type(opt), intent(in) :: pT1, pT2
    type(opt), intent(in) :: pR1, pR2
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact

    integer*8 :: i0, i1, maxCore
    type(opt) :: tmp

    i0 = istart
    maxCore = nCore - i0 + 1

    !--- c(a,i) -= 0.5*V(kl,cd)*r(cd,ki)*t(a,l)
    call tmp%init(iCore, popo, popo, zorb, zorb, 1_8, 1_8, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Fmi_r2abij(tmp, pR2, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
    call c1_r1Fmi(pC1, pT1, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- W(ak,ij) = V(ak,ic)*r(c,j)+0.5*V(ak,cd)*r(cd,ij)
    call tmp%init(iCore, popo, popo, popo, popv, iop, 1, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Wkbij_r1iabj(tmp, pR1, iabj, -1.0d0, iCore(i1), maxCore - tmp%len)
    call Wkbij_r2abci(tmp, pR2, abci, 1.0d0, iCore(i1), maxCore - tmp%len)
    !--- c(ab,ij) -= W(ak,ij)*t(b,k)
    call c2_r1ijka(pC2, pT1, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- c(ab,ij) += V(ab,cd)*t(c,i)*r(d,j)
    if (.not. vvvo_incore) then
        write (6, *) 'ERROR: current not support iostat=2'
        stop
    end if
    call tmp%init(iCore, popv, popv, popv, popo, iop, 1, 1, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Wabdj_r1abcd(tmp, pT1, abcd, 1.0d0, iCore(i1), maxCore - tmp%len)
    call c2_r1abci(pC2, pR1, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- W(mn,ij) = 0.5*V(mn,cd)*r(cd,ij)+V(mn,ic)*r(c,j)
    call tmp%init(iCore, popo, popo, popo, popo, iop, iop, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Wmnij_r2abij(tmp, pR2, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
    call Wmnij_r1ijka(tmp, pR1, ijka, 1.0d0, iCore(i1), maxCore - tmp%len)
    !--- c(ab,ij) += 0.5*W(mn,ij)*t(ab,mn)
    call c2_r2Wmnij(pC2, pT2, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- c(ab,ij) += V(ak,dc)*r(d,i)*t(cb,kj)
    call tmp%init(iCore, popv, popo, popv, popo, 1, 1, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Wmbej_r1abci(tmp, pR1, abci, 1.0d0, iCore(i1), maxCore - tmp%len)
    call c2_r2Wmbej(pC2, pT2, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- Gli = V(lk,ic)*r(c,k)+0.5*V(lk,dc)*r(dc,ik)
    call tmp%init(iCore, popo, popo, zorb, zorb, 1, 1, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Fmi_r1ijka(tmp, pR1, ijka, 1.0d0, iCore(i1), maxCore - tmp%len)
    call Fmi_r2abij(tmp, pR2, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
    !--- c(ab,ij) -= Gli*t(ab,lj)
    call c2_r2Fmi(pC2, pT2, tmp, fact, iCore(i1), maxCore - tmp%len)

    !--- Gad = V(ak,dc)*r(c,k)-0.5*V(lk,dc)*r(ca,kl)
    call tmp%init(iCore, popv, popv, zorb, zorb, 1, 1, symR, idx=i0)
    call tmp%arr_clear
    i1 = i0 + tmp%len
    call Fae_r1abci(tmp, pR1, abci, 1.0d0, iCore(i1), maxCore - tmp%len)
    call Fae_r2abij(tmp, pR2, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
    !----- c(ab,ij) += Gad*t(db,ij)
    call c2_r2Fae(pC2, pT2, tmp, fact, iCore(i1), maxCore - tmp%len)

end subroutine calC_VTR
