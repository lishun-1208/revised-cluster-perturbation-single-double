! 2026.10.27
!> calculate intermediate term of C += fact*J^(2)*R
subroutine calC_J2R(iCore, istart, pC1, pC2, pR1, pR2, fact)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC1, pC2
    type(opt), intent(in) :: pR1, pR2
    real*8, intent(in) :: fact

    integer*8 :: i0, maxCore

    i0 = istart
    maxCore = nCore - i0 + 1

    !--- c(ab,ij) -= W(ka,ji)*r(b,k)
    call c2_r1ijka(pC2, pR1, J2_Wovoo, fact, iCore(i0), maxCore)

    !--- c(a,i) += W(aj,ib)*r(b,j)
    call c1_r1iabj(pC1, pR1, J2_Wovvo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += W(ka,ci)*r(bc,jk)
    call c2_r2Wmbej(pC2, pR2, J2_Wovvo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += 0.5*V(mn,ij)*r(ab,mn)
    call c2_r2Wmnij(pC2, pR2, J2_Woooo, fact, iCore(i0), maxCore)

    !--- c(a,i) -= F(j,i)*r(a,j)
    call c1_r1Fmi(pC1, pR1, J2_Foo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) -= F(k,j)*r(ab,ik)
    call c2_r2Fmi(pC2, pR2, J2_Foo, fact, iCore(i0), maxCore)

    !--- c(a,i) += F(a,b)*r(b,i)
    call c1_r1Fae(pC1, pR1, J2_Fvv, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += F(a,c)*r(cb,ij)
    call c2_r2Fae(pC2, pR2, J2_Fvv, fact, iCore(i0), maxCore)

    !--- c(a,i) -= 0.5*W(jk,ic)*r(ca,kj)
    call c1_r2ijka(pC1, pR2, J2_Wooov, fact, iCore(i0), maxCore)

    !--- c(a,i) += F(j,b)*r(ab,ij)
    call c1_r2Fme(pC1, pR2, J2_Fov, fact, iCore(i0), maxCore)

end subroutine calC_J2R

! 2026.10.27
!> calculate intermediate term of C += fact*J^(3)*R
subroutine calC_J3R(iCore, istart, pC1, pC2, pR1, pR2, fact)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC1, pC2
    type(opt), intent(in) :: pR1, pR2
    real*8, intent(in) :: fact

    integer*8 :: i0, maxCore

    i0 = istart
    maxCore = nCore - i0 + 1

    !--- c(ab,ij) -= W(ka,ji)*r(b,k)
    call c2_r1ijka(pC2, pR1, J3_Wovoo, fact, iCore(i0), maxCore)

    !--- c(a,i) += W(aj,ib)*r(b,j)
    call c1_r1iabj(pC1, pR1, J3_Wovvo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += W(ka,ci)*r(bc,jk)
    call c2_r2Wmbej(pC2, pR2, J3_Wovvo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += 0.5*V(mn,ij)*r(ab,mn)
    call c2_r2Wmnij(pC2, pR2, J3_Woooo, fact, iCore(i0), maxCore)

    !--- c(a,i) -= F(j,i)*r(a,j)
    call c1_r1Fmi(pC1, pR1, J3_Foo, fact, iCore(i0), maxCore)

    !--- c(ab,ij) -= F(k,j)*r(ab,ik)
    call c2_r2Fmi(pC2, pR2, J3_Foo, fact, iCore(i0), maxCore)

    !--- c(a,i) += F(a,b)*r(b,i)
    call c1_r1Fae(pC1, pR1, J3_Fvv, fact, iCore(i0), maxCore)

    !--- c(ab,ij) += F(a,c)*r(cb,ij)
    call c2_r2Fae(pC2, pR2, J3_Fvv, fact, iCore(i0), maxCore)

    !--- c(a,i) -= 0.5*W(jk,ic)*r(ca,kj)
    call c1_r2ijka(pC1, pR2, J3_Wooov, fact, iCore(i0), maxCore)

    !--- c(a,i) += F(j,b)*r(ab,ij)
    call c1_r2Fme(pC1, pR2, J3_Fov, fact, iCore(i0), maxCore)

end subroutine calC_J3R