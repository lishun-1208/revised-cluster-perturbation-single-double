! write by lishun, 2024.10.21
! modified at 2024.11.25
!------------------------------------
! Wkaci = V(ka,ci)+V(kl,cd)*MP2(da,li), sort as (ck,ai)
! Fli = -0.5*v(kl,cd)*MP2(cd,ki)
! Fad = 0.5*v(kl,cd)*MP2(ca,kl)
!> c(a,i) -= fact*<a,i|[[V,MP2t],R1]|0>
subroutine calC1_VMP2R1(iCore, pC1, pR1, istart, fact)
    use mod_variable
    implicit none

    real*8, dimension(*), target, intent(inout) :: iCore !> full
    type(opt), intent(inout) :: pC1
    type(opt), intent(in) :: pR1
    integer*8, intent(in) :: istart ! free iCore
    real*8, intent(in) :: fact

    integer*8 :: i0, maxCore

    i0 = istart
    maxCore = nCore - i0 + 1

    !------ C(a,i) -= (Wkaci-V(ka,ci))*T(c,k)^(k)+Fli*T(a,l)^(k)-Fad*T(d,i)^(k)
    !--- -(Wkaci-V(ka,ci))*T(c,k)^(1)
    ! Wkaci-V(ka,ci), store in Wkaci
    call daxpy(Wkaci%len, -1.0d0, iabj%arr, 1, Wkaci%arr, 1)
    ! -Wkaci*T(c,k)^(1)
    call c1_r1iabj(pC1, pR1, Wkaci, -fact, iCore(i0), maxCore)
    ! return Wkaci
    call daxpy(Wkaci%len, 1.0d0, iabj%arr, 1, Wkaci%arr, 1)

    !--- calculate -Fli*t(a,l)
    call c1_r1Fmi(pC1, pR1, Fli, fact, iCore(i0), maxCore)
    !--- calculate Fad*t(d,i)
    call c1_r1Fae(pC1, pR1, Fad, fact, iCore(i0), maxCore)

end subroutine calC1_VMP2R1
