! write by lishun. 2024.10.15
! calculate Au1v^(0)*Tv += Wkaci*t(c,k)+Fli*t(a,l)-Fad*t(d,i)
!                        -0.5*V(kl,ci)*t(ca,kl)+0.5*V(ka,cd)*t(cd,ki)
! attention, not contain the OED term
!
subroutine calA0R1(iCore, istart)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore
    integer*8, intent(in) :: istart ! start pointer of free iCore

    integer*8 :: i0, nSize, maxCore
    integer*8, external :: irpdso
    real*8, external :: dnrm2

    i0 = istart
    maxCore = nCore - i0 + 1
    !---------- calculate Wkaci*t(c,k)
    call c1_r1iabj(R1, T1, Wkaci, 1.0d0, iCore(i0), maxCore)
    !--------- calculate Fli*t(a,l)
    call c1_r1Fmi(R1, T1, Fli, -1.0d0, iCore(i0), maxCore)
    !--------- calculate -Fad*t(d,i)
    call c1_r1Fae(R1, T1, Fad, -1.0d0, iCore(i0), maxCore)
    !-------- calculate -0.5*V(kl,ci)*t(ca,kl)
    call c1_r2ijka(R1, T2, ijka, 1.0d0, iCore(i0), maxCore)
    !------- calculate +0.5*V(ka,cd)*t(cd,ki)
    call c1_r2abci(R1, T2, abci, 1.0d0, iCore(i0), maxCore)

end subroutine calA0R1
