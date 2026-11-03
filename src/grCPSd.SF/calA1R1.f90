! write by lishun, 2024.10.21
! calculate Au1v^(1)*Tv^(n-1) term
! Au1v^(2)*Tv^(n-2) = -<a,i|[V*MP2T,R1^(n-2)]|0>
! ATTENTATION: MP2t with order 1
!
subroutine calA2R1(iCore, C1, istart, order)
    use mod_variable
    use mod_operator
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    type(opt), intent(inout) :: C1 ! constant singles operator
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> n>=3, the order of deltaT need calculate

    type(opt) :: tmp
    integer*8 :: i
    integer*8 :: i0, i1, nSize, idx, maxCore
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart
    call copy_opt(T1, tmp, T, iCore, i0)
    i0 = i0 + tmp%len
    maxCore = nCore - i0 + 1

    !---- read in T(a,i)^(n-2)
    idx = (T1%len + T2%len)*(order - 3) + 1
    call IO_stream(0.0d0, tmp%arr, idx, T1%len, f_T, '13', 'old')

    !-------- -<a,i|[V*MP2T,R1^(k)|0>=-[(Wkaci-V(ka,ci))*T(c,k)^(k)+Fli*T(a,l)^(k)-Fad*T(d,i)^(k)]

    !----- -(Wkaci-V(ka,ci))*T(c,k)^(1)
    ! Wkaci-V(ka,ci), store in Wkaci
    call daxpy(Wkaci%len, -1.0d0, iabj%arr, 1, Wkaci%arr, 1)
    ! -Wkaci*T(c,k)^(1)
    call c1_r1iabj(C1, tmp, Wkaci, -1.0d0, iCore(i0), maxCore)
    ! return Wkaci
    call daxpy(Wkaci%len, 1.0d0, iabj%arr, 1, Wkaci%arr, 1)

    !--------- calculate -Fli*t(a,l)
    call c1_r1Fmi(C1, tmp, Fli, 1.0d0, iCore(i0), maxCore)
    !--------- calculate Fad*t(d,i)
    call c1_r1Fae(C1, tmp, Fad, 1.0d0, iCore(i0), maxCore)

end subroutine calA2R1
