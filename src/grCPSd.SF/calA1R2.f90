! write by lishun, 2024.10.21
! calculate Au2v^(1)*Tv^(n-1) term
! Au2v^(1)*Tv^(n-1) = <ab,ij|[V,R2^(n-1)]|0>, total order is n
!
subroutine calA1R2(iCore, C2, istart, order)
    use mod_variable
    use mod_operator
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    type(opt), intent(inout) :: C2 !> C^(n) constant doubles operator
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order ! n, the order of deltaT need calculate

    integer*8 :: idx
    integer*8 :: i0, i1, nSize, maxCore
    type(opt) :: tmpT2
    integer*8, external :: irpdso, idsymsoc

    i0 = istart
    call copy_opt(T2, tmpT2, T, iCore, i0)
    i0 = i0 + T2%len
    maxCore = nCore - i0 + 1

    !---- read in t(ab,ij)^(k)
    idx = (T1%len + T2%len)*(order - 2) + T1%len + 1
    call IO_stream(0.0d0, tmpT2%arr, idx, tmpT2%len, f_T, '13', 'old')

    !------------ <ab,ij|[V,R2^(k)]|0> = PijPab{V(kb,cj)*t(ac,ik)^(k)}
    !                    +0.5*V(kl,ij)*t(ab,kl)^(k)+0.5*V(ab,cd)*t(cd,ij)^(k)
    !---- PijPab{V(kb,cj)*t(ac,ik)^(1)}, o3v3
    call c2_r2Wmbej(C2, tmpT2, iabj, 1.0d0, iCore(i0), maxCore)

    !---- 0.5*V(kl,ij)*t(ab,kl)^(1), o4v2
    call c2_r2Wmnij(C2, tmpT2, ijkl, 1.0d0, iCore(i0), maxCore)

    !---- 0.5*V(ab,cd)*t(cd,ij)^(1), o2v4
    call c2_r2abcd(C2, tmpT2, abcd, 1.0d0, iCore(i0), maxCore)

end subroutine calA1R2
