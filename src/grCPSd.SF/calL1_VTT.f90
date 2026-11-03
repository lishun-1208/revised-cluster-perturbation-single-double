! write by lishun, 2024.10.21
! calculate L(a,i) += 0.5*<a,i|[[V,T^(k)],T^(k)]|0>, total order is 2k+1
! L(a,i) += Gli*t(a,l)-Gad*t(d,i)+Gld*t(da,li)
! Gli = -V(kl,ci)*t(c,k)-0.5*V(kl,cd)*t(cd,ki)
! Gad = -V(ka,cd)*t(c,k)+0.5*V(kl,cd)*t(ca,kl)
! Gld = V(kl,cd)*t(c,k), store as (d,l)
subroutine calL1_VTT(iCore, L1, istart, order)
    use mod_variable
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    type(opt), intent(inout) :: L1
    integer*8, intent(in) :: istart ! free iCore
    integer*8, intent(in) :: order !> k for T^(k)

    integer*8 :: i0, i1, i2, nSize, idx
    integer*8 :: maxCore
    type(opt) :: tmpT1, tmpT2
    type(opt) :: Gli, Gad, Gld
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    i0 = istart
    call copy_opt(T1, tmpT1, T, iCore, i0)
    i0 = i0 + T1%len
    call copy_opt(T2, tmpT2, T, iCore, i0)
    i0 = i0 + T2%len
    call Gli%init(iCore, popo, popo, zorb, zorb, 1_8, 1_8, 1_8, idx=i0)
    i0 = i0 + Gli%len
    call Gad%init(iCore, popv, popv, zorb, zorb, 1_8, 1_8, 1_8, idx=i0)
    i0 = i0 + Gad%len
    call Gld%init(iCore, popv, popo, zorb, zorb, 1_8, 1_8, 1_8, idx=i0)
    i0 = i0 + Gld%len
    maxCore = nCore - i0 + 1

    !-------- read in T1^(k) and T2^(k)
    idx = (T1%len + T2%len)*(order - 1) + 1
    call IO_stream(0.0d0, tmpT1%arr, idx, T1%len + T2%len, f_T, '13', 'old')

    !------ calculate Gli = -V(kl,ci)*t(c,k)-0.5*V(kl,cd)*t(cd,ki)
    call Gli%arr_clear
    ! -V(kl,ci)*t(c,k) is -Fmi_t1ijka
    call Fmi_r1ijka(Gli, tmpT1, ijka, -1.0d0, iCore(i0), maxCore)
    ! -0.5*V(kl,cd)*t(cd,ki) is -Fmi_tau2abij
    call Fmi_r2abij(Gli, tmpT2, abij, -1.0d0, iCore(i0), maxCore)

    !----- calculate Gad = -V(ka,cd)*t(c,k)+0.5*V(kl,cd)*t(ca,kl)
    call Gad%arr_clear
    ! -V(ka,cd)*t(c,k) is -Fae_t1abci
    call Fae_r1abci(Gad, tmpT1, abci, -1.0d0, iCore(i0), maxCore)
    ! 0.5*V(kl,cd)*t(ca,kl) is -Fae_tau2abij
    call Fae_r2abij(Gad, tmpT2, abij, -1.0d0, iCore(i0), maxCore)

    !---- calculate Gld = V(kl,cd)*t(c,k), is Fme_t1abij
    call Gld%arr_clear
    call Fme_r1abij(Gld, tmpT1, abij, 1.0d0, iCore(i0), maxCore)

    !------ calculate Gli*t(a,l), is -r1_t1Fmi
    call c1_r1Fmi(L1, tmpT1, Gli, -1.0d0, iCore(i0), maxCore)

    !----- calculate -Gad*t(d,i), is -r1_t1Fae
    call c1_r1Fae(L1, tmpT1, Gad, -1.0d0, iCore(i0), maxCore)

    !---- calculate Gld*t(da,li), is same to r1_t2Fme
    call c1_r2Fme(L1, tmpT2, Gld, 1.0d0, iCore(i0), maxCore)

end subroutine calL1_VTT
