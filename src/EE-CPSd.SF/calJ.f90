! calculate intermediates Foo, Fvv, Fov, Woooo, Wovoo, Wooov, Wovvo
! not calculate Wvvvo, Wvovv, Wvvvv and higher particle rank intermediates
!> calculate J^(2) = <u|[[V,T^(1)],tau_v]|0> intermediates
! T1^(1)=0, R2^(0)=0, R1^(1)=0 cause:
!   for CPS(D-3), only need intermediates can contract with R1
subroutine calJ2(iCore, istart)
    use mod_variable
    use mod_orbit, only: zorb
    use mod_common, only: iop
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart

    type(opt) :: T1, T2
    integer*8 :: i0, maxCore

    i0 = istart
    call T1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, idx=i0) ! 0
    i0 = i0 + T1%len
    call T2%init(iCore, popv, popv, popo, popo, iop, iop, 1, idx=i0)
    i0 = i0 + T2%len
    maxCore = nCore - i0 + 1

    !--- read in T^(1)
    call IO_stream(0.0d0, T1%arr, 1, T1%len + T2%len, f_T, '13', 'old')

    !----------------------------
    ! F(j,i) = V(jk,ic)*t(c,k)+0.5*V(jk,dc)*t(dc,ik)
    !----------------------------
    call J2_Foo%arr_clear
!    call Fmi_r1ijka(J2_Foo, T1, ijka, 1.0d0, iCore(i0), maxCore)
    call Fmi_r2abij(J2_Foo, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! F(b,a) = V(bk,ac)*t(c,k)-0.5*V(lk,ac)*t(bc,lk)
    !----------------------------
    call J2_Fvv%arr_clear
!    call Fae_r1abci(J2_Fvv, T1, abci, 1.0d0, iCore(i0), maxCore)
    call Fae_r2abij(J2_Fvv, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! F(i,a) = V(ki,ca)*t(c,k)
    !----------------------------
!    call J2_Fov%arr_clear
!    call Fme_r1abij(J2_Fov, T1, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(jc,ik) = P(ik){V(cj,kd)*t(d,i)}-V(jl,ik)*t(c,l)
    !           +P(ik){V(jl,id)*t(dc,lk)}+0.5*V(cj,de)*t(de,ki)
    !----------------------------
    call J2_Wovoo%arr_clear
!    call Wkbij_r1iabj(J2_Wovoo, T1, iabj, -1.0d0, iCore(i0), maxCore)
!    call Wkbij_r1Uklij(J2_Wovoo, T1, ijkl, 1.0d0, iCore(i0), maxCore)
    call Wkbij_r2ijka(J2_Wovoo, T2, ijka, 1.0d0, iCore(i0), maxCore)
    call Wkbij_r2abci(J2_Wovoo, T2, abci, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(ja,bi) = -V(jk,bi)*t(a,k)+V(ja,bc)*t(c,i)
    !           +V(kj,cb)*t(ac,ik)
    !----------------------------
    call J2_Wovvo%arr_clear
!    call Wmbej_r1ijka(J2_Wovvo, T1, ijka, 1.0d0, iCore(i0), maxCore)
!    call Wmbej_r1abci(J2_Wovvo, T1, abci, 1.0d0, iCore(i0), maxCore)
    call Wmbej_r2abij(J2_Wovvo, T2, abij, 2.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(mn,ij) = V(mm,ia)*t(a,j)+0.5*V(mn,ab)*t(ab,ij)
    !----------------------------
!    call J2_Woooo%arr_clear
!    call Wmnij_r1ijka(J2_Woooo, T1, ijka, 1.0d0, iCore(i0), maxCore)
!    call Wmnij_r2abij(J2_Woooo, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(jk,ia) = V(jk,ca)*t(c,i)
    !----------------------------
!    call J2_Wooov%arr_clear
!    call Wlkjd_r1abij(J2_Wooov, T1, abij, 1.0d0, iCore(i0), maxCore)

    write (6, '(A,F)') 'norm of Foo^(2):', J2_Foo%norm()
    write (6, '(A,F)') 'norm of Fvv^(2):', J2_Fvv%norm()
!    write (6, '(A,F)') 'norm of Fov^(2):', J2_Fov%norm()
!    write (6, '(A,F)') 'norm of Woooo^(2):', J2_Woooo%norm()
!    write (6, '(A,F)') 'norm of Wooov^(2):', J2_Wooov%norm()
    write (6, '(A,F)') 'norm of Wovoo^(2):', J2_Wovoo%norm()
    write (6, '(A,F)') 'norm of Wovvo^(2):', J2_Wovvo%norm()

end subroutine calJ2

!> J^(3) = <u|[[V,T^(2)],tau_v]|0>+0.5*<u|[[[V,T^(1)],T^(1)],tau_v]|0>
! T1^(1)=0, R2^(0)=0, R1^(1)=0, cause
!   [V,T^(1)],T^(1)] term be 0
!   for CPS(D-3), only need intermediates can contract with R1 and result C1
subroutine calJ3(iCore, istart)
    use mod_variable
    use mod_orbit, only: zorb
    use mod_common, only: iop
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart

    type(opt) :: T1, T2, tmp
    integer*8 :: i0, i1, maxCore, idx

    i0 = istart
    call T1%init(iCore, popv, popo, zorb, zorb, 1, 1, 1, idx=i0)
    i0 = i0 + T1%len
    call T2%init(iCore, popv, popv, popo, popo, iop, iop, 1, idx=i0)
    i0 = i0 + T2%len
    maxCore = nCore - i0 + 1

    !------ [V,T] term
    !--- read in T^(2)
    idx = T1%len + T2%len + 1
    call IO_stream(0.0d0, T1%arr, idx, T1%len + T2%len, f_T, '13', 'old')

    !----------------------------
    ! F(j,i) = V(jk,ic)*t(c,k)+0.5*V(jk,dc)*t(dc,ik)
    !----------------------------
    call J3_Foo%arr_clear
    call Fmi_r1ijka(J3_Foo, T1, ijka, 1.0d0, iCore(i0), maxCore)
    call Fmi_r2abij(J3_Foo, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! F(b,a) = V(bk,ac)*t(c,k)-0.5*V(lk,ac)*t(bc,lk)
    !----------------------------
    call J3_Fvv%arr_clear
    call Fae_r1abci(J3_Fvv, T1, abci, 1.0d0, iCore(i0), maxCore)
    call Fae_r2abij(J3_Fvv, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! F(i,a) = V(ki,ca)*t(c,k)
    !----------------------------
!    call J3_Fov%arr_clear
!    call Fme_r1abij(J3_Fov, T1, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(jc,ik) = P(ik){V(cj,kd)*t(d,i)}-V(jl,ik)*t(c,l)
    !           +P(ik){V(jl,id)*t(dc,lk)}+0.5*V(cj,de)*t(de,ki)
    !----------------------------
!    call J3_Wovoo%arr_clear
!    call Wkbij_r1iabj(J3_Wovoo, T1, iabj, -1.0d0, iCore(i0), maxCore)
!    call Wkbij_r1Uklij(J3_Wovoo, T1, ijkl, 1.0d0, iCore(i0), maxCore)
!    call Wkbij_r2ijka(J3_Wovoo, T2, ijka, 1.0d0, iCore(i0), maxCore)
!    call Wkbij_r2abci(J3_Wovoo, T2, abci, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(ja,bi) = -V(jk,bi)*t(a,k)+V(ja,bc)*t(c,i)
    !           +V(kj,cb)*t(ac,ik)
    !----------------------------
    call J3_Wovvo%arr_clear
    call Wmbej_r1ijka(J3_Wovvo, T1, ijka, 1.0d0, iCore(i0), maxCore)
    call Wmbej_r1abci(J3_Wovvo, T1, abci, 1.0d0, iCore(i0), maxCore)
    call Wmbej_r2abij(J3_Wovvo, T2, abij, 2.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(mn,ij) = V(mm,ia)*t(a,j)+0.5*V(mn,ab)*t(ab,ij)
    !----------------------------
!    call J3_Woooo%arr_clear
!    call Wmnij_r1ijka(J3_Woooo, T1, ijka, 1.0d0, iCore(i0), maxCore)
!    call Wmnij_r2abij(J3_Woooo, T2, abij, 1.0d0, iCore(i0), maxCore)

    !----------------------------
    ! W(jk,ia) = V(jk,ca)*t(c,i)
    !----------------------------
!    call J3_Wooov%arr_clear
!    call Wlkjd_r1abij(J3_Wooov, T1, abij, 1.0d0, iCore(i0), maxCore)

    !------ [[V,T],T] term
    !--- read in T^(1)
!    call IO_stream(0.0d0, T1%arr, 1, T1%len + T2%len, f_T, '13', 'old')

    !----------------------------
    ! F(j,i) += V(jk,dc)*t(d,i)*t(c,k)
    ! F(b,a) -= V(lk,ac)*t(d,l)*t(c,k)
    ! W(jc,ik) += V(lj,de)*t(d,l)*t(ec,ik)
    !----------------------------
!    call copy_opt(fia, tmp, T, iCore, i0)
!    call tmp%arr_clear
!    i1 = i0 + tmp%len
!    call Fme_r1abij(tmp, T1, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
!
!    call Fmi_r1fia(J3_Foo, T1, tmp, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Fae_r1fia(J3_Fvv, T1, tmp, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Wkbij_r2Fkc(J3_Wovoo, T2, tmp, -1.0d0, iCore(i1), maxCore - tmp%len)

    !----------------------------
    ! W(aj,ib) -= V(kj,cb)*t(c,i)*t(a,k)
    ! W(mn,ij) += 0.5*V(mn,cd)*t(c,i)*t(d,j)
    ! W(jc,ik) += V(jl,ed)*t(e,i)*t(dc,lk)
    !----------------------------
!    call copy_opt(ijka, tmp, T, iCore, i0)
!    call tmp%arr_clear
!    i1 = i0 + tmp%len
!    call Wlkjd_r1abij(tmp, T1, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
!
!    call Wmbej_r1ijka(J3_Wovvo, T1, tmp, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Wmnij_r1ijka(J3_Woooo, T1, tmp, 0.5d0, iCore(i1), maxCore - tmp%len)
!    call Wkbij_r2ijka(J3_Wovoo, T2, tmp, 1.0d0, iCore(i1), maxCore - tmp%len)

    !----------------------------
    ! W(jc,ik) += -V(lj,kd)*t(c,l)*t(d,i)
    !            -0.5*V(jl,de)*t(de,ik)*t(c,l)
    !            +0.5*V(cj,de)*t(d,k)*t(e,i)
    !----------------------------
!    call copy_opt(ijkl, tmp, T, iCore, i0)
!    call tmp%arr_clear
!    i1 = i0 + tmp%len
!    call Wmnij_r1ijka(tmp, T1, ijka, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Wmnij_r2abij(tmp, T2, abij, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Wkbij_r1Uklij(J3_Wovoo, T1, tmp, 1.0d0, iCore(i1), maxCore - tmp%len)
!
!    call copy_opt(iabj, tmp, T, iCore, i0)
!    call tmp%arr_clear
!    i1 = i0 + tmp%len
!    call Wmbej_r1abci(tmp, T1, abci, 1.0d0, iCore(i1), maxCore - tmp%len)
!    call Wkbij_r1iabj(J3_Wovoo, T1, tmp, -0.5d0, iCore(i1), maxCore - tmp%len)

    write (6, '(A,F)') 'norm of Foo^(3):', J3_Foo%norm()
    write (6, '(A,F)') 'norm of Fvv^(3):', J3_Fvv%norm()
!    write (6, '(A,F)') 'norm of Fov^(3):', J3_Fov%norm()
!    write (6, '(A,F)') 'norm of Woooo^(3):', J3_Woooo%norm()
!    write (6, '(A,F)') 'norm of Wooov^(3):', J3_Wooov%norm()
!    write (6, '(A,F)') 'norm of Wovoo^(3):', J3_Wovoo%norm()
    write (6, '(A,F)') 'norm of Wovvo^(3):', J3_Wovvo%norm()

end subroutine calJ3
