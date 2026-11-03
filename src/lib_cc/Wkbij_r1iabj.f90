! 2026.7.6
!> W(ij,kb) -= fact*P(ij){Vkbcj(ck,bj)*r(c,i)}
subroutine Wkbij_r1iabj(Wkbij, r1, iabj, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wkbij
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: iabj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(iabj, r1, Wkbij, -fact, 'ckbj', 'ci', 'ijkb', iCore, maxCore, 1)

end subroutine Wkbij_r1iabj

! 2026.7.6
!> W(ij,kb) -= fact*P(ij){Vkbcj(ck,bj)*r(c,i)}
subroutine Wkbij_r1iabj_R(WkbijABAB, r1AA, iabjABAB, iabjBAAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: iabjABAB, iabjBAAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(iJ,kB) = -VkBcJ(ck,BJ)*r(c,i)
    call vv2v_424(iabjABAB, r1AA, WkbijABAB, -fact, 'ckBJ', 'ci', 'iJkB', iCore, maxCore, 0)
    ! W(iJ,kB) = VkBCi(Ck,Bi)*r(C,J) -> VKbcI(cK,bI)*r(c,j)
    call vv2v_424(iabjBAAB, r1AA, WkbijABAB, fact, 'ckbi', 'cj', 'ijkb', iCore, maxCore, 0)

end subroutine Wkbij_r1iabj_R

! 2026.7.6
!> W(ij,kb) -= fact*P(ij){Vkbcj(ck,bj)*r(c,i)}, need V(kb,cj)=V(jc,bk)
subroutine Wkbij_r1iabj_U(WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA, spin, &
                          r1AA, r1BB, iabjAAAA, iabjABAB, iabjBAAB, iabjBBBB, iabjABBA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA
    character(len=4) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: iabjAAAA, iabjABAB, iabjBAAB, iabjBBBB, iabjABBA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ij,kb) -= fact*P(ij){Vkbcj(ck,bj)*r(c,i)}
        call Wkbij_r1iabj(WkbijAAAA, r1AA, iabjAAAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(iJ,kB) = -VkBcJ(ck,BJ)*r(c,i)
        call vv2v_424(iabjABAB, r1AA, WkbijABAB, -fact, 'ckBJ', 'ci', 'iJkB', iCore, maxCore, 0)
        ! W(iJ,kB) = VkBCi(Ck,Bi)*r(C,J)
        call vv2v_424(iabjABBA, r1BB, WkbijABAB, fact, 'CkBi', 'CJ', 'iJkB', iCore, maxCore, 0)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(IJ,KB) -= fact*P(IJ){VKBCJ(CK,BJ)*r(C,I)}
        call Wkbij_r1iabj(WkbijBBBB, r1BB, iabjBBBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(4:4) == 'T') then
        ! W(Ij,Kb) = -VKbCj(CK,bj)*r(C,I) -> -VjCbK(bj,CK)*r(C,I)
        call vv2v_424(iabjABAB, r1BB, WkbijBABA, -fact, 'bjCK', 'CI', 'IjKb', iCore, maxCore, 0)
        ! W(Ij,Kb) = VKbcI(cK,bI)*r(c,j)
        call vv2v_424(iabjBAAB, r1AA, WkbijBABA, fact, 'cKbI', 'cj', 'IjKb', iCore, maxCore, 0)
    end if

end subroutine Wkbij_r1iabj_U
