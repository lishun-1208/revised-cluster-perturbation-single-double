! 2026.7.6
!> F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
subroutine Fmi_r2abij(Fmi, r2, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fmi
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (r2%ioppq == 1) then
        call vv2v_442(abij, r2, Fmi, 0.5d0*fact, 'efmn', 'efin', 'mi', iCore, maxCore)
    else
        call vv2v_442(abij, r2, Fmi, fact, 'efmn', 'efin', 'mi', iCore, maxCore)
    end if

end subroutine Fmi_r2abij

! 2026.7.6
!> F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
subroutine Fmi_r2abij_R(FmiAA, r2AAAA, r2ABAB, abijAAAA, abijABAB, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: FmiAA
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: abijAAAA, abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
    call Fmi_r2abij(FmiAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
    ! F(m,i) += fact*r(eF,iN)*V(eF,mN)
    call Fmi_r2abij(FmiAA, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

end subroutine Fmi_r2abij_R

! 2026.7.6
!> F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
subroutine Fmi_r2abij_U(FmiAA, FmiBB, r2AAAA, r2ABAB, r2BBBB, abijAAAA, abijABAB, abijBBBB, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmiAA, FmiBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*0.5*r(ef,in)*V(ef,mn)
    call Fmi_r2abij(FmiAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
    ! F(m,i) += fact*r(eF,iN)*V(eF,mN)
    call Fmi_r2abij(FmiAA, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

    ! F(M,I) += fact*0.5*r(EF,IN)*V(EF,MN)
    call Fmi_r2abij(FmiBB, r2BBBB, abijBBBB, fact, iCore, maxCore)
    ! F(M,I) += fact*r(fE,nI)*V(fE,nM)
    call vv2v_442(abijABAB, r2ABAB, FmiBB, fact, 'fEnM', 'fEnI', 'MI', iCore, maxCore)

end subroutine Fmi_r2abij_U
