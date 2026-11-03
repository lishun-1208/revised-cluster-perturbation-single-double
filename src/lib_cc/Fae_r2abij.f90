! 2026.7.6
!> F(a,e) -= fact*0.5*r(af,mn)*V(ef,mn)
subroutine Fae_r2abij(Fae, r2, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fae
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (r2%ioprs == 0) then
        call vv2v_442(r2, abij, Fae, -1.0d0*fact, 'afmn', 'efmn', 'ae', iCore, maxCore)
    else
        call vv2v_442(r2, abij, Fae, -0.5d0*fact, 'afmn', 'efmn', 'ae', iCore, maxCore)
    end if

end subroutine Fae_r2abij

! 2026.7.6
!> F(a,e) -= fact*0.5*r(af,mn)*V(ef,mn)
subroutine Fae_r2abij_R(FaeAA, r2AAAA, r2ABAB, abijAAAA, abijABAB, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: FaeAA
    type(opt), intent(in) :: r2AAAA, r2ABAB, abijAAAA, abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(a,e) = -0.5*r(af,mn)*V(ef,mn)
    call Fae_r2abij(FaeAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
    ! F(a,e) = -0.5*r(aF,mN)*V(eF,mN)
    !          -0.5*r(aF,nM)*V(eF,nM)
    !        = -r(aF,mN)*V(eF,mN)
    call Fae_r2abij(FaeAA, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

end subroutine Fae_r2abij_R

! 2026.7.6
!> F(a,e) -= fact*0.5*r(af,mn)*V(ef,mn)
subroutine Fae_r2abij_U(FaeAA, FaeBB, r2AAAA, r2ABAB, r2BBBB, abijAAAA, abijABAB, abijBBBB, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FaeAA, FaeBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB, abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(a,e) = -0.5*r(af,mn)*V(ef,mn)
    call Fae_r2abij(FaeAA, r2AAAA, abijAAAA, fact, iCore, maxCore)
    ! F(a,e) = -r(aF,mN)*V(eF,mN)
    call Fae_r2abij(FaeAA, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

    ! F(A,E) = -r(fA,mN)*V(fE,mN)
    call vv2v_442(r2ABAB, abijABAB, FaeBB, -fact, 'fAmN', 'fEmN', 'AE', iCore, maxCore)
    ! F(A,E) = -0.5*r(AF,MN)*V(EF,MN)
    call Fae_r2abij(FaeBB, r2BBBB, abijBBBB, fact, iCore, maxCore)

end subroutine Fae_r2abij_U
