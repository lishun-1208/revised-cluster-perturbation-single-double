! 2026.7.4
!> F(e,m) += fact*r(f,n)*V(ef,mn)
subroutine Fme_r1abij(Fme, r1, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fme
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(abij, r1, Fme, fact, 'efmn', 'fn', 'em', iCore, maxCore)

end subroutine Fme_r1abij

! 2026.7.4
!> F(e,m) += fact*r(f,n)*V(ef,mn)
subroutine Fme_r1abij_R(FmeAA, r1AA, abijAAAA, abijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmeAA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: abijAAAA, abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(e,m) += fact*r(f,n)*V(ef,mn)
    call Fme_r1abij(FmeAA, r1AA, abijAAAA, fact, iCore, maxCore)
    ! F(e,m) += fact*r(F,N)*V(eF,mN)
    call Fme_r1abij(FmeAA, r1AA, abijABAB, fact, iCore, maxCore)

end subroutine Fme_r1abij_R

! 2026.7.4
!> F(e,m) += fact*r(f,n)*V(ef,mn)
subroutine Fme_r1abij_U(FmeAA, FmeBB, r1AA, r1BB, abijAAAA, abijABAB, abijBBBB, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmeAA, FmeBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(e,m) = r(f,n)*V(ef,mn)
    call Fme_r1abij(FmeAA, r1AA, abijAAAA, fact, iCore, maxCore)
    ! F(e,m) = r(F,N)*V(eF,mN)
    call Fme_r1abij(FmeAA, r1BB, abijABAB, fact, iCore, maxCore)

    ! F(E,M) = r(F,N)*V(EF,MN)
    call Fme_r1abij(FmeBB, r1BB, abijBBBB, fact, iCore, maxCore)
    ! F(E,M) = r(f,n)*V(Ef,Mn) = r(f,n)*V(fE,nM)
    call vv2v_422(abijABAB, r1AA, FmeBB, fact, 'fEnM', 'fn', 'EM', iCore, maxCore)

end subroutine Fme_r1abij_U
