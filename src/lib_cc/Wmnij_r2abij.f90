! 2026.7.6
!> W(mn,ij) += fact*0.5*r(ef,ij)*V(ef,mn)
subroutine Wmnij_r2abij(Wmnij, r2, abij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wmnij
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: abij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (r2%ioppq == 0) then
        call vv2v_444(abij, r2, Wmnij, fact, 'efmn', 'efij', 'mnij', iCore, maxCore, 0)
    else
        call vv2v_444(abij, r2, Wmnij, 0.5d0*fact, 'efmn', 'efij', 'mnij', iCore, maxCore, 0)
    end if

end subroutine Wmnij_r2abij

! 2026.7.6
!> W(mn,ij) += fact*0.5*r(ef,ij)*V(ef,mn)
subroutine Wmnij_r2abij_R(WmnijABAB, r2ABAB, abijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmnijABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: abijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(mN,iJ) += r(eF,iJ)*V(eF,mN)
    call Wmnij_r2abij(WmnijABAB, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

end subroutine Wmnij_r2abij_R

! 2026.7.6
!> W(mn,ij) += fact*0.5*r(ef,ij)*V(ef,mn)
subroutine Wmnij_r2abij_U(WmnijAAAA, WmnijABAB, WmnijBBBB, &
                          r2AAAA, r2ABAB, r2BBBB, abijAAAA, abijABAB, abijBBBB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WmnijAAAA, WmnijABAB, WmnijBBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: abijAAAA, abijABAB, abijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(mn,ij) += 0.5*r(ef,ij)*V(ef,mn)
    call Wmnij_r2abij(WmnijAAAA, r2AAAA, abijAAAA, fact, iCore, maxCore)

    ! W(mN,iJ) += r(eF,iJ)*V(eF,mN)
    call Wmnij_r2abij(WmnijABAB, r2ABAB, abijABAB, 2.0d0*fact, iCore, maxCore)

    ! W(MN,IJ) += 0.5*r(EF,IJ)*V(EF,MN)
    call Wmnij_r2abij(WmnijBBBB, r2BBBB, abijBBBB, fact, iCore, maxCore)

end subroutine Wmnij_r2abij_U
