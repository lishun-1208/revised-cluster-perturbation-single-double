! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ab,mn)*W(mn,ij)
subroutine c2_r2Wmnij(c2, r2, Wmnij, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Wmnij
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (r2%ioppq == 0) then
        call vv2v_444(r2, Wmnij, c2, fact, 'abmn', 'mnij', 'abij', iCore, maxCore, 0)
    else
        call vv2v_444(r2, Wmnij, c2, 0.5d0*fact, 'abmn', 'mnij', 'abij', iCore, maxCore, 0)
    end if

end subroutine c2_r2Wmnij

! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ab,mn)*W(mn,ij)
subroutine c2_r2Wmnij_R(c2ABAB, r2ABAB, WmnijABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: WmnijABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) += r(aB,mN)*W(mN,iJ)
    call c2_r2Wmnij(c2ABAB, r2ABAB, WmnijABAB, 2.0d0*fact, iCore, maxCore)

end subroutine c2_r2Wmnij_R

! 2026.7.6
!> c(ab,ij) += 0.5*fact*r(ab,mn)*W(mn,ij)
subroutine c2_r2Wmnij_U(c2AAAA, c2ABAB, c2BBBB, r2AAAA, r2ABAB, r2BBBB, &
                        WmnijAAAA, WmnijABAB, WmnijBBBB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: WmnijAAAA, WmnijABAB, WmnijBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) += 0.5*fact*r(ab,mn)*W(mn,ij)
    call c2_r2Wmnij(c2AAAA, r2AAAA, WmnijAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) += r(aB,mN)*W(mN,iJ)
    call c2_r2Wmnij(c2ABAB, r2ABAB, WmnijABAB, 2.0d0*fact, iCore, maxCore)

    ! c(AB,IJ) += 0.5*fact*r(AB,MN)*W(MN,IJ)
    call c2_r2Wmnij(c2BBBB, r2BBBB, WmnijBBBB, fact, iCore, maxCore)

end subroutine c2_r2Wmnij_U
