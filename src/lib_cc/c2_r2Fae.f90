! 2026.7.6
!> c(ab,ij) += fact*P(ab){r(ae,ij)*F(b,e)}
subroutine c2_r2Fae(c2, r2, Fae, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Fae
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(r2, Fae, c2, fact, 'aeij', 'be', 'abij', iCore, maxCore, 1)

end subroutine c2_r2Fae

! 2026.7.6
!> c(ab,ij) += fact*P(ab){r(ae,ij)*F(b,e)}
subroutine c2_r2Fae_R(c2ABAB, r2ABAB, FaeAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: FaeAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) += fact*r(aE,iJ)*F(B,E)
    call vv2v_424(r2ABAB, FaeAA, c2ABAB, fact, 'aeij', 'be', 'abij', iCore, maxCore, 0)
    ! c(aB,iJ) -= fact*r(Be,iJ)*F(a,e) = r(eB,iJ)*F(a,e)
    call vv2v_424(r2ABAB, FaeAA, c2ABAB, fact, 'eBiJ', 'ae', 'aBiJ', iCore, maxCore, 0)

end subroutine c2_r2Fae_R

! 2026.7.6
!> c(ab,ij) += fact*P(ab){r(ae,ij)*F(b,e)}
subroutine c2_r2Fae_U(c2AAAA, c2ABAB, c2BBBB, r2AAAA, r2ABAB, r2BBBB, FaeAA, FaeBB, &
                      fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: FaeAA, FaeBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) += fact*P(ab){r(ae,ij)*F(b,e)}
    call c2_r2Fae(c2AAAA, r2AAAA, FaeAA, fact, iCore, maxCore)

    ! c(aB,iJ) += fact*r(aE,iJ)*F(B,E)
    call vv2v_424(r2ABAB, FaeBB, c2ABAB, fact, 'aEiJ', 'BE', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) -= fact*r(Be,iJ)*F(a,e) = r(eB,iJ)*F(a,e)
    call vv2v_424(r2ABAB, FaeAA, c2ABAB, fact, 'eBiJ', 'ae', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) += fact*P(AB){r(AE,IJ)*F(B,E)}
    call c2_r2Fae(c2BBBB, r2BBBB, FaeBB, fact, iCore, maxCore)

end subroutine c2_r2Fae_U
