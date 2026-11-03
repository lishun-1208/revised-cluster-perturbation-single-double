! 2026.7.6
!> c(a,i) += fact*r(ae,im)*F(e,m)
subroutine c1_r2Fme(c1, r2, Fme, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Fme
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(r2, Fme, c1, fact, 'aeim', 'em', 'ai', iCore, maxCore)

end subroutine c1_r2Fme

! 2026.7.6
!> c(a,i) += fact*r(ae,im)*F(e,m)
subroutine c1_r2Fme_R(c1AA, r2AAAA, r2ABAB, FmeAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: FmeAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) += fact*r(ae,im)*F(e,m)
    call c1_r2Fme(c1AA, r2AAAA, FmeAA, fact, iCore, maxCore)
    ! c(a,i) += fact*r(aE,iM)*F(E,M)
    call c1_r2Fme(c1AA, r2ABAB, FmeAA, fact, iCore, maxCore)

end subroutine c1_r2Fme_R

! 2026.7.6
!> c(a,i) += fact*r(ae,im)*F(e,m)
subroutine c1_r2Fme_U(c1AA, c1BB, r2AAAA, r2ABAB, r2BBBB, FmeAA, FmeBB, &
                      fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: FmeAA, FmeBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = r(ae,im)*F(e,m)
    call c1_r2Fme(c1AA, r2AAAA, FmeAA, fact, iCore, maxCore)
    ! c(a,i) = r(aE,iM)*F(E,M)
    call c1_r2Fme(c1AA, r2ABAB, FmeBB, fact, iCore, maxCore)

    ! c(A,I) = r(Ae,Im)*F(e,m) = r(eA,mI)*F(e,m)
    call vv2v_422(r2ABAB, FmeAA, c1BB, fact, 'eAmI', 'em', 'AI', iCore, maxCore)
    ! c(A,I) = r(AE,IM)*F(E,M)
    call c1_r2Fme(c1BB, r2BBBB, FmeBB, fact, iCore, maxCore)

end subroutine c1_r2Fme_U
