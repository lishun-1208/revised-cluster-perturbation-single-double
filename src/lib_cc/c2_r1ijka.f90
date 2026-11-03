! 2026.7.6
!> c(ab,ij) -= fact*P(ab){r(a,m)*V(ij,mb)}
subroutine c2_r1ijka(c2, r1, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(ijka, r1, c2, -fact, 'ijmb', 'am', 'abij', iCore, maxCore, 1)

end subroutine c2_r1ijka

! 2026.7.6
!> c(ab,ij) -= fact*P(ab){r(a,m)*V(ij,mb)}
subroutine c2_r1ijka_R(c2ABAB, r1AA, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) -= fact*r(a,m)*V(iJ,mB)
    call vv2v_424(ijkaABAB, r1AA, c2ABAB, -fact, 'iJmB', 'am', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = fact*r(B,M)*V(iJ,Ma) = -r(B,M)*V(Ji,Ma) -> -r(b,m)*V(jI,mA)
    call vv2v_424(ijkaABAB, r1AA, c2ABAB, -fact, 'jima', 'bm', 'abij', iCore, maxCore, 0)

end subroutine c2_r1ijka_R

! 2026.7.6
!> c(ab,ij) -= fact*P(ab){r(a,m)*V(ij,mb)}
subroutine c2_r1ijka_U(c2AAAA, c2ABAB, c2BBBB, r1AA, r1BB, &
                       ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) -= fact*P(ab){r(a,m)*V(ij,mb)}
    call c2_r1ijka(c2AAAA, r1AA, ijkaAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) -= fact*r(a,m)*V(iJ,mB)
    call vv2v_424(ijkaABAB, r1AA, c2ABAB, -fact, 'iJmB', 'am', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = fact*r(B,M)*V(iJ,Ma) = -r(B,M)*V(Ji,Ma)
    call vv2v_424(ijkaBABA, r1BB, c2ABAB, -fact, 'JiMa', 'BM', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) -= fact*P(AB){r(A,M)*V(IJ,MB)}
    call c2_r1ijka(c2BBBB, r1BB, ijkaBBBB, fact, iCore, maxCore)

end subroutine c2_r1ijka_U
