! 2026.7.13
!> c(a,i) -= fact*W(li,ka)*F(k,l)
subroutine c1d_FklWlika(c1d, Fkl, Wlika, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1d
    type(opt), intent(in) :: Fkl
    type(opt), intent(in) :: Wlika
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(Wlika, Fkl, c1d, -fact, 'lika', 'kl', 'ai', iCore, maxCore)

end subroutine c1d_FklWlika

! 2026.7.13
!> c(a,i) -= fact*W(li,ka)*F(k,l)
subroutine c1d_FklWlika_R(c1dAA, FklAA, WlikaAAAA, WlikaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA
    type(opt), intent(in) :: FklAA
    type(opt), intent(in) :: WlikaAAAA, WlikaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -W(li,ka)*F(k,l)
    !          -W(Li,Ka)*F(K,L) -> -W(lI,kA)*F(k,l)
    call c1d_FklWlika(c1dAA, FklAA, WlikaAAAA, fact, iCore, maxCore)
    call c1d_FklWlika(c1dAA, FklAA, WlikaABAB, fact, iCore, maxCore)

end subroutine c1d_FklWlika_R

! 2026.7.13
!> c(a,i) -= fact*W(li,ka)*F(k,l)
subroutine c1d_FklWlika_U(c1dAA, c1dBB, FklAA, FklBB, &
                          WlikaAAAA, WlikaABAB, WlikaBBBB, WlikaBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1dAA, c1dBB
    type(opt), intent(in) :: FklAA, FklBB
    type(opt), intent(in) :: WlikaAAAA, WlikaABAB, WlikaBBBB, WlikaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = -W(li,ka)*F(k,l)
    !          -W(Li,Ka)*F(K,L)
    call c1d_FklWlika(c1dAA, FklAA, WlikaAAAA, fact, iCore, maxCore)
    call c1d_FklWlika(c1dAA, FklBB, WlikaBABA, fact, iCore, maxCore)

    ! c(A,I) = -W(lI,kA)*F(k,l)
    !          -W(LI,KA)*F(K,L)
    call c1d_FklWlika(c1dBB, FklAA, WlikaABAB, fact, iCore, maxCore)
    call c1d_FklWlika(c1dBB, FklBB, WlikaBBBB, fact, iCore, maxCore)

end subroutine c1d_FklWlika_U
