! 2026.7.6
!> F(m,i) += fact*r(e,n)*V(mn,ie)
subroutine Fmi_r1ijka(Fmi, r1, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Fmi
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(ijka, r1, Fmi, fact, 'mnie', 'en', 'mi', iCore, maxCore)

end subroutine Fmi_r1ijka

! 2026.7.6
!> F(m,i) += fact*r(e,n)*V(mn,ie)
subroutine Fmi_r1ijka_R(FmiAA, r1AA, ijkaAAAA, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmiAA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*r(e,n)*V(mn,ie)
    call Fmi_r1ijka(FmiAA, r1AA, ijkaAAAA, fact, iCore, maxCore)
    ! F(m,i) += fact*r(E,N)*V(mN,iE)
    call Fmi_r1ijka(FmiAA, r1AA, ijkaABAB, fact, iCore, maxCore)

end subroutine Fmi_r1ijka_R

! 2026.7.6
!> F(m,i) += fact*r(e,n)*V(mn,ie)
subroutine Fmi_r1ijka_U(FmiAA, FmiBB, r1AA, r1BB, ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: FmiAA, FmiBB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! F(m,i) += fact*r(e,n)*V(mn,ie)
    call Fmi_r1ijka(FmiAA, r1AA, ijkaAAAA, fact, iCore, maxCore)
    ! F(m,i) += fact*r(E,N)*V(mN,iE)
    call Fmi_r1ijka(FmiAA, r1BB, ijkaABAB, fact, iCore, maxCore)

    ! F(M,I) += fact*r(e,n)*V(Mn,Ie)
    call Fmi_r1ijka(FmiBB, r1AA, ijkaBABA, fact, iCore, maxCore)
    ! F(M,I) += fact*r(E,N)*V(MN,IE)
    call Fmi_r1ijka(FmiBB, r1BB, ijkaBBBB, fact, iCore, maxCore)

end subroutine Fmi_r1ijka_U
