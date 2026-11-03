! 2026.7.6
!> c(ab,ij) -= fact*P(ij){r(ab,im)*F(m,j)}
subroutine c2_r2Fmi(c2, r2, Fmi, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Fmi
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(r2, Fmi, c2, -fact, 'abim', 'mj', 'abij', iCore, maxCore, 2)

end subroutine c2_r2Fmi

! 2026.7.6
!> c(ab,ij) -= fact*P(ij){r(ab,im)*F(m,j)}
subroutine c2_r2Fmi_R(c2ABAB, r2ABAB, FmiAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: FmiAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(aB,iJ) -= fact*r(aB,iM)*F(M,J)
    call vv2v_424(r2ABAB, FmiAA, c2ABAB, -fact, 'abim', 'mj', 'abij', iCore, maxCore, 0)
    ! c(aB,iJ) += fact*r(aB,Jm)*F(m,i) = -r(aB,mJ)*F(m,i)
    call vv2v_424(r2ABAB, FmiAA, c2ABAB, -fact, 'aBmJ', 'mi', 'aBiJ', iCore, maxCore, 0)

end subroutine c2_r2Fmi_R

! 2026.7.6
!> c(ab,ij) -= fact*P(ij){r(ab,im)*F(m,j)}
subroutine c2_r2Fmi_U(c2AAAA, c2ABAB, c2BBBB, r2AAAA, r2ABAB, r2BBBB, FmiAA, FmiBB, &
                      fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: FmiAA, FmiBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) -= fact*P(ij){r(ab,im)*F(m,j)}
    call c2_r2Fmi(c2AAAA, r2AAAA, FmiAA, fact, iCore, maxCore)

    ! c(aB,iJ) -= fact*r(aB,iM)*F(M,J)
    call vv2v_424(r2ABAB, FmiBB, c2ABAB, -fact, 'aBiM', 'MJ', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) += fact*r(aB,Jm)*F(m,i) = -r(aB,mJ)*F(m,i)
    call vv2v_424(r2ABAB, FmiAA, c2ABAB, -fact, 'aBmJ', 'mi', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) -= fact*P(IJ){r(AB,IM)*F(M,J)}
    call c2_r2Fmi(c2BBBB, r2BBBB, FmiBB, fact, iCore, maxCore)

end subroutine c2_r2Fmi_U
