! 2026.7.10
!> W(ij,kb) -= fact*F(c,k)*r(cb,ij)
subroutine Wkbij_r2Fkc(Wkbij, r2, Fkc, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wkbij
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Fkc
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_424(r2, Fkc, Wkbij, -fact, 'cbij', 'ck', 'ijkb', iCore, maxCore, 0)

end subroutine Wkbij_r2Fkc

! 2026.7.10
!> W(ij,kb) -= fact*F(c,k)*r(cb,ij)
subroutine Wkbij_r2Fkc_R(WkbijABAB, r2ABAB, FkcAA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijABAB
    type(opt), intent(in) :: r2ABAB
    type(opt), intent(in) :: FkcAA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(iJ,kB) = -F(c,k)*r(cB,iJ)
    call Wkbij_r2Fkc(WkbijABAB, r2ABAB, FkcAA, fact, iCore, maxCore)

end subroutine Wkbij_r2Fkc_R

! 2026.7.10
!> W(ij,kb) -= fact*F(c,k)*r(cb,ij)
subroutine Wkbij_r2Fkc_U(WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA, spin, &
                         r2AAAA, r2ABAB, r2BBBB, FkcAA, FkcBB, &
                         fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA
    character(len=4) :: spin !> 'T'/'F', whether calculate corresponding spin term
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: FkcAA, FkcBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ij,kb) = -F(c,k)*r(cb,ij)
        call Wkbij_r2Fkc(WkbijAAAA, r2AAAA, FkcAA, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(iJ,kB) = -F(c,k)*r(cB,iJ)
        call Wkbij_r2Fkc(WkbijABAB, r2ABAB, FkcAA, fact, iCore, maxCore)
    end if

    ! BBBB
    if (spin(3:3) == 'T') then
        ! W(IJ,KB) = -F(C,K)*r(CB,IJ)
        call Wkbij_r2Fkc(WkbijBBBB, r2BBBB, FkcBB, fact, iCore, maxCore)
    end if

    ! BABA
    if (spin(1:1) == 'T') then
        ! W(Ij,Kb) = -F(C,K)*r(Cb,Ij) = -F(C,K)*r(bC,jI)
        call vv2v_424(r2ABAB, FkcBB, WkbijBABA, -fact, 'bCjI', 'CK', 'IjKb', iCore, maxCore, 0)
    end if

end subroutine Wkbij_r2Fkc_U
