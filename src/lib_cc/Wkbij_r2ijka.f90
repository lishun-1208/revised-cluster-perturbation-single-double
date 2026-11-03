! 2026.7.10
!> W(ij,kb) += fact*P(ij){V(kl,ic)*r(bc,jl)}
subroutine Wkbij_r2ijka(Wkbij, r2, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: Wkbij
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_444(ijka, r2, Wkbij, fact, 'klic', 'bcjl', 'ijkb', iCore, maxCore, 1)

end subroutine Wkbij_r2ijka

! 2026.7.10
!> W(ij,kb) += fact*P(ij){V(kl,ic)*r(bc,jl)}
subroutine Wkbij_r2ijka_R(WkbijABAB, r2AAAA, r2ABAB, ijkaAAAA, ijkaABAB, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijABAB
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! W(iJ,kB) = V(kl,ic)*r(Bc,Jl) =  V(kl,ic)*r(cB,lJ)
    !           +V(kL,iC)*r(BC,JL) -> V(kL,iC)*r(bc,jl)
    !           -V(kL,Jc)*r(Bc,iL) = -V(Lk,Jc)*r(cB,iL) -> -V(lK,jC)*r(cB,iL)
    call vv2v_444(ijkaAAAA, r2ABAB, WkbijABAB, fact, 'klic', 'cblj', 'ijkb', iCore, maxCore, 0)
    call vv2v_444(ijkaABAB, r2AAAA, WkbijABAB, fact, 'klic', 'bcjl', 'ijkb', iCore, maxCore, 0)
    call vv2v_444(ijkaABAB, r2ABAB, WkbijABAB, -fact, 'lkjc', 'cbil', 'ijkb', iCore, maxCore, 0)

end subroutine Wkbij_r2ijka_R

! 2026.7.10
!> W(ij,kb) += fact*P(ij){V(kl,ic)*r(bc,jl)}
subroutine Wkbij_r2ijka_U(WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA, spin, &
                          r2AAAA, r2ABAB, r2BBBB, ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, &
                          fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: WkbijAAAA, WkbijABAB, WkbijBBBB, WkbijBABA
    character(len=4), intent(in) :: spin
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! AAAA
    if (spin(1:1) == 'T') then
        ! W(ij,kb) = P(ij){V(kl,ic)*r(bc,jl)}
        !           +P(ij){V(kL,iC)*r(bC,jL)}
        call Wkbij_r2ijka(WkbijAAAA, r2AAAA, ijkaAAAA, fact, iCore, maxCore)
        call Wkbij_r2ijka(WkbijAAAA, r2ABAB, ijkaABAB, fact, iCore, maxCore)
    end if

    ! ABAB
    if (spin(2:2) == 'T') then
        ! W(iJ,kB) = V(kl,ic)*r(Bc,Jl) =  V(kl,ic)*r(cB,lJ)
        !           +V(kL,iC)*r(BC,JL)
        !           -V(kL,Jc)*r(Bc,iL) = -V(Lk,Jc)*r(cB,iL)
        call vv2v_444(ijkaAAAA, r2ABAB, WkbijABAB, fact, 'klic', 'cBlJ', 'iJkB', iCore, maxCore, 0)
        call vv2v_444(ijkaABAB, r2BBBB, WkbijABAB, fact, 'kLiC', 'BCJL', 'iJkB', iCore, maxCore, 0)
        call vv2v_444(ijkaBABA, r2ABAB, WkbijABAB, -fact, 'LkJc', 'cBiL', 'iJkB', iCore, maxCore, 0)
    end if

    ! BBBB
    if (spin(1:1) == 'T') then
        ! W(IJ,KB) = P(IJ){V(KL,IC)*r(BC,JL)}
        !           +P(IJ){V(Kl,Ic)*r(Bc,Jl)} = P(IJ){V(Kl,Ic)*r(cB,lJ)}
        call Wkbij_r2ijka(WkbijBBBB, r2BBBB, ijkaBBBB, fact, iCore, maxCore)
        call vv2v_444(ijkaBABA, r2ABAB, WkbijBBBB, fact, 'KlIc', 'cBlJ', 'IJKB', &
                      iCore, maxCore, 1)
    end if

    ! BABA
    if (spin(1:1) == 'T') then
        ! W(Ij,Kb) = V(Kl,Ic)*r(bc,jl)
        !           +V(KL,IC)*r(bC,jL)
        !           -V(Kl,jC)*r(bC,Il) = -V(lK,jC)*r(bC,lI)
        call vv2v_444(ijkaBABA, r2AAAA, WkbijBABA, fact, 'KlIc', 'bcjl', 'IjKb', iCore, maxCore, 0)
        call vv2v_444(ijkaBBBB, r2ABAB, WkbijBABA, fact, 'KLIC', 'bCjL', 'IjKb', iCore, maxCore, 0)
        call vv2v_444(ijkaABAB, r2ABAB, WkbijBABA, -fact, 'lKjC', 'bClI', 'IjKb', iCore, maxCore, 0)
    end if

end subroutine Wkbij_r2ijka_U
