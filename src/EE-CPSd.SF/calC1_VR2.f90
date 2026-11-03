! write by lishun, 2024.12.5
! c(a,i) += fact*<a,i|[V,pR2]|0>
!        += 0.5*<ja||bc>*r(bc,ji)-0.5*<jk||bi>*r(ba,jk)
!-------------------------------------
! modified at 2025.1.13
subroutine calC1_VR2(iCore, istart, pC1, pR2, symR, fact)
    use mod_variable
    use mod_orbit, only: zorb
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC1
    type(opt), intent(in) :: pR2
    integer*8, intent(in) :: symR ! symmetry of pR2
    real*8, intent(in) :: fact

    integer*8 :: nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    ! c(a,i) += 0.5*<ja||bc>*r(bc,ji)
    call c1_r2abci(pC1, pR2, abci, fact, iCore(istart), nCore - istart + 1)

    ! c(a,i) -= 0.5*<jk||bi>*r(ba,jk)
    call c1_r2ijka(pC1, pR2, ijka, fact, iCore(istart), nCore - istart + 1)

end subroutine calC1_VR2
