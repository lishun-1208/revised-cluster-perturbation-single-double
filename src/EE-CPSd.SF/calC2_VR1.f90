! write by lishun, 2024.12.6
! c(ab,ij) += fact*<ab,ij|[V,R1]|>
!           = -Pab{<bk||ji>*t(a,k)}+Pij{<ba||jc)*t(c,i)}
!-----------------------------------------
! modified at 2025.1.13
subroutine calC2_VR1(iCore, istart, pC2, pR1, symR, fact)
    use mod_variable
    implicit none

    real*8, dimension(*), intent(inout) :: iCore !> full
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: pC2
    type(opt), intent(in) :: pR1
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact

    integer*8 :: nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    !------ c(ab,ij) -= fact*Pab{<bk||ji>*t(a,k)}
    call c2_r1ijka(pC2, pR1, ijka, fact, iCore(istart), nCore - istart + 1)

    !----- c(ab,ij) += fact*Pij{<ba||jc)*t(c,i)}
    call c2_r1abci(pC2, pR1, abci, fact, iCore(istart), nCore - istart + 1)

end subroutine calC2_VR1
