! write by lishun, 2024.12.5
! c(a,i) += fact*<a,i|[V,R2]|0>
!        += 0.5*<ja||bc>*r(bc,ji)-0.5*<jk||bi>*r(ba,jk)
subroutine calC1_VR2(iCore, maxCore, istart, C1, R2, symR, fact)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: istart
    type(opt), intent(inout) :: C1
    type(opt), intent(in) :: R2
    integer*8, intent(in) :: symR ! symmetry of R2
    real*8, intent(in) :: fact

    integer*8 :: i0

    i0 = istart

    if (debug) then
        write (6, '(A)') '********** start calC1_VR2 **********'
        write (6, *) 'norm of R2:', R2%norm()
    end if

    ! c(a,i) += 0.5*<ja||bc>*r(bc,ji)
    call c1_r2abci(C1, R2, abci, fact, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, *) 'norm of C1 1:', C1%norm()
    end if

    ! c(a,i) -= 0.5*<jk||bi>*r(ba,jk)
    call c1_r2ijka(C1, R2, ijka, fact, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, *) 'norm of C1 2:', C1%norm()
        write (6, '(A)') '*************************************'
    end if

end subroutine calC1_VR2
