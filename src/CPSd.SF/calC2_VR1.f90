! write by lishun, 2024.12.6
! c(ab,ij) += fact*<ab,ij|[V,R1]|>
!           = -Pab{<bk||ji>*t(a,k)}+Pij{<ba||jc)*t(c,i)}
subroutine calC2_VR1(iCore, maxCore, istart, C2, R1, symR, fact)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    type(opt), intent(inout) :: C2
    type(opt), intent(in) :: R1
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact

    integer*8 :: nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    if (debug) then
        write (6, '(A)') '********** start calC2_VR1 **********'
        write (6, *) 'norm of R1:', R1%norm()
        write (6, *) 'norm of C2 0:', C2%norm()
    end if

    !------ c(ab,ij) -= fact*Pab{<bk||ji>*t(a,k)}
    call c2_r1ijka(C2, R1, ijka, fact, iCore(istart), nCore - istart + 1)
    if (debug) then
        write (6, *) 'norm of C2 1:', C2%norm()
    end if

    !----- c(ab,ij) += fact*Pij{<ba||jc)*t(c,i)}
    call c2_r1abci(C2, R1, abci, fact, iCore(istart), nCore - istart + 1)
    if (debug) then
        write (6, *) 'norm of C2 2:', C2%norm()
        write (6, '(A)') '*********** end calC2_VR1 ***********'
    end if

end subroutine calC2_VR1
