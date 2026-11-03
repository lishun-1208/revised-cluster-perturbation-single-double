! write by lishun, 2024.12.6
! c(ab,ij) += fact*<ab,ij|[V,R2]|0>
!           = PijPab{<ak||ic>*r(cb,kj)}+0.5*<kl||ij>*r(ab,kl)+0.5*<ab||cd>*r(cd,ij)
subroutine calC2_VR2(iCore, maxCore, istart, C2, R2, symR, fact)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore, istart
    type(opt), intent(inout) :: C2
    type(opt), intent(in) :: R2
    integer*8, intent(in) :: symR
    real*8, intent(in) :: fact

    integer*8 :: nSize
    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    if (debug) then
        write (6, '(A)') '********** start calC2_VR2 **********'
        write (6, *) 'norm of R2:', R2%norm()
    end if

    !---- c(ab,ij) += fact*PijPab{<ak||ic>*r(cb,kj)}
    call c2_r2Wmbej(C2, R2, iabj, fact, iCore(istart), nCore - istart + 1)
    if (debug) then
        write (6, *) 'norm of C2 1:', C2%norm()
    end if

    !---- c(ab,ij) += fact*0.5*<kl||ij>*r(ab,kl)
    call c2_r2Wmnij(C2, R2, ijkl, fact, iCore(istart), nCore - istart + 1)
    if (debug) then
        write (6, *) 'norm of C2 2:', C2%norm()
    end if

    !---- c(ab,ij) += fact*0.5*<ab||cd>*r(cd,ij)
    call c2_r2abcd(C2, R2, abcd, fact, iCore(istart), nCore - istart + 1)
    if (debug) then
        write (6, *) 'norm of C2 3:', C2%norm()
        write (6, '(A)') '*********** end calC2_VR2 ***********'
    end if

end subroutine calC2_VR2
