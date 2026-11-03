! write by lishun, 2024.11.11
! calculate c(ab,ij) += fact*<ab,ij|[V,R2]|0>
! <ab,ij|[V,R2]|0> = PijPab{v(kb,cj)*r(ac,ik)}+0.5*v(ab,cd)*t(cd,ij)+0.5*v(kl,ij)*t(ab,kl)
subroutine calC2_VR2(iCore, pC2, pR2, symm, istart, fact)
    use mod_common
    use mod_variable
    use mod_operator
    implicit none

    real*8, dimension(*), target, intent(inout) :: iCore
    type(opt), intent(inout) :: pC2
    type(opt), intent(in) :: pR2
    integer*8, intent(in) :: symm ! symmetry of C2 and R2
    integer*8, intent(in) :: istart ! pointer of free iCore
    real*8, intent(in) :: fact

    integer*8 :: i0
    integer*8, external :: idsymsoc
    real*8, external :: dnrm2

    i0 = istart

    if (debug) then
        write (6, '(A)') '---------- start C2_VR2 ----------'
        write (6, '(A,F)') 'norm of C2 0:', pC2%norm()
    end if

    ! PijPab{v(kb,cj)*r(ac,ik)}
    call c2_r2Wmbej(pC2, pR2, iabj, fact, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of C2 1:', pC2%norm()
    end if

    ! 0.5*v(kl,ij)*t(ab,kl)
    call c2_r2Wmnij(pC2, pR2, ijkl, fact, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of C2 2:', pC2%norm()
    end if

    ! 0.5*v(ab,cd)*t(cd,ij)
    call c2_r2abcd(pC2, pR2, abcd, fact, iCore(i0), nCore - i0 + 1)
    if (debug) then
        write (6, '(A,F)') 'norm of C2 3:', pC2%norm()
        write (6, '(A)') '----------- end C2_VR2 -----------'
    end if

end subroutine calC2_VR2
