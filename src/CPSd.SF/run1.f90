! write by lishun, 2024.12.5
! calculate first order perturbation T
! A^(0)*T^(1)+<u|V|0>=0, as
! <u_1|[H,T1^(1)]|0>=0, so T1^(1)=0
! <u_2|[F,T2^(1)]|0>+<u_2|V|0>=0, so T2^(1)=MP2t
! seperate OED term:
! -e*T2^(1)=<u_2|V|0>
! E^(2) = 0.25*<ij||ab>*t(ab,ij)^(1)
subroutine run1(iCore, maxCore)
    use mod_variable
    implicit none

    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    real*8 :: tmp
    integer*8 :: i, idx, nSize

    integer*8, external :: irpdso, idsymsoc
    real*8, external :: dnrm2

    !---------- T1=0
    call dscal(T1%len, 0.0d0, T1%arr, 1)

    !---------- T2=MP2t
    ! t(ab,ij) = -<ab||ij>/e(ab,ij)
    call denomT2(T2, abij, iCore(iEorb))
    write (6, '(A,F)') 'norm of T2^(1) (MP2t):', T2%norm()
    call IO_stream(0.0d0, T1%arr, 1, T1%len + T2%len, f_T, '23', 'replace')

    !------ 2nd order correlation energy, same to MP2 energy
    call calEcorr(tmp, iCore, maxCore, iFree, T2, 1)
    write (6, '(A,F)') '2nd order correlation energy (MP2 energy):', tmp
    Ecorr = tmp

    write (6, '(A,F)') 'norm of T2:', T2%norm()
    call printEXinfo2_uhf(T2, 10)
end subroutine run1
