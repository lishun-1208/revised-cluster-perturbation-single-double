! write at 2024.12.30
!> calculate norm of V(pq,rs) with large size
!--------------------------------
! modified at 2025.10.11
subroutine calLargeNorm(norm, fileV, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, symV, &
                        iCore, maxCore)
    use mod_symm
    implicit none

    real*8, intent(out) :: norm
    character(len=8), intent(in) :: fileV
    integer*8, dimension(8), intent(in) :: pop_p, pop_q, pop_r, pop_s
    integer*8, intent(in) :: ioppq, ioprs
    integer*8, intent(in) :: symV
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: irreprs, numrs, rs
    integer*8 :: irreppq, numpq
    integer*8 :: nSize
    integer*8 :: idx
    integer*8, external :: idsymsoc, irpdso
    real*8, external :: dnrm2

    norm = 0.0d0
    ! calculate norm of V(pq,rs)
    do irreprs = 1, nirrep
        numrs = irpdso(irreprs, pop_r, pop_s, ioprs)
        irreppq = dirprd(symV, irreprs)
        do rs = 1, numrs
            ! get V(pq,**)
            nSize = irpdso(irreppq, pop_p, pop_q, ioppq)
            call getlstso(iCore, rs, 1, irreprs, 1, fileV)
            ! norm of V(pq,**)
            norm = norm + dnrm2(nSize, iCore, 1)**2
        end do
    end do
    norm = sqrt(norm)

end subroutine calLargeNorm
