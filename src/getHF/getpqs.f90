! write by lishun, 2024.10.30
! get <pq||rs> for one r from hard disk
! will set iop of vecout the same to file
! ATTENTION: not be test for ioprs=0
subroutine getpqs(vecout, r, irrepr, popp, popq, popr, pops, ioppq, ioprs, &
                  file, symV)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: vecout
    integer*8, intent(in) :: r, irrepr
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    character(len=8), intent(in) :: file
    integer*8, intent(in) :: symV

    integer*8 :: irreppq, numpq
    integer*8 :: irreprs, numrs
    integer*8, dimension(2) :: rs
    integer*8 :: irreps, nums, s
    integer*8 :: numr
    integer*8 :: idxout ! V(pq,*s) belong to irreps and s
    integer*8 :: i0
    integer*8, external :: irpdso, isymoffso

    i0 = 1
    idxout = 1
    do irreps = 1, nirrep
        irreprs = dirprd(irrepr, irreps)
        irreppq = dirprd(irreprs, symV)
        numpq = irpdso(irreppq, popp, popq, ioppq)
        numr = popr(irrepr)
        nums = pops(irreps)
        if (r > numr) then
            write (6, *) 'getpqs error', r, numr
            stop
        end if

        rs(1) = isymoffso(irreps, irreprs, popr, pops, ioprs) - 1
        do s = 1, nums
            if (ioprs == 1) then
                rs(2) = rs(1) + numr*(s - 1) + r
                call getlstso(vecout(idxout), rs(2), 1, irreprs, 1, file)
            else if (ioprs == 0) then
                if (r < s) then
                    rs(2) = rs(1) + (s - 1)*(s - 2)/2 + r
                    call getlstso(vecout(idxout), rs(2), 1, irreprs, 1, file)
                else if (r == s) then
                    call dscal(numpq, 0.0d0, vecout(idxout), 1)
                else if (r > s) then
                    rs(2) = rs(1) + (r - 1)*(r - 2)/2 + s
                    call getlstso(vecout(idxout), rs(2), 1, irreprs, 1, file)
                    call dscal(numpq, -1.0d0, vecout(idxout), 1)
                end if
            end if
            idxout = idxout + numpq
        end do
    end do

end subroutine getpqs
