! write by LiShun, 2024.9.20
!> calculate index of V(pq,rs) for specified (irrepp,iirepq,irrepr,irreps,p,q,r,s)
! if one index is none, ie. s, need set pops=zorb,irreps=1,s=1
!----------------------------------
! modified at 2024.11.13
subroutine putidx(idxout, irrepp, irrepq, irrepr, irreps, p, q, r, s, &
                  popp, popq, ioppq, popr, pops, ioprs, symV)
    use mod_symm
    implicit none

    integer*8, intent(out) :: idxout !> the idx of V(pq,rs) in 1d array
    integer*8, intent(in) :: irrepp, irrepq, irrepr, irreps
    integer*8, intent(in) :: p, q, r, s !> symmerty adapted p,q,r,s
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    integer*8, intent(in) :: symV !> irrep of V(pq,rs)

    integer*8 :: irreppq, irreprs
    integer*8 :: numpq, numrs
    integer*8 :: numr, nump
    integer*8, dimension(8) :: ioff ! idx for V belong to irreprs
    integer*8, dimension(7) :: idx ! idx for V belong to irreprs,irreps,s,r,irrepq,q,p
    integer*8, external :: irpdso, isymoffso

    !----------- give ioff
    ioff(1) = 0
    do irreprs = 1, nirrep - 1
        irreppq = dirprd(irreprs, symV)
        numpq = irpdso(irreppq, popp, popq, ioppq)
        numrs = irpdso(irreprs, popr, pops, ioprs)
        ioff(irreprs + 1) = ioff(irreprs) + numpq*numrs
    end do

    !----------- give idx
    irreprs = dirprd(irrepr, irreps)
    irreppq = dirprd(irrepp, irrepq)
    if (dirprd(irreprs, irreppq) /= symV) then
        write (6, *) 'symmetry error in putidx', irrepp, irrepq, irrepr, irreps, symV
        stop
    end if
    numrs = irpdso(irreprs, popr, pops, ioprs)
    numpq = irpdso(irreppq, popp, popq, ioppq)
    numr = popr(irrepr)
    nump = popp(irrepp)

    idx(1) = ioff(irreprs)
    idx(2) = idx(1) + numpq*(isymoffso(irreps, irreprs, popr, pops, ioprs) - 1)
    ! r,s
    if (ioprs == 1) then
        idx(3) = idx(2) + numpq*numr*(s - 1)
        idx(4) = idx(3) + numpq*(r - 1)
    else if (ioprs == 0) then
        if (irrepr > irreps) then ! ioprs=0 need irrepr<irreps
            write (6, *) 'symmetry of rs error in putidx', irrepr, irreps
            stop
        else if (irrepr < irreps) then
            idx(3) = idx(2) + numpq*numr*(s - 1)
            idx(4) = idx(3) + numpq*(r - 1)
        else if (irrepr == irreps) then ! need r<s
            if (r >= s) then
                write (6, *) 'rs error in putidx', r, s
                stop
            else if (r < s) then
                idx(3) = idx(2) + numpq*((s - 1)*(s - 2)/2)
                idx(4) = idx(3) + numpq*(r - 1)
            end if
        end if
    end if
    idx(5) = idx(4) + isymoffso(irrepq, irreppq, popp, popq, ioppq) - 1
    ! p,q
    if (ioppq == 1) then
        idx(6) = idx(5) + nump*(q - 1)
        idx(7) = idx(6) + p
    else if (ioppq == 0) then
        if (irrepp > irrepq) then ! ioppq=0 need irrepp<irrepq
            write (6, *) 'symmetry of pq error in putidx', irrepp, irrepq
            stop
        else if (irrepp < irrepq) then
            idx(6) = idx(5) + nump*(q - 1)
            idx(7) = idx(6) + p
        else if (irrepp == irrepq) then ! need p<q
            if (p >= q) then
                write (6, *) 'pq error in putidx', p, q
                stop
            else if (p < q) then
                idx(6) = idx(5) + (q - 1)*(q - 2)/2
                idx(7) = idx(6) + p
            end if
        end if
    end if

    idxout = idx(7)

end subroutine putidx
