! write at 2025.6.4, lishun
!> calculate the (irrepP,irrepQ,irrepR,irrepS,p,q,r,s) for a specified V(pq,rs)
!> the p,q,r,s is symmerty adapted
!-------------------------
! modified at 2026.1.28
subroutine getidx(idxin, idxout, popp, popq, ioppq, popr, pops, ioprs, sym)
    use mod_symm
    implicit none

    integer, intent(in) :: idxin !> the idx of speciffied V(pq,rs) in 1d array
    integer, dimension(8), intent(out) :: idxout !> (irrepp,irrepq,irrepr,irreps,p,q,r,s)
    integer, dimension(8), intent(in) :: popp, popq
    integer, intent(in) :: ioppq
    integer, dimension(8), intent(in) :: popr, pops
    integer, intent(in) :: ioprs
    integer, intent(in) :: sym !> irrep of V(pq,rs)

    integer :: irrep
    integer :: irreprs, irreppq
    integer :: numrs, numpq
    integer :: irrepp, irrepq, irrepr, irreps
    integer :: nump, numq, numr, nums
    integer :: p, q, r, s
    integer, dimension(8) :: ioff ! pointer of V(pq,rs) belong to irreprs
    integer, dimension(6) :: idxV ! pointer of V belong to irreprs,irreps,s,r,irrepq,q
    integer :: idx, i
    integer, external :: irpdso, isymoffso, idsymsoc

    i = idsymsoc(sym, popp, popq, ioppq, popr, pops, ioprs)
    if (idxin > i) then
        write (6, *) 'WARN: getidx have some error', idxin, i
    end if

    !--- determine irreprs
    ioff = 1
    do irreprs = 1, nirrep - 1
        irreppq = dirprd(irreprs, sym)
        numrs = irpdso(irreprs, popr, pops, ioprs)
        numpq = irpdso(irreppq, popp, popq, ioppq)
        ioff(irreprs + 1) = ioff(irreprs) + numpq*numrs
    end do
    do irrep = 1, nirrep
        if (idxin >= ioff(irrep)) then
            irreprs = irrep
        else
            exit
        end if
    end do
    idxV(1) = ioff(irreprs)
    irreppq = dirprd(irreprs, sym)
    numrs = irpdso(irreprs, popr, pops, ioprs)
    numpq = irpdso(irreppq, popp, popq, ioppq)

    !--- determine irreps
    do irrep = 1, nirrep
        idx = idxV(1) + numpq*(isymoffso(irrep, irreprs, popr, pops, ioprs) - 1)
        if (idxin >= idx) then
            irreps = irrep
        else
            exit
        end if
    end do
    idxV(2) = idxV(1) + numpq*(isymoffso(irreps, irreprs, popr, pops, ioprs) - 1)
    irrepr = dirprd(irreps, irreprs)
    numr = popr(irrepr)
    nums = pops(irreps)

    ! determine s and r
    if (ioprs == 0) then
        if (irrepr == irreps) then
            !--- s
            do i = 2, nums
                idx = idxV(2) + numpq*((i - 1)*(i - 2)/2)
                if (idxin >= idx) then
                    s = i
                else
                    exit
                end if
            end do
            idxV(3) = idxV(2) + numpq*((s - 1)*(s - 2)/2)
            !--- r
            do i = 1, s - 1
                idx = idxV(3) + numpq*(i - 1)
                if (idxin >= idx) then
                    r = i
                else
                    exit
                end if
            end do
            idxV(4) = idxV(3) + numpq*(r - 1)
        else ! irrepr<irreps
            !--- s
            do i = 1, nums
                idx = idxV(2) + numpq*numr*(i - 1)
                if (idxin >= idx) then
                    s = i
                else
                    exit
                end if
            end do
            idxV(3) = idxV(2) + numpq*numr*(s - 1)
            !--- r
            do i = 1, numr
                idx = idxV(3) + numpq*(i - 1)
                if (idxin >= idx) then
                    r = i
                else
                    exit
                end if
            end do
            idxV(4) = idxV(3) + numpq*(r - 1)
        end if
    else if (ioprs == 1) then
        !--- s
        do i = 1, nums
            idx = idxV(2) + numpq*numr*(i - 1)
            if (idxin >= idx) then
                s = i
            else
                exit
            end if
        end do
        idxV(3) = idxV(2) + numpq*numr*(s - 1)
        !--- r
        do i = 1, numr
            idx = idxV(3) + numpq*(i - 1)
            if (idxin >= idx) then
                r = i
            else
                exit
            end if
        end do
        idxV(4) = idxV(3) + numpq*(r - 1)
    else if (ioprs == 2) then
        write (6, *) "current not support ioprs=2"
        stop
    end if

    !--- determine irrepq
    do irrep = 1, nirrep
        idx = idxV(4) + isymoffso(irrep, irreppq, popp, popq, ioppq) - 1
        if (idxin >= idx) then
            irrepq = irrep
        else
            exit
        end if
    end do
    idxV(5) = idxV(4) + isymoffso(irrepq, irreppq, popp, popq, ioppq) - 1
    irrepp = dirprd(irrepq, irreppq)
    nump = popp(irrepp)
    numq = popq(irrepq)

    ! determine q and p
    if (ioppq == 0) then
        if (irrepp == irrepq) then
            !--- q
            do i = 2, numq
                idx = idxV(5) + (i - 1)*(i - 2)/2
                if (idxin >= idx) then
                    q = i
                else
                    exit
                end if
            end do
            idxV(6) = idxV(5) + (q - 1)*(q - 2)/2
            !--- p
            p = idxin - idxV(6) + 1
        else ! irrepp<irrepq
            !--- q
            do i = 1, numq
                idx = idxV(5) + nump*(i - 1)
                if (idxin >= idx) then
                    q = i
                else
                    exit
                end if
            end do
            idxV(6) = idxV(5) + nump*(q - 1)
            !--- p
            p = idxin - idxV(6) + 1
        end if
    else if (ioppq == 1) then
        !--- q
        do i = 1, numq
            idx = idxV(5) + nump*(i - 1)
            if (idxin >= idx) then
                q = i
            else
                exit
            end if
        end do
        idxV(6) = idxV(5) + nump*(q - 1)
        !--- p
        p = idxin - idxV(6) + 1
    else if (ioppq == 2) then
        write (6, *) "current not support ioppq=2"
        stop
    end if

    idxout(1) = irrepp
    idxout(2) = irrepq
    idxout(3) = irrepr
    idxout(4) = irreps
    idxout(5) = p
    idxout(6) = q
    idxout(7) = r
    idxout(8) = s

end subroutine getidx
