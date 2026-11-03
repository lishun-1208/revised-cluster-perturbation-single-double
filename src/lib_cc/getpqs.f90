! write by lishun, 2024.10.30
! modified at 2026.7.14
!------------------------------------------
!> Get I(pq,rs) for given r from hard disk
! will set iop of vecout the same to file
subroutine getpqs(vecout, r, irrepr, popp, popq, ioppq, popr, pops, ioprs, &
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
    integer*8 :: i0, nSize
    integer*8, external :: irpdso, isymoffso, idsymsoc

    i0 = 1

    if (ioprs == 1) then
        idxout = 1
        do irreps = 1, nirrep
            irreprs = dirprd(irrepr, irreps)
            irreppq = dirprd(irreprs, symV)
            numpq = irpdso(irreppq, popp, popq, ioppq)
            numr = popr(irrepr)
            nums = pops(irreps)
            if (r > numr) then
                write (6, '(A,2I4)') 'ERRPR: r in getpqs error.', r, numr
                stop
            end if

            rs(1) = isymoffso(irreps, irreprs, popr, pops, ioprs) - 1
            do s = 1, nums
                rs(2) = rs(1) + numr*(s - 1) + r
                call getlstso(vecout(idxout), rs(2), 1, irreprs, 1, file)
                idxout = idxout + numpq
            end do
        end do
    else if ((ioprs == 0)) then
        ! I(pq,rs) = -I(pq,sr)
        call getpqrso(vecout, r, irrepr, 1, file)
        nSize = idsymsoc(dirprd(symV, irrepr), popp, popq, ioppq, pops, zorb, 1)
        call dscal(nSize, -1.0d0, vecout, 1)
    else
        write (6, '(A.I2)') 'ERROR: getpqs current only support ioprs=0/1.', ioprs
        stop
    end if

end subroutine getpqs

! 2026.7.10
! modified at 2026.7.14
!> Get I(pq,rs) for a given r in memory
subroutine getpqs_incore(vecout, r, irrepr, popp, popq, ioppq, popr, pops, ioprs, &
                         vecin, symV)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(out) :: vecout
    integer*8, intent(in) :: r !> r belong to irrepr
    integer*8, intent(in) :: irrepr
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    real*8, dimension(*), intent(in) :: vecin
    integer*8, intent(in) :: symV

    integer*8 :: irreppq, numpq
    integer*8 :: irreprs, numrs
    integer*8, dimension(2) :: rs
    integer*8 :: irreps, nums, s
    integer*8 :: numr
    integer*8 :: idxout ! V(pq,*s) belong to irreps and s
    integer*8 :: i0, nSize
    integer*8, external :: irpdso, isymoffso, idsymsoc

    i0 = 1
    if (ioprs == 1) then
        idxout = 1
        do irreps = 1, nirrep
            irreprs = dirprd(irrepr, irreps)
            irreppq = dirprd(irreprs, symV)
            numpq = irpdso(irreppq, popp, popq, ioppq)
            numr = popr(irrepr)
            nums = pops(irreps)
            if (r > numr) then
                write (6, '(A,2I4)') 'ERRPR: r in getpqs_incore error.', r, numr
                stop
            end if

            rs(1) = isymoffso(irreps, irreprs, popr, pops, ioprs) - 1
            do s = 1, nums
                rs(2) = rs(1) + numr*(s - 1) + r
                call getlstso_incore(vecout(idxout), rs(2), 1, irreprs, &
                                     popp, popq, ioppq, popr, pops, ioprs, &
                                     vecin, symV)
                idxout = idxout + numpq
            end do
        end do
    else if (ioprs == 0) then
        ! I(pq,rs) = -I(pq,sr)
        call getpqrso_incore(vecout, r, irrepr, &
                             popp, popq, ioppq, popr, pops, ioprs, &
                             vecin, symV)
        nSize = idsymsoc(dirprd(symV, irrepr), popp, popq, ioppq, pops, zorb, 1)
        call dscal(nSize, -1.0d0, vecout, 1)
    else
        write (6, '(A.I2)') 'ERROR: getpqs_incore current only support ioprs=0/1.', ioprs
        stop
    end if

end subroutine getpqs_incore

! 2026.7.10
! modified at 2026.7.14
!> Put I(pq,rs) for a given r in memory
subroutine putpqs_incore(vecin, r, irrepr, popp, popq, ioppq, popr, pops, ioprs, &
                         vecout, symV)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: vecin
    integer*8, intent(in) :: r !> r belong to irrepr
    integer*8, intent(in) :: irrepr
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    real*8, dimension(*), intent(inout) :: vecout
    integer*8, intent(in) :: symV

    integer*8 :: irreppq, numpq
    integer*8 :: irreprs, numrs
    integer*8, dimension(2) :: rs
    integer*8 :: irreps, nums, s
    integer*8 :: numr
    integer*8 :: idxin ! V(pq,*s) belong to irreps and s
    integer*8 :: i0, nSize
    integer*8, external :: irpdso, isymoffso, idsymsoc

    i0 = 1
    if (ioprs == 1) then
        idxin = 1
        do irreps = 1, nirrep
            irreprs = dirprd(irrepr, irreps)
            irreppq = dirprd(irreprs, symV)
            numpq = irpdso(irreppq, popp, popq, ioppq)
            numr = popr(irrepr)
            nums = pops(irreps)
            if (r > numr) then
                write (6, '(A,2I4)') 'ERRPR: r in putpqs_incore error.', r, numr
                stop
            end if

            rs(1) = isymoffso(irreps, irreprs, popr, pops, ioprs) - 1
            do s = 1, nums
                rs(2) = rs(1) + numr*(s - 1) + r
                call putlstso_incore(vecin(idxin), rs(2), 1, irreprs, &
                                     popp, popq, ioppq, popr, pops, ioprs, &
                                     vecout, symV)
                idxin = idxin + numpq
            end do
        end do
    else if (ioprs == 0) then
        ! I(pq,rs)=-I(pq,sr)
        nSize = idsymsoc(dirprd(symV, irrepr), popp, popq, ioppq, pops, zorb, 1)
        call dscal(nSize, -1.0d0, vecin, 1)
        call putpqrso_incore(vecin, r, irrepr, &
                             popp, popq, ioppq, popr, pops, ioprs, &
                             vecout, symV)
        call dscal(nSize, -1.0d0, vecin, 1)
    else
        write (6, '(A.I2)') 'ERROR: putpqs_incore current only support ioprs=0/1.', ioprs
        stop
    end if

end subroutine putpqs_incore

! 2026.7.10
! modified at 2026.7.14
!> Put I(pq,rs) for a given r in hard disk
subroutine putpqs(vecin, r, irrepr, popp, popq, ioppq, popr, pops, ioprs, &
                  file, symV)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: vecin
    integer*8, intent(in) :: r !> r belong to irrepr
    integer*8, intent(in) :: irrepr
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    character(len=8), intent(in) :: file
    integer*8, intent(in) :: symV

    integer*8 :: irreppq, numpq
    integer*8 :: irreprs, numrs
    integer*8, dimension(2) :: rs
    integer*8 :: irreps, nums, s
    integer*8 :: numr
    integer*8 :: idxin ! V(pq,*s) belong to irreps and s
    integer*8 :: i0, nSize
    integer*8, external :: irpdso, isymoffso, idsymsoc

    i0 = 1
    if (ioprs == 1) then
        idxin = 1
        do irreps = 1, nirrep
            irreprs = dirprd(irrepr, irreps)
            irreppq = dirprd(irreprs, symV)
            numpq = irpdso(irreppq, popp, popq, ioppq)
            numr = popr(irrepr)
            nums = pops(irreps)
            if (r > numr) then
                write (6, '(A,2I4)') 'ERRPR: r in putpqs error.', r, numr
                stop
            end if

            rs(1) = isymoffso(irreps, irreprs, popr, pops, ioprs) - 1
            do s = 1, nums
                rs(2) = rs(1) + numr*(s - 1) + r
                call putlstso(vecin(idxin), rs(2), 1, irreprs, 1, file)
                idxin = idxin + numpq
            end do
        end do
    else if (ioprs == 0) then
        ! I(pq,rs)=-I(pq,sr)
        nSize = idsymsoc(dirprd(symV, irrepr), popp, popq, ioppq, pops, zorb, 1)
        call dscal(nSize, -1.0d0, vecin, 1)
        call putpqrso(vecin, r, irrepr, 1, file)
        call dscal(nSize, -1.0d0, vecin, 1)
    else
        write (6, '(A.I2)') 'ERROR: putpqs current only support ioprs=0/1.', ioprs
        stop
    end if

end subroutine putpqs
