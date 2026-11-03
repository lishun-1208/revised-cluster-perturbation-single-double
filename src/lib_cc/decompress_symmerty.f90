! 2026.6.4
!> decompress a DPD formate vector to a no symmerty vector
subroutine decompress_DPD(in, out, len_out, popp, popq, ioppq, &
                          popr, pops, ioprs, sym, iCore, maxCore)
    use mod_symm
    implicit none

    real*8, dimension(*), intent(in) :: in !> DPD tensor
    real*8, dimension(len_out), intent(out):: out !> flat tensor
    integer*8, intent(in) :: len_out
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    integer*8, intent(in) :: sym
    real*8, dimension(maxCore), intent(inout) :: iCore !> free space
    integer*8, intent(in) :: maxCore

    target :: iCore, in
    real*8, dimension(:), pointer :: ptr
    integer*8 :: i0, i1, i2, nSize
    integer*8 :: irrepPQ, irrepRS
    integer*8 :: numpq, numrs
    integer*8, dimension(8) :: ioff = 0
    integer*8 :: irrepP, irrepQ, irrepR, irrepS
    integer*8 :: nump, numq, numr, nums ! total number for a given symmerty
    integer*8 :: p, q, r, s ! for a given symmerty
    integer*8, dimension(6) :: idx !> (irrepRS,irrepS,s,r,irrepQ,q)
    integer*8 :: idx_ns
    integer*8 :: p_ns, q_ns, r_ns, s_ns ! without out symmerty
    integer*8 :: np, nq, nr, ns ! total number without symmerty
    integer*8, external :: irpdso, idsymsoc, isymoffso

    np = sum(popp(1:nirrep))
    nq = sum(popq(1:nirrep))
    nr = sum(popr(1:nirrep))
    ns = sum(pops(1:nirrep))

    if (np*nq*nr*ns /= len_out) then
        write (6, *) 'ERROR: decompress_DPD error.', np*nq*nr*ns, len_out
        stop
    end if

    i0 = 1
    call dscal(len_out, 0.0d0, out, 1)

    !--- expand
    nSize = idsymsoc(sym, popp, popq, 1, popr, pops, 1)
    if (ioppq == 0) then
        i1 = i0 + nSize
        call iexpso(in, iCore(i0), popp, popq, 0, 1, popr, pops, ioprs, 0, sym)
        ptr => iCore(i0:i0 - 1 + nSize)
    else
        i1 = i0
        ptr => in(1:nSize)
    end if
    if (ioprs == 0) then
        call iexpso(ptr, iCore(i1), popp, popq, 1, 0, popr, pops, 0, 1, sym)
        ptr => iCore(i1:i1 - 1 + nSize)
    end if

    ioff(1) = 1
    do irrepRS = 1, nirrep - 1
        irrepPQ = dirprd(irrepRS, sym)
        numrs = irpdso(irrepRS, popr, pops, 1)
        numpq = irpdso(irrepPQ, popp, popq, 1)
        ioff(irrepRS + 1) = ioff(irrepRS) + numpq*numrs
    end do

    !------ decompress
    ! cycle for irrepRS
    do irrepRS = 1, nirrep
        irrepPQ = dirprd(irrepRS, sym)
        numpq = irpdso(irrepPQ, popp, popq, 1)
        numrs = irpdso(irrepRS, popr, pops, 1)
        if (numpq == 0 .or. numrs == 0) then
            cycle
        end if

        idx(1) = ioff(irrepRS)
        ! cycle for irrepS
        do irrepS = 1, nirrep
            irrepR = dirprd(irrepRS, irrepS)
            numr = popr(irrepR)
            nums = pops(irrepS)
            if (numr == 0 .or. nums == 0) then
                cycle
            end if

            idx(2) = idx(1) + numpq*(isymoffso(irrepS, irrepRS, popr, pops, 1) - 1)
            ! cycle for s
            do s = 1, nums
                s_ns = s + sum(pops(1:irrepS - 1))

                idx(3) = idx(2) + numpq*numr*(s - 1)
                ! cycle for r
                do r = 1, numr
                    r_ns = r + sum(popr(1:irrepR - 1))

                    idx(4) = idx(3) + numpq*(r - 1)
                    ! cycle for irrepQ
                    do irrepQ = 1, nirrep
                        irrepP = dirprd(irrepPQ, irrepQ)
                        nump = popp(irrepP)
                        numq = popq(irrepQ)
                        if (nump == 0 .or. numq == 0) then
                            cycle
                        end if

                        idx(5) = idx(4) + isymoffso(irrepQ, irrepPQ, popp, popq, 1) - 1
                        ! cycle for q
                        do q = 1, numq
                            q_ns = q + sum(popq(1:irrepQ - 1))

                            idx(6) = idx(5) + nump*(q - 1)
                            ! for p
                            p_ns = 1 + sum(popp(1:irrepP - 1))
                            idx_ns = p_ns + np*(q_ns - 1 + nq*(r_ns - 1 + nr*(s_ns - 1)))
                            call dcopy(nump, ptr(idx(6)), 1, out(idx_ns), 1)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine decompress_DPD

! 2026.7.4
!> compress a flat formate tensor to a DPD tensor
subroutine compress_DPD(in, out, len_in, popp, popq, ioppq, &
                        popr, pops, ioprs, sym, iCore, maxCore)
    use mod_symm
    implicit none

    real*8, dimension(len_in), intent(in) :: in !> flat tensor
    real*8, dimension(*), intent(out):: out !> DPD tensor
    integer*8, intent(in) :: len_in
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    integer*8, intent(in) :: sym
    real*8, dimension(maxCore), intent(inout) :: iCore !> free space
    integer*8, intent(in) :: maxCore

    target :: iCore, out
    real*8, dimension(:), pointer :: ptr
    integer*8 :: i0, i1, i2, nSize
    integer*8 :: irrepPQ, irrepRS
    integer*8 :: numpq, numrs
    integer*8, dimension(8) :: ioff = 0
    integer*8 :: irrepP, irrepQ, irrepR, irrepS
    integer*8 :: nump, numq, numr, nums ! total number for a given symmerty
    integer*8 :: p, q, r, s ! for a given symmerty
    integer*8, dimension(6) :: idx !> (irrepRS,irrepS,s,r,irrepQ,q)
    integer*8 :: idx_ns
    integer*8 :: p_ns, q_ns, r_ns, s_ns ! without out symmerty
    integer*8 :: np, nq, nr, ns ! total number without symmerty
    integer*8, external :: irpdso, idsymsoc, isymoffso

    np = sum(popp(1:nirrep))
    nq = sum(popq(1:nirrep))
    nr = sum(popr(1:nirrep))
    ns = sum(pops(1:nirrep))

    if (np*nq*nr*ns /= len_in) then
        write (6, *) 'ERROR: compress_DPD error.', np*nq*nr*ns, len_in
        stop
    end if

    i0 = 1

    if ((ioppq == 1) .and. (ioprs == 1)) then
        i1 = i0
        ptr(1:1) => out(1:1)
    else
        i1 = i0 + nSize
        ptr(1:1) => iCore(i0:i0)
    end if
    nSize = idsymsoc(sym, popp, popq, 1, popr, pops, 1)
    call dscal(nSize, 0.0d0, ptr, 1)

    ioff(1) = 1
    do irrepRS = 1, nirrep - 1
        irrepPQ = dirprd(irrepRS, sym)
        numrs = irpdso(irrepRS, popr, pops, 1)
        numpq = irpdso(irrepPQ, popp, popq, 1)
        ioff(irrepRS + 1) = ioff(irrepRS) + numpq*numrs
    end do

    !------ compress
    ! cycle for irrepRS
    do irrepRS = 1, nirrep
        irrepPQ = dirprd(irrepRS, sym)
        numpq = irpdso(irrepPQ, popp, popq, 1)
        numrs = irpdso(irrepRS, popr, pops, 1)
        if (numpq == 0 .or. numrs == 0) then
            cycle
        end if

        idx(1) = ioff(irrepRS)
        ! cycle for irrepS
        do irrepS = 1, nirrep
            irrepR = dirprd(irrepRS, irrepS)
            numr = popr(irrepR)
            nums = pops(irrepS)
            if (numr == 0 .or. nums == 0) then
                cycle
            end if

            idx(2) = idx(1) + numpq*(isymoffso(irrepS, irrepRS, popr, pops, 1) - 1)
            ! cycle for s
            do s = 1, nums
                s_ns = s + sum(pops(1:irrepS - 1))

                idx(3) = idx(2) + numpq*numr*(s - 1)
                ! cycle for r
                do r = 1, numr
                    r_ns = r + sum(popr(1:irrepR - 1))

                    idx(4) = idx(3) + numpq*(r - 1)
                    ! cycle for irrepQ
                    do irrepQ = 1, nirrep
                        irrepP = dirprd(irrepPQ, irrepQ)
                        nump = popp(irrepP)
                        numq = popq(irrepQ)
                        if (nump == 0 .or. numq == 0) then
                            cycle
                        end if

                        idx(5) = idx(4) + isymoffso(irrepQ, irrepPQ, popp, popq, 1) - 1
                        ! cycle for q
                        do q = 1, numq
                            q_ns = q + sum(popq(1:irrepQ - 1))

                            idx(6) = idx(5) + nump*(q - 1)
                            ! for p
                            p_ns = 1 + sum(popp(1:irrepP - 1))
                            idx_ns = p_ns + np*(q_ns - 1 + nq*(r_ns - 1 + nr*(s_ns - 1)))
                            call dcopy(nump, in(idx_ns), 1, ptr(idx(6)), 1)
                        end do
                    end do
                end do
            end do
        end do
    end do

    !--- compress
    if ((ioppq == 0) .and. (ioprs == 0)) then
        call icompso(ptr, iCore(i1), popp, popq, 1, 0, popr, pops, 1, 1, sym)
        call icompso(iCore(i1), out, popp, popq, 0, 1, popr, pops, 1, 0, sym)
    else if (ioppq == 0) then
        call icompso(ptr, out, popp, popq, 1, 0, popr, pops, 1, 1, sym)
    else if (ioprs == 0) then
        call icompso(ptr, out, popp, popq, 1, 1, popr, pops, 1, 0, sym)
    end if

end subroutine compress_DPD
