! copy form libsrint and lib_mod
module mod_py2cc
    implicit none
    ! mod_orb
!    integer*8 :: pop(8, 2) ! occ orbital for symmetry and spin
!    integer*8 :: vrt(8, 2) ! vrt orbital for symmetry and spin
!    integer*8 :: vrta(8), vrtb(8), popa(8), popb(8) ! virtual and occupied orbital for alpha and beta orbital
!    integer*8 :: nv(8), no(8) ! total virtual and occupied orbitals for different symmetry
    integer*8, parameter, dimension(8) :: zorb = (/1, 0, 0, 0, 0, 0, 0, 0/) ! zorb(1)=1, zorb(2:8)=0

    ! mod_size
!    integer*8 :: ilnbuf
    integer*8, parameter :: iprcln = 32768
    integer*8, parameter :: iprcwd = 4096
!    integer*8 :: icrsiz !icrsiz is the memory size for cpu, ichcsz=4kib=4096
    integer*8, parameter :: ichcsz = 4096
!    integer*8 :: igpusize, igpuleft !memory size of gpu and memory left in gpu for dynamic calculation

    ! mod_symm
!    integer*8 :: symRx, symRy, symRz ! symmetry of Rx, Ry, and Rz
    integer*8 :: dirprd(8, 8) ! direct product table
    integer*8 :: nirrep ! number of irreducible representation

    ! mod_iop
!    integer*8, parameter :: prec=8
    integer*8, parameter :: isd = 2
!    integer*8, parameter :: igpu
!    integer*8, parameter :: totgpu

    character(len=200) :: path ! path to tep file

    private inipckso, updmoiso, zero, wrdirso, rddirso, myicopy, sqsymso, sq2symso

contains
    integer*8 function irpdso(sym, popp, popq, ioppq)
        ! calculate the total number of pq belong to sym and with iop
        ! ioppq=0: p<q
        ! ioppq=1: p,q
        ! ioppq=2: p<=q
        implicit none

        integer*8, intent(in) :: sym
        integer*8, intent(in), dimension(8) :: popp, popq
        integer*8, intent(in) :: ioppq

        integer*8 :: irrep_q, irrep_p

        if (sym > 8) then
            write (6, *) 'ERROR: sym in irpdso error', sym
            stop
        end if

        irpdso = 0

        if (ioppq == 1) then
            do irrep_q = 1, nirrep
                irrep_p = dirprd(irrep_q, sym)
                irpdso = irpdso + popp(irrep_p)*popq(irrep_q)
            end do
        else if (ioppq == 0) then
            do irrep_q = 1, nirrep
                irrep_p = dirprd(irrep_q, sym)
                if (irrep_p > irrep_q) then
                    cycle
                else if (irrep_p < irrep_q) then
                    irpdso = irpdso + popp(irrep_p)*popq(irrep_q)
                else if (irrep_p == irrep_q) then
                    irpdso = irpdso + popp(irrep_p)*(popq(irrep_q) - 1)/2
                end if
            end do
        else if (ioppq == 2) then
            do irrep_q = 1, nirrep
                irrep_p = dirprd(irrep_q, sym)
                if (irrep_p > irrep_q) then
                    cycle
                else if (irrep_p < irrep_q) then
                    irpdso = irpdso + popp(irrep_p)*popq(irrep_q)
                else if (irrep_p == irrep_q) then
                    irpdso = irpdso + popp(irrep_p)*(popq(irrep_q) + 1)/2
                end if
            end do
        else
            write (6, *) 'ERROR: ioppq in irpdso error', ioppq
            stop
        end if

        return
    end function irpdso

    integer*8 function idsymsoc(irrepx, popp, popq, ioppq, &
                                popr, pops, ioprs)
        ! calculate total number of (pq,rs) belong to irrepx and with iop
        implicit none

        integer*8, intent(in) :: irrepx
        integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
        integer*8, intent(in) :: ioppq, ioprs

        integer*8 :: irreppq, irreprs

        if (irrepx > 8) then
            write (6, *) 'ERROR: irrepx in idsymsoc', irrepx
            stop
        end if

        idsymsoc = 0

        do irreprs = 1, nirrep
            irreppq = dirprd(irreprs, irrepx)
            idsymsoc = idsymsoc + &
                       irpdso(irreppq, popp, popq, ioppq)*irpdso(irreprs, popr, pops, ioprs)
        end do

        return
    end function idsymsoc

    subroutine clmoioso0(irrepx, icore, &
                         popp, popq, ioppqin, popr, pops, ioprsin, filename)
        ! create integral file
        implicit none

        integer*8, intent(in) :: irrepx, ioppqin, ioprsin
        integer*8, intent(in), dimension(8) :: popp, popq, popr, pops
        real*8, intent(inout) :: icore(*)
        character(len=8), intent(in) :: filename

        integer*8 :: moionuml, moionumr, rectmp
        COMMON/MOIOSO/MOIONUML(10, 20), MOIONUMR(10, 20), RECTMP(10, 20)
        integer*8 :: irrepr(8), irrepi(8)
        integer*8 :: ioppq(8), ioprs(8), i, j
        character(len=10) :: filename2

        moionuml = 0
        moionumr = 0
        rectmp = 0

        open (unit=412, file=trim(path)//trim(filename), &
              form='unformatted', &
              access='direct', recl=iprcln, status='replace')
        irrepr(1) = irrepx
        ioppq(1) = ioppqin
        ioprs(1) = ioprsin
        call inipckso(irrepx, icore, popp, popq, ioppq(1), &
                      popr, pops, ioprs(1), 1)
        close (412)

        filename2(1:8) = filename(1:8)
        filename2(9:10) = '.c'
        open (unit=413, file=trim(path)//trim(filename2), &
              form='unformatted', status='replace')
        rewind (413)
        write (413) irrepx
        write (413) (popp(i), i=1, 8)
        write (413) (popq(i), i=1, 8)
        write (413) (popr(i), i=1, 8)
        write (413) (pops(i), i=1, 8)
        write (413) (ioppq(i), i=1, 8)
        write (413) (ioprs(i), i=1, 8)
        write (413) (irrepr(i), i=1, 8)
        write (413) (irrepi(i), i=1, 8)
        write (413) ((moionuml(i, j), i=1, 10), j=1, 20)
        write (413) ((moionumr(i, j), i=1, 10), j=1, 20)
        write (413) ((rectmp(i, j), i=1, 10), j=1, 20)
        close (413)
        return
    end subroutine clmoioso0

    subroutine inipckso(irrepx, icore, popp, popq, ioppq, &
                        popr, pops, ioprs, list)
        implicit none
        integer*8, intent(in) :: ioppq, ioprs
        integer*8, intent(in) :: popp(8), popq(8)
        integer*8, intent(in) :: popr(8), pops(8)
        integer*8, intent(in) :: list, irrepx
        real*8, intent(inout) :: icore(*)
        integer*8 irrepr, irrepl, numl, numr

        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, irrepx)
            numl = irpdso(irrepl, popp, popq, ioppq)
            numr = irpdso(irrepr, popr, pops, ioprs)
            call updmoiso(numl, numr, irrepr, list, icore)
        end do
        return
    end subroutine inipckso

    subroutine updmoiso(numl, numr, irrepx, list, icore)
        implicit none
        integer*8 :: moionuml, moionumr, rectmp
        common/moioso/moionuml(10, 20), moionumr(10, 20), rectmp(10, 20)
        integer*8, intent(in) :: numl, numr, list, irrepx
        real*8, intent(inout) :: icore(*)
        real*8 itmp(ichcsz)
        integer*8 ifile
        integer*8 nwords
        integer*8 nfull, npart, totrec, totwrd, listrec
        integer*8 iwrite, iout, irecord

        ifile = 412
        moionuml(irrepx, list) = numl
        moionumr(irrepx, list) = numr

        nwords = numl*numr

        nfull = nwords/iprcwd
        npart = mod(nwords, iprcwd)

        totrec = nfull
        if (npart /= 0) totrec = totrec + 1

        totwrd = npart
        if (irrepx == nirrep) then
            rectmp(1, list + 1) = rectmp(irrepx, list) &
                                  + totrec
        else
            rectmp(irrepx + 1, list) = rectmp(irrepx, list) &
                                       + totrec
        end if
        listrec = rectmp(irrepx, list)

        call zero(itmp, ichcsz)
        iwrite = listrec + totrec
        do iout = 1, totrec
            irecord = iwrite - iout + 1
            call wrdirso(ifile, irecord, itmp, ichcsz)
        end do
        return
    end subroutine updmoiso

    subroutine zero(a, len)
        implicit none

        integer*8, intent(in) :: len
        real*8, dimension(len), intent(inout) :: a

        integer*8 :: i
        !$OMP PARALLEL DO
        do i = 1, len
            a(i) = 0.0d0
        end do
        !$OMP END PARALLEL DO
        return
    end subroutine zero

    subroutine wrdirso(iunit, irec, ivec, length)
        implicit none
        integer*8, intent(in) :: length, irec, iunit
        real*8, intent(inout) :: ivec(length)

        integer*8 :: irealunit, irealrec, ier

        if (irec .le. 0) return
        irealunit = iunit
        irealrec = irec
        write (irealunit, rec=irealrec, err=555, iostat=ier) ivec
        return
555     write (6, *) 'write wrdir error'
    end subroutine wrdirso

    subroutine putlstso(vec, rs, nrs, symrs, list, file)
        implicit none

        real*8, dimension(*), intent(in) :: vec ! vec use to store data read from file
        integer*8, intent(in) :: rs ! the first rs belong to symrs need read
        integer*8, intent(in) :: nrs ! total number of rs need read
        integer*8, intent(in) :: symrs ! symmetry of rs
        integer*8, intent(in) :: list ! file list, need be 1
        character(len=8), intent(in) :: file

        character(len=10) :: file_c
        integer*8 :: symI ! symmetry of I(pq,rs)
        integer*8, dimension(8) :: popp, popq, popr, pops ! population of p,q,r,s for symmetry
        integer*8, dimension(8) :: ioppq_8, ioprs_8
        integer*8 :: ioppq, ioprs
        integer*8 :: i
        integer*8, dimension(8) :: ioff ! position of I belong to irreprs
        integer*8, dimension(2) :: idx ! position of I belong to irreprs,rs
        integer*8 :: irrep
        integer*8 :: irreppq, numpq, pq
        integer*8 :: irreprs, numrs
        integer*8 :: totn
        integer*8 :: position ! position of bit in file

        ! read information of I(pq,rs)
        file_c = file//'.c'
        open (unit=413, file=trim(path)//trim(file_c), form='unformatted', status='old', action='read')
        rewind (413)
        read (413) symI
        read (413) (popp(i), i=1, 8)
        read (413) (popq(i), i=1, 8)
        read (413) (popr(i), i=1, 8)
        read (413) (pops(i), i=1, 8)
        read (413) (ioppq_8(i), i=1, 8)
        read (413) (ioprs_8(i), i=1, 8)
        close (413)
        ioppq = ioppq_8(1)
        ioprs = ioprs_8(1)

        ! give ioff
        ioff(1) = 0
        ! in file generate by libsrint, 4096s real*8 as a record
        do irrep = 1, nirrep - 1
            irreprs = irrep
            irreppq = dirprd(irreprs, symI)
            numrs = irpdso(irreprs, popr, pops, ioprs)
            numpq = irpdso(irreppq, popp, popq, ioppq)
            totn = (numpq*numrs - 1 + 4096)/4096
            ioff(irrep + 1) = ioff(irrep) + totn*4096
        end do

        ! read in data
        irreppq = dirprd(symrs, symI)
        numpq = irpdso(irreppq, popp, popq, ioppq)
        idx(1) = ioff(symrs)
        idx(2) = idx(1) + numpq*(rs - 1)
        open (unit=412, file=trim(path)//trim(file), form='unformatted', status='old', action='write', access='stream')
        position = idx(2)*8 + 1 ! real*8 need 8 byte
        write (412, pos=position) (vec(i), i=1, numpq*nrs)
        close (412)

    end subroutine putlstso

    subroutine rddirso(iunit, irec, ivec, length)
        implicit none
        integer*8, intent(in) :: length, irec, iunit
        real*8, intent(inout) :: ivec(length)

        integer*8 irealunit, irealrec, ier
        irealunit = iunit
        irealrec = irec
        read (irealunit, rec=irealrec, err=555, iostat=ier) ivec
        return
555     write (6, *) 'error in rddriso'
        stop
    end subroutine rddirso

    subroutine myicopy(i1, i2, len)
        implicit none
        real*8, intent(in) :: i1(len)
        real*8, intent(out) :: i2(len)
        integer*8, intent(in) :: len

        call dcopy(len, i1, 1, i2, 1)
        return
    end subroutine myicopy

    integer*8 function isymoffso(irrepq1, irrep, popp, popq, iop)
        implicit none

        integer*8, intent(in) :: irrepq1, irrep, iop
        integer*8, dimension(8), intent(in) :: popp, popq
        integer*8 :: irrepp, irrepq
        isymoffso = 1
        if (iop == 1) then
        do irrepq = 1, irrepq1 - 1
            irrepp = dirprd(irrepq, irrep)
            isymoffso = isymoffso + &
                        popp(irrepp)*popq(irrepq)
        end do
        elseif (iop == 0) then
        do irrepq = 1, irrepq1 - 1
            irrepp = dirprd(irrepq, irrep)
            if (irrepp < irrepq) then
                isymoffso = isymoffso + &
                            popp(irrepp)*popq(irrepq)
            elseif (irrepp == irrepq) then
                isymoffso = isymoffso + &
                            popp(irrepp)*(popp(irrepq) - 1)/2
            end if
        end do
        else if (iop == 2) then
        do irrepq = 1, irrepq1 - 1
            irrepp = dirprd(irrepq, irrep)
            if (irrepp < irrepq) then
                isymoffso = isymoffso + &
                            popp(irrepp)*popq(irrepq)
            else if (irrepp == irrepq) then
                isymoffso = isymoffso + &
                            popp(irrepp)*(popp(irrepq) + 1)/2
            end if
        end do
        end if
        return
    end function isymoffso

    subroutine icompso(iin, iout, popp, popq, ioppq, icpq, &
                       popr, pops, ioprs, icrs, sym)
        implicit none

        integer*8, intent(in) :: popp(8), popq(8), popr(8), pops(8)
        integer*8, intent(in) :: ioppq, icpq, ioprs, icrs
        integer*8, intent(in) :: sym
        real*8, intent(in) :: iin(*)
        !real*4, intent(in) :: iin(*)
        real*8, intent(out) :: iout(*)
        !real*4, intent(out) :: iout(*)
        real*8, allocatable :: itmp(:)
        !real*4, allocatable :: itmp(:)
        integer*8 irrepr, irrepl, i, numl, numr, numl2, numr2, popl
        integer*8 nsize1, nsize2, ioff, ioff1, nsize
        if (ioppq /= 0 .and. ioppq /= 1 .and. ioppq /= 2) then
            write (6, *) 'error of ioppq in isymtrso'
            stop
        end if
        if (ioprs /= 0 .and. ioprs /= 1 .and. ioprs /= 2) then
            write (6, *) 'error of ioprs in isymtrso'
            stop
        end if
        if (icpq == 0 .and. icrs == 1) then
            ioff = 1
            ioff1 = 1
            do irrepr = 1, nirrep
                irrepl = dirprd(irrepr, sym)
                numr = irpdso(irrepr, popr, pops, ioprs)
                numl = irpdso(irrepl, popp, popq, ioppq)
                numl2 = irpdso(irrepl, popp, popq, icpq)
                call sqsymso(irrepl, popp, numl2, numl, numr, &
                             iout(ioff1), iin(ioff))
                ioff = ioff + numl*numr
                !ioff = ioff + numl*numr*isd
                ioff1 = ioff1 + numl2*numr
                !ioff1 = ioff1 + numl2*numr*isd
            end do
        elseif (icpq == 1 .and. icrs == 0) then
            ioff = 1
            ioff1 = 1
            do irrepr = 1, nirrep
                irrepl = dirprd(irrepr, sym)
                numl = irpdso(irrepl, popp, popq, ioppq)
                numr = irpdso(irrepr, popr, pops, ioprs)
                numr2 = irpdso(irrepr, popr, pops, icrs)
                call sq2symso(irrepr, popr, numl, iout(ioff1), iin(ioff))
                ioff = ioff + numl*numr
                !ioff = ioff + numl*numr*isd
                ioff1 = ioff1 + numl*numr2
                !ioff1 = ioff1 + numl*numr2*isd
            end do
        else
            write (6, *) 'attention!!!, error in icompso'
            stop
        end if
        return
    end subroutine icompso

    subroutine sqsymso(IRREPPQ, POPP, NUMPQ2, &
                       NUMPQ, NUMRS, IOUT, IIN)
        use omp_lib
        implicit none

        integer*8, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ
        integer*8, intent(in) :: NUMRS
        real*8, intent(in) :: iin(*)
        real*8, intent(out) :: iout(*)
        integer*8 nsize1, nsize2
        integer*8 ioffo(8), ioffn(8), vstarto, vstartn
        integer*8 vstartn2
        integer*8 irrepp, irrepq, irrepr1, irreps
        integer*8 irreprs
        integer*8 nump, numq, numr, nums
        integer*8 p, q, r, s, pq, rs
        integer*8 vstartso, vstartqo
        integer*8 vstarts, vstartr, vstartp, vstartq
        integer*8 idxin, idxout, idxout2
        integer*8 numout, numout2, numin, numout3, numout4
        ioffo(1) = 0; ioffn(1) = 0
        do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, irreppq)
            if (irrepp > irrepq) cycle
            if (irrepp == irrepq) then
                numq = popp(irrepq)
                nump = popp(irrepp)
                vstarto = isymoffso(irrepq, irreppq, popp, popp, 1)
                vstartn = isymoffso(irrepq, irreppq, popp, popp, 0)
                !$OMP PARALLEL DO

                do rs = 1, numrs
                do q = 2, numq
                do p = 1, q - 1
                    numin = (rs - 1)*numpq + vstarto - 1 + (q - 1)*nump + p
                    numout = (rs - 1)*numpq2 + vstartn - 1 + (q - 1)*(q - 2)/2 + p
                    iout(numout) = iin(numin)
                end do
                end do
                end do
                !$OMP END PARALLEL DO
            elseif (irrepp < irrepq) then
                numq = popp(irrepq)
                nump = popp(irrepp)
                vstarto = isymoffso(irrepq, irreppq, popp, popp, 1)
                vstartn = isymoffso(irrepq, irreppq, popp, popp, 0)
                !$OMP PARALLEL DO
                do rs = 1, numrs
                do q = 1, numq
                do p = 1, nump
                    numin = (rs - 1)*numpq + vstarto - 1 + (q - 1)*nump + p
                    numout = (rs - 1)*numpq2 + vstartn - 1 + (q - 1)*nump + p
                    iout(numout) = iin(numin)
                end do
                end do
                end do
                !$OMP END PARALLEL DO
            end if
        end do
        return
    end subroutine sqsymso

    subroutine sq2symso(sym, POPR, NUMPQ, IOUT, IIN)
        use omp_lib
        implicit none

        integer*8, intent(in) :: popr(8)
        integer*8, intent(in) :: sym, numpq
        real*8, intent(in) :: iin(*)
        real*8, intent(out) :: iout(*)
        integer*8 irrepr, irrepl, i, numl, numl2, numr2, popl
        integer*8 nsize1, nsize2
        integer*8 ioffo(8), ioffn(8), vstarto, vstartn
        integer*8 vstartn2
        integer*8 irrepp, irrepq, irrepr1, irreps
        integer*8 irreprs, irreppq, numrs, numrs2, numpq2
        integer*8 nump, numq, numr, nums
        integer*8 p, q, r, s, pq, rs
        integer*8 vstartso, vstartqo
        integer*8 vstarts, vstartr, vstartp, vstartq
        integer*8 idxin, idxout, idxout2
        integer*8 numout, numout2, numin, numout3, numout4
        irreprs = sym
        do irreps = 1, nirrep
            irrepr = dirprd(irreprs, irreps)
            if (irrepr > irreps) cycle
            if (irrepr == irreps) then
                numr = popr(irrepr)
                nums = popr(irreps)
                vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
                vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
                !$OMP PARALLEL DO
                do s = 2, nums
                do r = 1, s - 1
                do pq = 1, numpq
                    numin = (vstarto - 1 + (s - 1)*numr + r - 1)*numpq + pq
                    numout = (vstartn - 1 + (s - 1)*(s - 2)/2 + r - 1)*numpq + pq
                    iout(numout) = iin(numin)
                end do
                end do
                end do
                !$OMP END PARALLEL DO
            elseif (irrepr < irreps) then
                numr = popr(irrepr)
                nums = popr(irreps)
                vstarto = isymoffso(irreps, irreprs, popr, popr, 1)
                vstartn = isymoffso(irreps, irreprs, popr, popr, 0)
                !$OMP PARALLEL DO
                do s = 1, nums
                do r = 1, numr
                do pq = 1, numpq
                    numin = (vstarto - 1 + (s - 1)*numr + r - 1)*numpq + pq
                    numout = (vstartn - 1 + (s - 1)*numr + r - 1)*numpq + pq
                    iout(numout) = iin(numin)
                end do
                end do
                end do
                !$OMP END PARALLEL DO
            end if
        end do
        return
    end subroutine sq2symso

    !> calculate the (irrepP,irrepQ,irrepR,irrepS,p,q,r,s) for a specified V(pq,rs)
    subroutine getidx(idxin, idxout, popp, popq, ioppq, popr, pops, ioprs, sym)
        implicit none

        integer*8, intent(in) :: idxin ! the idx of speciffied V(pq,rs) in 1d array
        integer*8, dimension(8), intent(out) :: idxout !> (irrepp,irrepq,irrepr,irreps,p,q,r,s)
        integer*8, dimension(8), intent(in) :: popp, popq
        integer*8, intent(in) :: ioppq
        integer*8, dimension(8), intent(in) :: popr, pops
        integer*8, intent(in) :: ioprs
        integer*8, intent(in) :: sym ! irrep of V(pq,rs)

        integer*8 :: irrep
        integer*8 :: irreprs, irreppq
        integer*8 :: numrs, numpq
        integer*8 :: irrepp, irrepq, irrepr, irreps
        integer*8 :: nump, numq, numr, nums
        integer*8 :: p, q, r, s
        integer*8, dimension(8) :: ioff ! pointer of V(pq,rs) belong to irreprs
        integer*8, dimension(6) :: idxV ! pointer of V belong to irreprs,irreps,s,r,irrepq,q
        integer*8 :: idx, i

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

    subroutine dcopy(n, a, ia, b, ib)
        implicit none
        integer*8, intent(in) :: n, ia, ib
        real*8, dimension(n), intent(in) :: a
        real*8, dimension(n), intent(out) :: b

        integer*8 :: i
        !$OMP PARALLEL DO
        do i = 1, n
            b(i) = a(i)
        end do
        !$OMP END PARALLEL DO
    end subroutine dcopy

!    function dnrm2(n, a, inca) result(norm)
!        integer*8, intent(in) :: n, inca
!        real*8, dimension(n), intent(in) :: a
!        real*8 :: norm
!
!        integer*8 :: i
!
!        norm = 0.0d0
!        do i = 1, n
!            norm = norm + a(i)**2
!        end do
!        norm = sqrt(norm)
!    end function dnrm2

    ! write by LiShun, 2024.9.20
    !> calculate index of V(pq,rs) for specified (irrepp,iirepq,irrepr,irreps,p,q,r,s)
    ! if one index is none, ie. s, need set pops=zorb,irreps=1,s=1
    !----------------------------------
    ! modified at 2024.11.13
    subroutine putidx(idxout, irrepp, irrepq, irrepr, irreps, p, q, r, s, &
                      popp, popq, ioppq, popr, pops, ioprs, symV)
        implicit none

        integer*8, intent(out) :: idxout !> the idx of V(pq,rs) in 1d array
        integer*8, intent(in) :: irrepp, irrepq, irrepr, irreps
        integer*8, intent(in) :: p, q, r, s
        integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
        integer*8, intent(in) :: ioppq, ioprs
        integer*8, intent(in) :: symV !> irrep of V(pq,rs)

        integer*8 :: irreppq, irreprs
        integer*8 :: numpq, numrs
        integer*8 :: numr, nump
        integer*8, dimension(8) :: ioff ! idx for V belong to irreprs
        integer*8, dimension(7) :: idx ! idx for V belong to irreprs,irreps,s,r,irrepq,q,p

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

end module mod_py2cc
