! 2026.4.28
!> I(pq,rs)=I(pQ,rS)-I(Qp,Rs)=I(pQ,rS)-I(qP,rS)
!> need abab and aaaa with the same IO state
subroutine abab2aaaa(Iabab, Iaaaa, type, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(in) :: Iabab
    type(opt), intent(inout) :: Iaaaa
    character(len=2), intent(in) :: type !> 'pq' or 'rs', exchange pq or rs
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (Iabab%incore) then
        call p_abab2aaaa(Iabab%arr, Iaaaa%arr, Iaaaa%popp, Iaaaa%popq, Iaaaa%ioppq, &
                         Iaaaa%popr, Iaaaa%pops, Iaaaa%ioprs, Iaaaa%sym, &
                         iCore, type)
    else
        call p_abab2aaaa_h(Iabab%f, Iaaaa%f, Iaaaa%popp, Iaaaa%popq, Iaaaa%ioppq, &
                           Iaaaa%popr, Iaaaa%pops, Iaaaa%ioprs, Iaaaa%sym, &
                           iCore, type)
    end if

end subroutine abab2aaaa

! write by LiShun, 2024.5.20
! for closed shell:
! I(pq,rs)=I(pQ,rS)-I(Qp,Rs)=I(pQ,rS)-I(qP,rS)
! need alpha and beta have the same pop
!---------------------------------
! modified at 2025.11.16
subroutine p_abab2aaaa(Iabab, Iaaaa, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, sym, iCore, type)
    implicit none

    real*8, dimension(*), intent(in) :: Iabab
    real*8, dimension(*), intent(out) :: Iaaaa
    integer*8, dimension(8), intent(in) :: pop_p, pop_q, pop_r, pop_s
    integer*8, intent(in) :: ioppq, ioprs !> for I(pq,rs)
    integer*8, intent(in) :: sym
    real*8, dimension(*), intent(in) :: iCore
    character(len=2), intent(in) :: type !> 'pq'/'rs', exchange pq or rs

    integer*8 :: i0, i1, i2, nSize, iI
    integer*8, external :: idsymsoc

    i0 = 1
    ! I(pq,rs)=I(pQ,rS), i0
    nSize = idsymsoc(sym, pop_p, pop_q, 1, pop_r, pop_s, 1)
    i1 = i0 + nSize
    call dcopy(nSize, Iabab, 1, iCore(i0), 1)
    ! exchange pq or rs
    if (type == 'pq') then
        ! I(pQ,rS)->I(qP,rS), i1
        call isymtrso(Iabab, iCore(i1), pop_p, pop_q, 1, 1, pop_r, pop_s, 1, 0, sym)
    else if (type == 'rs') then
        ! I(pQ,rS)->I(pQ,sR), i1
        call isymtrso(Iabab, iCore(i1), pop_p, pop_q, 1, 0, pop_r, pop_s, 1, 1, sym)
    end if
    ! I(pq,rs)=I(pQ,rS)-I(qP,rS), i0
    call daxpy(nSize, -1.0d0, iCore(i1), 1, iCore(i0), 1)
    iI = i0
    ! compress pq and rs
    if (ioppq == 0 .and. ioprs == 0) then
        call icompso(iCore(i0), iCore(i1), pop_p, pop_q, 1, 1, pop_r, pop_s, 1, 0, sym)
        call icompso(iCore(i1), iCore(i0), pop_p, pop_q, 1, 0, pop_r, pop_s, 0, 1, sym)
        iI = i0
    else if (ioppq == 0 .and. ioprs == 1) then
        call icompso(iCore(i0), iCore(i1), pop_p, pop_q, 1, 0, pop_r, pop_s, 1, 1, sym)
        iI = i1
    else if (ioppq == 1 .and. ioprs == 0) then
        call icompso(iCore(i0), iCore(i1), pop_p, pop_q, 1, 1, pop_r, pop_s, 1, 0, sym)
        iI = i1
    end if
    nSize = idsymsoc(sym, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs)
    call dcopy(nSize, iCore(iI), 1, Iaaaa, 1)

end subroutine p_abab2aaaa

! 2025.12.16
subroutine p_abab2aaaa_h(f_abab, f_aaaa, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, sym, iCore, type)
    use mod_symm
    use mod_orbit
    implicit none

    character(len=8), intent(in) :: f_abab, f_aaaa
    integer*8, dimension(8), intent(in) :: pop_p, pop_q, pop_r, pop_s
    integer*8, intent(in) :: ioppq, ioprs !> for I(pq,rs)
    integer*8, intent(in) :: sym
    real*8, dimension(*), intent(in) :: iCore
    character(len=2), intent(in):: type !> 'pq'/'rs', exchange pq or rs

    integer*8 :: i0, i1, i2, nSize
    integer*8 :: irreprs, numrs, rs
    integer*8 :: irreppq
    integer*8 :: rsAB
    integer*8, dimension(8) :: idx
    integer*8, external :: idsymsoc, irpdso

    if (type == 'rs') then
        write (6, *) 'ERROR: p_abab2aaaa_h current only support type=pq'
        stop
    end if

    i0 = 1
    call clmoioso0(sym, iCore(i0), pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, f_aaaa)
    if (type == 'pq') then
        !--- cycle for rsAA
        do irreprs = 1, nirrep
            irreppq = dirprd(sym, irreprs)
            numrs = irpdso(irreprs, pop_r, pop_s, ioprs)
            nSize = irpdso(irreppq, pop_p, pop_q, 1)
            do rs = 1, numrs
                call getidx(rs, idx, zorb, zorb, 1, pop_r, pop_s, ioprs, irreprs)
                call putidx(rsAB, 1, 1, idx(3), idx(4), 1, 1, idx(7), idx(8), &
                            zorb, zorb, 1, pop_r, pop_s, 1, irreprs)
                ! read V(pQ,rS), i0
                i1 = i0 + nSize
                call getlstso(iCore(i0), rsAB, 1, irreprs, 1, f_abab)
                ! V(pQ,rS)->V(qP,rS), i1
                i2 = i1 + nSize
                call isymtrso(iCore(i0), iCore(i1), pop_p, pop_q, 1, 1, zorb, zorb, 1, 0, irreppq)
                ! V(pq,rs)=V(pQ,rS)-V(qP,rS), i0
                call daxpy(nSize, -1.0d0, iCore(i1), 1, iCore(i0), 1)

                if (ioppq == 0) then
                    ! compress pq
                    call icompso(iCore(i0), iCore(i1), pop_p, pop_q, 1, 0, zorb, zorb, 1, 1, irreppq)
                    call putlstso(iCore(i1), rs, 1, irreprs, 1, f_aaaa)
                else if (ioppq == 1) then
                    call putlstso(iCore(i0), rs, 1, irreprs, 1, f_aaaa)
                end if
            end do
        end do
    end if

end subroutine p_abab2aaaa_h
