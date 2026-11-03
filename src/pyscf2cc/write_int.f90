! lish, 2026.4.5
subroutine init(dirprd_in, nirrep_in, path_in)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: dirprd_in(8, 8), nirrep_in
    character(len=200), intent(in) :: path_in

    dirprd = dirprd_in
    nirrep = nirrep_in
    path = path_in

end subroutine init

subroutine creat_file(symI, popp, popq, ioppq, popr, pops, ioprs, file, icore)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: symI ! symmetry I(pq,rs)
    integer*8, dimension(8), intent(in) :: popp, popq, popr, pops
    integer*8, intent(in) :: ioppq, ioprs
    character(len=8), intent(in) :: file
    real*8, dimension(*), intent(inout) :: icore

    ! create file
    call clmoioso0(symI, icore, popp, popq, ioppq, popr, pops, ioprs, file)

end subroutine creat_file

!> Write I(pq,rs) for a sepcified rs,
!> the I(pq,rs) need be adapted symmerty, and ioppq=1
subroutine write_block(len_v, vec, irrep_pq, symI, popp, popq, ioppq, file, rs, icore)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: len_v ! len of such block
    real*8, dimension(*), intent(in) :: vec ! I(pq,*) ?????????????????
    integer*8, intent(in) :: irrep_pq, symI ! symmetry of pq and I(pq,rs)
    integer*8, dimension(8), intent(in) :: popp, popq
    integer*8, intent(in) :: ioppq
    character(len=8), intent(in) :: file
    integer*8, intent(in) :: rs
    real*8, dimension(*), intent(inout) :: icore

    integer*8 :: i0, nsize
    integer*8 :: irrep_rs
    real*8 :: norm

    i0 = 1
    irrep_rs = dirprd(irrep_pq, symI)
    if (ioppq == 0) then
        ! compress pq
        call icompso(vec, icore(i0), popp, popq, 1, 0, zorb, zorb, 1, 1, irrep_pq)
        nsize = irpdso(irrep_pq, popp, popq, 0)
        call putlstso(icore(i0), rs, 1, irrep_rs, 1, file)
    else
        call putlstso(vec, rs, 1, irrep_rs, 1, file)
    end if

end subroutine write_block

!> Compress integer which ignore symmerty to symmetry adapted.
!> For example, sym_pq=1,
!> the original I(pq,*) be:
!>                 irrep_q
!>            1    2    3    4
!> irrep_p  -------------------
!>   1     | I11   0    0    0
!>   2     |  0   I22   0    0
!>   3     |  0    0   I33   0
!>   4     |  0    0    0   I44
!> the compressed integral be:
!>  (I11,I22,I33,I44)
subroutine compress_by_symmerty(n_p, n_q, n_pq, pop_p, pop_q, sym, v_in, v_out)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: n_p, n_q, n_pq ! total number of p,q and symmetry adapted I(pq,*)
    integer*8, dimension(8), intent(in) :: pop_p, pop_q
    integer*8, intent(in) :: sym !> symmetry of pq
    real*8, dimension(n_p*n_q), intent(in) :: v_in ! vector in
    real*8, dimension(n_pq), intent(out) :: v_out !> vector out

    integer*8 :: irrep_p, irrep_q, num_p, num_q, q
    integer*8 :: irrep
    integer*8 :: idx_in, idx_out

    if (n_p /= sum(pop_p)) then
        write (6, *) 'ERROR: n_p in compress_by_symmerty error'
        stop
    end if
    if (n_q /= sum(pop_q)) then
        write (6, *) 'ERROR: n_q in compress_by_symmerty error'
        stop
    end if
    if (n_pq /= irpdso(sym, pop_p, pop_q, 1)) then
        write (6, *) 'ERROR: n_pq in compress_by_symmerty error'
        stop
    end if

    ! cycle for irrep_q
    idx_in = 1 ! position of v_in for irrep_q
    do irrep_q = 1, nirrep
        irrep_p = dirprd(sym, irrep_q)
        num_p = pop_p(irrep_p)
        num_q = pop_q(irrep_q)
        idx_out = isymoffso(irrep_q, sym, pop_p, pop_q, 1) ! position of v_out for irrep_q
        ! cycle for q
        do q = 1, num_q
            idx_in = idx_in + sum(pop_p(1:irrep_p - 1))
            call dcopy(num_p, v_in(idx_in), 1, v_out(idx_out), 1)
            idx_in = idx_in + num_p + sum(pop_p(irrep_p + 1:nirrep))
            idx_out = idx_out + num_p
        end do
    end do

end subroutine compress_by_symmerty

!> calculate total number of pq adapted symmerty
function calc_num_of_pq(pop_p, pop_q, sym, iop_pq) result(n_pq)
    use mod_py2cc
    implicit none

    integer*8, dimension(8), intent(in) :: pop_p, pop_q
    integer*8, intent(in) :: sym !> symmetry of pq with iop=1
    integer*8, intent(in) :: iop_pq
    integer*8 :: n_pq

    n_pq = irpdso(sym, pop_p, pop_q, iop_pq)

end function calc_num_of_pq

!> calculate r and s without symmerty for a specified rs in I(pq,rs)
subroutine get_rs(rs, irrep_rs, pop_r, pop_s, ioprs, r, s)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: rs, irrep_rs, ioprs
    integer*8, dimension(8), intent(in) :: pop_r, pop_s
    integer*8, intent(out) :: r, s

    integer*8, dimension(8) :: idx
    integer*8 :: irrep_r, irrep_s, sym_s, sym_r

    call getidx(rs, idx, zorb, zorb, 1, pop_r, pop_s, ioprs, irrep_rs)
    irrep_r = idx(3)
    irrep_s = idx(4)
    sym_r = idx(7) ! r adapted symmerty
    sym_s = idx(8) ! s adapted symmerty
    r = sum(pop_r(1:irrep_r - 1)) + sym_r
    s = sum(pop_s(1:irrep_s - 1)) + sym_s

end subroutine get_rs

!> calculate symmetry adapted rs for a specified r and s in I(pq,rs)
subroutine put_rs(irrep_rs, pop_r, pop_s, irrep_r, irrep_s, ioprs, r, s, rs)
    use mod_py2cc
    implicit none

    integer*8, intent(in) :: irrep_rs, irrep_r, irrep_s
    integer*8, intent(in) :: ioprs
    integer*8, dimension(8), intent(in) :: pop_r, pop_s
    integer*8, intent(in) :: r !> symmerty adapted r
    integer*8, intent(in) :: s !> symmerty adapted s
    integer*8, intent(out) :: rs

    call putidx(rs, 1, 1, irrep_r, irrep_s, 1, 1, r, s, &
                zorb, zorb, 1, pop_r, pop_s, ioprs, irrep_rs)

end subroutine put_rs
