! 2026.4.23
!> combine two 1e operator to 2e operator
!> V(pq,rs) += fact*A(pq)*B(rs)
subroutine combine_1e_to_2e(out, in1, in2, fact, iCore, maxCore)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: out !> V
    type(opt), intent(in) :: in1 !> A
    type(opt), intent(in) :: in2 !> B
    real*8, intent(in) :: fact !> A*B
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if ((out%ioppq /= 1) .or. (out%ioprs /= 1)) then
        write (6, *) 'ERROR: combine_1e_to_2e need ioppq=ioprs=1'
        stop
    end if

    call p_combine_1e_to_2e(out%arr, in1%arr, in2%arr, out%popp, out%popq, out%popr, out%pops, &
                            in1%sym, in2%sym, fact)

end subroutine combine_1e_to_2e

subroutine p_combine_1e_to_2e(out, in1, in2, pop_p, pop_q, pop_r, pop_s, sym1, sym2, fact)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: out
    real*8, dimension(*), intent(in) :: in1, in2
    integer*8, dimension(8), intent(in) :: pop_p, pop_q, pop_r, pop_s
    integer*8, intent(in) :: sym1, sym2
    real*8, intent(in) :: fact

    integer*8 :: sym

    sym = dirprd(sym1, sym2) ! symmetry of V=A*B
    ! V(pq,rs) = A(pq,**)*B(rs,**)
    call VpqmnVrsmn_to_Vpqrs(in1, in2, sym1, sym2, zorb, zorb, 1, &
                             pop_p, pop_q, 1, pop_r, pop_s, 1, out, sym, &
                             fact, 1.0d0)

end subroutine p_combine_1e_to_2e
