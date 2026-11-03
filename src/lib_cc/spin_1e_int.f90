!> Read spin adapted 1e integral form file, and combin to no spin form
subroutine read_combine_1e_int(Ins, f_IAA, f_IBB, pop_p, pop_q, pop_pa, pop_pb, pop_qa, pop_qb, &
                               sym, iCore, maxCore)
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(out) :: Ins
    character(len=8), intent(in) :: f_IAA, f_IBB
    integer*8, dimension(8), intent(in) :: pop_p, pop_q, pop_pa, pop_pb, pop_qa, pop_qb
    integer*8, intent(in) :: sym
    real*8, dimension(maxCore), intent(inout) :: iCore !> full work space
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, nsize
    integer, external :: irpdso

    i0 = 1

    nsize = irpdso(sym, pop_p, pop_q, 1)
    call dscal(nsize, 0.0d0, Ins, 1)

    call getallso(iCore(i0), 1, f_IAA)
    call combintegral(Ins, iCore(i0), pop_pa, pop_pb, pop_qa, pop_qb, 1, &
                      zorb, zorb, zorb, zorb, 1, sym, &
                      pop_p, pop_q, 1, zorb, zorb, 1, sym, 'AAAA')

    call getallso(iCore(i0), 1, f_IBB)
    call combintegral(Ins, iCore(i0), pop_pa, pop_pb, pop_qa, pop_qb, 1, &
                      zorb, zorb, zorb, zorb, 1, sym, &
                      pop_p, pop_q, 1, zorb, zorb, 1, sym, 'BBAA')

end subroutine read_combine_1e_int
