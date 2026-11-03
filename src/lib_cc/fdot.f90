! 2026.4.26
! modified at 2026.4.28
!---------------------------------------
!> calculate val=fact*f1(pq,rs)*f2(pq,rs)
function fdot(f1, f2, fact, iCore, maxCore) result(val)
    use mod_orbit
    use mod_symm
    use mod_operator
    implicit none

    type(opt), intent(in) :: f1, f2
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    real*8 :: val

    integer*8 :: irrep, numrs, rs
    integer*8 :: sym
    integer*8 :: i0, i1, nSize
    integer*8, external :: idsymsoc, irpdso
    real*8, external :: ddot

    if ((f1%len /= f2%len) .or. (f1%sym /= f2%sym)) then
        write (6, *) 'ERROR: f1 and f2 with different structure'
        write (6, *) 'f1:'
        call f1%print
        write (6, *) 'f2:'
        call f1%print
        stop
    end if

    i0 = 1
    val = 0.0d0
    ! cycle for symmetry of rs
    do irrep = 1, nirrep
        sym = dirprd(irrep, f1%sym) ! symmerty of pq
        numrs = irpdso(irrep, f1%popr, f1%pops, f1%ioprs)

        nSize = irpdso(sym, f1%popp, f1%popq, f1%ioppq)
        i1 = i0 + nSize
        do rs = 1, numrs
            call getlstso(iCore(i0), rs, 1, irrep, 1, f1%f)
            call getlstso(iCore(i1), rs, 1, irrep, 1, f2%f)
            val = val + fact*ddot(nSize, iCore(i0), 1, iCore(i1), 1)
        end do
    end do

end function fdot
