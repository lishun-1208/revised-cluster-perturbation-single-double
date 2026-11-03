! 2026.7.13
!> c(ab,ij) -= fact*W(ij,kd)*U(ba,dk)
subroutine c2d_WijkdUdkba(c2d, Wijkd, Udkba, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: Wijkd
    type(opt), intent(in) :: Udkba
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp, tmp2
    integer*8 :: i0, i1, i2

    i0 = 1
    if (Udkba%incore) then
        ! W(ij,kd)->W(ji,dk), i0
        call tmp%init(iCore, Wijkd%popq, Wijkd%popp, Wijkd%pops, Wijkd%popr, &
                      Wijkd%ioppq, 1, Wijkd%sym, T, i0)
        i1 = i0 + tmp%len
        if (Wijkd%ioppq == 1) then
            call isymtrso(Wijkd%arr, tmp%arr, Wijkd%popp, Wijkd%popq, 1, 1, &
                          Wijkd%popr, Wijkd%pops, 1, 1, Wijkd%sym)
        else if (Wijkd%ioppq == 0) then
            call isymtrso(Wijkd%arr, tmp%arr, Wijkd%popp, Wijkd%popq, 0, 0, &
                          Wijkd%popr, Wijkd%pops, 1, 1, Wijkd%sym)
            call dscal(tmp%len, -1.0d0, tmp%arr, 1)
        end if
        ! tmp2(ba,ji) = -fact*tmp(ji,dk)*U(ba,dk)
        call tmp2%init(iCore, c2d%popq, c2d%popp, c2d%pops, c2d%popr, &
                       c2d%ioppq, c2d%ioprs, c2d%sym, T, i1)
        i2 = i1 + tmp2%len
        call dscal(tmp2%len, 0.0d0, tmp2%arr, 1)
        call vv2v_444(Udkba, tmp, tmp2, -fact, 'badk', 'jidk', 'baji', iCore, maxCore, 0)
        ! tmp2(ba,ji) -> tmp2(ab,ij)
        if ((tmp2%ioppq == 0) .and. (tmp2%ioprs == 0)) then
            call daxpy(tmp2%len, 1.0d0, tmp2%arr, 1, c2d%arr, 1)
        else if ((tmp2%ioppq == 1) .and. (tmp2%ioprs == 1)) then
            call isymtrso(tmp2%arr, iCore(i2), tmp2%popp, tmp2%popq, 1, 1, &
                          tmp2%popr, tmp2%pops, 1, 1, tmp2%sym)
            call daxpy(tmp2%len, 1.0d0, iCore(i2), 1, c2d%arr, 1)
        else
            write (6, *) 'ERROR: c2d_WijkdUdkba need ioppq=ioprs.'
            stop
        end if
    else
        call p_c2d_WijkdUdkba_h(c2d%arr, Wijkd%arr, Udkba%f, c2d%popr, c2d%pops, c2d%popp, c2d%popq, &
                                Wijkd%pops, Wijkd%popr, c2d%ioppq, Wijkd%sym, Udkba%sym, &
                                fact, iCore, maxCore)
    end if

end subroutine c2d_WijkdUdkba

! write at 2025.7.4, lihsun
! calculate Cijab -= fact*Wijkd*Udkba
! Cijab and Udkba sort as (ab,ij) and (ba,dk)
subroutine p_c2d_WijkdUdkba_h(c2d, Wijkd, file_Udkba, pop_i, pop_j, pop_a, pop_b, &
                              pop_d, pop_k, iop, symW, symU, fact, iCore, maxCore)
    use mod_symm
    use mod_orbit
    implicit none

    real*8, dimension(*), intent(inout) :: c2d
    real*8, dimension(*), intent(in) :: Wijkd
    character(len=8), intent(in) :: file_Udkba
    integer*8, dimension(8), intent(in) :: pop_i, pop_j, pop_a, pop_b, pop_d, pop_k
    integer*8, intent(in) :: iop
    integer*8, intent(in) :: symW, symU
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore
    integer*8, intent(in) :: maxCore

    integer*8 :: symC, sym1, sym2
    integer*8 :: irrepk, numk, k
    integer*8 :: i0, i1, i2, i3, nSize, iEnd
    integer*8, external :: irpdso, idsymsoc

    i0 = 1
    iEnd = 0
    symC = dirprd(symW, symU)
    ! W(ij,kd)->W(ij,dk), i0
    nSize = idsymsoc(symW, pop_i, pop_j, iop, pop_k, pop_d, 1)
    i1 = i0 + nSize
    call isymtrso(Wijkd, iCore(i0), pop_i, pop_j, iop, 0, pop_k, pop_d, 1, 1, symW)
    ! cycle for k
    do irrepk = 1, nirrep
        numk = pop_k(irrepk)
        sym1 = dirprd(irrepk, symW) ! symmetry of W(ij,d*)
        sym2 = dirprd(irrepk, symU) ! symmerty of U(ab,d*)
        do k = 1, numk
            ! get W(ij,d*), i1
            nSize = idsymsoc(sym1, pop_i, pop_j, iop, pop_d, zorb, 1)
            i2 = i1 + nSize
            call getpqrso_incore(iCore(i1), k, irrepk, &
                                 pop_i, pop_j, iop, pop_d, pop_k, 1, iCore(i0), symW)
            ! get U(ab,d*), i2
            nSize = idsymsoc(sym2, pop_a, pop_b, iop, pop_d, zorb, 1)
            i3 = i2 + nSize
            call getpqrso(iCore(i2), k, irrepk, 1, file_Udkba)
            ! Cabij += fact*U(ab,d*)*W(ij,d*)
            call VpqmnVrsmn_to_Vpqrs(iCore(i2), iCore(i1), sym2, sym1, pop_d, zorb, 1, &
                                     pop_a, pop_b, iop, pop_i, pop_j, iop, c2d, symC, fact, 1.0d0)
        end do
        iEnd = max(i3, iEnd)
    end do

    if (iEnd > maxCore) then
        write (6, *) 'ERROR: memory not enough for p_c2d_WijkdUdkba_h', iEnd, maxCore
        stop
    end if
end subroutine p_c2d_WijkdUdkba_h

! 2026.7.13
!> c(ab,ij) -= fact*W(ij,kd)*U(ba,dk)
subroutine c2d_WijkdUdkba_R(c2dABAB, WijkdABAB, UdkbaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dABAB
    type(opt), intent(in) :: WijkdABAB
    type(opt), intent(in) :: UdkbaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    i0 = 1
    ! c(aB,iJ) = -W(iJ,kD)*U(Ba,Dk) -> -W(iJ,kD)*U(bA,dK)
    !            -W(iJ,Kd)*U(Ba,dK) = -W(Ji,Kd)*U(aB,dK) -> -W(jI,kD)*U(aB,dK)
    call tmp%init(iCore, c2dABAB%popp, c2dABAB%popq, c2dABAB%popr, c2dABAB%pops, &
                  1, 1, c2dABAB%sym, T, i0)
    i1 = i0 + tmp%len
    call dscal(tmp%len, 0.0d0, tmp%arr, 1)
    ! tmp(aB,iJ) = -W(iJ,kD)*U(bA,dK)
    call c2d_WijkdUdkba(tmp, WijkdABAB, UdkbaABAB, fact, iCore(i1), maxCore - tmp%len)
    ! tmp2(aB,iJ) = tmp(aB,iJ)+tmp(bA,jI), i1
    call iassymso(tmp%arr, iCore(i1), tmp%sym, tmp%popp, tmp%popq, tmp%popr, tmp%pops, 3)
    call daxpy(tmp%len, 1.0d0, iCore(i1), 1, c2dABAB%arr, 1)

end subroutine c2d_WijkdUdkba_R

! 2026.7.13
!> c(ab,ij) -= fact*W(ij,kd)*U(ba,dk)
subroutine c2d_WijkdUdkba_U(c2dAAAA, c2dABAB, c2dBBBB, WijkdAAAA, WijkdABAB, WijkdBBBB, WijkdBABA, &
                            UdkbaAAAA, UdkbaABAB, UdkbaBBBB, UdkbaBABA, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dAAAA, c2dABAB, c2dBBBB
    type(opt), intent(in) :: WijkdAAAA, WijkdABAB, WijkdBBBB, WijkdBABA
    type(opt), intent(in) :: UdkbaAAAA, UdkbaABAB, UdkbaBBBB, UdkbaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    i0 = 1
    ! c(ab,ij) = -W(ij,kd)*U(ba,dk)
    call c2d_WijkdUdkba(c2dAAAA, WijkdAAAA, UdkbaAAAA, fact, iCore, maxCore)

    ! c(aB,iJ) = -W(iJ,kD)*U(Ba,Dk)
    call c2d_WijkdUdkba(c2dABAB, WijkdABAB, UdkbaBABA, fact, iCore, maxCore)
    ! c(aB,iJ) = -W(iJ,Kd)*U(Ba,dK) = -W(Ji,Kd)*U(aB,dK)
    call tmp%init(iCore, c2dABAB%popq, c2dABAB%popp, c2dABAB%pops, c2dABAB%popr, &
                  1, 1, c2dABAB%sym, T, i0)
    i1 = i0 + tmp%len
    call dscal(tmp%len, 0.0d0, tmp%arr, 1)
    ! tmp(Ba,Ji) = -W(Ji,Kd)*U(aB,dK)
    call c2d_WijkdUdkba(tmp, WijkdBABA, UdkbaABAB, fact, iCore(i1), maxCore - tmp%len)
    call isymtrso(tmp%arr, iCore(i1), tmp%popp, tmp%popq, 1, 1, &
                  tmp%popr, tmp%pops, 1, 1, tmp%sym)
    call daxpy(tmp%len, 1.0d0, iCore(i1), 1, c2dABAB%arr, 1)

    ! c(AB,IJ) = -W(IJ,KD)*U(BA,DK)
    call c2d_WijkdUdkba(c2dBBBB, WijkdBBBB, UdkbaBBBB, fact, iCore, maxCore)

end subroutine c2d_WijkdUdkba_U
