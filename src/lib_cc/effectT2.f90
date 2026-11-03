! 2026.4.15
! calculate antisymmetrization of t(a,i)*f(b,j)
! type=1: t(ab,ij) += fact*Pab{t(a,i)*t(b,j)} = fact*Pij{t(a,i)*t(b,j)}, ab,ij is antisym
!                  += fact*[t(a,i)*t(b,j)-t(a,j)*t(b,i)]
! type=2: t(ab,ij) += fact*PabPij{t(a,i)*f(b,j)}, ab,ij is antisym
!                  += fact*[t(a,i)*f(b,j)-t(a,j)*f(b,i)-t(b,i)*f(a,j)+t(b,j)*f(a,i)]
! type=3: t(ab,ij) += fact*Pab{t(a,i)*f(b,j)}, only ab is antisym
!                  += fact*[t(a,i)*f(b,j)-t(b,i)*f(a,j)]
! type=4: t(ab,ij) += fact*Pij{t(a,i)*f(b,j)}, only ij is antisym
!                  += fact*[t(a,i)*f(b,j)-t(a,j)*f(b,i)]
subroutine effectT2(T2, T1_1, T1_2, fact, iCore, maxCore, type)
    use mod_operator
    implicit none

    type(opt), intent(inout) :: T2
    type(opt), intent(in) :: T1_1
    type(opt), intent(in) :: T1_2
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: type

    if (T2%incore) then
        call p_effectT2(T2%arr, T1_1%arr, T1_2%arr, iCore, fact, maxCore, &
                        T2%popp, T2%popq, T2%popr, T2%pops, T2%ioppq, T2%ioprs, &
                        T1_1%sym, T1_2%sym, T2%sym, type)
    else
        call p_effectT2_h(T2%f, T1_1%arr, T1_2%arr, iCore, fact, maxCore, &
                          T2%popp, T2%popq, T2%popr, T2%pops, T2%ioppq, T2%ioprs, &
                          T1_1%sym, T1_2%sym, T2%sym, type)
    end if

end subroutine effectT2

! write by LiShun, 2024.5.14
! calculate antisymmetrization of t(a,i)*f(b,j)
! type=1: t(ab,ij) += fact*Pab{t(a,i)*t(b,j)} = fact*Pij{t(a,i)*t(b,j)}, ab,ij is antisym
!                  += fact*[t(a,i)*t(b,j)-t(a,j)*t(b,i)]
! type=2: t(ab,ij) += fact*PabPij{t(a,i)*f(b,j)}, ab,ij is antisym
!                  += fact*[t(a,i)*f(b,j)-t(a,j)*f(b,i)-t(b,i)*f(a,j)+t(b,j)*f(a,i)]
! type=3: t(ab,ij) += fact*Pab{t(a,i)*f(b,j)}, only ab is antisym
!                  += fact*[t(a,i)*f(b,j)-t(b,i)*f(a,j)]
! type=4: t(ab,ij) += fact*Pij{t(a,i)*f(b,j)}, only ij is antisym
!                  += fact*[t(a,i)*f(b,j)-t(a,j)*f(b,i)]
!----------------------------
! modified at 2024.12.4
subroutine p_effectT2(t2, t1_1, t1_2, iCore, fact, maxCore, &
                      pop_a, pop_b, pop_i, pop_j, iopab, iopij, &
                      symT1_1, symT1_2, symT2, type)
    use mod_orbit
    use mod_symm
    implicit none

    real*8, dimension(*), intent(inout) :: t2
    real*8, dimension(*), intent(in) :: t1_1, t1_2
    real*8, dimension(maxCore), intent(inout) :: iCore
    real*8, intent(in) :: fact
    integer*8, intent(in) :: maxCore
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: iopab, iopij
    integer*8, intent(in) :: symT1_1, symT1_2, symT2
    integer*8, intent(in) :: type

    integer*8 :: symm
    integer*8 :: i0, i1, i2, iEnd, nSize, itmp
    integer*8, external :: idsymsoc

    symm = dirprd(symT1_1, symT1_2)
    if (symm /= symT2) then
        write (6, *) 'ERROR: symmerty of effectT error', symT1_1, symT1_2, symT2
        stop
    end if
    i0 = 1
    if (type == 1) then ! iop can be 0 or 1
        if (iopab /= iopij) then
            write (6, *) 'ERROR: p_effectT2 type=1 only support iopab=iopij'
            stop
        end if
        !---------- t1*t1
        ! tmp(ai,bj) = t1(ai,**)*t1(**,bj), i0
        nSize = idsymsoc(symT2, pop_a, pop_i, 1, pop_b, pop_j, 1)
        i1 = i0 + nSize
        call VpqmnVmnrs_to_Vpqrs(t1_1, t1_2, symT1_1, symT1_2, zorb, zorb, 1, &
                                 pop_a, pop_i, 1, pop_b, pop_j, 1, iCore(i0), symT2, 1.0D0, 0.0D0)
        ! tmp(ai,bj) -> tmp(ab,ij), i1
        i2 = i1 + nSize
        call sstgenso(iCore(i0), iCore(i1), nSize, pop_a, pop_i, pop_b, pop_j, iCore(i2), symT2, '1324')
        iEnd = i2 + nSize
        ! tmp(ab,ij) -> tmp(ab,ji), i0
        call isymtrso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_i, pop_j, 1, 1, symT2)
        ! tmp2(ab,ij) = tmp(ab,ij)-tmp(ab,ji), is antisymm, i1
        call daxpy(nSize, -1.0d0, iCore(i0), 1, iCore(i1), 1)
        if (iopab == 0) then
            ! compress (ab,ij)
            call icompso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_i, pop_j, 1, 1, symT2)
            call icompso(iCore(i0), iCore(i1), pop_a, pop_b, 0, 1, pop_i, pop_j, 1, 0, symT2)
        end if
        ! t(ab,ij) += fact*tmp2(ab,ij)
        nSize = idsymsoc(symT2, pop_a, pop_b, iopab, pop_i, pop_j, iopij)
        call daxpy(nSize, fact, iCore(i1), 1, t2, 1)
    else if (type == 2) then ! iop can be 0 or 1
        !---------- t1*t1
        ! tmp(ai,bj) = t1_1(ai,**)*t1_2(**,bj), i0
        nSize = idsymsoc(symT2, pop_a, pop_i, 1, pop_b, pop_j, 1)
        i1 = i0 + nSize
        call VpqmnVmnrs_to_Vpqrs(t1_1, t1_2, symT1_1, symT1_2, zorb, zorb, 1, &
                                 pop_a, pop_i, 1, pop_b, pop_j, 1, iCore(i0), symT2, 1.0D0, 0.0D0)
        ! tmp(ai,bj) -> tmp(ab,ij), i1
        i2 = i1 + nSize
        call sstgenso(iCore(i0), iCore(i1), nSize, pop_a, pop_i, pop_b, pop_j, iCore(i2), symT2, '1324')
        ! tmp(ab,ij) -> tmp(ab,ji), i0
        call isymtrso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_i, pop_j, 1, 1, symT2)
        ! tmp2(ab,ij) = tmp(ab,ij)-tmp(ab,ji), ij is antisymm, i1
        call daxpy(nSize, -1.0d0, iCore(i0), 1, iCore(i1), 1)
        itmp = i1
        if (iopij == 0) then
            ! compress tmp2(ab,ij), i0
            call icompso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 0, symT2)
            itmp = i0
        end if
        ! tmp2(ab,ij)->tmp2(ba,ij), i2
        nSize = idsymsoc(symT2, pop_a, pop_b, 1, pop_i, pop_j, iopij)
        iEnd = i2 + nSize
        call isymtrso(iCore(itmp), iCore(i2), pop_a, pop_b, 1, 1, pop_i, pop_j, iopij, 0, symT2)
        ! tmp3(ab,ij)=tmp2(ab,ij)-tmp2(ba,ij), itmp
        call daxpy(nSize, -1.0d0, iCore(i2), 1, iCore(itmp), 1)
        if (iopab == 0) then
            ! compress tmp3(ab,ij), i2
            call icompso(iCore(itmp), iCore(i2), pop_a, pop_b, 1, 0, pop_i, pop_j, 0, 1, symT2)
            itmp = i2
        end if
        ! t(ab,ij) += fact*tmp3(ab,ij)
        nSize = idsymsoc(symT2, pop_a, pop_b, iopab, pop_i, pop_j, iopij)
        call daxpy(nSize, fact, iCore(itmp), 1, t2, 1)
    else if (type == 3 .or. type == 4) then ! iop of one side can be 0
        ! tmp1(ai,bj)=t(**,ai)*r(**,bj), i0
        nSize = idsymsoc(symT2, pop_a, pop_i, 1, pop_b, pop_j, 1)
        i1 = i0 + nSize
        call VmnpqVmnrs_to_Vpqrs(t1_1, t1_2, symT1_1, symT1_2, zorb, zorb, 1, &
                                 pop_a, pop_i, 1, pop_b, pop_j, 1, iCore(i0), symT2, 1.0d0, 0.0d0)
        ! tmp1(ai,bj)->tmp1(ab,ij), i1
        i2 = i1 + nSize
        call sstgenso(iCore(i0), iCore(i1), nSize, pop_a, pop_i, pop_b, pop_j, &
                      iCore(i2), symT2, '1324')
        if (type == 3) then
            ! tmp1(ab,ij)->tmp1(ba,ij), i0
            call isymtrso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 0, symT2)
        else if (type == 4) then
            ! tmp1(ab,ij)->tmp1(ab,ji), i0
            call isymtrso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_i, pop_j, 1, 1, symT2)
        end if
        ! tmp2(ab,ij)=tmp1(ab,ij)-tmp1(ba,ij) or tmp1(ab,ji), i1
        call daxpy(nSize, -1.0d0, iCore(i0), 1, iCore(i1), 1)
        itmp = i1
        if ((iopab == 0) .or. (iopij == 0)) then
            ! compress ab or ij
            if ((type == 3) .and. (iopab == 0)) then
                call icompso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 0, pop_i, pop_j, 1, 1, symT2)
            else if ((type == 4) .and. (iopij == 0)) then
                call icompso(iCore(i1), iCore(i0), pop_a, pop_b, 1, 1, pop_i, pop_j, 1, 0, symT2)
            end if
            itmp = i0
        end if
        ! t2(ab,ij) += fact*tmp2(ab,ij)
        nSize = idsymsoc(symT2, pop_a, pop_b, iopab, pop_i, pop_j, iopij)
        call daxpy(nSize, fact, iCore(itmp), 1, t2, 1)
        iEnd = i2
    end if

    if (iEnd >= maxCore) then
        write (6, *) 'iCore not enough for p_effectT2'
        stop
    end if

end subroutine p_effectT2

! 2026.4.23
subroutine p_effectT2_h(f_t2, t1_1, t1_2, iCore, fact, maxCore, &
                        pop_a, pop_b, pop_i, pop_j, iopab, iopij, &
                        symT1_1, symT1_2, symT2, type)
    use mod_orbit
    use mod_symm
    implicit none

    character(len=8), intent(in) :: f_t2
    real*8, dimension(*), intent(in) :: t1_1, t1_2
    real*8, dimension(maxCore), intent(inout) :: iCore
    real*8, intent(in) :: fact
    integer*8, intent(in) :: maxCore
    integer*8, dimension(8), intent(in) :: pop_a, pop_b, pop_i, pop_j ! the occupation number for each irrep
    integer*8, intent(in) :: iopab, iopij
    integer*8, intent(in) :: symT1_1, symT1_2, symT2
    integer*8, intent(in) :: type

    integer*8 :: irrep, numj, j
    integer*8 :: sym1, sym2
    integer*8 :: i0, i1, i2, i3, i4, iT, nSize
    integer*8, external :: idsymsoc

    i0 = 1
    if (type == 3) then
        ! cycle for symmetry of j
        do irrep = 1, nirrep
            numj = pop_j(irrep)
            sym1 = dirprd(symT2, irrep) ! symmetry of T2(ab,i*)
            sym2 = dirprd(symT1_2, irrep) ! symmetry of T1_2(b,*)

            nSize = idsymsoc(sym1, pop_a, pop_b, iopab, pop_i, zorb, 1)
            i1 = i0 + nSize
            i2 = i1 + pop_b(sym2)
            ! cycle for j
            do j = 1, numj
                ! get T2(ab,i*), i0
                call getpqrso(iCore(i0), j, irrep, 1, f_t2)
                ! get T1_2(b,*), i1
                call getpqrso_incore(iCore(i1), j, irrep, zorb, zorb, 1, &
                                     pop_b, pop_j, 1, t1_2, symT1_2)
                ! tmp(ai,b*) = T1_1(a,i)*T1_2(b,*), i2
                nSize = idsymsoc(sym1, pop_a, pop_i, 1, pop_b, zorb, 1)
                i3 = i2 + nSize
                call VmnpqVmnrs_to_Vpqrs(t1_1, iCore(i1), symT1_1, sym2, zorb, zorb, 1, &
                                         pop_a, pop_i, 1, pop_b, zorb, 1, iCore(i2), sym1, &
                                         1.0d0, 0.0d0)
                ! tmp(ai,b*)->tmp(ab,i*), i3
                i4 = i3 + nSize
                call sstgenso(iCore(i2), iCore(i3), nSize, pop_a, pop_i, &
                              pop_b, zorb, iCore(i4), sym1, '1324')
                ! tmp(ab,i*)->tmp(ba,i*), i2
                call isymtrso(iCore(i3), iCore(i2), pop_a, pop_b, 1, 1, &
                              pop_i, zorb, 1, 0, sym1)
                ! tmp2(ab,i*) = tmp(ab,i*)-tmp(ba,i*), i3
                call daxpy(nSize, -1.0d0, iCore(i2), 1, iCore(i3), 1)
                iT = i3
                if (iopab == 0) then
                    ! compress ab, i2
                    call icompso(iCore(i3), iCore(i2), pop_a, pop_b, 1, 0, &
                                 pop_i, zorb, 1, 1, sym1)
                    iT = i2
                end if
                ! T2(ab,i*) += fact*tmp(ab,i*)
                nSize = idsymsoc(sym1, pop_a, pop_b, iopab, pop_i, zorb, 1)
                call daxpy(nSize, fact, iCore(iT), 1, iCore(i0), 1)
                ! put T2(ab,i*)
                call putpqrso(iCore(i0), j, irrep, 1, f_t2)
            end do
        end do
    else
        write (6, *) 'ERROR: p_effectT2_h only support type=3'
        stop
    end if

end subroutine p_effectT2_h

! 2026.7.6
! type=1: t(ab,ij) += fact*P(ab){t(a,i)*t(b,j)}
! type=2: t(ab,ij) += fact*P(ab|ij){t(a,i)*f(b,j)}
! type=3: t(ab,ij) += fact*P(ab){t(a,i)*f(b,j)}
! type=4: t(ab,ij) += fact*P(ij){t(a,i)*f(b,j)}
subroutine effectT2_R(T2AAAA, T2ABAB, T1AA, R1AA, fact, iCore, maxCore, type)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: T2AAAA, T2ABAB
    type(opt), intent(in) :: T1AA, R1AA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: type

    if (type == 1) then
        ! t(ab,ij) += P(ab){t(a,i)*t(b,j)}
        call vv2v_224(T1AA, T1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 1)
        ! t(aB,iJ) += t(a,i)*t(B,J)
        call vv2v_224(T1AA, T1AA, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
    else if (type == 2) then
        ! t(ab,ij) += P(ab|ij){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 3)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        !            +t(B,J)*r(a,i)
        call vv2v_224(T1AA, R1AA, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        call vv2v_224(R1AA, T1AA, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
    else if (type == 3) then
        ! t(ab,ij) += P(ab){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 1)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        call vv2v_224(T1AA, R1AA, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
    else if (type == 4) then
        ! t(ab,ij) += P(ij){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 2)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        call vv2v_224(T1AA, R1AA, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
    end if

end subroutine effectT2_R

! 2026.7.14
! type=1: t(ab,ij) += fact*P(ab){t(a,i)*t(b,j)}
! type=2: t(ab,ij) += fact*P(ab|ij){t(a,i)*f(b,j)}
! type=3: t(ab,ij) += fact*P(ab){t(a,i)*f(b,j)}
! type=4: t(ab,ij) += fact*P(ij){t(a,i)*f(b,j)}
subroutine effectT2_U(T2AAAA, T2ABAB, T2BBBB, T2BABA, T1AA, T1BB, R1AA, R1BB, &
                      fact, iCore, maxCore, type)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: T2AAAA, T2ABAB, T2BBBB, T2BABA
    type(opt), intent(in) :: T1AA, T1BB, R1AA, R1BB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore
    integer*8, intent(in) :: type

    if (type == 1) then
        ! t(ab,ij) += P(ab){t(a,i)*t(b,j)}
        call vv2v_224(T1AA, T1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 1)
        ! t(aB,iJ) += t(a,i)*t(B,J)
        call vv2v_224(T1AA, T1BB, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        ! t(AB,IJ) += P(AB){t(A,I)*t(B,J)}
        call vv2v_224(T1BB, T1BB, T2BBBB, fact, 'AI', 'BJ', 'ABIJ', iCore, maxCore, 1)
    else if (type == 2) then
        ! t(ab,ij) += P(ab|ij){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 3)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        !            +t(B,J)*r(a,i)
        call vv2v_224(T1AA, R1BB, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        call vv2v_224(R1AA, T1BB, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        ! t(AB,IJ) += P(AB|IJ){t(A,I)*r(B,J)}
        call vv2v_224(T1BB, R1BB, T2BBBB, fact, 'AI', 'BJ', 'ABIJ', iCore, maxCore, 3)
    else if (type == 3) then
        ! t(ab,ij) += P(ab){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 1)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        call vv2v_224(T1AA, R1BB, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        ! t(Ab,Ij) += t(A,I)*r(b,j)
        call vv2v_224(T1BB, R1AA, T2BABA, fact, 'AI', 'bj', 'AbIj', iCore, maxCore, 0)
        ! t(AB,IJ) += P(AB){t(A,I)*r(B,J)}
        call vv2v_224(T1BB, R1BB, T2BBBB, fact, 'AI', 'BJ', 'ABIJ', iCore, maxCore, 1)
    else if (type == 4) then
        ! t(ab,ij) += P(ij){t(a,i)*r(b,j)}
        call vv2v_224(T1AA, R1AA, T2AAAA, fact, 'ai', 'bj', 'abij', iCore, maxCore, 2)
        ! t(aB,iJ) += t(a,i)*r(B,J)
        call vv2v_224(T1AA, R1BB, T2ABAB, fact, 'ai', 'BJ', 'aBiJ', iCore, maxCore, 0)
        ! t(Ab,Ij) += t(A,I)*r(b,j)
        call vv2v_224(T1BB, R1AA, T2BABA, fact, 'AI', 'bj', 'AbIj', iCore, maxCore, 0)
        ! t(AB,IJ) += P(IJ){t(A,I)*r(B,J)}
        call vv2v_224(T1BB, R1BB, T2BBBB, fact, 'AI', 'BJ', 'ABIJ', iCore, maxCore, 2)
    end if

end subroutine effectT2_U
