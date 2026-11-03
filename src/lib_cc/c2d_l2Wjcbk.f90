! 2026.7.10
!> c(ab,ij) += fact*P(ab|ij){Wjcbk(bj,ck)*l(ac,ik)}
subroutine c2d_l2Wjcbk(c2d, l2, Wjcbk, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2d
    type(opt), intent(in) :: l2
    type(opt), intent(in) :: Wjcbk
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_444(l2, Wjcbk, c2d, fact, 'acik', 'bjck', 'abij', iCore, maxCore, 3)

end subroutine c2d_l2Wjcbk

! 2026.7.10
!> c(ab,ij) += fact*P(ab|ij){Wjcbk(bj,ck)*l(ac,ik)}
subroutine c2d_l2Wjcbk_R(c2dABAB, l2AAAA, l2ABAB, WjcbkAAAA, WjcbkABAB, WjcbkBAAB, &
                         fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dABAB
    type(opt), intent(in) :: l2AAAA, l2ABAB
    type(opt), intent(in) :: WjcbkAAAA, WjcbkABAB, WjcbkBAAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    i0 = 1

    ! c(aB,iJ) = WJcBk(BJ,ck)*l(ac,ik) -> WjCbK(bj,CK)*l(ac,ik)
    !           +WJCBK(BJ,CK)*l(aC,iK) -> Wjcbk(bj,ck)*l(aC,iK)
    !           +Wicak(ai,ck)*l(Bc,Jk) -> Wicak(ai,ck)*l(bC,jK)
    !           +WiCaK(ai,CK)*l(BC,JK) -> WiCaK(ai,CK)*l(bc,jk)
    !           -WJcaK(aJ,cK)*l(Bc,iK) = WJcaK(aJ,cK)*l(cB,iK)
    !           -WiCBk(Bi,Ck)*l(aC,Jk) = WiCBk(Bi,Ck)*l(aC,kJ) -> WIcbK(bI,cK)*l(aC,kJ)
    call copy_opt(c2dABAB, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call dscal(tmp%len, 0.0d0, tmp%arr, 1)
    ! c(aB,iJ) = WjCbK(bj,CK)*l(ac,ik)+WiCaK(ai,CK)*l(bc,jk)
    call vv2v_444(l2AAAA, WjcbkABAB, tmp, fact, 'acik', 'bjck', 'abij', &
                  iCore(i1), maxCore - tmp%len, 0)
    ! c(aB,iJ) = Wjcbk(bj,ck)*l(aC,iK)+Wicak(ai,ck)*l(bC,jK)
    call vv2v_444(l2ABAB, WjcbkAAAA, tmp, fact, 'acik', 'bjck', 'abij', &
                  iCore(i1), maxCore - tmp%len, 0)
    call iassymso(iCore(i0), iCore(i1), tmp%sym, tmp%popp, tmp%popq, tmp%popr, tmp%pops, 3)
    call daxpy(tmp%len, 1.0d0, iCore(i1), 1, c2dABAB%arr, 1)
    ! c(aB,iJ) = WJcaK(aJ,cK)*l(cB,iK)
    call vv2v_444(l2ABAB, WjcbkBAAB, c2dABAB, fact, 'cbik', 'ajck', 'abij', &
                  iCore, maxCore, 0)
    ! c(aB,iJ) = WIcbK(bI,cK)*l(aC,kJ)
    call vv2v_444(l2ABAB, WjcbkBAAB, c2dABAB, fact, 'ackj', 'bick', 'abij', &
                  iCore, maxCore, 0)

end subroutine c2d_l2Wjcbk_R

! 2026.7.10
!> c(ab,ij) += fact*P(ab|ij){Wjcbk(bj,ck)*l(ac,ik)}
subroutine c2d_l2Wjcbk_U(c2dAAAA, c2dABAB, c2dBBBB, l2AAAA, l2ABAB, l2BBBB, &
                         WjcbkAAAA, WjcbkABAB, WjcbkBAAB, WjcbkBBBB, WjcbkBABA, WjcbkABBA, &
                         fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2dAAAA, c2dABAB, c2dBBBB
    type(opt), intent(in) :: l2AAAA, l2ABAB, l2BBBB
    type(opt), intent(in) :: WjcbkAAAA, WjcbkABAB, WjcbkBAAB, WjcbkBBBB, WjcbkBABA, WjcbkABBA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    type(opt) :: tmp
    integer*8 :: i0, i1

    i0 = 1

    ! c(ab,ij) = P(ab|ij){Wjcbk(bj,ck)*l(ac,ik)}
    !           +P(ab|ij){WjCbK(bj,CK)*l(aC,iK)}
    call c2d_l2Wjcbk(c2dAAAA, l2AAAA, WjcbkAAAA, fact, iCore, maxCore)
    call c2d_l2Wjcbk(c2dAAAA, l2ABAB, WjcbkABAB, fact, iCore, maxCore)

    ! c(aB,iJ) = WJcBk(BJ,ck)*l(ac,ik)
    call vv2v_444(l2AAAA, WjcbkBABA, c2dABAB, fact, 'acik', 'BJck', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = WJCBK(BJ,CK)*l(aC,iK)
    call vv2v_444(l2ABAB, WjcbkBBBB, c2dABAB, fact, 'aCiK', 'BJCK', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = Wicak(ai,ck)*l(Bc,Jk) = Wicak(ai,ck)*l(cB,kJ)
    call vv2v_444(l2ABAB, WjcbkAAAA, c2dABAB, fact, 'cBkJ', 'aick', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = WiCaK(ai,CK)*l(BC,JK)
    call vv2v_444(l2BBBB, WjcbkABAB, c2dABAB, fact, 'BCJK', 'aiCK', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) =-WJcaK(aJ,cK)*l(Bc,iK) = WJcaK(aJ,cK)*l(cB,iK)
    call vv2v_444(l2ABAB, WjcbkBAAB, c2dABAB, fact, 'cBiK', 'aJcK', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) =-WiCBk(Bi,Ck)*l(aC,Jk) = WiCBk(Bi,Ck)*l(aC,kJ)
    call vv2v_444(l2ABAB, WjcbkABBA, c2dABAB, fact, 'aCkJ', 'BiCk', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) = P(AB|IJ){WJCBK(BJ,CK)*l(AC,IK)}
    !           +P(AB|IJ){WJcBk(BJ,ck)*l(Ac,Ik)} = P(AB|IJ){WJcBk(BJ,ck)*l(cA,kI)}
    call c2d_l2Wjcbk(c2dBBBB, l2BBBB, WjcbkBBBB, fact, iCore, maxCore)
    call vv2v_444(l2ABAB, WjcbkBABA, c2dBBBB, fact, 'cAkI', 'BJck', 'ABIJ', iCore, maxCore, 3)

end subroutine c2d_l2Wjcbk_U
