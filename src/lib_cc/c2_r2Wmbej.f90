! 2026.7.6
!> c(ab,ij) += fact*P(ab|ij){r(ae,im)*Wmbej(em,bj)}
subroutine c2_r2Wmbej(c2, r2, Wmbej, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: Wmbej
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_444(r2, Wmbej, c2, fact, 'aeim', 'embj', 'abij', iCore, maxCore, 3)

end subroutine c2_r2Wmbej

! 2026.7.6
!> c(ab,ij) += fact*P(ab|ij){r(ae,im)*Wmbej(em,bj)}
subroutine c2_r2Wmbej_R(c2ABAB, r2AAAA, r2ABAB, WmbejAAAA, WmbejABAB, WmbejBAAB, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2ABAB
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: WmbejAAAA, WmbejABAB, WmbejBAAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1
    type(opt) :: tmp

    i0 = 1

    ! c(aB,iJ) = r(ae,im)*WmBeJ(em,BJ)
    !            r(aE,iM)*WMBEJ(EM,BJ) -> r(aE,iM)*Wmbej(em,bj)
    !           -r(aE,Jm)*WmBEi(Em,Bi) = r(aE,mJ)*WmBEi(Em,Bi) -> r(aE,mJ)*WMbeI(eM,bI)
    !           -r(Be,iM)*WMaeJ(eM,aJ) = r(eB,iM)*WMaeJ(eM,aJ)
    !            r(Be,Jm)*Wmaei(em,ai) = r(eB,mJ)*Wmaei(em,ai)
    !            r(BE,JM)*WMaEi(EM,ai) -> r(BE,JM)*WmAeI(em,AI)
    !--- c(aB,iJ) = r(ae,im)*WmBeJ(em,BJ)+r(BE,JM)*WmAeI(em,AI)
    !--- c(aB,iJ) = r(aE,iM)*Wmbej(em,bj)+r(eB,mJ)*Wmaei(em,ai)
    call copy_opt(c2ABAB, tmp, T, iCore, i0)
    i1 = i0 + tmp%len
    call dscal(tmp%len, 0.0d0, tmp%arr, 1)
    call vv2v_444(r2AAAA, WmbejABAB, tmp, fact, 'aeim', 'embj', 'abij', &
                  iCore(i1), maxCore - tmp%len, 0)
    call vv2v_444(r2ABAB, WmbejAAAA, tmp, fact, 'aeim', 'embj', 'abij', &
                  iCore(i1), maxCore - tmp%len, 0)
    call iassymso(iCore(i0), iCore(i1), tmp%sym, tmp%popp, tmp%popq, tmp%popr, tmp%pops, 3)
    call daxpy(tmp%len, 1.0d0, iCore(i1), 1, c2ABAB%arr, 1)
    !--- c(aB,iJ) = r(aE,mJ)*WMbeI(eM,bI)
    call vv2v_444(r2ABAB, WmbejBAAB, c2ABAB, fact, 'aemj', 'embi', 'abij', &
                  iCore, maxCore, 0)
    !--- c(aB,iJ) = r(eB,iM)*WMaeJ(eM,aJ)
    call vv2v_444(r2ABAB, WmbejBAAB, c2ABAB, fact, 'ebim', 'emaj', 'abij', &
                  iCore, maxCore, 0)

end subroutine c2_r2Wmbej_R

! 2026.7.6
!> c(ab,ij) += fact*P(ab|ij){r(ae,im)*Wmbej(em,bj)}
subroutine c2_r2Wmbej_U(c2AAAA, c2ABAB, c2BBBB, r2AAAA, r2ABAB, r2BBBB, &
                        WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA, &
                        fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c2AAAA, c2ABAB, c2BBBB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: WmbejAAAA, WmbejABAB, WmbejBAAB, WmbejBBBB, WmbejBABA, WmbejABBA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(ab,ij) = P(ab|ij){r(ae,im)*Wmbej(em,bj)}
    !           +P(ab|ij){r(aE,iM)*WMbEj(EM,bj)}
    call c2_r2Wmbej(c2AAAA, r2AAAA, WmbejAAAA, fact, iCore, maxCore)
    call c2_r2Wmbej(c2AAAA, r2ABAB, WmbejBABA, fact, iCore, maxCore)

    ! c(aB,iJ) = r(ae,im)*WmBeJ(em,BJ)
    call vv2v_444(r2AAAA, WmbejABAB, c2ABAB, fact, 'aeim', 'emBJ', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = r(aE,iM)*WMBEJ(EM,BJ)
    call vv2v_444(r2ABAB, WmbejBBBB, c2ABAB, fact, 'aEiM', 'EMBJ', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = -r(aE,Jm)*WmBEi(Em,Bi) = r(aE,mJ)*WmBEi(Em,Bi)
    call vv2v_444(r2ABAB, WmbejABBA, c2ABAB, fact, 'aEmJ', 'EmBi', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = -r(Be,iM)*WMaeJ(eM,aJ) = r(eB,iM)*WMaeJ(eM,aJ)
    call vv2v_444(r2ABAB, WmbejBAAB, c2ABAB, fact, 'eBiM', 'eMaJ', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = r(Be,Jm)*Wmaei(em,ai) = r(eB,mJ)*Wmaei(em,ai)
    call vv2v_444(r2ABAB, WmbejAAAA, c2ABAB, fact, 'eBmJ', 'emai', 'aBiJ', iCore, maxCore, 0)
    ! c(aB,iJ) = r(BE,JM)*WMaEi(EM,ai)
    call vv2v_444(r2BBBB, WmbejBABA, c2ABAB, fact, 'BEJM', 'EMai', 'aBiJ', iCore, maxCore, 0)

    ! c(AB,IJ) = P(AB|IJ){r(AE,IM)*WMBEJ(EM,BJ)}
    !           +P(AB|IJ){r(Ae,Im)*WmBeJ(em,BJ)} = P(AB|IJ){r(eA,mI)*WmBeJ(em,BJ)}
    call c2_r2Wmbej(c2BBBB, r2BBBB, WmbejBBBB, fact, iCore, maxCore)
    call vv2v_444(r2ABAB, WmbejABAB, c2BBBB, fact, 'eAmI', 'emBJ', 'ABIJ', iCore, maxCore, 3)

end subroutine c2_r2Wmbej_U
