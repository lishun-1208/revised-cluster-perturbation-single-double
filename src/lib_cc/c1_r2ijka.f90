! 2026.7.6
!> c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
subroutine c1_r2ijka(c1, r2, ijka, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r2
    type(opt), intent(in) :: ijka
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    if (r2%ioprs == 0) then
        call vv2v_442(r2, ijka, c1, -fact, 'aemn', 'mnie', 'ai', iCore, maxCore)
    else
        call vv2v_442(r2, ijka, c1, -0.5d0*fact, 'aemn', 'mnie', 'ai', iCore, maxCore)
    end if

end subroutine c1_r2ijka

! 2026.7.6
!> c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
subroutine c1_r2ijka_R(c1AA, r2AAAA, r2ABAB, ijkaAAAA, ijkaABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r2AAAA, r2ABAB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
    call c1_r2ijka(c1AA, r2AAAA, ijkaAAAA, fact, iCore, maxCore)
    ! c(a,i) -= fact*r(aE,mN)*V(mN,iE)
    call c1_r2ijka(c1AA, r2ABAB, ijkaABAB, 2.0d0*fact, iCore, maxCore)

end subroutine c1_r2ijka_R

! 2026.7.6
!> c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
subroutine c1_r2ijka_U(c1AA, c1BB, r2AAAA, r2ABAB, r2BBBB, ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA, &
                       fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r2AAAA, r2ABAB, r2BBBB
    type(opt), intent(in) :: ijkaAAAA, ijkaABAB, ijkaBBBB, ijkaBABA
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    integer*8 :: i0, i1
    type(opt) :: tmp

    ! c(a,i) -= fact*0.5*r(ae,mn)*V(mn,ie)
    call c1_r2ijka(c1AA, r2AAAA, ijkaAAAA, fact, iCore, maxCore)
    ! c(a,i) -= fact*r(aE,mN)*V(mN,iE)
    call c1_r2ijka(c1AA, r2ABAB, ijkaABAB, 2.0d0*fact, iCore, maxCore)

    ! c(A,I) -= fact*0.5*r(AE,MN)*V(MN,IE)
    call c1_r2ijka(c1BB, r2BBBB, ijkaBBBB, fact, iCore, maxCore)

    !--- c(A,I) -= fact*r(Ae,Mn)*V(Mn,Ie)
    ! r(eA,nM)->r(Ae,Mn)
    i0 = 1
    call tmp%init(iCore, r2ABAB%popq, r2ABAB%popp, r2ABAB%pops, r2ABAB%popr, &
                  1, 1, r2ABAB%sym, T, i0)
    i1 = i0 + tmp%len
    call isymtrso(r2ABAB%arr, tmp%arr, r2ABAB%popp, r2ABAB%popq, 1, 1, &
                  r2ABAB%popr, r2ABAB%pops, 1, 1, tmp%sym)
    call c1_r2ijka(c1BB, tmp, ijkaBABA, 2.0d0*fact, iCore(i1), maxCore - tmp%len)

end subroutine c1_r2ijka_U
