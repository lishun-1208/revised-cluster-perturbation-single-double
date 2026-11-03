! 2026.7.6
!> c(a,i) += fact*r(f,n)*Vnafi(fn,ai)
subroutine c1_r1iabj(c1, r1, iabj, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1
    type(opt), intent(in) :: r1
    type(opt), intent(in) :: iabj
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    call vv2v_422(iabj, r1, c1, fact, 'fnai', 'fn', 'ai', iCore, maxCore)

end subroutine c1_r1iabj

! 2026.7.6
!> c(a,i) += fact*r(f,n)*Vnafi(fn,ai), need V(ia,bj)=V(jb,ai)
subroutine c1_r1iabj_R(c1AA, r1AA, iabjAAAA, iabjABAB, fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA
    type(opt), intent(in) :: r1AA
    type(opt), intent(in) :: iabjAAAA, iabjABAB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = r(f,n)*Vnafi(fn,ai)
    call c1_r1iabj(c1AA, r1AA, iabjAAAA, fact, iCore, maxCore)
    ! c(a,i) = r(F,N)*VNaFi(FN,ai)-> r(F,N)*ViFaN(ai,FN)
    call vv2v_422(iabjABAB, r1AA, c1AA, fact, 'aiFN', 'FN', 'ai', iCore, maxCore)

end subroutine c1_r1iabj_R

! 2026.7.6
!> c(a,i) += fact*r(f,n)*Vnafi(fn,ai), need V(ia,bj)=V(jb,ai)
subroutine c1_r1iabj_U(c1AA, c1BB, r1AA, r1BB, iabjAAAA, iabjABAB, iabjBBBB, &
                       fact, iCore, maxCore)
    use mod_operator
    use mod_contract
    implicit none

    type(opt), intent(inout) :: c1AA, c1BB
    type(opt), intent(in) :: r1AA, r1BB
    type(opt), intent(in) :: iabjAAAA, iabjABAB, iabjBBBB
    real*8, intent(in) :: fact
    real*8, dimension(maxCore), intent(inout) :: iCore !> free work space
    integer*8, intent(in) :: maxCore

    ! c(a,i) = r(f,n)*Vnafi(fn,ai)
    call c1_r1iabj(c1AA, r1AA, iabjAAAA, fact, iCore, maxCore)
    ! c(a,i) = r(F,N)*VNaFi(FN,ai) -> r(F,N)*ViFaN(ai,FN)
    call vv2v_422(iabjABAB, r1BB, c1AA, fact, 'aiFN', 'FN', 'ai', iCore, maxCore)

    ! c(A,I) = r(f,n)*VnAfI(fn,AI)
    call c1_r1iabj(c1BB, r1AA, iabjABAB, fact, iCore, maxCore)
    ! c(A,I) = r(F,N)*VNAFI(FN,AI)
    call c1_r1iabj(c1BB, r1BB, iabjBBBB, fact, iCore, maxCore)

end subroutine c1_r1iabj_U
