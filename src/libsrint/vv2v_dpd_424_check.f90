!=======================================================================
!  subroutine: vv2v_dpd_424_check
!
!  功能说明:
!    总入口/调度器。实际测试逻辑拆分在四个独立文件中：
!      vv2v_dpd_424_check_d.f90   : CPU double precision
!      vv2v_dpd_424_check_s.f90   : CPU single precision
!      vv2v_dpd_424_check_dg.f90  : GPU double precision
!      vv2v_dpd_424_check_sg.f90  : GPU single precision
!
!  如果只想验证某一路径，可临时注释掉下面对应 call。
!=======================================================================
subroutine vv2v_dpd_424_check()
  implicit none

  call vv2v_dpd_424_check_d()
  call vv2v_dpd_424_check_s()
  call vv2v_dpd_424_check_dg()
  call vv2v_dpd_424_check_sg()

end subroutine vv2v_dpd_424_check
