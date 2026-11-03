!=======================================================================
!  subroutine: vv2v_dpd_424_check_sg
!
!  功能说明:
!    针对 4指标 * 2指标 -> 4指标 (收缩 1 个指标) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    调用参考: VpqrmVms_to_Vpqrs / VpqrmVsm_to_Vpqrs
!=======================================================================
subroutine vv2v_dpd_424_check_sg() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npqrm, nms, npqrs
  integer :: npqrs_iter, nwork  
  integer :: symvpqrm, symvms, symvpqrs
  integer :: popm(8), popp(8), popq(8), popr(8), pops(8), popdum(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_ioprm, def_ioprs

  ! --- 外部函数 ---
  integer, external :: idsymsoc, irpdso
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpqrm(:), vms(:), vpqrs_core(:), vpqrs_ref(:)
  real*8, allocatable :: icore(:)
  real*8, allocatable :: v1_tr(:), vms_tr(:), v3_tmp(:)
  real*4, allocatable :: vpqrm_s(:), vms_s(:), vpqrs_core_s(:), vpqrs_ref_s(:)
  real*4, allocatable :: icore_s(:)
  real*4, allocatable :: v1_tr_s(:), vms_tr_s(:), v3_tmp_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpqrm_g(:), vms_g(:), vpqrs_core_g(:), vpqrs_ref_g(:)
  real*8, allocatable, device :: icore_g(:)
  real*4, allocatable, device :: vpqrm_sg(:), vms_sg(:), vpqrs_core_sg(:), vpqrs_ref_sg(:)
  real*4, allocatable, device :: icore_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpqrm, symvms, perm_in)
  
  call init_parameters(npqrm, nms, npqrs, &
                       symvpqrm, symvms, symvpqrs, &
                       popm, popdum, popp, popq, popr, pops, &
                       def_ioprm, def_ioppq, def_ioprs, &
                       nirrep_in, nvrt_in, npop_in) 

  ! =====================================================================
  ! 严谨的内存管理：利用 DPD 物理规则精确计算各个张量的大小
  ! =====================================================================
  ! Base + special-test maximum sizes.
  ! The special 424 tests intentionally use different index spaces for
  ! a/b/e/i/j/m.  Therefore the test arrays must be sized by the maximum
  ! of every layout used below, not only by the default pqrm/ms/pqrs layout.
  npqrm = max( &
       idsymsoc(symvpqrm, popp, popq, 1, popr, popm, 1), & ! pqrm / aeij full
       idsymsoc(symvpqrm, popp, popq, 0, popr, popm, 1), & ! ijbm
       idsymsoc(symvpqrm, popp, popq, 0, popm, popr, 1), & ! ijmb
       idsymsoc(symvpqrm, popp, popq, 1, popr, popm, 0), & ! aeij, ij compressed
       idsymsoc(symvpqrm, pops, popq, 1, popr, popm, 0), & ! beij, ij compressed
       idsymsoc(symvpqrm, popq, pops, 1, popr, popm, 0), & ! ebij, ij compressed
       idsymsoc(symvpqrm, pops, popq, 1, popm, popr, 1), & ! beji full
       idsymsoc(symvpqrm, pops, popq, 1, popr, popm, 1), & ! beij full
       idsymsoc(symvpqrm, popq, pops, 1, popr, popm, 1)  & ! ebij full
       )

  nms = max( &
       irpdso(symvms, popm, pops, 1), & ! ms
       irpdso(symvms, pops, popm, 1), & ! sm / am in ijmb tests
       irpdso(symvms, popm, pops, 1), & ! ma
       irpdso(symvms, pops, popq, 1), & ! be
       irpdso(symvms, popp, popq, 1), & ! ae
       irpdso(symvms, popq, popp, 1), & ! ea
       irpdso(symvms, popq, pops, 1)  & ! eb
       )

  npqrs = max( &
       idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1), & ! pqrs
       idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1), & ! ijba
       idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 0), & ! abij, ij compressed
       idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 1), & ! abij full
       idsymsoc(symvpqrs, popp, pops, 1, popm, popr, 1), & ! abji full
       idsymsoc(symvpqrs, pops, popp, 1, popr, popm, 0), & ! baij compressed
       idsymsoc(symvpqrs, pops, popp, 1, popr, popm, 1)  & ! baij full
       )

  nwork = max(npqrm, npqrs, nms)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_424_check_sg] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  ! --- 打印矩阵空间大小信息 ---
  write(6,*) " [ Memory Allocation Info (Base DP / 8 Bytes) ] "
  write(6,"(A15, I15, A, F10.2, A)") "    V1(pqrm) : ", npqrm, " elements | ", real(npqrm, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2(ms)   : ", nms,   " elements | ", real(nms, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V3(pqrs) : ", npqrs, " elements | ", real(npqrs, 8)*8.0d0/1048576.0d0, " MB"
  write(6,*) "=========================================================="

  ! --- 精确分配缓冲维度 ---
  allocate(vpqrm(npqrm), vms(nms), vpqrs_core(npqrs), vpqrs_ref(npqrs))
  allocate(icore(nwork))
  allocate(v1_tr(npqrm), vms_tr(nms), v3_tmp(npqrs))
  allocate(vpqrm_s(npqrm), vms_s(nms), vpqrs_core_s(npqrs), vpqrs_ref_s(npqrs))
  allocate(icore_s(nwork))
  allocate(v1_tr_s(npqrm), vms_tr_s(nms), v3_tmp_s(npqrs))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpqrm, vms, npqrm, nms)
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpqrm, vpqrm_s, npqrm)
  call dtos(vms, vms_s, nms)

  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqrm_sg(npqrm), vms_sg(nms), vpqrs_core_sg(npqrs), vpqrs_ref_sg(npqrs))
  allocate(icore_sg(nwork))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpqrm_s, vpqrm_sg, npqrm)
  call myicopy_g(vms_s, vms_sg, nms)

  ! --- GPU SP: Standard 424 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqrm * ms -> pqrs, IOP=1) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_sg, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Swapped 424 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqrm * sm -> pqrs, IOP=1) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVsm_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_sg, pops, popm, 1, symvms, 'sm', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Extreme Scramble (rmpq * ms -> prsq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * ms -> prsq, IOP=1) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, pops, popq, 1) 
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_core_sg, symvpqrs, fact, fact2)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call sstgenso1(vpqrs_core_sg, vpqrs_ref_sg, popp, popq, popr, pops, symvpqrs, '1342') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqrm_sg, icore_sg, popp, popq, popr, popm, symvpqrm, '3412')

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(icore_sg, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       vms_sg, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_sg, popp, popr, 1, pops, popq, 1, symvpqrs, 'prsq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_sg(1:npqrs_iter) = 888.0e0
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_sg, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_sg(1:npqrs_iter) = 888.0e0  
     call vv2v_dpd_424(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_sg, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * ms -> pqrs, IOP_pq=0) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 0, popr, pops, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_sg, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_sg, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_sg, popp, popq, 0, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * ms -> pqrs, IOP_pq=2) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_sg, vms_sg, symvpqrm, symvms, popm, popp, popq, 2, popr, pops, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_sg, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_sg, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_sg, popp, popq, 2, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))


end subroutine vv2v_dpd_424_check_sg

