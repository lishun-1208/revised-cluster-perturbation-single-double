!=======================================================================
!  subroutine: vv2v_dpd_424_check_s
!
!  功能说明:
!    针对 4指标 * 2指标 -> 4指标 (收缩 1 个指标) 的统一基准测试套件。
!    执行顺序: CPU-SP
!    调用参考: VpqrmVms_to_Vpqrs / VpqrmVsm_to_Vpqrs
!=======================================================================
subroutine vv2v_dpd_424_check_s() 
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
       idsymsoc(symvpqrm, popp, popq, 1, popr, popm, 1), & ! pqrm / all full
       idsymsoc(symvpqrm, popp, popq, 0, popr, popm, 1), & ! abim, a<b
       idsymsoc(symvpqrm, popp, popq, 2, popr, popm, 1), & ! abim, a<=b
       idsymsoc(symvpqrm, popp, popq, 0, popm, popr, 1), & ! abmi, a<b
       idsymsoc(symvpqrm, popp, popq, 2, popm, popr, 1), & ! abmi, a<=b
       idsymsoc(symvpqrm, popr, popm, 1, popp, popq, 0), & ! imab, a<b in pair2
       idsymsoc(symvpqrm, popr, popm, 1, popp, popq, 2), & ! imab, a<=b in pair2
       idsymsoc(symvpqrm, popm, popr, 1, popp, popq, 0), & ! miab, a<b in pair2
       idsymsoc(symvpqrm, popm, popr, 1, popp, popq, 2), & ! miab, a<=b in pair2
       idsymsoc(symvpqrm, popp, popq, 0, popm, popr, 1), & ! ijmb legacy
       idsymsoc(symvpqrm, popp, popq, 1, popr, popm, 0), & ! aeij legacy
       idsymsoc(symvpqrm, pops, popq, 1, popr, popm, 0), & ! beij legacy
       idsymsoc(symvpqrm, popq, pops, 1, popr, popm, 0), & ! ebij legacy
       idsymsoc(symvpqrm, pops, popq, 1, popm, popr, 1), & ! beji legacy
       idsymsoc(symvpqrm, pops, popq, 1, popr, popm, 1), & ! beij legacy full
       idsymsoc(symvpqrm, popq, pops, 1, popr, popm, 1)  & ! ebij legacy full
       )

  nms = max( &
       irpdso(symvms, popm, pops, 1), & ! mj/ms
       irpdso(symvms, pops, popm, 1), & ! jm/sm/am
       irpdso(symvms, pops, popq, 1), & ! be legacy
       irpdso(symvms, popp, popq, 1), & ! ae legacy
       irpdso(symvms, popq, popp, 1), & ! ea legacy
       irpdso(symvms, popq, pops, 1)  & ! eb legacy
       )

  npqrs = max( &
       idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1), & ! pqrs / all full
       idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1), & ! abij, a<b
       idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1), & ! abij, a<=b
       idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1), & ! abji, a<b
       idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1), & ! abji, a<=b
       idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0), & ! ijab, a<b in pair2
       idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2), & ! ijab, a<=b in pair2
       idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0), & ! jiab, a<b in pair2
       idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2), & ! jiab, a<=b in pair2
       idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1), & ! ijba legacy
       idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 0), & ! abij legacy, ij compressed
       idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 1), & ! abij legacy full
       idsymsoc(symvpqrs, popp, pops, 1, popm, popr, 1), & ! abji legacy full
       idsymsoc(symvpqrs, pops, popp, 1, popr, popm, 0), & ! baij legacy compressed
       idsymsoc(symvpqrs, pops, popp, 1, popr, popm, 1)  & ! baij legacy full
       )

  nwork = max(npqrm, npqrs, nms)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_424_check_s] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  ! --- 打印矩阵空间大小信息 ---
  write(6,*) " [ Memory Allocation Info (Base SP / 4 Bytes) ] "
  write(6,"(A15, I15, A, F10.2, A)") "    V1(pqrm) : ", npqrm, " elements | ", real(npqrm, 8)*4.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2(ms)   : ", nms,   " elements | ", real(nms, 8)*4.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V3(pqrs) : ", npqrs, " elements | ", real(npqrs, 8)*4.0d0/1048576.0d0, " MB"
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
  ! [ 阶段 1：CPU 单精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  ! --- CPU SP: Standard 424 (pqrm * ms -> pqrs) ---
  write(6,*) '>>> TEST 1: Standard Mode (pqrm * ms -> pqrs, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_s(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_s, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Swapped 424 (pqrm * sm -> pqrs) ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqrm * sm -> pqrs, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVsm_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_s(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_s, pops, popm, 1, symvms, 'sm', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Extreme Scramble (rmpq * ms -> prsq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * ms -> prsq, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, pops, popq, 1) 
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_core_s, symvpqrs, fact, fact2)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call sstgenso1(vpqrs_core_s, vpqrs_ref_s, popp, popq, popr, pops, symvpqrs, '1342') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqrm_s, icore_s, popp, popq, popr, popm, symvpqrm, '3412')

  vpqrs_core_s(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(icore_s, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       vms_s, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_s, popp, popr, 1, pops, popq, 1, symvpqrs, 'prsq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_s(1:npqrs_iter) = 888.0d0
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 1, popr, pops, vpqrs_ref_s, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_s(1:npqrs_iter) = 888.0d0  
     call vv2v_dpd_424(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_s, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * ms -> pqrs, IOP_pq=0) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 0, popr, pops, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_s(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_s, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * ms -> pqrs, IOP_pq=2) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VpqrmVms_to_Vpqrs(vpqrm_s, vms_s, symvpqrm, symvms, popm, popp, popq, 2, popr, pops, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_s(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                       vms_s, popm, pops, 1, symvms, 'ms', &
                       vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  !======================================================================
  ! [ 重构重点：CPU SP compressed virtual-pair 424 tests ]
  !
  ! 目标:
  !   1. 优先覆盖 a<b,im * mj -> a<b,ij
  !   2. 同时覆盖 a<=b,im * mj -> a<=b,ij
  !   3. 覆盖 im/mi, mj/jm, ij/ji 的任意重排
  !   4. 覆盖 a<b / a<=b pair 在 L/W 的左右 pair 之间整体移动
  !
  ! 合法性:
  !   * a,b 使用同一虚轨道空间；iopab=0/2 时必须 popp=popq。
  !   * i,m,j 使用占据轨道空间；本组测试中 im/mj/ij 均保持 full iop=1。
  !   * 不测试非法的 a,i 或 a,m 压缩 pair。
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) ' [ CPU SP PRIORITY 424: a<b/a<=b virtual-pair tests ]'
  write(6,*) '----------------------------------------------------------'

  write(6,*) ''
  write(6,*) '--- P1: a<b/a<=b, im/mi * mj/jm -> a<b/a<=b, ij/ji ---'

  ! --- P1-01: abim * mj -> abij, iopab=0 ---
  write(6,*) '>>> P1-01: abim * mj -> abij, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-02: abim * mj -> abji, iopab=0 ---
  write(6,*) '>>> P1-02: abim * mj -> abji, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-03: abim * jm -> abij, iopab=0 ---
  write(6,*) '>>> P1-03: abim * jm -> abij, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-04: abim * jm -> abji, iopab=0 ---
  write(6,*) '>>> P1-04: abim * jm -> abji, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-05: abmi * mj -> abij, iopab=0 ---
  write(6,*) '>>> P1-05: abmi * mj -> abij, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 0, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-06: abmi * mj -> abji, iopab=0 ---
  write(6,*) '>>> P1-06: abmi * mj -> abji, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 0, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-07: abmi * jm -> abij, iopab=0 ---
  write(6,*) '>>> P1-07: abmi * jm -> abij, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 0, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-08: abmi * jm -> abji, iopab=0 ---
  write(6,*) '>>> P1-08: abmi * jm -> abji, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 0, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-09: abim * mj -> abij, iopab=2 ---
  write(6,*) '>>> P1-09: abim * mj -> abij, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-10: abim * mj -> abji, iopab=2 ---
  write(6,*) '>>> P1-10: abim * mj -> abji, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-11: abim * jm -> abij, iopab=2 ---
  write(6,*) '>>> P1-11: abim * jm -> abij, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-12: abim * jm -> abji, iopab=2 ---
  write(6,*) '>>> P1-12: abim * jm -> abji, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-13: abmi * mj -> abij, iopab=2 ---
  write(6,*) '>>> P1-13: abmi * mj -> abij, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 2, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-14: abmi * mj -> abji, iopab=2 ---
  write(6,*) '>>> P1-14: abmi * mj -> abji, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 2, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-15: abmi * jm -> abij, iopab=2 ---
  write(6,*) '>>> P1-15: abmi * jm -> abij, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 2, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P1-16: abmi * jm -> abji, iopab=2 ---
  write(6,*) '>>> P1-16: abmi * jm -> abji, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 2, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  write(6,*) ''
  write(6,*) '--- P2: im/mi,a<b/a<=b * mj/jm -> ij/ji,a<b/a<=b ---'

  ! --- P2-01: imab * mj -> ijab, iopab=0 ---
  write(6,*) '>>> P2-01: imab * mj -> ijab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 0, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-02: imab * mj -> jiab, iopab=0 ---
  write(6,*) '>>> P2-02: imab * mj -> jiab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 0, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 0, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-03: imab * jm -> ijab, iopab=0 ---
  write(6,*) '>>> P2-03: imab * jm -> ijab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'imab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 0, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-04: imab * jm -> jiab, iopab=0 ---
  write(6,*) '>>> P2-04: imab * jm -> jiab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 0, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'imab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 0, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-05: miab * mj -> ijab, iopab=0 ---
  write(6,*) '>>> P2-05: miab * mj -> ijab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 0, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 0, symvpqrm, 'miab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 0, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-06: miab * mj -> jiab, iopab=0 ---
  write(6,*) '>>> P2-06: miab * mj -> jiab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 0, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 0, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 0, symvpqrm, 'miab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 0, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-07: miab * jm -> ijab, iopab=0 ---
  write(6,*) '>>> P2-07: miab * jm -> ijab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 0, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 0, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 0, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-08: miab * jm -> jiab, iopab=0 ---
  write(6,*) '>>> P2-08: miab * jm -> jiab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 0, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 0, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 0, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 0, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-09: imab * mj -> ijab, iopab=2 ---
  write(6,*) '>>> P2-09: imab * mj -> ijab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 2, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 2, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-10: imab * mj -> jiab, iopab=2 ---
  write(6,*) '>>> P2-10: imab * mj -> jiab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 2, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 2, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 2, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-11: imab * jm -> ijab, iopab=2 ---
  write(6,*) '>>> P2-11: imab * jm -> ijab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 2, symvpqrm, 'imab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 2, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-12: imab * jm -> jiab, iopab=2 ---
  write(6,*) '>>> P2-12: imab * jm -> jiab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 2, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 2, symvpqrm, 'imab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 2, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-13: miab * mj -> ijab, iopab=2 ---
  write(6,*) '>>> P2-13: miab * mj -> ijab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 2, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 2, symvpqrm, 'miab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 2, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-14: miab * mj -> jiab, iopab=2 ---
  write(6,*) '>>> P2-14: miab * mj -> jiab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 2, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 2, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 2, symvpqrm, 'miab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 2, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-15: miab * jm -> ijab, iopab=2 ---
  write(6,*) '>>> P2-15: miab * jm -> ijab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 2, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 2, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 2, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P2-16: miab * jm -> jiab, iopab=2 ---
  write(6,*) '>>> P2-16: miab * jm -> jiab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 2, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 2, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 2, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 2, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  write(6,*) ''
  write(6,*) '--- P3: cross-output movement of intact a<b/a<=b pair ---'

  ! --- P3-01: abim * mj -> ijab, iopab=0 ---
  write(6,*) '>>> P3-01: abim * mj -> ijab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 0, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-02: abim * jm -> jiab, iopab=0 ---
  write(6,*) '>>> P3-02: abim * jm -> jiab, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 0, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 0, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 0, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-03: imab * mj -> abij, iopab=0 ---
  write(6,*) '>>> P3-03: imab * mj -> abij, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-04: miab * jm -> abji, iopab=0 ---
  write(6,*) '>>> P3-04: miab * jm -> abji, iopab=0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 0, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 0, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 0, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-05: abim * mj -> ijab, iopab=2 ---
  write(6,*) '>>> P3-05: abim * mj -> ijab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popr, pops, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call itranspso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popr, pops, 1, popp, popq, 2, symvpqrs, 'ijab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-06: abim * jm -> jiab, iopab=2 ---
  write(6,*) '>>> P3-06: abim * jm -> jiab, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, pops, popr, 1, popp, popq, 2)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call zeroso(icore_s, nwork)
  call itranspso(v3_tmp_s, icore_s, popp, popq, 2, popr, pops, 1, symvpqrs)
  call isymtrso(icore_s, vpqrs_ref_s, popr, pops, 1, 1, popp, popq, 2, 0, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, pops, popr, 1, popp, popq, 2, symvpqrs, 'jiab', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-07: imab * mj -> abij, iopab=2 ---
  write(6,*) '>>> P3-07: imab * mj -> abij, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call itranspso(vpqrm_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  vpqrs_ref_s(1:npqrs_iter) = v3_tmp_s(1:npqrs_iter)
  call vv2v_dpd_424(vpqrm_s, popr, popm, 1, popp, popq, 2, symvpqrm, 'imab', &
                    vms_s, popm, pops, 1,                 symvms,   'mj',   &
                    vpqrs_core_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P3-08: miab * jm -> abji, iopab=2 ---
  write(6,*) '>>> P3-08: miab * jm -> abji, iopab=2 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(icore_s, nwork)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, icore_s, popm, popr, 1, 1, popp, popq, 2, 0, symvpqrm)
  call itranspso(icore_s, v1_tr_s, popr, popm, 1, popp, popq, 2, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popm, popr, 1, popp, popq, 2, symvpqrm, 'miab', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  write(6,*) ''
  write(6,*) '--- P4: sign / accumulation checks on priority a<b path ---'

  ! --- P4-01: abmi * jm -> abji, iopab=0, fact=-1.0 ---
  write(6,*) '>>> P4-01: abmi * jm -> abji, iopab=0, fact=-1.0 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, pops, popr, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 0, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  ! IMPORTANT: the reference route must use the same fact/fact2 as the optimized route.
  ! Otherwise a negative-sign test reports equal norms but relative error = 2.0.
  call vv2v_dpd_424(v1_tr_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'abij', -1.0d0, 0.0d0)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 0, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s, pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'abji', -1.0d0, 0.0d0)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- P4-02: accumulation on transformed priority case ---
  write(6,*) '>>> P4-02: abmi * jm -> abji, iopab=2, fact=1.0, fact2=0.5 [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, pops, popr, 1)
  vpqrs_ref_s(1:npqrs_iter)  = 777.0d0
  vpqrs_core_s(1:npqrs_iter) = 777.0d0
  call zeroso(v1_tr_s, npqrm)
  call isymtrso(vpqrm_s, v1_tr_s, popp, popq, 2, 0, popm, popr, 1, 1, symvpqrm)
  call zeroso(vms_tr_s, nms)
  call tr2dpd_s(vms_s, vms_tr_s, pops, popm, symvms)
  call zeroso(v3_tmp_s, npqrs)
  v3_tmp_s(1:idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 1)) = 777.0d0
  call vv2v_dpd_424(v1_tr_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'abim', &
                    vms_tr_s, popm, pops, 1,                 symvms,   'mj',   &
                    v3_tmp_s, popp, popq, 2, popr, pops, 1, symvpqrs, 'abij', 1.0d0, 0.5d0)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, popq, 2, 0, popr, pops, 1, 1, symvpqrs)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 2, popm, popr, 1, symvpqrm, 'abmi', &
                    vms_s,   pops, popm, 1,                 symvms,   'jm',   &
                    vpqrs_core_s, popp, popq, 2, pops, popr, 1, symvpqrs, 'abji', 1.0d0, 0.5d0)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))


  !======================================================================
  ! [ 阶段 1B：CPU SP 额外排列/压缩 pair 回归测试 ]
  !
  ! 注意：
  !   * 如果某一侧 pair 的 iop < 1（例如 ij 为 i<j 压缩 pair），
  !     不能用 sstgenso/sstgenso1 把它当作普通四维指标重排。
  !     只允许用 isymtrso 在保持压缩 pair 完整的前提下交换另一侧 pair。
  !   * 下面的 TEST 7/8/14 专门验证 ij 压缩 pair 下的 ijmb 路径。
  !   * TEST 9/10/11 专门验证输出右 pair 为 i<j 压缩的 aeij/beij/ebij 路径。
  !   * TEST 12/13 验证 all-iop=1 时右 pair 反序 beji/ea 的路径。
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1B: CPU SP SPECIAL 424 PERMUTATION TESTS ]'
  write(6,*) '----------------------------------------------------------'

  ! --- CPU SP SPECIAL 7: ijmb * am -> ijba, compare with ijbm * am -> ijba ---
  write(6,*) '>>> TEST 7: ijmb * am -> ijba, ij compressed, compare with ijbm route [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)

  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'ijbm', &
                    vms_s,   pops, popm, 1,                 symvms,   'am',   &
                    vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', fact, fact2)

  call zeroso(icore_s, nwork)
  call isymtrso(vpqrm_s, icore_s, popp, popq, 0, 0, popr, popm, 1, 1, symvpqrm)
  call vv2v_dpd_424(icore_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'ijmb', &
                    vms_s,   pops, popm, 1,                 symvms,   'am',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 8: ijmb * ma -> ijba, swapped 2-index input ---
  write(6,*) '>>> TEST 8: ijmb * ma -> ijba, ij compressed, swapped T1 pair [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)

  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'ijbm', &
                    vms_s,   popm, pops, 1,                 symvms,   'ma',   &
                    vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', fact, fact2)

  call zeroso(icore_s, nwork)
  call isymtrso(vpqrm_s, icore_s, popp, popq, 0, 0, popr, popm, 1, 1, symvpqrm)
  call vv2v_dpd_424(icore_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'ijmb', &
                    vms_s,   popm, pops, 1,                 symvms,   'ma',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))
  call vv2v_dpd_424(icore_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'ijmb', &
                    vms_s,   popm, pops, 1,                 symvms,   'ma',   &
                    vpqrs_core_s, popp, popq, 0, pops, popr, 1, symvpqrs, 'ijab', fact, fact2)
  call isymtrso(vpqrs_core_s, icore_s, popp, popq, 0, 0, pops, popr, 1, 1, symvpqrs)
  call diffv(vpqrs_ref_s, icore_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 9: aeij * be -> abij, ij compressed ---
  ! This is the canonical path that should work even when the output right pair ij is compressed.
  write(6,*) '>>> TEST 9: aeij * be -> abij, ij compressed smoke/regression [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)

  call vv2v_dpd_424(vpqrm_s, popp, popq, 1, popr, popm, 0, symvpqrm, 'aeij', &
                    vms_s,   pops, popq, 1,                 symvms,   'be',   &
                    vpqrs_ref_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)
  call vv2v_dpd_424(vpqrm_s, popp, popq, 1, popr, popm, 0, symvpqrm, 'aeij', &
                    vms_s,   pops, popq, 1,                 symvms,   'be',   &
                    vpqrs_core_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)
  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 10: beij * ae -> abij, ij compressed ---
  ! Reference uses the mathematically equivalent path:
  !   beij * ae -> abij
  !   ebij * ea -> abij
  ! where ebij is obtained by isymtrso on the left pair, and ea is a 2-index transpose of ae.
  write(6,*) '>>> TEST 10: beij * ae -> abij, ij compressed, compare with ebij * ea [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call zeroso(vms_tr_s, nms)

  call vv2v_dpd_424(vpqrm_s, pops, popq, 1, popr, popm, 0, symvpqrm, 'beij', &
                    vms_s,   popp, popq, 1,                 symvms,   'ae',   &
                    vpqrs_core_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)

  call isymtrso(vpqrm_s, v1_tr_s, pops, popq, 1, 1, popr, popm, 0, 0, symvpqrm)
  call tr2dpd_s(vms_s, vms_tr_s, popp, popq, symvms)
  call vv2v_dpd_424(v1_tr_s, popq, pops, 1, popr, popm, 0, symvpqrm, 'ebij', &
                    vms_tr_s, popq, popp, 1,                symvms,   'ea',   &
                    vpqrs_ref_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)

  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 11: ebij * ae -> abij, ij compressed ---
  ! Reference uses beij * ea after transposing the first pair of V1 and the 2-index T.
  write(6,*) '>>> TEST 11: ebij * ae -> abij, ij compressed, compare with beij * ea [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 0)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call zeroso(vms_tr_s, nms)

  call vv2v_dpd_424(vpqrm_s, popq, pops, 1, popr, popm, 0, symvpqrm, 'ebij', &
                    vms_s,   popp, popq, 1,                 symvms,   'ae',   &
                    vpqrs_core_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)

  call isymtrso(vpqrm_s, v1_tr_s, popq, pops, 1, 1, popr, popm, 0, 0, symvpqrm)
  call tr2dpd_s(vms_s, vms_tr_s, popp, popq, symvms)
  call vv2v_dpd_424(v1_tr_s, pops, popq, 1, popr, popm, 0, symvpqrm, 'beij', &
                    vms_tr_s, popq, popp, 1,                symvms,   'ea',   &
                    vpqrs_ref_s, popp, pops, 1, popr, popm, 0, symvpqrs, 'abij', fact, fact2)

  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 12: beji * ae -> abij, all iop = 1 ---
  ! Reference:
  !   beji --isymtrso(right pair)--> beij
  !   beij * ae -> abji
  !   abji --isymtrso(right pair)--> abij
  write(6,*) '>>> TEST 12: beji * ae -> abij, all iop=1, compare with beij route [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call zeroso(v3_tmp_s, npqrs)

  call vv2v_dpd_424(vpqrm_s, pops, popq, 1, popm, popr, 1, symvpqrm, 'beji', &
                    vms_s,   popp, popq, 1,                 symvms,   'ae',   &
                    vpqrs_core_s, popp, pops, 1, popr, popm, 1, symvpqrs, 'abij', fact, fact2)

  call isymtrso(vpqrm_s, v1_tr_s, pops, popq, 1, 0, popm, popr, 1, 1, symvpqrm)
  call vv2v_dpd_424(v1_tr_s, pops, popq, 1, popr, popm, 1, symvpqrm, 'beij', &
                    vms_s,   popp, popq, 1,                 symvms,   'ae',   &
                    v3_tmp_s, popp, pops, 1, popm, popr, 1, symvpqrs, 'abji', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, pops, 1, 0, popm, popr, 1, 1, symvpqrs)

  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 13: beji * ea -> abij, all iop = 1, swapped 2-index ---
  write(6,*) '>>> TEST 13: beji * ea -> abij, all iop=1, swapped 2-index [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, pops, 1, popr, popm, 1)
  call zeroso(vpqrs_ref_s,  npqrs_iter)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call zeroso(v1_tr_s, npqrm)
  call zeroso(v3_tmp_s, npqrs)

  call vv2v_dpd_424(vpqrm_s, pops, popq, 1, popm, popr, 1, symvpqrm, 'beji', &
                    vms_s,   popq, popp, 1,                 symvms,   'ea',   &
                    vpqrs_core_s, popp, pops, 1, popr, popm, 1, symvpqrs, 'abij', fact, fact2)

  call isymtrso(vpqrm_s, v1_tr_s, pops, popq, 1, 0, popm, popr, 1, 1, symvpqrm)
  call vv2v_dpd_424(v1_tr_s, pops, popq, 1, popr, popm, 1, symvpqrm, 'beij', &
                    vms_s,   popq, popp, 1,                 symvms,   'ea',   &
                    v3_tmp_s, popp, pops, 1, popm, popr, 1, symvpqrs, 'abji', fact, fact2)
  call isymtrso(v3_tmp_s, vpqrs_ref_s, popp, pops, 1, 0, popm, popr, 1, 1, symvpqrs)

  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP SPECIAL 14: fact2 accumulation on problematic ijmb path ---
  write(6,*) '>>> TEST 14: ijmb * am -> ijba with fact2=0.5 accumulation [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  vpqrs_ref_s(1:npqrs_iter)  = 777.0d0
  vpqrs_core_s(1:npqrs_iter) = 777.0d0

  call vv2v_dpd_424(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'ijbm', &
                    vms_s,   pops, popm, 1,                 symvms,   'am',   &
                    vpqrs_ref_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', 1.0d0, 0.5d0)

  call zeroso(icore_s, nwork)
  call isymtrso(vpqrm_s, icore_s, popp, popq, 0, 0, popr, popm, 1, 1, symvpqrm)
  call vv2v_dpd_424(icore_s, popp, popq, 0, popm, popr, 1, symvpqrm, 'ijmb', &
                    vms_s,   pops, popm, 1,                 symvms,   'am',   &
                    vpqrs_core_s, popp, popq, 0, popr, pops, 1, symvpqrs, 'ijba', 1.0d0, 0.5d0)

  write(6,"(A25, E15.8)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1)
  write(6,"(A25, E15.8)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1)
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))


end subroutine vv2v_dpd_424_check_s

!=======================================================================
! 2-index DPD transpose helper for CPU double check only:
!   vin(p,q) -> vout(q,p), iop = 1
! This is used only to build reference tensors in vv2v_dpd_424_check.
!=======================================================================
subroutine tr2dpd_s(vin, vout, popp, popq, sym)
  use mod_iop
  use mod_symm
  implicit none

  integer, intent(in) :: popp(8), popq(8), sym
  real*4, intent(in)  :: vin(*)
  real*4, intent(out) :: vout(*)

  integer :: irrepq, irrepp, nump, numq
  integer :: p, q
  integer :: vstarto, vstartn
  integer :: base_in, base_out
  integer, external :: isymoffso

  do irrepq = 1, nirrep
     irrepp = dirprd(irrepq, sym)
     numq = popq(irrepq)
     nump = popp(irrepp)
     if (nump == 0 .or. numq == 0) cycle

     vstarto = isymoffso(irrepq, sym, popp, popq, 1)
     vstartn = isymoffso(irrepp, sym, popq, popp, 1)

     do q = 1, numq
        base_in  = vstarto - 1 + (q - 1) * nump
        base_out = vstartn - 1 + q
        do p = 1, nump
           vout(base_out + (p - 1) * numq) = vin(base_in + p)
        end do
     end do
  end do
end subroutine tr2dpd_s


