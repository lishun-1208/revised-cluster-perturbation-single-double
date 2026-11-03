!=======================================================================
!  subroutine: vv2v_dpd_422_check
!
!  功能说明:
!    针对 4指标 * 2指标 -> 2指标 (收缩 2 个指标) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    【核心机制】: 利用 VmnpqVmnrs_to_Vpqrs 配合 pdum 策略作为绝对 Reference。
!=======================================================================
subroutine vv2v_dpd_422_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npqmn, nmn, npq, nqn, npm
  integer :: npq_iter, npm_iter  
  integer :: symvpqmn, symvmn, symvpq
  integer :: popm(8), popn(8), popp(8), popq(8), popdum(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter, i
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_iopmn

  ! --- 外部函数 ---
  integer, external :: idsymsoc, irpdso
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpqmn(:), vmn(:), vpq_core(:), vpq_ref(:)
  real*8, allocatable :: vqn(:), vpm_core(:), vpm_ref(:)
  real*8, allocatable :: icore(:)
  
  real*4, allocatable :: vpqmn_s(:), vmn_s(:), vpq_core_s(:), vpq_ref_s(:)
  real*4, allocatable :: vqn_s(:), vpm_core_s(:), vpm_ref_s(:)
  real*4, allocatable :: icore_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpqmn_g(:), vmn_g(:), vpq_core_g(:), vpq_ref_g(:)
  real*8, allocatable, device :: vqn_g(:), vpm_core_g(:), vpm_ref_g(:)
  real*8, allocatable, device :: icore_g(:)
  
  real*4, allocatable, device :: vpqmn_sg(:), vmn_sg(:), vpq_core_sg(:), vpq_ref_sg(:)
  real*4, allocatable, device :: vqn_sg(:), vpm_core_sg(:), vpm_ref_sg(:)
  real*4, allocatable, device :: icore_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpqmn, symvmn, perm_in)
  
  call init_parameters(npqmn, nmn, npq, &
                       symvpqmn, symvmn, symvpq, &
                       popm, popn, popp, popq, popdum, popdum, &
                       def_iopmn, def_ioppq, def_iopmn, &
                       nirrep_in, nvrt_in, npop_in) 

  ! 【关键】：重置并构建恒为 1 的 Dummy 维度
  popdum = 0
  popdum(1) = 1

  ! =====================================================================
  ! 严谨的内存管理：利用 DPD 物理规则精确计算各个张量的大小
  ! =====================================================================
  npqmn = idsymsoc(symvpqmn, popp, popq, 1, popm, popn, 1)
  nmn   = irpdso(symvmn, popm, popn, 1)
  npq   = irpdso(symvpq, popp, popq, 1)
  
  ! 为 TEST 7 (pqmn * qn = pm) 准备专用的不规则大小张量
  nqn   = irpdso(symvmn, popq, popn, 1)
  npm   = irpdso(symvpq, popp, popm, 1)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_422_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  write(6,*) " [ Memory Allocation Info (Base DP / 8 Bytes) ] "
  write(6,"(A15, I15, A, F10.2, A)") "    V1(pqmn) : ", npqmn, " elements | ", real(npqmn, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2(mn)   : ", nmn,   " elements | ", real(nmn, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V3(pq)   : ", npq,   " elements | ", real(npq, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V_ext(qn): ", nqn,   " elements | ", real(nqn, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V_ext(pm): ", npm,   " elements | ", real(npm, 8)*8.0d0/1048576.0d0, " MB"
  write(6,*) "=========================================================="

  ! --- 精确分配缓冲维度 ---
  allocate(vpqmn(npqmn), vmn(nmn), vpq_core(npq), vpq_ref(npq))
  allocate(vqn(nqn), vpm_core(npm), vpm_ref(npm))
  allocate(icore(npqmn))
  
  allocate(vpqmn_s(npqmn), vmn_s(nmn), vpq_core_s(npq), vpq_ref_s(npq))
  allocate(vqn_s(nqn), vpm_core_s(npm), vpm_ref_s(npm))
  allocate(icore_s(npqmn))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpqmn, vmn, npqmn, nmn)
  
  ! 初始化 vqn 以供 Test 7 独立使用
  do i = 1, nqn
     vqn(i) = sin(real(i, 8)) * 100.0d0
  end do
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpqmn, vpqmn_s, npqmn)
  call dtos(vmn, vmn_s, nmn)
  call dtos(vqn, vqn_s, nqn)

  !======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  ! --- CPU DP: Standard 422 (pqmn * mn -> pq) ---
  write(6,*) '>>> TEST 1: Standard Mode (pqmn * mn -> pq, IOP=1) <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn, icore, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore, vmn, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn, popm, popn, 1, symvmn, 'mn', &
                       vpq_core, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Swapped 422 (pqmn * nm -> pq) ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmn * nm -> pq, IOP=1) <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn, icore, popp, popq, popm, popn, symvpqmn, '4312')
  call VmnpqVmnrs_to_Vpqrs(icore, vmn, symvpqmn, symvmn, popn, popm, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn, popn, popm, 1, symvmn, 'nm', &
                       vpq_core, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Extreme Scramble (mnpq * mn -> pq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (mnpq * mn -> pq, IOP=1) <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1) 
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn, icore, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore, vmn, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(icore, popm, popn, 1, popp, popq, 1, symvpqmn, 'mnpq', &
                       vmn, popm, popn, 1, symvmn, 'mn', &
                       vpq_core, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  vpq_ref(1:npq_iter) = 888.0d0
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn, icore, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore, vmn, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref, symvpq, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpq_core(1:npq_iter) = 888.0d0  
     call vv2v_dpd_422(vpqmn, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn, popm, popn, 1, symvmn, 'mn', &
                       vpq_core, popp, popq, 1, symvpq, 'pq', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (mnpq * mn -> pq, IOP_pq=0) [DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 0)
  call zeroso(vpq_ref, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn, vmn, symvpqmn, symvmn, popm, popn, 1, popp, popq, 0, popdum, popdum, 1, vpq_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn, popm, popn, 1, popp, popq, 0, symvpqmn, 'mnpq', &
                       vmn, popm, popn, 1, symvmn, 'mn', &
                       vpq_core, popp, popq, 0, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (mnpq * mn -> pq, IOP_pq=2) [DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 2)
  call zeroso(vpq_ref, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn, vmn, symvpqmn, symvmn, popm, popn, 1, popp, popq, 2, popdum, popdum, 1, vpq_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npq_iter, vpq_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn, popm, popn, 1, popp, popq, 2, symvpqmn, 'mnpq', &
                       vmn, popm, popn, 1, symvmn, 'mn', &
                       vpq_core, popp, popq, 2, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npq_iter, vpq_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref, vpq_core, int(npq_iter, 8))

  ! --- CPU DP: Non-Adjacent Contraction Test ---
  write(6,*) '>>> TEST 7: Non-Adjacent Contraction (pqmn * qn -> pm, IOP=1) [DP] <<<'
  npm_iter = irpdso(symvpq, popp, popm, 1)
  call zeroso(vpm_ref, npm_iter)

  call walltimeso(t_s)
  ! 将 pqmn 物理转置为 qnpm (即 2, 4, 1, 3 顺序)，以满足参考程序只能前部收缩的死板要求
  call sstgenso1(vpqmn, icore, popp, popq, popm, popn, symvpqmn, '2413')
  call VmnpqVmnrs_to_Vpqrs(icore, vqn, symvpqmn, symvmn, popq, popn, 1, popp, popm, 1, popdum, popdum, 1, vpm_ref, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npm_iter, vpm_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpm_core(1:npm_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     ! 新引擎发威：完全无需内存转置，仅靠命名映射直接算出结果
     call vv2v_dpd_422(vpqmn, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vqn, popq, popn, 1, symvmn, 'qn', &
                       vpm_core, popp, popm, 1, symvpq, 'pm', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npm_iter, vpm_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpm_ref, vpm_core, int(npm_iter, 8))


  !======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  ! --- CPU SP: Standard 422 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqmn * mn -> pq, IOP=1) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_s, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_s, icore_s, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_s, vmn_s, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_s(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_s, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_s, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_s, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Swapped 422 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmn * nm -> pq, IOP=1) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_s, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_s, icore_s, popp, popq, popm, popn, symvpqmn, '4312')
  call VmnpqVmnrs_to_Vpqrs(icore_s, vmn_s, symvpqmn, symvmn, popn, popm, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_s(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_s, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_s, popn, popm, 1, symvmn, 'nm', &
                       vpq_core_s, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Extreme Scramble (mnpq * mn -> pq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (mnpq * mn -> pq, IOP=1) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_s, icore_s, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_s, vmn_s, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_s(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(icore_s, popm, popn, 1, popp, popq, 1, symvpqmn, 'mnpq', &
                       vmn_s, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_s, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  vpq_ref_s(1:npq_iter) = 888.0e0
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_s, icore_s, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_s, vmn_s, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_s, symvpq, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpq_core_s(1:npq_iter) = 888.0e0  
     call vv2v_dpd_422(vpqmn_s, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_s, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_s, popp, popq, 1, symvpq, 'pq', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (mnpq * mn -> pq, IOP_pq=0) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 0)
  call zeroso(vpq_ref_s, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_s, vmn_s, symvpqmn, symvmn, popm, popn, 1, popp, popq, 0, popdum, popdum, 1, vpq_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_s(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_s, popm, popn, 1, popp, popq, 0, symvpqmn, 'mnpq', &
                       vmn_s, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_s, popp, popq, 0, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (mnpq * mn -> pq, IOP_pq=2) [SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 2)
  call zeroso(vpq_ref_s, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_s, vmn_s, symvpqmn, symvmn, popm, popn, 1, popp, popq, 2, popdum, popdum, 1, vpq_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_s(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_s, popm, popn, 1, popp, popq, 2, symvpqmn, 'mnpq', &
                       vmn_s, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_s, popp, popq, 2, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_s, vpq_core_s, int(npq_iter, 8))

  ! --- CPU SP: Non-Adjacent Contraction Test ---
  write(6,*) '>>> TEST 7: Non-Adjacent Contraction (pqmn * qn -> pm, IOP=1) [SP] <<<'
  npm_iter = irpdso(symvpq, popp, popm, 1)
  call zeroso(vpm_ref_s, npm_iter)

  call walltimeso(t_s)
  call sstgenso1(vpqmn_s, icore_s, popp, popq, popm, popn, symvpqmn, '2413')
  call VmnpqVmnrs_to_Vpqrs(icore_s, vqn_s, symvpqmn, symvmn, popq, popn, 1, popp, popm, 1, popdum, popdum, 1, vpm_ref_s, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npm_iter, vpm_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpm_core_s(1:npm_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_s, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vqn_s, popq, popn, 1, symvmn, 'qn', &
                       vpm_core_s, popp, popm, 1, symvpq, 'pm', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npm_iter, vpm_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpm_ref_s, vpm_core_s, int(npm_iter, 8))


  !======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqmn_g(npqmn), vmn_g(nmn), vpq_core_g(npq), vpq_ref_g(npq))
  allocate(vqn_g(nqn), vpm_core_g(npm), vpm_ref_g(npm))
  allocate(icore_g(npqmn))
  
  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vpqmn, vpqmn_g, npqmn)
  call myicopy_g(vmn, vmn_g, nmn)
  call myicopy_g(vqn, vqn_g, nqn)

  ! --- GPU DP: Standard 422 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqmn * mn -> pq, IOP=1) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_g, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_g, icore_g, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_g, vmn_g, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_g, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_g(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_g, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_g, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_g, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Swapped 422 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmn * nm -> pq, IOP=1) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_g, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_g, icore_g, popp, popq, popm, popn, symvpqmn, '4312')
  call VmnpqVmnrs_to_Vpqrs(icore_g, vmn_g, symvpqmn, symvmn, popn, popm, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_g, symvpq, fact, fact2) 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_g(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_g, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_g, popn, popm, 1, symvmn, 'nm', &
                       vpq_core_g, popp, popq, 1, symvpq, 'pq', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Extreme Scramble (mnpq * mn -> pq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (mnpq * mn -> pq, IOP=1) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_g, icore_g, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_g, vmn_g, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_g, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_g(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(icore_g, popm, popn, 1, popp, popq, 1, symvpqmn, 'mnpq', &
                       vmn_g, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_g, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  vpq_ref_g(1:npq_iter) = 888.0d0
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_g, icore_g, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_g, vmn_g, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_g, symvpq, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpq_core_g(1:npq_iter) = 888.0d0  
     call vv2v_dpd_422(vpqmn_g, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_g, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_g, popp, popq, 1, symvpq, 'pq', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (mnpq * mn -> pq, IOP_pq=0) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 0)
  call zeroso(vpq_ref_g, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_g, vmn_g, symvpqmn, symvmn, popm, popn, 1, popp, popq, 0, popdum, popdum, 1, vpq_ref_g, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_g(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_g, popm, popn, 1, popp, popq, 0, symvpqmn, 'mnpq', &
                       vmn_g, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_g, popp, popq, 0, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (mnpq * mn -> pq, IOP_pq=2) [GPU DP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 2)
  call zeroso(vpq_ref_g, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_g, vmn_g, symvpqmn, symvmn, popm, popn, 1, popp, popq, 2, popdum, popdum, 1, vpq_ref_g, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_g(1:npq_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_g, popm, popn, 1, popp, popq, 2, symvpqmn, 'mnpq', &
                       vmn_g, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_g, popp, popq, 2, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_g, vpq_core_g, int(npq_iter, 8))

  ! --- GPU DP: Non-Adjacent Contraction Test ---
  write(6,*) '>>> TEST 7: Non-Adjacent Contraction (pqmn * qn -> pm, IOP=1) [GPU DP] <<<'
  npm_iter = irpdso(symvpq, popp, popm, 1)
  call zeroso(vpm_ref_g, npm_iter)

  call walltimeso(t_s)
  call sstgenso1(vpqmn_g, icore_g, popp, popq, popm, popn, symvpqmn, '2413')
  call VmnpqVmnrs_to_Vpqrs(icore_g, vqn_g, symvpqmn, symvmn, popq, popn, 1, popp, popm, 1, popdum, popdum, 1, vpm_ref_g, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npm_iter, vpm_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpm_core_g(1:npm_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_g, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vqn_g, popq, popn, 1, symvmn, 'qn', &
                       vpm_core_g, popp, popm, 1, symvpq, 'pm', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npm_iter, vpm_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpm_ref_g, vpm_core_g, int(npm_iter, 8))


  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqmn_sg(npqmn), vmn_sg(nmn), vpq_core_sg(npq), vpq_ref_sg(npq))
  allocate(vqn_sg(nqn), vpm_core_sg(npm), vpm_ref_sg(npm))
  allocate(icore_sg(npqmn))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpqmn_s, vpqmn_sg, npqmn)
  call myicopy_g(vmn_s, vmn_sg, nmn)
  call myicopy_g(vqn_s, vqn_sg, nqn)

  ! --- GPU SP: Standard 422 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqmn * mn -> pq, IOP=1) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_sg, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_sg, icore_sg, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_sg, vmn_sg, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_sg(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_sg, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_sg, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_sg, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Swapped 422 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmn * nm -> pq, IOP=1) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  call zeroso(vpq_ref_sg, npq_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_sg, icore_sg, popp, popq, popm, popn, symvpqmn, '4312')
  call VmnpqVmnrs_to_Vpqrs(icore_sg, vmn_sg, symvpqmn, symvmn, popn, popm, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_sg(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_sg, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_sg, popn, popm, 1, symvmn, 'nm', &
                       vpq_core_sg, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Extreme Scramble (mnpq * mn -> pq) ---
  write(6,*) '>>> TEST 3: Extreme Scramble (mnpq * mn -> pq, IOP=1) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1) 
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_sg, icore_sg, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_sg, vmn_sg, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpq_core_sg(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(icore_sg, popm, popn, 1, popp, popq, 1, symvpqmn, 'mnpq', &
                       vmn_sg, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_sg, popp, popq, 1, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 1)
  vpq_ref_sg(1:npq_iter) = 888.0e0
  
  call walltimeso(t_s)
  call sstgenso1(vpqmn_sg, icore_sg, popp, popq, popm, popn, symvpqmn, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore_sg, vmn_sg, symvpqmn, symvmn, popm, popn, 1, popp, popq, 1, popdum, popdum, 1, vpq_ref_sg, symvpq, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpq_core_sg(1:npq_iter) = 888.0e0  
     call vv2v_dpd_422(vpqmn_sg, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vmn_sg, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_sg, popp, popq, 1, symvpq, 'pq', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Intact Pair IOP=0 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (mnpq * mn -> pq, IOP_pq=0) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 0)
  call zeroso(vpq_ref_sg, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_sg, vmn_sg, symvpqmn, symvmn, popm, popn, 1, popp, popq, 0, popdum, popdum, 1, vpq_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_sg(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_sg, popm, popn, 1, popp, popq, 0, symvpqmn, 'mnpq', &
                       vmn_sg, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_sg, popp, popq, 0, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Intact Pair IOP=2 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (mnpq * mn -> pq, IOP_pq=2) [GPU SP] <<<'
  npq_iter = irpdso(symvpq, popp, popq, 2)
  call zeroso(vpq_ref_sg, npq_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpqmn_sg, vmn_sg, symvpqmn, symvmn, popm, popn, 1, popp, popq, 2, popdum, popdum, 1, vpq_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npq_iter, vpq_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpq_core_sg(1:npq_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_sg, popm, popn, 1, popp, popq, 2, symvpqmn, 'mnpq', &
                       vmn_sg, popm, popn, 1, symvmn, 'mn', &
                       vpq_core_sg, popp, popq, 2, symvpq, 'pq', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npq_iter, vpq_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpq_ref_sg, vpq_core_sg, int(npq_iter, 8))

  ! --- GPU SP: Non-Adjacent Contraction Test ---
  write(6,*) '>>> TEST 7: Non-Adjacent Contraction (pqmn * qn -> pm, IOP=1) [GPU SP] <<<'
  npm_iter = irpdso(symvpq, popp, popm, 1)
  call zeroso(vpm_ref_sg, npm_iter)

  call walltimeso(t_s)
  call sstgenso1(vpqmn_sg, icore_sg, popp, popq, popm, popn, symvpqmn, '2413')
  call VmnpqVmnrs_to_Vpqrs(icore_sg, vqn_sg, symvpqmn, symvmn, popq, popn, 1, popp, popm, 1, popdum, popdum, 1, vpm_ref_sg, symvpq, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npm_iter, vpm_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpm_core_sg(1:npm_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_422(vpqmn_sg, popp, popq, 1, popm, popn, 1, symvpqmn, 'pqmn', &
                       vqn_sg, popq, popn, 1, symvmn, 'qn', &
                       vpm_core_sg, popp, popm, 1, symvpq, 'pm', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npm_iter, vpm_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpm_ref_sg, vpm_core_sg, int(npm_iter, 8))

end subroutine vv2v_dpd_422_check
