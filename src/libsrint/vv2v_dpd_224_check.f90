!=======================================================================
!  subroutine: vv2v_dpd_224_check
!
!  功能说明:
!    针对 2指标 * 2指标 -> 4指标 (外积张量生成) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    【核心机制】: 
!       1. 验证 224 Wrapper 虚拟扩展调用 444 引擎的绝对准确性。
!       2. 利用 VmnpqVmnrs_to_Vpqrs 配合 popm=popn=popdum 策略作为绝对 Reference。
!       3. 重点测试 V1 和 V2 在各种 IOP 压缩态 (0, 1, 2) 下的正确映射。
!=======================================================================
subroutine vv2v_dpd_224_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npq, nrs, nqp, nsr
  integer :: npqrs, npqrs_iter  
  integer :: symvpq, symvrs, symvpqrs
  integer :: popp(8), popq(8), popr(8), pops(8), popdum(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter, i
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_ioprs

  ! --- 外部函数 ---
  integer, external :: irpdso, idsymsoc
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpq(:), vrs(:), vpqrs_core(:), vpqrs_ref(:)
  real*8, allocatable :: vqp(:), vsr(:), icore_pq(:), icore_rs(:)
  
  real*4, allocatable :: vpq_s(:), vrs_s(:), vpqrs_core_s(:), vpqrs_ref_s(:)
  real*4, allocatable :: vqp_s(:), vsr_s(:), icore_pq_s(:), icore_rs_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpq_g(:), vrs_g(:), vpqrs_core_g(:), vpqrs_ref_g(:)
  real*8, allocatable, device :: vqp_g(:), vsr_g(:), icore_pq_g(:), icore_rs_g(:)
  
  real*4, allocatable, device :: vpq_sg(:), vrs_sg(:), vpqrs_core_sg(:), vpqrs_ref_sg(:)
  real*4, allocatable, device :: vqp_sg(:), vsr_sg(:), icore_pq_sg(:), icore_rs_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpq, symvrs, perm_in)
  
  call init_parameters(npq, nrs, npqrs, &
                       symvpq, symvrs, symvpqrs, &
                       popp, popq, popr, pops, popdum, popdum, &
                       def_ioppq, def_ioprs, def_ioppq, &
                       nirrep_in, nvrt_in, npop_in) 

  ! 【关键】：重置并构建恒为 1 的 Dummy 维度，用于 Reference 欺骗
  popdum = 0
  popdum(1) = 1

  ! =====================================================================
  ! 严谨的内存管理：利用 DPD 物理规则精确计算各个张量的大小
  ! =====================================================================
  npq = irpdso(symvpq, popp, popq, 1)
  nrs = irpdso(symvrs, popr, pops, 1)
  nqp = irpdso(symvpq, popq, popp, 1)
  nsr = irpdso(symvrs, pops, popr, 1)
  npqrs = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_224_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  write(6,*) " [ Memory Allocation Info (Base DP / 8 Bytes) ] "
  write(6,"(A15, I15, A, F10.2, A)") "    V1(pq)   : ", npq, " elements | ", real(npq, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2(rs)   : ", nrs, " elements | ", real(nrs, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V3(pqrs) : ", npqrs, " elements | ", real(npqrs, 8)*8.0d0/1048576.0d0, " MB"
  write(6,*) "=========================================================="

  ! --- 分配缓冲 ---
  allocate(vpq(npq), vrs(nrs), vpqrs_core(npqrs), vpqrs_ref(npqrs))
  allocate(vqp(nqp), vsr(nsr), icore_pq(npq), icore_rs(nrs))
  
  allocate(vpq_s(npq), vrs_s(nrs), vpqrs_core_s(npqrs), vpqrs_ref_s(npqrs))
  allocate(vqp_s(nqp), vsr_s(nsr), icore_pq_s(npq), icore_rs_s(nrs))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpq, vrs, npq, nrs)
  
  ! 初始化独立的乱序矩阵
  do i = 1, nqp; vqp(i) = sin(real(i, 8)) * 100.0d0; end do
  do i = 1, nsr; vsr(i) = cos(real(i, 8)) * 100.0d0; end do
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpq, vpq_s, npq)
  call dtos(vrs, vrs_s, nrs)
  call dtos(vqp, vqp_s, nqp)
  call dtos(vsr, vsr_s, nsr)


  !======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  ! --- CPU DP: Standard 224 (pq * rs -> pqrs) ---
  write(6,*) '>>> TEST 1: Standard Mode (pq * rs -> pqrs, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  
  call walltimeso(t_s)
  ! Ref: 利用 popdum 作为收缩维度 m,n 进行无开销外积
  call VmnpqVmnrs_to_Vpqrs(vpq, vrs, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq, popp, popq, 1, symvpq, 'pq', &
                       vrs, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Swapped Inputs (qp * sr -> pqrs) ---
  write(6,*) '>>> TEST 2: Swapped Inputs (qp * sr -> pqrs, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  
  call walltimeso(t_s)
  ! 将 (q,p,dum,dum) 转置为 (p,q,dum,dum) 送入 Ref -> 对应 2134 映射
  call sstgenso1(vqp, icore_pq, popq, popp, popdum, popdum, symvpq, '2134')
  call sstgenso1(vsr, icore_rs, pops, popr, popdum, popdum, symvrs, '2134')
  call VmnpqVmnrs_to_Vpqrs(icore_pq, icore_rs, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vqp, popq, popp, 1, symvpq, 'qp', &
                       vsr, pops, popr, 1, symvrs, 'sr', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref(1:npqrs_iter) = 888.0d0
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq, vrs, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core(1:npqrs_iter) = 888.0d0  
     call vv2v_dpd_224(vpq, popp, popq, 1, symvpq, 'pq', &
                       vrs, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Compression Test 1 (pq:0, rs:1) ---
  write(6,*) '>>> TEST 4: Compression Mix 1 (pq * rs -> pqrs, IOP_pq=0, IOP_rs=1) [DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq, vrs, symvpq, symvrs, popdum, popdum, 1, popp, popq, 0, popr, pops, 1, vpqrs_ref, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq, popp, popq, 0, symvpq, 'pq', &
                       vrs, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core, popp, popq, 0, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Compression Test 2 (pq:2, rs:2) ---
  write(6,*) '>>> TEST 5: Full Compression (pq * rs -> pqrs, IOP_pq=2, IOP_rs=2) [DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 2)
  call zeroso(vpqrs_ref, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq, vrs, symvpq, symvrs, popdum, popdum, 1, popp, popq, 2, popr, pops, 2, vpqrs_ref, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq, popp, popq, 2, symvpq, 'pq', &
                       vrs, popr, pops, 2, symvrs, 'rs', &
                       vpqrs_core, popp, popq, 2, popr, pops, 2, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))


  !======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  ! --- CPU SP: Standard 224 ---
  write(6,*) '>>> TEST 1: Standard Mode (pq * rs -> pqrs, IOP=1) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_s, vrs_s, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_s, popp, popq, 1, symvpq, 'pq', &
                       vrs_s, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Swapped Inputs ---
  write(6,*) '>>> TEST 2: Swapped Inputs (qp * sr -> pqrs, IOP=1) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vqp_s, icore_pq_s, popq, popp, popdum, popdum, symvpq, '2134')
  call sstgenso1(vsr_s, icore_rs_s, pops, popr, popdum, popdum, symvrs, '2134')
  call VmnpqVmnrs_to_Vpqrs(icore_pq_s, icore_rs_s, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vqp_s, popq, popp, 1, symvpq, 'qp', &
                       vsr_s, pops, popr, 1, symvrs, 'sr', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_s(1:npqrs_iter) = 888.0e0
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_s, vrs_s, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_s(1:npqrs_iter) = 888.0e0  
     call vv2v_dpd_224(vpq_s, popp, popq, 1, symvpq, 'pq', &
                       vrs_s, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Compression Test 2 ---
  write(6,*) '>>> TEST 5: Full Compression (pq * rs -> pqrs, IOP_pq=2, IOP_rs=2) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 2)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_s, vrs_s, symvpq, symvrs, popdum, popdum, 1, popp, popq, 2, popr, pops, 2, vpqrs_ref_s, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_s, popp, popq, 2, symvpq, 'pq', &
                       vrs_s, popr, pops, 2, symvrs, 'rs', &
                       vpqrs_core_s, popp, popq, 2, popr, pops, 2, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))


  !======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpq_g(npq), vrs_g(nrs), vpqrs_core_g(npqrs), vpqrs_ref_g(npqrs))
  allocate(vqp_g(nqp), vsr_g(nsr), icore_pq_g(npq), icore_rs_g(nrs))
  
  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vpq, vpq_g, npq)
  call myicopy_g(vrs, vrs_g, nrs)
  call myicopy_g(vqp, vqp_g, nqp)
  call myicopy_g(vsr, vsr_g, nsr)

  ! --- GPU DP: Standard 224 ---
  write(6,*) '>>> TEST 1: Standard Mode (pq * rs -> pqrs, IOP=1) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_g, vrs_g, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_g, popp, popq, 1, symvpq, 'pq', &
                       vrs_g, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Swapped Inputs ---
  write(6,*) '>>> TEST 2: Swapped Inputs (qp * sr -> pqrs, IOP=1) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  
  call walltimeso(t_s)
  ! 【严格传入 _g 显存指针，杜绝 Error 700】
  call sstgenso1(vqp_g, icore_pq_g, popq, popp, popdum, popdum, symvpq, '2134')
  call sstgenso1(vsr_g, icore_rs_g, pops, popr, popdum, popdum, symvrs, '2134')
  call VmnpqVmnrs_to_Vpqrs(icore_pq_g, icore_rs_g, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2) 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vqp_g, popq, popp, 1, symvpq, 'qp', &
                       vsr_g, pops, popr, 1, symvrs, 'sr', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_g(1:npqrs_iter) = 888.0d0
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_g, vrs_g, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_g(1:npqrs_iter) = 888.0d0  
     call vv2v_dpd_224(vpq_g, popp, popq, 1, symvpq, 'pq', &
                       vrs_g, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Compression Test 1 ---
  write(6,*) '>>> TEST 4: Compression Mix 1 (pq * rs -> pqrs, IOP_pq=0, IOP_rs=1) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_g, vrs_g, symvpq, symvrs, popdum, popdum, 1, popp, popq, 0, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_g, popp, popq, 0, symvpq, 'pq', &
                       vrs_g, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_g, popp, popq, 0, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))


  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpq_sg(npq), vrs_sg(nrs), vpqrs_core_sg(npqrs), vpqrs_ref_sg(npqrs))
  allocate(vqp_sg(nqp), vsr_sg(nsr), icore_pq_sg(npq), icore_rs_sg(nrs))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpq_s, vpq_sg, npq)
  call myicopy_g(vrs_s, vrs_sg, nrs)
  call myicopy_g(vqp_s, vqp_sg, nqp)
  call myicopy_g(vsr_s, vsr_sg, nsr)

  ! --- GPU SP: Standard 224 ---
  write(6,*) '>>> TEST 1: Standard Mode (pq * rs -> pqrs, IOP=1) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_sg, vrs_sg, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_sg, popp, popq, 1, symvpq, 'pq', &
                       vrs_sg, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Swapped Inputs ---
  write(6,*) '>>> TEST 2: Swapped Inputs (qp * sr -> pqrs, IOP=1) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vqp_sg, icore_pq_sg, popq, popp, popdum, popdum, symvpq, '2134')
  call sstgenso1(vsr_sg, icore_rs_sg, pops, popr, popdum, popdum, symvrs, '2134')
  call VmnpqVmnrs_to_Vpqrs(icore_pq_sg, icore_rs_sg, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vqp_sg, popq, popp, 1, symvpq, 'qp', &
                       vsr_sg, pops, popr, 1, symvrs, 'sr', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_sg(1:npqrs_iter) = 888.0e0
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_sg, vrs_sg, symvpq, symvrs, popdum, popdum, 1, popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_sg(1:npqrs_iter) = 888.0e0  
     call vv2v_dpd_224(vpq_sg, popp, popq, 1, symvpq, 'pq', &
                       vrs_sg, popr, pops, 1, symvrs, 'rs', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, 'pqrs', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Compression Test 2 ---
  write(6,*) '>>> TEST 5: Full Compression (pq * rs -> pqrs, IOP_pq=2, IOP_rs=2) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 2, popr, pops, 2)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(vpq_sg, vrs_sg, symvpq, symvrs, popdum, popdum, 1, popp, popq, 2, popr, pops, 2, vpqrs_ref_sg, symvpqrs, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_224(vpq_sg, popp, popq, 2, symvpq, 'pq', &
                       vrs_sg, popr, pops, 2, symvrs, 'rs', &
                       vpqrs_core_sg, popp, popq, 2, popr, pops, 2, symvpqrs, 'pqrs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

end subroutine vv2v_dpd_224_check
