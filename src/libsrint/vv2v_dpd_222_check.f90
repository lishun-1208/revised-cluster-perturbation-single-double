!=======================================================================
!  subroutine: vv2v_dpd_222_check
!
!  功能说明:
!    针对 2指标 * 2指标 -> 2指标 (矩阵乘法) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    【核心机制】: 
!       1. 验证 222 Wrapper 调用 444 引擎的绝对准确性。
!       2. 依旧利用 VmnpqVmnrs_to_Vpqrs 配合 pdum(1)=1 策略作为绝对 Reference。
!       3. 严禁 IOP=0/2 压缩态 (222 模式强制要求 IOP=1 完全展开)。
!=======================================================================
subroutine vv2v_dpd_222_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npm, nmr, nrm, npr
  integer :: npr_iter  
  integer :: symvpm, symvmr, symvpr
  integer :: popp(8), popm(8), popr(8), popdum(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter, i
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_iopmn

  ! --- 外部函数 ---
  integer, external :: irpdso
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpm(:), vmr(:), vrm(:), vpr_core(:), vpr_ref(:)
  real*8, allocatable :: icore1(:), icore2(:)
  
  real*4, allocatable :: vpm_s(:), vmr_s(:), vrm_s(:), vpr_core_s(:), vpr_ref_s(:)
  real*4, allocatable :: icore1_s(:), icore2_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpm_g(:), vmr_g(:), vrm_g(:), vpr_core_g(:), vpr_ref_g(:)
  real*8, allocatable, device :: icore1_g(:), icore2_g(:)
  
  real*4, allocatable, device :: vpm_sg(:), vmr_sg(:), vrm_sg(:), vpr_core_sg(:), vpr_ref_sg(:)
  real*4, allocatable, device :: icore1_sg(:), icore2_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpm, symvmr, perm_in)
  
  call init_parameters(npm, nmr, npr, &
                       symvpm, symvmr, symvpr, &
                       popp, popm, popr, popdum, popdum, popdum, &
                       def_iopmn, def_ioppq, def_iopmn, &
                       nirrep_in, nvrt_in, npop_in) 

  ! 【关键】：重置并构建恒为 1 的 Dummy 维度
  popdum = 0
  popdum(1) = 1

  ! =====================================================================
  ! 严谨的内存管理：利用 DPD 物理规则精确计算各个张量的大小
  ! =====================================================================
  npm = irpdso(symvpm, popp, popm, 1)
  nmr = irpdso(symvmr, popm, popr, 1)
  npr = irpdso(symvpr, popp, popr, 1)
  nrm = irpdso(symvmr, popr, popm, 1)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_222_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  write(6,*) " [ Memory Allocation Info (Base DP / 8 Bytes) ] "
  write(6,"(A15, I15, A, F10.2, A)") "    V1(pm)   : ", npm, " elements | ", real(npm, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2(mr)   : ", nmr, " elements | ", real(nmr, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V2_sw(rm): ", nrm, " elements | ", real(nrm, 8)*8.0d0/1048576.0d0, " MB"
  write(6,"(A15, I15, A, F10.2, A)") "    V3(pr)   : ", npr, " elements | ", real(npr, 8)*8.0d0/1048576.0d0, " MB"
  write(6,*) "=========================================================="

  ! --- 分配缓冲 ---
  allocate(vpm(npm), vmr(nmr), vrm(nrm), vpr_core(npr), vpr_ref(npr))
  allocate(icore1(npm), icore2(nmr))
  
  allocate(vpm_s(npm), vmr_s(nmr), vrm_s(nrm), vpr_core_s(npr), vpr_ref_s(npr))
  allocate(icore1_s(npm), icore2_s(nmr))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpm, vmr, npm, nmr)
  
  ! 初始化独立的 vrm 矩阵用于测试乱序收缩
  do i = 1, nrm
     vrm(i) = sin(real(i, 8)) * 100.0d0
  end do
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpm, vpm_s, npm)
  call dtos(vmr, vmr_s, nmr)
  call dtos(vrm, vrm_s, nrm)


  !======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  ! --- CPU DP: Standard 222 (pm * mr -> pr) ---
  write(6,*) '>>> TEST 1: Standard Mode (pm * mr -> pr, IOP=1) <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref, npr_iter)
  
  call walltimeso(t_s)
  ! V1: 物理内存是 (p, dum, m, dum)，要送进 mnpq 接口，需要转置成 (m, dum, p, dum)，对应 3412
  call sstgenso1(vpm, icore1, popp, popdum, popm, popdum, symvpm, '3412')
  ! V2: 物理内存是 (m, dum, r, dum)，刚好匹配参考程序的 mnrs，直接传入 vmr
  call VmnpqVmnrs_to_Vpqrs(icore1, vmr, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npr_iter, vpr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpr_core(1:npr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm, popp, popm, 1, symvpm, 'pm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpr_core, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npr_iter, vpr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref, vpr_core, int(npr_iter, 8))

  ! --- CPU DP: Swapped 222 (pm * rm -> pr) ---
  write(6,*) '>>> TEST 2: Swapped Contraction (pm * rm -> pr, IOP=1) <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm, icore1, popp, popdum, popm, popdum, symvpm, '3412')
  ! vrm 物理内存是 (r, dum, m, dum)，参考程序 V2 必须是 (m, dum, r, dum)，转置 3412
  call sstgenso1(vrm, icore2, popr, popdum, popm, popdum, symvmr, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1, icore2, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npr_iter, vpr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpr_core(1:npr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm, popp, popm, 1, symvpm, 'pm', &
                       vrm, popr, popm, 1, symvmr, 'rm', &
                       vpr_core, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npr_iter, vpr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref, vpr_core, int(npr_iter, 8))

  ! --- CPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [DP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  vpr_ref(1:npr_iter) = 888.0d0
  
  call walltimeso(t_s)
  call sstgenso1(vpm, icore1, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1, vmr, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref, symvpr, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npr_iter, vpr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpr_core(1:npr_iter) = 888.0d0  
     call vv2v_dpd_222(vpm, popp, popm, 1, symvpm, 'pm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpr_core, popp, popr, 1, symvpr, 'pr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npr_iter, vpr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref, vpr_core, int(npr_iter, 8))


  !======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  ! --- CPU SP: Standard 222 ---
  write(6,*) '>>> TEST 1: Standard Mode (pm * mr -> pr, IOP=1) [SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_s, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm_s, icore1_s, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_s, vmr_s, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_s, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpr_core_s(1:npr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_s, popp, popm, 1, symvpm, 'pm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_s, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_s, vpr_core_s, int(npr_iter, 8))

  ! --- CPU SP: Swapped 222 ---
  write(6,*) '>>> TEST 2: Swapped Contraction (pm * rm -> pr, IOP=1) [SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_s, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm_s, icore1_s, popp, popdum, popm, popdum, symvpm, '3412')
  call sstgenso1(vrm_s, icore2_s, popr, popdum, popm, popdum, symvmr, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_s, icore2_s, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_s, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpr_core_s(1:npr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_s, popp, popm, 1, symvpm, 'pm', &
                       vrm_s, popr, popm, 1, symvmr, 'rm', &
                       vpr_core_s, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_s, vpr_core_s, int(npr_iter, 8))

  ! --- CPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  vpr_ref_s(1:npr_iter) = 888.0e0
  
  call walltimeso(t_s)
  call sstgenso1(vpm_s, icore1_s, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_s, vmr_s, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_s, symvpr, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpr_core_s(1:npr_iter) = 888.0e0  
     call vv2v_dpd_222(vpm_s, popp, popm, 1, symvpm, 'pm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_s, popp, popr, 1, symvpr, 'pr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_s, vpr_core_s, int(npr_iter, 8))


  !======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpm_g(npm), vmr_g(nmr), vrm_g(nrm), vpr_core_g(npr), vpr_ref_g(npr))
  allocate(icore1_g(npm), icore2_g(nmr))
  
  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vpm, vpm_g, npm)
  call myicopy_g(vmr, vmr_g, nmr)
  call myicopy_g(vrm, vrm_g, nrm)

  ! --- GPU DP: Standard 222 ---
  write(6,*) '>>> TEST 1: Standard Mode (pm * mr -> pr, IOP=1) [GPU DP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_g, npr_iter)
  
  call walltimeso(t_s)
  ! 【严格传入 _g 显存指针，杜绝 Error 700】
  call sstgenso1(vpm_g, icore1_g, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_g, vmr_g, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_g, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpr_core_g(1:npr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_g, popp, popm, 1, symvpm, 'pm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_g, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_g, vpr_core_g, int(npr_iter, 8))

  ! --- GPU DP: Swapped 222 ---
  write(6,*) '>>> TEST 2: Swapped Contraction (pm * rm -> pr, IOP=1) [GPU DP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_g, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm_g, icore1_g, popp, popdum, popm, popdum, symvpm, '3412')
  call sstgenso1(vrm_g, icore2_g, popr, popdum, popm, popdum, symvmr, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_g, icore2_g, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_g, symvpr, fact, fact2) 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpr_core_g(1:npr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_g, popp, popm, 1, symvpm, 'pm', &
                       vrm_g, popr, popm, 1, symvmr, 'rm', &
                       vpr_core_g, popp, popr, 1, symvpr, 'pr', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_g, vpr_core_g, int(npr_iter, 8))

  ! --- GPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [GPU DP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  vpr_ref_g(1:npr_iter) = 888.0d0
  
  call walltimeso(t_s)
  call sstgenso1(vpm_g, icore1_g, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_g, vmr_g, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_g, symvpr, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpr_core_g(1:npr_iter) = 888.0d0  
     call vv2v_dpd_222(vpm_g, popp, popm, 1, symvpm, 'pm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_g, popp, popr, 1, symvpr, 'pr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_g, vpr_core_g, int(npr_iter, 8))


  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpm_sg(npm), vmr_sg(nmr), vrm_sg(nrm), vpr_core_sg(npr), vpr_ref_sg(npr))
  allocate(icore1_sg(npm), icore2_sg(nmr))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpm_s, vpm_sg, npm)
  call myicopy_g(vmr_s, vmr_sg, nmr)
  call myicopy_g(vrm_s, vrm_sg, nrm)

  ! --- GPU SP: Standard 222 ---
  write(6,*) '>>> TEST 1: Standard Mode (pm * mr -> pr, IOP=1) [GPU SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_sg, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm_sg, icore1_sg, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_sg, vmr_sg, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_sg, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpr_core_sg(1:npr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_sg, popp, popm, 1, symvpm, 'pm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_sg, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_sg, vpr_core_sg, int(npr_iter, 8))

  ! --- GPU SP: Swapped 222 ---
  write(6,*) '>>> TEST 2: Swapped Contraction (pm * rm -> pr, IOP=1) [GPU SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  call zeroso(vpr_ref_sg, npr_iter)
  
  call walltimeso(t_s)
  call sstgenso1(vpm_sg, icore1_sg, popp, popdum, popm, popdum, symvpm, '3412')
  call sstgenso1(vrm_sg, icore2_sg, popr, popdum, popm, popdum, symvmr, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_sg, icore2_sg, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_sg, symvpr, fact, fact2)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpr_core_sg(1:npr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_222(vpm_sg, popp, popm, 1, symvpm, 'pm', &
                       vrm_sg, popr, popm, 1, symvmr, 'rm', &
                       vpr_core_sg, popp, popr, 1, symvpr, 'pr', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_sg, vpr_core_sg, int(npr_iter, 8))

  ! --- GPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 3: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  npr_iter = irpdso(symvpr, popp, popr, 1)
  vpr_ref_sg(1:npr_iter) = 888.0e0
  
  call walltimeso(t_s)
  call sstgenso1(vpm_sg, icore1_sg, popp, popdum, popm, popdum, symvpm, '3412')
  call VmnpqVmnrs_to_Vpqrs(icore1_sg, vmr_sg, symvpm, symvmr, popm, popdum, 1, popp, popdum, 1, popr, popdum, 1, vpr_ref_sg, symvpr, 1.0d0, 0.5d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npr_iter, vpr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpr_core_sg(1:npr_iter) = 888.0e0  
     call vv2v_dpd_222(vpm_sg, popp, popm, 1, symvpm, 'pm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpr_core_sg, popp, popr, 1, symvpr, 'pr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npr_iter, vpr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpr_ref_sg, vpr_core_sg, int(npr_iter, 8))

end subroutine vv2v_dpd_222_check
