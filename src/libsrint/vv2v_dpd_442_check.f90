!=======================================================================
!  subroutine: vv2v_dpd_442_check
!
!  功能说明:
!    针对 4指标 * 4指标 -> 2指标 (收缩 3 个指标) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    重构说明:
!      1. 严格遵循 DPD 物理约束，去除了所有违背空间同构性的伪测试。
!      2. 新增 TEST 10: 专门验证 Macro-RL 交叉宏维在压缩态 (IOP=0) 下的收缩正确性，
!         通过 iexpso 全展开降维 0.5 倍法生成绝对可靠的 Reference。
!=======================================================================
subroutine vv2v_dpd_442_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npqrm, npqrn, nmn
  integer :: nmn_iter  
  integer :: symvpqrm, symvpqrn, symvmn
  integer :: popp(8), popq(8), popr(8), popm(8), popn(8), popdum(8)
  
  ! 用于 3指标收缩 的维度变量
  integer :: pops(8)
  integer :: npqms, nrs, nrs_iter
  integer :: symvpqms, symvrs

  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_ioprm, def_ioprn
  
  ! 大内存评估 (用于重排缓存与解压)
  integer :: maxcore_exp

  ! --- 外部函数 ---
  integer, external :: idsymsoc, irpdso
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpqrm(:), vpqrn(:), vmn_core(:), vmn_ref(:)
  real*8, allocatable :: icore1(:), icore2(:)
  real*4, allocatable :: vpqrm_s(:), vpqrn_s(:), vmn_core_s(:), vmn_ref_s(:)
  real*4, allocatable :: icore1_s(:), icore2_s(:)

  real*8, allocatable :: vpqms(:), vrs_core(:), vrs_ref(:)
  real*4, allocatable :: vpqms_s(:), vrs_core_s(:), vrs_ref_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpqrm_g(:), vpqrn_g(:), vmn_core_g(:), vmn_ref_g(:)
  real*8, allocatable, device :: icore1_g(:), icore2_g(:)
  real*4, allocatable, device :: vpqrm_sg(:), vpqrn_sg(:), vmn_core_sg(:), vmn_ref_sg(:)
  real*4, allocatable, device :: icore1_sg(:), icore2_sg(:)

  real*8, allocatable, device :: vpqms_g(:), vrs_core_g(:), vrs_ref_g(:)
  real*4, allocatable, device :: vpqms_sg(:), vrs_core_sg(:), vrs_ref_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpqrm, symvpqrn, perm_in)
  
  call init_parameters(npqrm, npqrn, nmn, &
                       symvpqrm, symvpqrn, symvmn, &
                       popp, popq, popr, popm, popn, popdum, &
                       def_ioppq, def_ioprm, def_ioprn, &
                       nirrep_in, nvrt_in, npop_in) 

  ! 计算分配大小 (默认全解压 IOP=1 态)
  write(6,*)'pop', popp
  write(6,*)'pop', popq
  popq = popp
  npqrm = idsymsoc(symvpqrm, popp, popq, 1, popr, popm, 1)
  npqrn = idsymsoc(symvpqrn, popp, popq, 1, popr, popn, 1)
  nmn   = irpdso(symvmn, popm, popn, 1)
  maxcore_exp = npqrm + npqrn

  ! 初始化 3指标收缩 相关维度
  pops = popn
  symvpqms = symvpqrn
  symvrs = dirprd(symvpqrm, symvpqms)
  npqms = idsymsoc(symvpqms, popp, popq, 1, popm, pops, 1)
  nrs   = irpdso(symvrs, popr, pops, 1)
  nrs_iter = nrs

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_442_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="

  allocate(vpqrm(npqrm), vpqrn(npqrn), vmn_core(nmn), vmn_ref(nmn))
  allocate(icore1(maxcore_exp), icore2(maxcore_exp))
  allocate(vpqrm_s(npqrm), vpqrn_s(npqrn), vmn_core_s(nmn), vmn_ref_s(nmn))
  allocate(icore1_s(maxcore_exp), icore2_s(maxcore_exp))

  allocate(vpqms(npqms), vrs_core(nrs), vrs_ref(nrs))
  allocate(vpqms_s(npqms), vrs_core_s(nrs), vrs_ref_s(nrs))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpqrm, vpqrn, npqrm, npqrn)
  call init_integrals_dp(vpqms, vpqms, npqms, npqms)
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpqrm, vpqrm_s, npqrm)
  call dtos(vpqrn, vpqrn_s, npqrn)
  call dtos(vpqms, vpqms_s, npqms)


  !======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  write(6,*) '>>> TEST 1: Standard Mode (pqrm * pqrn -> mn, IOP=1) <<<'
  nmn_iter = irpdso(symvmn, popm, popn, 1)
  call zeroso(vmn_ref, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  vmn_core(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))

  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmr * pqnr -> mn, IOP=1) <<<'
  call zeroso(vmn_ref, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm, icore1, popp, popq, popr, popm, symvpqrm, '1243')
  call sstgenso1(vpqrn, icore2, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1, popp, popq, 1, popm, popr, 1, symvpqrm, 'pqmr', &
                       icore2, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))

  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * rnpq -> mn, IOP=1) <<<'
  call zeroso(vmn_ref, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm, icore1, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn, icore2, popp, popq, popr, popn, symvpqrn, '3412')
  vmn_core(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2, popr, popn, 1, popp, popq, 1, symvpqrn, 'rnpq', &
                       vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))

  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [DP] <<<'
  vmn_ref(1:nmn_iter) = 444.0d0
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref, symvmn, 1.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  call walltimeso(t_s)
  do iter = 1, ncycle
     vmn_core(1:nmn_iter) = 888.0d0  
     call vv2v_dpd_442(vpqrm, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core, popm, popn, 1, symvmn, 'mn', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))

  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * pqrn -> mn, IOP_pq=0) [DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(vmn_ref, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 0, popr, popm, popn, vmn_ref, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn, popp, popq, 0, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * pqrn -> mn, IOP_pq=2) [DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=2 for different orbital spaces.'
  else
      call zeroso(vmn_ref, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 2, popr, popm, popn, vmn_ref, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn, popp, popq, 2, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 7: 3-Index Contraction (pqrm * pqms -> rs, IOP=1) [DP] <<<'
  call zeroso(vrs_ref, nrs_iter)
  call sstgenso1(vpqrm, icore1, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1, vpqms, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nrs_iter, vrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  vrs_core(1:nrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqms, popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nrs_iter, vrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref, vrs_core, int(nrs_iter, 8))

  write(6,*) '>>> TEST 8: Scramble Contract (pmrq * pqms -> rs, IOP=1) [DP] <<<'
  call zeroso(vrs_ref, nrs_iter)
  call sstgenso1(vpqrm, icore1, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1, vpqms, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nrs_iter, vrs_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm, icore2, popp, popq, popr, popm, symvpqrm, '1432')
  vrs_core(1:nrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore2, popp, popm, 1, popr, popq, 1, symvpqrm, 'pmrq', &
                       vpqms,  popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nrs_iter, vrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref, vrs_core, int(nrs_iter, 8))

  write(6,*) '>>> TEST 9: Macro-RL Routing (rmpq * pqnr -> mn, IOP=1) [DP] <<<'
  call zeroso(vmn_ref, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm, vpqrn, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm, icore1, popp, popq, popr, popm, symvpqrm, '3412') 
  call sstgenso1(vpqrn, icore2, popp, popq, popr, popn, symvpqrn, '1243') 
  vmn_core(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))

  ! --- 【新增】CPU DP: Macro-RL Routing with IOP=0 ---
  write(6,*) '>>> TEST 10: Macro-RL Routing (rmpq * pqnr -> mn, IOP_pq=0) [DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
      write(6,*)
  else
      call zeroso(icore1, maxcore_exp)
      call zeroso(icore2, maxcore_exp)
      ! 将随机数缓冲当作 IOP=0 数组，全展开为 IOP=1 用于生成绝对 Reference
      call iexpso(vpqrm, icore1, popr, popm, 1, 0, popp, popq, 0, 1, symvpqrm)
      call iexpso(vpqrn, icore2, popp, popq, 0, 1, popn, popr, 1, 0, symvpqrn)
      
      call zeroso(vmn_ref, nmn_iter)
      call walltimeso(t_s)
      call vv2v_dpd_442(icore1, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                        icore2, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                        vmn_ref, popm, popn, 1, symvmn, 'mn', fact, fact2)
      call walltimeso(t_e)
      ! 原理: 全展开(IOP=1)的反对称收缩恰好是压缩(IOP=0)收缩的 2 倍
      vmn_ref(1:nmn_iter) = vmn_ref(1:nmn_iter) * 0.5d0
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(nmn_iter, vmn_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         ! 直接测试底层引擎的原生压缩路由
         call vv2v_dpd_442(vpqrm, popr, popm, 1, popp, popq, 0, symvpqrm, 'rmpq', &
                           vpqrn, popp, popq, 0, popn, popr, 1, symvpqrn, 'pqnr', &
                           vmn_core, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(nmn_iter, vmn_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref, vmn_core, int(nmn_iter, 8))
  end if


  !======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  write(6,*) '>>> TEST 1: Standard Mode (pqrm * pqrn -> mn, IOP=1) [SP] <<<'
  call zeroso(vmn_ref_s, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_s, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  vmn_core_s(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_s, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))

  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmr * pqnr -> mn, IOP=1) [SP] <<<'
  call zeroso(vmn_ref_s, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_s, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_s, icore1_s, popp, popq, popr, popm, symvpqrm, '1243')
  call sstgenso1(vpqrn_s, icore2_s, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_s(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_s, popp, popq, 1, popm, popr, 1, symvpqrm, 'pqmr', &
                       icore2_s, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))

  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * rnpq -> mn, IOP=1) [SP] <<<'
  call zeroso(vmn_ref_s, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_s, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_s, icore1_s, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_s, icore2_s, popp, popq, popr, popn, symvpqrn, '3412')
  vmn_core_s(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_s, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_s, popr, popn, 1, popp, popq, 1, symvpqrn, 'rnpq', &
                       vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))

  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [SP] <<<'
  vmn_ref_s(1:nmn_iter) = 444.0e0
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_s, symvmn, 1.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  call walltimeso(t_s)
  do iter = 1, ncycle
     vmn_core_s(1:nmn_iter) = 888.0e0  
     call vv2v_dpd_442(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_s, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_s, popm, popn, 1, symvmn, 'mn', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))

  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * pqrn -> mn, IOP_pq=0) [SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(vmn_ref_s, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 0, popr, popm, popn, vmn_ref_s, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
      vmn_core_s(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_s, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_s, popp, popq, 0, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * pqrn -> mn, IOP_pq=2) [SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=2 for different orbital spaces.'
  else
      call zeroso(vmn_ref_s, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 2, popr, popm, popn, vmn_ref_s, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
      vmn_core_s(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_s, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_s, popp, popq, 2, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 7: 3-Index Contraction (pqrm * pqms -> rs, IOP=1) [SP] <<<'
  call zeroso(vrs_ref_s, nrs_iter)
  call sstgenso1(vpqrm_s, icore1_s, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_s, vpqms_s, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_s, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  vrs_core_s(1:nrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_s, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqms_s, popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_s, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_s, vrs_core_s, int(nrs_iter, 8))

  write(6,*) '>>> TEST 8: Scramble Contract (pmrq * pqms -> rs, IOP=1) [SP] <<<'
  call zeroso(vrs_ref_s, nrs_iter)
  call sstgenso1(vpqrm_s, icore1_s, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_s, vpqms_s, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_s, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_s, icore2_s, popp, popq, popr, popm, symvpqrm, '1432')
  vrs_core_s(1:nrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore2_s, popp, popm, 1, popr, popq, 1, symvpqrm, 'pmrq', &
                       vpqms_s,  popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_s, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_s, vrs_core_s, int(nrs_iter, 8))

  write(6,*) '>>> TEST 9: Macro-RL Routing (rmpq * pqnr -> mn, IOP=1) [SP] <<<'
  call zeroso(vmn_ref_s, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_s, vpqrn_s, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_s, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_s, icore1_s, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_s, icore2_s, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_s(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_s, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_s, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))

  ! --- 【新增】CPU SP: Macro-RL Routing with IOP=0 ---
  write(6,*) '>>> TEST 10: Macro-RL Routing (rmpq * pqnr -> mn, IOP_pq=0) [SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(icore1_s, maxcore_exp)
      call zeroso(icore2_s, maxcore_exp)
      call iexpso(vpqrm_s, icore1_s, popr, popm, 1, 0, popp, popq, 0, 1, symvpqrm)
      call iexpso(vpqrn_s, icore2_s, popp, popq, 0, 1, popn, popr, 1, 0, symvpqrn)
      
      call zeroso(vmn_ref_s, nmn_iter)
      call walltimeso(t_s)
      call vv2v_dpd_442(icore1_s, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                        icore2_s, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                        vmn_ref_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
      call walltimeso(t_e)
      vmn_ref_s(1:nmn_iter) = vmn_ref_s(1:nmn_iter) * 0.5e0
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_s(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_s, popr, popm, 1, popp, popq, 0, symvpqrm, 'rmpq', &
                           vpqrn_s, popp, popq, 0, popn, popr, 1, symvpqrn, 'pqnr', &
                           vmn_core_s, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_s, vmn_core_s, int(nmn_iter, 8))
  end if


  !======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqrm_g(npqrm), vpqrn_g(npqrn), vmn_core_g(nmn), vmn_ref_g(nmn))
  allocate(icore1_g(maxcore_exp), icore2_g(maxcore_exp))
  
  allocate(vpqms_g(npqms), vrs_core_g(nrs), vrs_ref_g(nrs))

  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vpqrm, vpqrm_g, npqrm)
  call myicopy_g(vpqrn, vpqrn_g, npqrn)
  call myicopy_g(vpqms, vpqms_g, npqms)

  write(6,*) '>>> TEST 1: Standard Mode (pqrm * pqrn -> mn, IOP=1) [GPU DP] <<<'
  call zeroso(vmn_ref_g, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_g, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  vmn_core_g(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_g, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_g, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))

  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmr * pqnr -> mn, IOP=1) [GPU DP] <<<'
  call zeroso(vmn_ref_g, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_g, symvmn, fact) 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_g, icore1_g, popp, popq, popr, popm, symvpqrm, '1243')
  call sstgenso1(vpqrn_g, icore2_g, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_g(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_g, popp, popq, 1, popm, popr, 1, symvpqrm, 'pqmr', &
                       icore2_g, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))

  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * rnpq -> mn, IOP=1) [GPU DP] <<<'
  call zeroso(vmn_ref_g, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_g, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_g, icore1_g, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_g, icore2_g, popp, popq, popr, popn, symvpqrn, '3412')
  vmn_core_g(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_g, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_g, popr, popn, 1, popp, popq, 1, symvpqrn, 'rnpq', &
                       vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))

  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [GPU DP] <<<'
  vmn_ref_g(1:nmn_iter) = 444.0d0
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_g, symvmn, 1.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  call walltimeso(t_s)
  do iter = 1, ncycle
     vmn_core_g(1:nmn_iter) = 888.0d0  
     call vv2v_dpd_442(vpqrm_g, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_g, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_g, popm, popn, 1, symvmn, 'mn', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))

  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * pqrn -> mn, IOP_pq=0) [GPU DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(vmn_ref_g, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 0, popr, popm, popn, vmn_ref_g, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_g(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_g, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_g, popp, popq, 0, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * pqrn -> mn, IOP_pq=2) [GPU DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=2 for different orbital spaces.'
  else
      call zeroso(vmn_ref_g, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 2, popr, popm, popn, vmn_ref_g, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_g(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_g, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_g, popp, popq, 2, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 7: 3-Index Contraction (pqrm * pqms -> rs, IOP=1) [GPU DP] <<<'
  call zeroso(vrs_ref_g, nrs_iter)
  call sstgenso1(vpqrm_g, icore1_g, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_g, vpqms_g, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_g, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  vrs_core_g(1:nrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_g, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqms_g, popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_g, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_g, vrs_core_g, int(nrs_iter, 8))

  write(6,*) '>>> TEST 8: Scramble Contract (pmrq * pqms -> rs, IOP=1) [GPU DP] <<<'
  call zeroso(vrs_ref_g, nrs_iter)
  call sstgenso1(vpqrm_g, icore1_g, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_g, vpqms_g, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_g, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_g, icore2_g, popp, popq, popr, popm, symvpqrm, '1432')
  vrs_core_g(1:nrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore2_g, popp, popm, 1, popr, popq, 1, symvpqrm, 'pmrq', &
                       vpqms_g,  popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_g, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_g, vrs_core_g, int(nrs_iter, 8))

  write(6,*) '>>> TEST 9: Macro-RL Routing (rmpq * pqnr -> mn, IOP=1) [GPU DP] <<<'
  call zeroso(vmn_ref_g, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_g, vpqrn_g, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_g, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_g, icore1_g, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_g, icore2_g, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_g(1:nmn_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_g, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_g, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))

  ! --- 【新增】GPU DP: Macro-RL Routing with IOP=0 ---
  write(6,*) '>>> TEST 10: Macro-RL Routing (rmpq * pqnr -> mn, IOP_pq=0) [GPU DP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(icore1_g, maxcore_exp)
      call zeroso(icore2_g, maxcore_exp)
      call iexpso(vpqrm_g, icore1_g, popr, popm, 1, 1, popp, popq, 0, 1, symvpqrm)
      call iexpso(vpqrn_g, icore2_g, popp, popq, 0, 1, popn, popr, 1, 1, symvpqrn)
      
      call zeroso(vmn_ref_g, nmn_iter)
      call walltimeso(t_s)
      call vv2v_dpd_442(icore1_g, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                        icore2_g, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                        vmn_ref_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
      call walltimeso(t_e)
      vmn_ref(1:nmn_iter) = vmn_ref_g(1:nmn_iter)
      vmn_ref(1:nmn_iter) = vmn_ref(1:nmn_iter) * 0.5d0
      vmn_ref_g(1:nmn_iter) = vmn_ref(1:nmn_iter)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_g(1:nmn_iter) = 999.0d0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_g, popr, popm, 1, popp, popq, 0, symvpqrm, 'rmpq', &
                           vpqrn_g, popp, popq, 0, popn, popr, 1, symvpqrn, 'pqnr', &
                           vmn_core_g, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_g, vmn_core_g, int(nmn_iter, 8))
  end if


  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqrm_sg(npqrm), vpqrn_sg(npqrn), vmn_core_sg(nmn), vmn_ref_sg(nmn))
  allocate(icore1_sg(maxcore_exp), icore2_sg(maxcore_exp))
  allocate(vpqms_sg(npqms), vrs_core_sg(nrs), vrs_ref_sg(nrs))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpqrm_s, vpqrm_sg, npqrm)
  call myicopy_g(vpqrn_s, vpqrn_sg, npqrn)
  call myicopy_g(vpqms_s, vpqms_sg, npqms)

  write(6,*) '>>> TEST 1: Standard Mode (pqrm * pqrn -> mn, IOP=1) [GPU SP] <<<'
  call zeroso(vmn_ref_sg, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_sg, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  vmn_core_sg(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_sg, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))

  write(6,*) '>>> TEST 2: Swapped 2-Index (pqmr * pqnr -> mn, IOP=1) [GPU SP] <<<'
  call zeroso(vmn_ref_sg, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_sg, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_sg, icore1_sg, popp, popq, popr, popm, symvpqrm, '1243')
  call sstgenso1(vpqrn_sg, icore2_sg, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_sg(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_sg, popp, popq, 1, popm, popr, 1, symvpqrm, 'pqmr', &
                       icore2_sg, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))

  write(6,*) '>>> TEST 3: Extreme Scramble (rmpq * rnpq -> mn, IOP=1) [GPU SP] <<<'
  call zeroso(vmn_ref_sg, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_sg, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_sg, icore1_sg, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_sg, icore2_sg, popp, popq, popr, popn, symvpqrn, '3412')
  vmn_core_sg(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_sg, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_sg, popr, popn, 1, popp, popq, 1, symvpqrn, 'rnpq', &
                       vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))

  write(6,*) '>>> TEST 4: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  vmn_ref_sg(1:nmn_iter) = 444.0e0
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_sg, symvmn, 1.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  call walltimeso(t_s)
  do iter = 1, ncycle
     vmn_core_sg(1:nmn_iter) = 888.0e0  
     call vv2v_dpd_442(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqrn_sg, popp, popq, 1, popr, popn, 1, symvpqrn, 'pqrn', &
                       vmn_core_sg, popm, popn, 1, symvmn, 'mn', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))

  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqrm * pqrn -> mn, IOP_pq=0) [GPU SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(vmn_ref_sg, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 0, popr, popm, popn, vmn_ref_sg, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_sg(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_sg, popp, popq, 0, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_sg, popp, popq, 0, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqrm * pqrn -> mn, IOP_pq=2) [GPU SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=2 for different orbital spaces.'
  else
      call zeroso(vmn_ref_sg, nmn_iter)
      call walltimeso(t_s)
      call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 2, popr, popm, popn, vmn_ref_sg, symvmn, fact)
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_sg(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_sg, popp, popq, 2, popr, popm, 1, symvpqrm, 'pqrm', &
                           vpqrn_sg, popp, popq, 2, popr, popn, 1, symvpqrn, 'pqrn', &
                           vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))
  end if

  write(6,*) '>>> TEST 7: 3-Index Contraction (pqrm * pqms -> rs, IOP=1) [GPU SP] <<<'
  call zeroso(vrs_ref_sg, nrs_iter)
  call sstgenso1(vpqrm_sg, icore1_sg, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_sg, vpqms_sg, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_sg, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  vrs_core_sg(1:nrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(vpqrm_sg, popp, popq, 1, popr, popm, 1, symvpqrm, 'pqrm', &
                       vpqms_sg, popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_sg, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_sg, vrs_core_sg, int(nrs_iter, 8))

  write(6,*) '>>> TEST 8: Scramble Contract (pmrq * pqms -> rs, IOP=1) [GPU SP] <<<'
  call zeroso(vrs_ref_sg, nrs_iter)
  call sstgenso1(vpqrm_sg, icore1_sg, popp, popq, popr, popm, symvpqrm, '1243')
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(icore1_sg, vpqms_sg, symvpqrm, symvpqms, popp, popq, 1, popm, popr, pops, vrs_ref_sg, symvrs, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nrs_iter, vrs_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_sg, icore2_sg, popp, popq, popr, popm, symvpqrm, '1432')
  vrs_core_sg(1:nrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore2_sg, popp, popm, 1, popr, popq, 1, symvpqrm, 'pmrq', &
                       vpqms_sg,  popp, popq, 1, popm, pops, 1, symvpqms, 'pqms', &
                       vrs_core_sg, popr, pops, 1, symvrs, 'rs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nrs_iter, vrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vrs_ref_sg, vrs_core_sg, int(nrs_iter, 8))

  write(6,*) '>>> TEST 9: Macro-RL Routing (rmpq * pqnr -> mn, IOP=1) [GPU SP] <<<'
  call zeroso(vmn_ref_sg, nmn_iter)
  call walltimeso(t_s)
  call VpqrmVpqrn_to_Vmn(vpqrm_sg, vpqrn_sg, symvpqrm, symvpqrn, popp, popq, 1, popr, popm, popn, vmn_ref_sg, symvmn, fact)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  call sstgenso1(vpqrm_sg, icore1_sg, popp, popq, popr, popm, symvpqrm, '3412')
  call sstgenso1(vpqrn_sg, icore2_sg, popp, popq, popr, popn, symvpqrn, '1243')
  vmn_core_sg(1:nmn_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_442(icore1_sg, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                       icore2_sg, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                       vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))

  ! --- 【新增】GPU SP: Macro-RL Routing with IOP=0 ---
  write(6,*) '>>> TEST 10: Macro-RL Routing (rmpq * pqnr -> mn, IOP_pq=0) [GPU SP] <<<'
  if (any(popp /= popq)) then
      write(6,*) ' -> SKIP: popp /= popq. DPD rules forbid IOP=0 for different orbital spaces.'
  else
      call zeroso(icore1_sg, maxcore_exp)
      call zeroso(icore2_sg, maxcore_exp)
      call iexpso(vpqrm_sg, icore1_sg, popr, popm, 1, 1, popp, popq, 0, 1, symvpqrm)
      call iexpso(vpqrn_sg, icore2_sg, popp, popq, 0, 1, popn, popr, 1, 1, symvpqrn)
      
      call zeroso(vmn_ref_sg, nmn_iter)
      call walltimeso(t_s)
      call vv2v_dpd_442(icore1_sg, popr, popm, 1, popp, popq, 1, symvpqrm, 'rmpq', &
                        icore2_sg, popp, popq, 1, popn, popr, 1, symvpqrn, 'pqnr', &
                        vmn_ref_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
      call walltimeso(t_e)
! 规避 NVFORTRAN S-0519 错误：Device -> Host -> Scale -> Device
      vmn_ref_s(1:nmn_iter) = vmn_ref_sg(1:nmn_iter)
      vmn_ref_s(1:nmn_iter) = vmn_ref_s(1:nmn_iter) * 0.5e0
      vmn_ref_sg(1:nmn_iter) = vmn_ref_s(1:nmn_iter)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(nmn_iter, vmn_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
      
      vmn_core_sg(1:nmn_iter) = 999.0e0
      call walltimeso(t_s)
      do iter = 1, ncycle
         call vv2v_dpd_442(vpqrm_sg, popr, popm, 1, popp, popq, 0, symvpqrm, 'rmpq', &
                           vpqrn_sg, popp, popq, 0, popn, popr, 1, symvpqrn, 'pqnr', &
                           vmn_core_sg, popm, popn, 1, symvmn, 'mn', fact, fact2)
      end do
      call walltimeso(t_e)
      write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(nmn_iter, vmn_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
      call diffv(vmn_ref_sg, vmn_core_sg, int(nmn_iter, 8))
  end if

end subroutine vv2v_dpd_442_check
