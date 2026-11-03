!=======================================================================
!  subroutine: vv2v_dpd_444_check
!
!  功能说明:
!    统一基准测试套件 (动态维度对齐版)。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    【核心升级】: 引入 npqrs_iter 动态计算每次收缩的实际有效空间，
!                  完美兼容 IOP=0 的压缩维度测试，彻底解决脏数据误报问题。
!=======================================================================
subroutine vv2v_dpd_444_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npq, nmn, nrs
  integer :: npqrs, nmnpq, nmnrs
  integer :: npqrs_iter  ! 动态记录当前轮次目标张量的实际大小
  integer :: symvmnpq, symvmnrs, symvpqrs
  integer :: popm(8), popn(8), popp(8), popq(8), popr(8), pops(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_iopmn, def_ioppq, def_ioprs

  ! --- 外部函数 ---
  integer, external :: idsymsoc
  real*8, external :: nrm2sol, nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vmnpq(:), vmnrs(:), vpqrs_core(:), vpqrs_ref(:)
  real*8, allocatable :: icore(:), icore2(:)
  real*4, allocatable :: vmnpq_s(:), vmnrs_s(:), vpqrs_core_s(:), vpqrs_ref_s(:)
  real*4, allocatable :: icore_s(:), icore2_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vmnpq_g(:), vmnrs_g(:), vpqrs_core_g(:), vpqrs_ref_g(:)
  real*8, allocatable, device :: icore_g(:), icore2_g(:)
  real*4, allocatable, device :: vmnpq_sg(:), vmnrs_sg(:), vpqrs_core_sg(:), vpqrs_ref_sg(:)
  real*4, allocatable, device :: icore_sg(:), icore2_sg(:)

  ! 1. 解析参数与初始化
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvmnpq, symvmnrs, perm_in)
  call init_parameters(nmnpq, nmnrs, npqrs, &
                       symvmnpq, symvmnrs, symvpqrs, &
                       popm, popn, popp, popq, popr, pops, &
                       def_iopmn, def_ioppq, def_ioprs, &
                       nirrep_in, nvrt_in, npop_in) 

  if (ncycle <= 0) ncycle = 1  
  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_444_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="

  ! 分配并生成全局基础数据 (保持最大尺寸分配)
  allocate(vmnpq(nmnpq), vmnrs(nmnrs), vpqrs_core(npqrs), vpqrs_ref(npqrs))
  allocate(icore(nmnpq), icore2(nmnrs))
  allocate(vmnpq_s(nmnpq), vmnrs_s(nmnrs), vpqrs_core_s(npqrs), vpqrs_ref_s(npqrs))
  allocate(icore_s(nmnpq), icore2_s(nmnrs))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vmnpq, vmnrs, nmnpq, nmnrs)
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vmnpq, vmnpq_s, nmnpq)
  call dtos(vmnrs, vmnrs_s, nmnrs)


  ! ======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  ! ======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  ! --- CPU DP: Fast Path ---
  write(6,*) '>>> TEST 1: Fast Path (3412 * 3456 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq, vmnrs, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, fact, fact2)
  
  vpqrs_core(1:npqrs_iter) = 999.0d0  ! 精准注入脏数据
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Fast Path AntiSym ---
  write(6,*) '>>> TEST 2: Fast Path AntiSym (3412 * 3456 -> 1256, IOP=0 for Inputs) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(Vmnpq, Vmnrs, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 1, popr, pops, 1, Vpqrs_ref, symvpqrs, fact, fact2) 
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popm, popn, 0, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs, popm, popn, 0, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Fast Path Smart Route ---
  write(6,*) '>>> TEST 3: Fast Path Smart Route (1234 * 3456 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref, npqrs_iter)
  call VpqmnVmnrs_to_Vpqrs(Vmnpq, Vmnrs, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, Vpqrs_ref, symvpqrs, fact, fact2) 
  
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popp, popq, 1, popm, popn, 1, symvmnpq, '1234', &
                       vmnrs, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Slow Path Fallback ---
  write(6,*) '>>> TEST 4: Slow Path Fallback (3142 * 3546 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call sstgenso1(vmnpq, icore, popm, popp, popn, popq, symvmnpq, '1324')
  call sstgenso1(vmnrs, icore2, popm, popr, popn, pops, symvmnrs, '1324')
  call zeroso(vpqrs_ref, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore, icore2, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, fact, fact2)

  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popm, popp, 1, popn, popq, 1, symvmnpq, '3142', &
                       vmnrs, popm, popr, 1, popn, pops, 1, symvmnrs, '3546', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Ultimate Scramble ---
  write(6,*) '>>> TEST 5: Ultimate Scramble (4132 * 5364 -> 2165, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popq, popp, 1, pops, popr, 1) ! 目标维度 qpsr
  call sstgenso1(vmnpq, icore, popn, popp, popm, popq, symvmnpq, '3142')
  call sstgenso1(vmnrs, icore2, popr, popm, pops, popn, symvmnrs, '2431')
  call zeroso(vpqrs_ref, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore, icore2, symvmnpq, symvmnrs, popm, popn, 1, &
                           popq, popp, 1, pops, popr, 1, vpqrs_ref, symvpqrs, fact, fact2)

  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popn, popp, 1, popm, popq, 1, symvmnpq, '4132', &
                       vmnrs, popr, popm, 1, pops, popn, 1, symvmnrs, '5364', &
                       vpqrs_core, popq, popp, 1, pops, popr, 1, symvpqrs, '2165', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Ultimate Scramble 2 (With Reference Transpose) ---
  write(6,*) '>>> TEST 6: Ultimate Scramble 2 (pmqn * rmsn -> prqs, IOP=1) <<<'
  call sstgenso1(vmnpq, icore, popp, popm, popq, popn, symvmnpq, '2413')
  call sstgenso1(vmnrs, icore2, popr, popm, pops, popn, symvmnrs, '2413')
  
  ! 步骤 3: 借用 vpqrs_core 作为临时数组，计算 pqrs
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1) 
  call zeroso(vpqrs_core, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore, icore2, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_core, symvpqrs, fact, fact2)
  
  ! 步骤 4: 转置为目标的 prqs
  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, popq, pops, 1) 
  call zeroso(vpqrs_ref, npqrs_iter)
  call sstgenso1(vpqrs_core, vpqrs_ref, popp, popq, popr, pops, symvpqrs, '1324')

  ! 步骤 5: 验证引擎
  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popp, popm, 1, popq, popn, 1, symvmnpq, 'pmqn', &
                       vmnrs, popr, popm, 1, pops, popn, 1, symvmnrs, 'rmsn', &
                       vpqrs_core, popp, popr, 1, popq, pops, 1, symvpqrs, 'prqs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Fast Path IOP=0 ---
  write(6,*) '>>> TEST 7: Fast Path (1234 * 5634 -> 1256, IOP=0) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 0)
  call zeroso(vpqrs_ref, npqrs_iter)
  call VpqmnVrsmn_to_Vpqrs(vmnpq, vmnrs, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 0, popr, pops, 0, vpqrs_ref, symvpqrs, fact, fact2)

  vpqrs_core(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq, popp, popq, 0, popm, popn, 0, symvmnpq, '1234', &
                       vmnrs, popr, pops, 0, popm, popn, 0, symvmnrs, '5634', &
                       vpqrs_core, popp, popq, 0, popr, pops, 0, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))

  ! --- CPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 8: Accumulation Test (fact2 = 0.5) [DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref(1:npqrs_iter) = 888.0d0
  call VmnpqVmnrs_to_Vpqrs(vmnpq, vmnrs, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref, symvpqrs, 1.0d0, 0.5d0)
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core(1:npqrs_iter) = 888.0d0  ! 循环内重置累加基础
     call vv2v_dpd_444(vmnpq, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2sol(npqrs_iter, vpqrs_ref, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2sol(npqrs_iter, vpqrs_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref, vpqrs_core, int(npqrs_iter, 8))


  ! ======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  ! ======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  ! --- CPU SP: Fast Path ---
  write(6,*) '>>> TEST 1: Fast Path (3412 * 3456 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_s, vmnrs_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_s, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Fast Path AntiSym ---
  write(6,*) '>>> TEST 2: Fast Path AntiSym (3412 * 3456 -> 1256, IOP=0 for Inputs) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_s, vmnrs_s, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popm, popn, 0, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_s, popm, popn, 0, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Fast Path Smart Route ---
  write(6,*) '>>> TEST 3: Fast Path Smart Route (1234 * 3456 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VpqmnVmnrs_to_Vpqrs(Vmnpq_s, Vmnrs_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, Vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popp, popq, 1, popm, popn, 1, symvmnpq, '1234', &
                       vmnrs_s, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Slow Path Fallback ---
  write(6,*) '>>> TEST 4: Slow Path Fallback (3142 * 3546 -> 1256, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call sstgenso1(vmnpq_s, icore_s, popm, popp, popn, popq, symvmnpq, '1324')
  call sstgenso1(vmnrs_s, icore2_s, popm, popr, popn, pops, symvmnrs, '1324')
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_s, icore2_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(Vmnpq_s, popm, popp, 1, popn, popq, 1, symvmnpq, '3142', &
                       Vmnrs_s, popm, popr, 1, popn, pops, 1, symvmnrs, '3546', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Ultimate Scramble ---
  write(6,*) '>>> TEST 5: Ultimate Scramble (4132 * 5364 -> 2165, IOP=1) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popq, popp, 1, pops, popr, 1)
  call sstgenso1(vmnpq_s, icore_s, popn, popp, popm, popq, symvmnpq, '3142')
  call sstgenso1(vmnrs_s, icore2_s, popr, popm, pops, popn, symvmnrs, '2431')
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_s, icore2_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popq, popp, 1, pops, popr, 1, vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popn, popp, 1, popm, popq, 1, symvmnpq, '4132', &
                       vmnrs_s, popr, popm, 1, pops, popn, 1, symvmnrs, '5364', &
                       vpqrs_core_s, popq, popp, 1, pops, popr, 1, symvpqrs, '2165', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Ultimate Scramble 2 ---
  write(6,*) '>>> TEST 6: Ultimate Scramble 2 (pmqn * rmsn -> prqs, IOP=1) [SP] <<<'
  call sstgenso1(vmnpq_s, icore_s, popp, popm, popq, popn, symvmnpq, '2413')
  call sstgenso1(vmnrs_s, icore2_s, popr, popm, pops, popn, symvmnrs, '2413')
  
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_core_s, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_s, icore2_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_core_s, symvpqrs, fact, fact2)
  
  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, popq, pops, 1)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call sstgenso1(vpqrs_core_s, vpqrs_ref_s, popp, popq, popr, pops, symvpqrs, '1324')

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popp, popm, 1, popq, popn, 1, symvmnpq, 'pmqn', &
                       vmnrs_s, popr, popm, 1, pops, popn, 1, symvmnrs, 'rmsn', &
                       vpqrs_core_s, popp, popr, 1, popq, pops, 1, symvpqrs, 'prqs', fact, fact2)
  end do
  call walltimeso(t_e)

  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Fast Path IOP=0 ---
  write(6,*) '>>> TEST 7: Fast Path (1234 * 5634 -> 1256, IOP=0) <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 0)
  call zeroso(vpqrs_ref_s, npqrs_iter)
  call VpqmnVrsmn_to_Vpqrs(vmnpq_s, vmnrs_s, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 0, popr, pops, 0, vpqrs_ref_s, symvpqrs, fact, fact2)

  vpqrs_core_s(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_s, popp, popq, 0, popm, popn, 0, symvmnpq, '1234', &
                       vmnrs_s, popr, pops, 0, popm, popn, 0, symvmnrs, '5634', &
                       vpqrs_core_s, popp, popq, 0, popr, pops, 0, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))

  ! --- CPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 8: Accumulation Test (fact2 = 0.5) [SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_s(1:npqrs_iter) = 888.0e0
  call VmnpqVmnrs_to_Vpqrs(vmnpq_s, vmnrs_s, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_s, symvpqrs, 1.0d0, 0.5d0)
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_s(1:npqrs_iter) = 888.0e0 
     call vv2v_dpd_444(vmnpq_s, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_s, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_s, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_s, vpqrs_core_s, int(npqrs_iter, 8))


  ! ======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  ! ======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vmnpq_g(nmnpq), vmnrs_g(nmnrs), vpqrs_core_g(npqrs), vpqrs_ref_g(npqrs))
  allocate(icore_g(nmnpq), icore2_g(nmnrs))
  
  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vmnpq, vmnpq_g,nmnpq )
  call myicopy_g(vmnrs, vmnrs_g,nmnrs )

  ! --- GPU DP: Fast Path ---
  write(6,*) '>>> TEST 1: Fast Path (3412 * 3456 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_g, vmnrs_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2)

  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_g, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Fast Path AntiSym ---
  write(6,*) '>>> TEST 2: Fast Path AntiSym (3412 * 3456 -> 1256, IOP=0 for Inputs) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_g, vmnrs_g, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2) 
  
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popm, popn, 0, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_g, popm, popn, 0, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Fast Path Smart Route ---
  write(6,*) '>>> TEST 3: Fast Path Smart Route (1234 * 3456 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VpqmnVmnrs_to_Vpqrs(vmnpq_g, vmnrs_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2) 
  
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popp, popq, 1, popm, popn, 1, symvmnpq, '1234', &
                       vmnrs_g, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2) 
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Slow Path Fallback ---
  write(6,*) '>>> TEST 4: Slow Path Fallback (3142 * 3546 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call sstgenso1(vmnpq_g, icore_g, popm, popp, popn, popq, symvmnpq, '1324')
  call sstgenso1(vmnrs_g, icore2_g, popm, popr, popn, pops, symvmnrs, '1324')
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_g, icore2_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, fact, fact2)

  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popm, popp, 1, popn, popq, 1, symvmnpq, '3142', &
                       vmnrs_g, popm, popr, 1, popn, pops, 1, symvmnrs, '3546', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Ultimate Scramble ---
  write(6,*) '>>> TEST 5: Ultimate Scramble (4132 * 5364 -> 2165, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popq, popp, 1, pops, popr, 1)
  call sstgenso1(vmnpq_g, icore_g, popn, popp, popm, popq, symvmnpq, '3142')
  call sstgenso1(vmnrs_g, icore2_g, popr, popm, pops, popn, symvmnrs, '2431')
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_g, icore2_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popq, popp, 1, pops, popr, 1, vpqrs_ref_g, symvpqrs, fact, fact2)
 
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popn, popp, 1, popm, popq, 1, symvmnpq, '4132', &
                       vmnrs_g, popr, popm, 1, pops, popn, 1, symvmnrs, '5364', &
                       vpqrs_core_g, popq, popp, 1, pops, popr, 1, symvpqrs, '2165', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Ultimate Scramble 2 ---
  write(6,*) '>>> TEST 6: Ultimate Scramble 2 (pmqn * rmsn -> prqs, IOP=1) [GPU] <<<'
  call sstgenso1(vmnpq_g, icore_g, popp, popm, popq, popn, symvmnpq, '2413')
  call sstgenso1(vmnrs_g, icore2_g, popr, popm, pops, popn, symvmnrs, '2413')
  
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_core_g, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_g, icore2_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_core_g, symvpqrs, fact, fact2)
  
  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, popq, pops, 1)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call sstgenso1(vpqrs_core_g, vpqrs_ref_g, popp, popq, popr, pops, symvpqrs, '1324')
  
  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popp, popm, 1, popq, popn, 1, symvmnpq, 'pmqn', &
                       vmnrs_g, popr, popm, 1, pops, popn, 1, symvmnrs, 'rmsn', &
                       vpqrs_core_g, popp, popr, 1, popq, pops, 1, symvpqrs, 'prqs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Fast Path IOP=0 ---
  write(6,*) '>>> TEST 7: Fast Path (1234 * 5634 -> 1256, IOP=0) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 0)
  call zeroso(vpqrs_ref_g, npqrs_iter)
  call VpqmnVrsmn_to_Vpqrs(vmnpq_g, vmnrs_g, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 0, popr, pops, 0, vpqrs_ref_g, symvpqrs, fact, fact2)

  vpqrs_core_g(1:npqrs_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_g, popp, popq, 0, popm, popn, 0, symvmnpq, '1234', &
                       vmnrs_g, popr, pops, 0, popm, popn, 0, symvmnrs, '5634', &
                       vpqrs_core_g, popp, popq, 0, popr, pops, 0, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))

  ! --- GPU DP: Accumulation Test ---
  write(6,*) '>>> TEST 8: Accumulation Test (fact2 = 0.5) [GPU DP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_g(1:npqrs_iter) = 888.0d0
  call VmnpqVmnrs_to_Vpqrs(vmnpq_g, vmnrs_g, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_g, symvpqrs, 1.0d0, 0.5d0)
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_g(1:npqrs_iter) = 888.0d0
     call vv2v_dpd_444(vmnpq_g, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_g, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_g, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_g, vpqrs_core_g, int(npqrs_iter, 8))


  ! ======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  ! ======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vmnpq_sg(nmnpq), vmnrs_sg(nmnrs), vpqrs_core_sg(npqrs), vpqrs_ref_sg(npqrs))
  allocate(icore_sg(nmnpq), icore2_sg(nmnrs))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vmnpq_s,vmnpq_sg,nmnpq)
  call myicopy_g(vmnrs_s,vmnrs_sg,nmnrs)

  ! --- GPU SP: Fast Path ---
  write(6,*) '>>> TEST 1: Fast Path (3412 * 3456 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_sg, vmnrs_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_sg, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Fast Path AntiSym ---
  write(6,*) '>>> TEST 2: Fast Path AntiSym (3412 * 3456 -> 1256, IOP=0 for Inputs) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(vmnpq_sg, vmnrs_sg, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popm, popn, 0, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_sg, popm, popn, 0, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Fast Path Smart Route ---
  write(6,*) '>>> TEST 3: Fast Path Smart Route (1234 * 3456 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VpqmnVmnrs_to_Vpqrs(vmnpq_sg, vmnrs_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popp, popq, 1, popm, popn, 1, symvmnpq, '1234', &
                       vmnrs_sg, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Slow Path Fallback ---
  write(6,*) '>>> TEST 4: Slow Path Fallback (3142 * 3546 -> 1256, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call sstgenso1(vmnpq_sg, icore_sg, popm, popp, popn, popq, symvmnpq, '1324')
  call sstgenso1(vmnrs_sg, icore2_sg, popm, popr, popn, pops, symvmnrs, '1324')
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_sg, icore2_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popm, popp, 1, popn, popq, 1, symvmnpq, '3142', &
                       vmnrs_sg, popm, popr, 1, popn, pops, 1, symvmnrs, '3546', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Ultimate Scramble ---
  write(6,*) '>>> TEST 5: Ultimate Scramble (4132 * 5364 -> 2165, IOP=1) [GPU] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popq, popp, 1, pops, popr, 1)
  call sstgenso1(vmnpq_sg, icore_sg, popn, popp, popm, popq, symvmnpq, '3142')
  call sstgenso1(vmnrs_sg, icore2_sg, popr, popm, pops, popn, symvmnrs, '2431')
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_sg, icore2_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popq, popp, 1, pops, popr, 1, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popn, popp, 1, popm, popq, 1, symvmnpq, '4132', &
                       vmnrs_sg, popr, popm, 1, pops, popn, 1, symvmnrs, '5364', &
                       vpqrs_core_sg, popq, popp, 1, pops, popr, 1, symvpqrs, '2165', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Ultimate Scramble 2 ---
  write(6,*) '>>> TEST 6: Ultimate Scramble 2 (pmqn * rmsn -> prqs, IOP=1) [GPU SP] <<<'
  call sstgenso1(vmnpq_sg, icore_sg, popp, popm, popq, popn, symvmnpq, '2413')
  call sstgenso1(vmnrs_sg, icore2_sg, popr, popm, pops, popn, symvmnrs, '2413')
  
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  call zeroso(vpqrs_core_sg, npqrs_iter)
  call VmnpqVmnrs_to_Vpqrs(icore_sg, icore2_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_core_sg, symvpqrs, fact, fact2)

  npqrs_iter = idsymsoc(symvpqrs, popp, popr, 1, popq, pops, 1)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call sstgenso1(vpqrs_core_sg, vpqrs_ref_sg, popp, popq, popr, pops, symvpqrs, '1324')

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popp, popm, 1, popq, popn, 1, symvmnpq, 'pmqn', &
                       vmnrs_sg, popr, popm, 1, pops, popn, 1, symvmnrs, 'rmsn', &
                       vpqrs_core_sg, popp, popr, 1, popq, pops, 1, symvpqrs, 'prqs', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Fast Path IOP=0 ---
  write(6,*) '>>> TEST 7: Fast Path (1234 * 5634 -> 1256, IOP=0) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 0, popr, pops, 0)
  call zeroso(vpqrs_ref_sg, npqrs_iter)
  call VpqmnVrsmn_to_Vpqrs(vmnpq_sg, vmnrs_sg, symvmnpq, symvmnrs, popm, popn, 0, &
                           popp, popq, 0, popr, pops, 0, vpqrs_ref_sg, symvpqrs, fact, fact2)

  vpqrs_core_sg(1:npqrs_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_444(vmnpq_sg, popp, popq, 0, popm, popn, 0, symvmnpq, '1234', &
                       vmnrs_sg, popr, pops, 0, popm, popn, 0, symvmnrs, '5634', &
                       vpqrs_core_sg, popp, popq, 0, popr, pops, 0, symvpqrs, '1256', fact, fact2)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

  ! --- GPU SP: Accumulation Test ---
  write(6,*) '>>> TEST 8: Accumulation Test (fact2 = 0.5) [GPU SP] <<<'
  npqrs_iter = idsymsoc(symvpqrs, popp, popq, 1, popr, pops, 1)
  vpqrs_ref_sg(1:npqrs_iter) = 888.0e0
  call VmnpqVmnrs_to_Vpqrs(vmnpq_sg, vmnrs_sg, symvmnpq, symvmnrs, popm, popn, 1, &
                           popp, popq, 1, popr, pops, 1, vpqrs_ref_sg, symvpqrs, 1.0d0, 0.5d0)
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqrs_core_sg(1:npqrs_iter) = 888.0e0
     call vv2v_dpd_444(vmnpq_sg, popm, popn, 1, popp, popq, 1, symvmnpq, '3412', &
                       vmnrs_sg, popm, popn, 1, popr, pops, 1, symvmnrs, '3456', &
                       vpqrs_core_sg, popp, popq, 1, popr, pops, 1, symvpqrs, '1256', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqrs_iter, vpqrs_ref_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqrs_iter, vpqrs_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqrs_ref_sg, vpqrs_core_sg, int(npqrs_iter, 8))

end subroutine vv2v_dpd_444_check
