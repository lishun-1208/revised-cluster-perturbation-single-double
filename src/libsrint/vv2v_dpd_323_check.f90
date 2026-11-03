!=======================================================================
!  subroutine: vv2v_dpd_323_check
!
!  功能说明:
!    针对 3指标 * 2指标 -> 3指标 (例如 Vpqm * Vmr -> Vpqr) 的统一基准测试套件。
!    执行顺序: CPU-DP -> CPU-SP -> GPU-DP -> GPU-SP
!    【终极防污染修复】: 彻底解决了 v1_ref 全局变量被跨测试污染的问题。
!                  在 GPU 拷贝之前，当场重新进行压缩态内存重排打包，
!                  确保传输到显存中的数据绝对纯净！
!=======================================================================
subroutine vv2v_dpd_323_check() 
  use mod_iop
  use mod_size
  use mod_orbit
  use mod_ioff 
  use mod_symm
  implicit none

  ! --- 体系与对称性变量 ---
  integer :: npqm, nmr, npqr
  integer :: npqr_iter, npqm_iter  
  integer :: symvpqm, symvmr, symvpqr
  integer :: popp(8), popq(8), popr(8), popm(8), popdum(8)
  
  integer :: nirrep_in, nvrt_in, npop_in, ncycle, iter
  character(len=4) :: perm_in
  real*8 :: tt0, tt1, t_s, t_e
  real*8 :: fact, fact2
  integer :: def_ioppq, def_ioprm, def_ioppr

  ! --- 内存重排辅助变量 ---
  integer :: idx, num_m, num_pq, m_loc, pq_loc, i

  ! --- 外部函数 ---
  integer, external :: idsymsoc, irpdso
  real*8, external :: nrm2so

  ! --- CPU 数组声明 ---
  real*8, allocatable :: vpqm(:), vmr(:), vpqr_core(:), vpqr_ref(:)
  real*8, allocatable :: icore1(:), icore2(:), v1_ref(:)
  real*4, allocatable :: vpqm_s(:), vmr_s(:), vpqr_core_s(:), vpqr_ref_s(:)
  real*4, allocatable :: icore1_s(:), icore2_s(:), v1_ref_s(:)

  ! --- GPU 数组声明 ---
  real*8, allocatable, device :: vpqm_g(:), vmr_g(:), vpqr_core_g(:), vpqr_ref_g(:)
  real*8, allocatable, device :: icore1_g(:), icore2_g(:), v1_ref_g(:)
  real*4, allocatable, device :: vpqm_sg(:), vmr_sg(:), vpqr_core_sg(:), vpqr_ref_sg(:)
  real*4, allocatable, device :: icore1_sg(:), icore2_sg(:), v1_ref_sg(:)

  ! 1. 解析参数与初始化 
  call parse_cmd_args(nirrep_in, nvrt_in, npop_in, ncycle, symvpqm, symvmr, perm_in)
  
  call init_parameters(npqm, nmr, npqr, &
                       symvpqm, symvmr, symvpqr, &
                       popp, popq, popr, popm, popdum, popdum, &
                       def_ioppq, def_ioprm, def_ioppr, &
                       nirrep_in, nvrt_in, npop_in) 

  ! 设置虚拟维度 Dummy Pop (大小恒为 1，用于 323 模拟 444)
  popdum = 0; popdum(1) = 1

  npqm = idsymsoc(symvpqm, popp, popq, 1, popm, popdum, 1)
  nmr  = irpdso(symvmr, popm, popr, 1)
  npqr = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)

  if (ncycle <= 0) ncycle = 1  
  fact = 1.0d0
  fact2 = 0.0d0

  write(6,*) "=========================================================="
  write(6,*) "       [vv2v_dpd_323_check] Dynamic Correctness Benchmark "
  write(6,*) "=========================================================="
  write(6,"(A, I10, A, I10, A, I10)") "  Matrix Size: Vpqm =", npqm, "  Vmr =", nmr, "  Vpqr =", npqr
  write(6,*) "=========================================================="

  allocate(vpqm(npqm), vmr(nmr), vpqr_core(npqr), vpqr_ref(npqr))
  allocate(icore1(npqm), icore2(nmr), v1_ref(npqm))
  
  allocate(vpqm_s(npqm), vmr_s(nmr), vpqr_core_s(npqr), vpqr_ref_s(npqr))
  allocate(icore1_s(npqm), icore2_s(nmr), v1_ref_s(npqm))

  igpu = 0; isd = 2; prec = 8
  call init_integrals_dp(vpqm, vmr, npqm, nmr)
  
  igpu = 0; isd = 1; prec = 4
  call dtos(vpqm, vpqm_s, npqm)
  call dtos(vmr, vmr_s, nmr)

  !======================================================================
  ! [ 阶段 1：CPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 1: CPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 2; prec = 8

  ! --- TEST 1 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqm * mr -> pqr, IOP=1) [DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  call sstgenso1(vpqm, v1_ref, popp, popq, popm, popdum, symvpqm, '3412')
  call zeroso(vpqr_ref, npqr_iter)
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 2 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqm * rm -> pqr, IOP=1) [DP] <<<'
  call zeroso(vpqr_ref, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vmr, icore2, popm, popr, popdum, popdum, symvmr, '2134')
  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm, popp, popq, 1, popm, symvpqm, 'pqm', &
                       icore2, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 3 ---
  write(6,*) '>>> TEST 3: Extreme Scramble (qpm * rm -> pqr, IOP=1) [DP] <<<'
  call zeroso(vpqr_ref, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqm, icore1, popp, popq, popm, popdum, symvpqm, '2134')
  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(icore1, popq, popp, 1, popm, symvpqm, 'qpm', &
                       icore2, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 4 ---
  write(6,*) '>>> TEST 4: Accumulation Test (pqm * mr -> pqr, fact2 = 0.5) [DP] <<<'
  call zeroso(vpqr_ref, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  vpqr_ref(1:npqr_iter) = 444.0d0 + vpqr_ref(1:npqr_iter)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqr_core(1:npqr_iter) = 888.0d0  
     call vv2v_dpd_323(vpqm, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 5 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqm * mr -> pqr, IOP_pq=0) [DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 0, popr, popdum, 1)
  
  ! 【当场打包】：绝对避免被上一个测试污染
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 0)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref(idx + (pq_loc-1)*num_m + m_loc) = vpqm(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do

  call zeroso(vpqr_ref, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 0, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm, popp, popq, 0, popm, symvpqm, 'pqm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core, popp, popq, 0, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 6 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqm * mr -> pqr, IOP_pq=2) [DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 2, popr, popdum, 1)
  
  ! 【当场打包】：绝对避免被污染
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 2)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref(idx + (pq_loc-1)*num_m + m_loc) = vpqm(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do

  call zeroso(vpqr_ref, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 2, popr, popdum, 1, vpqr_ref, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm, popp, popq, 2, popm, symvpqm, 'pqm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core, popp, popq, 2, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))

  ! --- TEST 7 ---
  write(6,*) '>>> TEST 7: Output Scramble (pqm * mr -> rpq, IOP=1) [DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  
  ! 恢复 v1_ref 的全尺寸转置供 TEST 7 使用
  call sstgenso1(vpqm, v1_ref, popp, popq, popm, popdum, symvpqm, '3412')
  
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref, vmr, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_core, symvpqr, 1.0d0, 0.0d0)
  call zeroso(vpqr_ref, npqr_iter)
  call sstgenso1(vpqr_core, vpqr_ref, popp, popq, popr, popdum, symvpqr, '3124') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core, popr, popp, 1, popq, symvpqr, 'rpq', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref, vpqr_core, int(npqr_iter, 8))


  !======================================================================
  ! [ 阶段 2：CPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 2: CPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  igpu = 0; isd = 1; prec = 4

  call sstgenso1(vpqm_s, v1_ref_s, popp, popq, popm, popdum, symvpqm, '3412')

  ! --- TEST 1 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqm * mr -> pqr, IOP=1) [SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_s, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_s, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 2 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqm * rm -> pqr, IOP=1) [SP] <<<'
  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vmr_s, icore2_s, popm, popr, popdum, popdum, symvmr, '2134')
  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_s, popp, popq, 1, popm, symvpqm, 'pqm', &
                       icore2_s, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_s, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 3 ---
  write(6,*) '>>> TEST 3: Extreme Scramble (qpm * rm -> pqr, IOP=1) [SP] <<<'
  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqm_s, icore1_s, popp, popq, popm, popdum, symvpqm, '2134')
  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(icore1_s, popq, popp, 1, popm, symvpqm, 'qpm', &
                       icore2_s, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_s, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 4 ---
  write(6,*) '>>> TEST 4: Accumulation Test (pqm * mr -> pqr, fact2 = 0.5) [SP] <<<'
  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  vpqr_ref_s(1:npqr_iter) = 444.0e0 + vpqr_ref_s(1:npqr_iter)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqr_core_s(1:npqr_iter) = 888.0e0  
     call vv2v_dpd_323(vpqm_s, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_s, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 5 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqm * mr -> pqr, IOP_pq=0) [SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 0, popr, popdum, 1)
  
  ! 【当场打包】：SP 压缩态即时重排
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 0)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref_s(idx + (pq_loc-1)*num_m + m_loc) = vpqm_s(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do

  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 0, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_s, popp, popq, 0, popm, symvpqm, 'pqm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_s, popp, popq, 0, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 6 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqm * mr -> pqr, IOP_pq=2) [SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 2, popr, popdum, 1)
  
  ! 【当场打包】：SP 压缩态即时重排
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 2)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref_s(idx + (pq_loc-1)*num_m + m_loc) = vpqm_s(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do

  call zeroso(vpqr_ref_s, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 2, popr, popdum, 1, vpqr_ref_s, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_s, popp, popq, 2, popm, symvpqm, 'pqm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_s, popp, popq, 2, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))

  ! --- TEST 7 ---
  write(6,*) '>>> TEST 7: Output Scramble (pqm * mr -> rpq, IOP=1) [SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  
  ! 恢复 v1_ref_s 的全尺寸转置供 TEST 7 使用
  call sstgenso1(vpqm_s, v1_ref_s, popp, popq, popm, popdum, symvpqm, '3412')

  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_s, vmr_s, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_core_s, symvpqr, 1.0d0, 0.0d0)
  call zeroso(vpqr_ref_s, npqr_iter)
  call sstgenso1(vpqr_core_s, vpqr_ref_s, popp, popq, popr, popdum, symvpqr, '3124') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_s, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_s(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_s, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_s, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_s, popr, popp, 1, popq, symvpqr, 'rpq', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_s, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_s, vpqr_core_s, int(npqr_iter, 8))


  !======================================================================
  ! [ 阶段 3：GPU 双精度 (DP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 3: GPU DOUBLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqm_g(npqm), vmr_g(nmr), vpqr_core_g(npqr), vpqr_ref_g(npqr))
  allocate(icore1_g(npqm), icore2_g(nmr), v1_ref_g(npqm))
  
  igpu = 1; isd = 2; prec = 8
  call myicopy_g(vpqm, vpqm_g, npqm)
  call myicopy_g(vmr, vmr_g, nmr)
  call sstgenso1(vpqm_g, v1_ref_g, popp, popq, popm, popdum, symvpqm, '3412')

  ! --- TEST 1 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqm * mr -> pqr, IOP=1) [GPU DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_g, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_g, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 2 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqm * rm -> pqr, IOP=1) [GPU DP] <<<'
  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vmr_g, icore2_g, popm, popr, popdum, popdum, symvmr, '2134')
  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_g, popp, popq, 1, popm, symvpqm, 'pqm', &
                       icore2_g, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_g, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 3 ---
  write(6,*) '>>> TEST 3: Extreme Scramble (qpm * rm -> pqr, IOP=1) [GPU DP] <<<'
  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqm_g, icore1_g, popp, popq, popm, popdum, symvpqm, '2134')
  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(icore1_g, popq, popp, 1, popm, symvpqm, 'qpm', &
                       icore2_g, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_g, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 4 ---
  write(6,*) '>>> TEST 4: Accumulation Test (pqm * mr -> pqr, fact2 = 0.5) [GPU DP] <<<'
  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
!$acc kernels present(vpqr_ref_g)
  vpqr_ref_g(1:npqr_iter) = 444.0d0 + vpqr_ref_g(1:npqr_iter)
!$acc end kernels
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqr_core_g(1:npqr_iter) = 888.0d0  
     call vv2v_dpd_323(vpqm_g, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_g, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 5 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqm * mr -> pqr, IOP_pq=0) [GPU DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 0, popr, popdum, 1)
  
  ! 【当场打包】：重新生成 IOP=0 对齐的张量并传入 GPU，彻底消灭残留脏数据！
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 0)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref(idx + (pq_loc-1)*num_m + m_loc) = vpqm(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do
  npqm_iter = idsymsoc(symvpqm, popp, popq, 0, popm, popdum, 1)
  call myicopy_g(v1_ref, v1_ref_g, npqm_iter)

  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 0, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_g, popp, popq, 0, popm, symvpqm, 'pqm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_g, popp, popq, 0, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 6 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqm * mr -> pqr, IOP_pq=2) [GPU DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 2, popr, popdum, 1)
  
  ! 【当场打包】：重新生成 IOP=2 对齐的张量并传入 GPU
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 2)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref(idx + (pq_loc-1)*num_m + m_loc) = vpqm(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do
  npqm_iter = idsymsoc(symvpqm, popp, popq, 2, popm, popdum, 1)
  call myicopy_g(v1_ref, v1_ref_g, npqm_iter)

  call zeroso(vpqr_ref_g, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 2, popr, popdum, 1, vpqr_ref_g, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_g, popp, popq, 2, popm, symvpqm, 'pqm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_g, popp, popq, 2, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))

  ! --- TEST 7 ---
  write(6,*) '>>> TEST 7: Output Scramble (pqm * mr -> rpq, IOP=1) [GPU DP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  
  ! 恢复全尺寸转置
  call sstgenso1(vpqm_g, v1_ref_g, popp, popq, popm, popdum, symvpqm, '3412')

  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_g, vmr_g, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_core_g, symvpqr, 1.0d0, 0.0d0)
  call zeroso(vpqr_ref_g, npqr_iter)
  call sstgenso1(vpqr_core_g, vpqr_ref_g, popp, popq, popr, popdum, symvpqr, '3124') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_g, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_g(1:npqr_iter) = 999.0d0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_g, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_g, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_g, popr, popp, 1, popq, symvpqr, 'rpq', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_g, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_g, vpqr_core_g, int(npqr_iter, 8))


  !======================================================================
  ! [ 阶段 4：GPU 单精度 (SP) ]
  !======================================================================
  write(6,*) ''
  write(6,*) '----------------------------------------------------------'
  write(6,*) '    [ PHASE 4: GPU SINGLE PRECISION ]'
  write(6,*) '----------------------------------------------------------'
  
  allocate(vpqm_sg(npqm), vmr_sg(nmr), vpqr_core_sg(npqr), vpqr_ref_sg(npqr))
  allocate(icore1_sg(npqm), icore2_sg(nmr), v1_ref_sg(npqm))
  
  igpu = 1; isd = 1; prec = 4
  call myicopy_g(vpqm_s, vpqm_sg, npqm)
  call myicopy_g(vmr_s, vmr_sg, nmr)
  call sstgenso1(vpqm_sg, v1_ref_sg, popp, popq, popm, popdum, symvpqm, '3412')

  ! --- TEST 1 ---
  write(6,*) '>>> TEST 1: Standard Mode (pqm * mr -> pqr, IOP=1) [GPU SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_sg, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_sg, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 2 ---
  write(6,*) '>>> TEST 2: Swapped 2-Index (pqm * rm -> pqr, IOP=1) [GPU SP] <<<'
  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vmr_sg, icore2_sg, popm, popr, popdum, popdum, symvmr, '2134')
  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_sg, popp, popq, 1, popm, symvpqm, 'pqm', &
                       icore2_sg, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_sg, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 3 ---
  write(6,*) '>>> TEST 3: Extreme Scramble (qpm * rm -> pqr, IOP=1) [GPU SP] <<<'
  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  call sstgenso1(vpqm_sg, icore1_sg, popp, popq, popm, popdum, symvpqm, '2134')
  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(icore1_sg, popq, popp, 1, popm, symvpqm, 'qpm', &
                       icore2_sg, popr, popm, 1, symvmr, 'rm', &
                       vpqr_core_sg, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 4 ---
  write(6,*) '>>> TEST 4: Accumulation Test (pqm * mr -> pqr, fact2 = 0.5) [GPU SP] <<<'
  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
!$acc kernels present(vpqr_ref_sg)
  vpqr_ref_sg(1:npqr_iter) = 444.0e0 + vpqr_ref_sg(1:npqr_iter)
!$acc end kernels
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  call walltimeso(t_s)
  do iter = 1, ncycle
     vpqr_core_sg(1:npqr_iter) = 888.0e0  
     call vv2v_dpd_323(vpqm_sg, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_sg, popp, popq, 1, popr, symvpqr, 'pqr', 1.0d0, 0.5d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 5 ---
  write(6,*) '>>> TEST 5: Intact Pair IOP=0 (pqm * mr -> pqr, IOP_pq=0) [GPU SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 0, popr, popdum, 1)
  
  ! 【当场打包】：重新生成 IOP=0 对齐的张量并传入 GPU，彻底消灭残留脏数据！
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 0)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref_s(idx + (pq_loc-1)*num_m + m_loc) = vpqm_s(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do
  npqm_iter = idsymsoc(symvpqm, popp, popq, 0, popm, popdum, 1)
  call myicopy_g(v1_ref_s, v1_ref_sg, npqm_iter)

  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 0, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_sg, popp, popq, 0, popm, symvpqm, 'pqm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_sg, popp, popq, 0, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 6 ---
  write(6,*) '>>> TEST 6: Intact Pair IOP=2 (pqm * mr -> pqr, IOP_pq=2) [GPU SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 2, popr, popdum, 1)
  
  ! 【当场打包】：重新生成 IOP=2 对齐的张量并传入 GPU
  idx = 0
  do i = 1, nirrep
     num_m = popm(i)
     num_pq = irpdso(dirprd(i, symvpqm), popp, popq, 2)
     do m_loc = 1, num_m
        do pq_loc = 1, num_pq
           v1_ref_s(idx + (pq_loc-1)*num_m + m_loc) = vpqm_s(idx + (m_loc-1)*num_pq + pq_loc)
        end do
     end do
     idx = idx + num_m * num_pq
  end do
  npqm_iter = idsymsoc(symvpqm, popp, popq, 2, popm, popdum, 1)
  call myicopy_g(v1_ref_s, v1_ref_sg, npqm_iter)

  call zeroso(vpqr_ref_sg, npqr_iter)
  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 2, popr, popdum, 1, vpqr_ref_sg, symvpqr, 1.0d0, 0.0d0)
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'
  
  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_sg, popp, popq, 2, popm, symvpqm, 'pqm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_sg, popp, popq, 2, popr, symvpqr, 'pqr', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

  ! --- TEST 7 ---
  write(6,*) '>>> TEST 7: Output Scramble (pqm * mr -> rpq, IOP=1) [GPU SP] <<<'
  npqr_iter = idsymsoc(symvpqr, popp, popq, 1, popr, popdum, 1)
  
  ! 恢复全尺寸转置
  call sstgenso1(vpqm_sg, v1_ref_sg, popp, popq, popm, popdum, symvpqm, '3412')

  call walltimeso(t_s)
  call VmnpqVmnrs_to_Vpqrs(v1_ref_sg, vmr_sg, symvpqm, symvmr, popm, popdum, 1, &
                           popp, popq, 1, popr, popdum, 1, vpqr_core_sg, symvpqr, 1.0d0, 0.0d0)
  call zeroso(vpqr_ref_sg, npqr_iter)
  call sstgenso1(vpqr_core_sg, vpqr_ref_sg, popp, popq, popr, popdum, symvpqr, '3124') 
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Ref Norm:', nrm2so(npqr_iter, vpqr_ref_sg, 1), 'Ref Time:', (t_e-t_s), ' s'

  vpqr_core_sg(1:npqr_iter) = 999.0e0
  call walltimeso(t_s)
  do iter = 1, ncycle
     call vv2v_dpd_323(vpqm_sg, popp, popq, 1, popm, symvpqm, 'pqm', &
                       vmr_sg, popm, popr, 1, symvmr, 'mr', &
                       vpqr_core_sg, popr, popp, 1, popq, symvpqr, 'rpq', 1.0d0, 0.0d0)
  end do
  call walltimeso(t_e)
  write(6,"(A25, E15.8, A15, F10.6, A2)") 'Opt Norm:', nrm2so(npqr_iter, vpqr_core_sg, 1), 'Avg Time:', (t_e-t_s)/ncycle, ' s'
  call diffv(vpqr_ref_sg, vpqr_core_sg, int(npqr_iter, 8))

end subroutine vv2v_dpd_323_check
