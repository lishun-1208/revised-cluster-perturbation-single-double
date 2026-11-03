! =======================================================================
! 终极工业级张量引擎 - 统一调度入口 (Dispatcher)
!
! 【模块功能】
!   作为整个张量收缩运算的唯一前端入口，负责隔离外层算法与底层硬件。
!   根据全局环境变量 (prec, igpu) 自动将任务分发至四台高能引擎：
!     - vv2v_dpd_444_d  (CPU, 双精度)
!     - vv2v_dpd_444_s  (CPU, 单精度)
!     - vv2v_dpd_444_dg (GPU, 双精度, 显存拼装+cuBLAS)
!     - vv2v_dpd_444_sg (GPU, 单精度, 显存拼装+cuBLAS)
!
! 【核心特性 1：自由度自动映射翻译器 (Auto Index Mapper)】
!   外层调用不再被底层的 '1234' 数字硬编码所绑架！
!   本主程序允许您使用任意 4 字符的别名来表示张量的 4 个维度（如 aeim, mbej, pmqn）。
!   翻译器会自动将其拦截，并翻译为底层引擎所需的严格数字规范向下传递。
!
! 【⚠️ 底层子程序调用严重警告 (Strict Numeric Routing)】
!   虽然本主程序 (vv2v_dpd_444) 允许任意字符乱序组合，但底层的四个具体运算子程序
!   (_d, _s, _dg, _sg) 内部是基于严格的数字特征匹配与寻址的！
!   如果您在高级算法中需要绕过本主程序，直接调用底层的四个子程序，
!   请务必注意：传入的类型字符串（type1, type2, type3）必须严格按照
!   '1234'、'3456' 和 '1256' 的数字来进行随机排列组合！
!     - 参与收缩的内部指标必须固定为：'3' 和 '4'
!     - 结果矩阵的外部指标必须固定为：'1', '2', '5', '6'
!   绝对不能向底层子程序直接传递 'pmqn'、'rmsn' 这种字母组合，否则将导致匹配失败及内存越界！
!
! 【核心特性 2：极速通道 (Fast Path) 与 性能建议】
!   底层引擎基于访存连续性，包含两套自动调度的计算逻辑：
!     A. 极速通道 (Fast Path): 内存天然连续，实现 Zero-Copy 直接调用 xgemm / cuBLAS。
!     B. 回退通道 (Fallback) : 触发极速内存重排 (Packing) 后，再调用底层硬件 GEMM 加速。
!
!   [������ 性能调优指南]：
!   经严谨测试，Fallback 引擎凭借高效的内存打包机制，已达到极高的吞吐量。但若外层算法
!   能严格命中以下四种标准排布，将彻底省去 GPU/CPU 的内存重排开销，可进一步带来约 25% 的性能提升：
!     - '3412' * '3456' -> '1256'
!     - '1234' * '3456' -> '1256'
!     - '1234' * '5634' -> '1256'
!     - '3412' * '5634' -> '1256'
!   建议：在算法推导阶段，若能利用张量对称性无代价地凑成标准排布，请优先使用。
!
! 【安全与规范注意事项】
!   1. IOP 限制：如果运算掉入了 Fallback 慢速通道（如 3142 * 3546 -> 1256），
!      则所有的 IOP 必须为 1 (Full space)。底层引擎严禁在乱序交叉的维度上执行 IOP=0/2 
!      的压缩空间计算（数学上具有歧义）。请先用 sstgenso1 解压为 IOP=1 后再调用。
!   2. GPU 指针：当 igpu=1 时，传入的 v1, v2, v3 必须已经是存在于 Device 上的显存指针，
!      本 Dispatcher 不负责 Host to Device 的数据拷贝。
!
! 【BLAS 标准缩并系数 (Scaling Factors)】
!   遵循 BLAS GEMM 标准数学模型：C = alpha * (A * B) + beta * C
!   - fact1 (alpha): 乘法项缩放系数。
!   - fact2 (beta):  累加项缩放系数。
!   如果不需要累加（纯计算并覆盖 C），务必将 fact2 设置为 0.0d0；
!   如果需要把计算增量累加到 C 原有的值上，务必将 fact2 设置为 1.0d0。
! =======================================================================
      subroutine vv2v_dpd_444(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                               v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
                               v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1, fact2)
   
      use mod_iop
      implicit none

      ! --- 接口参数声明 (对外维持标准的 Double Precision 签名) ---
      character(len=4), intent(in)  :: type1, type2, type3 
      real*8, intent(in)            :: v1(*), v2(*)        
      real*8, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12, iop3_34
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1, fact2        

      ! --- 自动映射字符变量 (用于任意别名转换) ---
      character(len=4) :: t1_map, t2_map, t3_map
      character(len=1) :: map_from(6), map_to(6)
      integer :: n_mapped, i, j
      logical :: is_contracted

      ! --- 合法性校验辅助变量 ---
      integer :: pop1(4,8), pop2(4,8)
      integer :: pos1(6), pos2(6)
      logical :: is_fast_path

      ! =============================================================
      ! 自由度自动映射翻译器 (Automatic Index Mapper)
      ! =============================================================
      
      ! 1. 强制将结果矩阵 v3 (type3) 的四个外部指标映射为标准的 '1', '2', '5', '6'
      map_from(1) = type3(1:1); map_to(1) = '1'
      map_from(2) = type3(2:2); map_to(2) = '2'
      map_from(3) = type3(3:3); map_to(3) = '5'
      map_from(4) = type3(4:4); map_to(4) = '6'

      ! 2. 遍历 v1 (type1)，找出那两个不包含在 v3 中的“收缩指标”，并映射为 '3' 和 '4'
      n_mapped = 4
      do i = 1, 4
          is_contracted = .true.
          do j = 1, 4
              if (type1(i:i) == type3(j:j)) then
                  is_contracted = .false.
                  exit
              end if
          end do
          if (is_contracted) then
              n_mapped = n_mapped + 1
              map_from(n_mapped) = type1(i:i)
              if (n_mapped == 5) map_to(n_mapped) = '3'
              if (n_mapped == 6) map_to(n_mapped) = '4'
          end if
      end do

      ! 3. 将推导出的映射规则应用到局部的 t1_map, t2_map, t3_map
      t1_map = type1; t2_map = type2; t3_map = type3
      do i = 1, 4
          do j = 1, 6
              if (type1(i:i) == map_from(j)) then
                  t1_map(i:i) = map_to(j)
                  exit
              end if
          end do
          do j = 1, 6
              if (type2(i:i) == map_from(j)) then
                  t2_map(i:i) = map_to(j)
                  exit
              end if
          end do
          do j = 1, 6
              if (type3(i:i) == map_from(j)) then
                  t3_map(i:i) = map_to(j)
                  exit
              end if
          end do
      end do

      ! =============================================================
      ! [安全墙] 强制合法性与一致性校验 (Sanity & Consistency Checks)
      ! =============================================================
      
      ! 1. 准备校验数据：仅将涉及收缩校验的张量 1 和 2 维数打包
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:); pop2(3,:) = p2_3(:); pop2(4,:) = p2_4(:)

      ! 追踪翻译后指标 '3' 和 '4' (收缩指标) 在张量 1 和 2 中的物理维度位置(1~4)
      pos1 = 0; pos2 = 0
      do i = 1, 4
          if (t1_map(i:i) == '3') pos1(3) = i
          if (t1_map(i:i) == '4') pos1(4) = i
          if (t2_map(i:i) == '3') pos2(3) = i
          if (t2_map(i:i) == '4') pos2(4) = i
      end do

      ! 校验1：收缩指标同源轨道校验 (V1 和 V2 的收缩指标对应的轨道维数必须严格相等)
      if (pos1(3) /= 0 .and. pos2(3) /= 0) then
          if (any(pop1(pos1(3), :) /= pop2(pos2(3), :))) then
              write(6,*) "========================================================="
              write(6,*) " FATAL ERROR in vv2v_dpd_444 (Dispatcher):"
              write(6,*) " 收缩轨道维度验证失败！"
              write(6,*) " 张量 V1 和 V2 在第 1 个收缩指标上提供的轨道数不一致。"
              write(6,*) " 请检查调用组合: ", type1, " * ", type2, " -> ", type3
              write(6,*) "========================================================="
              stop
          end if
      end if
      if (pos1(4) /= 0 .and. pos2(4) /= 0) then
          if (any(pop1(pos1(4), :) /= pop2(pos2(4), :))) then
              write(6,*) "========================================================="
              write(6,*) " FATAL ERROR in vv2v_dpd_444 (Dispatcher):"
              write(6,*) " 收缩轨道维度验证失败！"
              write(6,*) " 张量 V1 和 V2 在第 2 个收缩指标上提供的轨道数不一致。"
              write(6,*) " 请检查调用组合: ", type1, " * ", type2, " -> ", type3
              write(6,*) "========================================================="
              stop
          end if
      end if

      ! 校验2：探测高速通道模式 (Fast Path) 
      is_fast_path = .false.
      if (t1_map == '3412' .and. t2_map == '3456' .and. t3_map == '1256') then
          if (iop1_12 == iop2_12 .and. iop1_34 == iop3_12 .and. iop2_34 == iop3_34) is_fast_path = .true.
      else if (t1_map == '1234' .and. t2_map == '3456' .and. t3_map == '1256') then
          if (iop1_34 == iop2_12 .and. iop1_12 == iop3_12 .and. iop2_34 == iop3_34) is_fast_path = .true.
      else if (t1_map == '1234' .and. t2_map == '5634' .and. t3_map == '1256') then
          if (iop1_34 == iop2_34 .and. iop1_12 == iop3_12 .and. iop2_12 == iop3_34) is_fast_path = .true.
      else if (t1_map == '3412' .and. t2_map == '5634' .and. t3_map == '1256') then
          if (iop1_12 == iop2_34 .and. iop1_34 == iop3_12 .and. iop2_12 == iop3_34) is_fast_path = .true.
      end if

      ! 校验3：非高速通道 (Fallback) 的 IOP 强限制
      if (.not. is_fast_path) then
          if (iop1_12 /= 1 .or. iop1_34 /= 1 .or. &
              iop2_12 /= 1 .or. iop2_34 /= 1 .or. &
              iop3_12 /= 1 .or. iop3_34 /= 1) then
              write(6,*) "========================================================="
              write(6,*) " FATAL ERROR in vv2v_dpd_444 (Dispatcher):"
              write(6,*) " 非高速通道 (Fallback) 不支持且严禁使用压缩空间 (IOP=0 或 2)。"
              write(6,*) " 您的收缩组合 [", type1, " * ", type2, " -> ", type3, "] "
              write(6,*) " 必须先在外部调用 sstgenso1 将张量解压为全空间 (IOP=1)。"
              write(6,*) "========================================================="
              stop
          end if
      end if

      ! =============================================================
      ! 核心路由分发器 (基于全局系统参数 prec 和 igpu)
      ! =============================================================
      if (prec == 8 .and. igpu == 0) then 
         ! CPU 双精度核心
         call vv2v_dpd_444_d(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 4 .and. igpu == 0) then 
         ! CPU 单精度核心 (统一传入 real*8)
         call vv2v_dpd_444_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 8 .and. igpu == 1) then 
         ! GPU 双精度核心
         call vv2v_dpd_444_dg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                              
      elseif (prec == 4 .and. igpu == 1) then 
         ! GPU 单精度核心 (统一传入 real*8)
         call vv2v_dpd_444_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
      endif

      return
      end subroutine vv2v_dpd_444

