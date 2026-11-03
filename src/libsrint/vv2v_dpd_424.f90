! =======================================================================
! 终极工业级张量引擎 - 统一调度入口 (Dispatcher) [424 模式]
!
! 【模块功能】
!   作为 4-index * 2-index -> 4-index 收缩运算的唯一前端入口。
!   负责隔离外层算法与底层硬件，根据全局环境变量 (prec, igpu) 自动分发任务：
!     - vv2v_dpd_424_d  (CPU, 双精度)
!     - vv2v_dpd_424_s  (CPU, 单精度)
!     - vv2v_dpd_424_dg (GPU, 双精度)
!     - vv2v_dpd_424_sg (GPU, 单精度)
!
! 【核心特性 1：自由度自动映射翻译器 (Auto Index Mapper)】
!   允许外部传入任意别名组合（如 pqrm * ms -> pqrs 或 pmqn * rmsn -> prqs）。
!   引擎会自动将其提取，将结果张量映射为 '1234'，将收缩指标映射为 '5'，
!   并向下透传给底层数字敏感的核心。
!
! 【核心特性 2：智能 IOP 探测墙 (Smart IOP Barrier)】
!   完美遵循 DPD 物理特性：只要 Pair (如 pq) 在收缩前后保持整体完整，
!   则放行 IOP=0 或 2。涉及收缩或拆分的维度强制要求 IOP=1。
! =======================================================================
      subroutine vv2v_dpd_424(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1, fact2)
   
      use mod_size  ! 导入全局硬件与精度配置变量 (prec, igpu)
      use mod_symm  ! 对称性校验
      use mod_iop
      implicit none

      ! --- 接口参数声明 (对外维持标准的 Double Precision 签名) ---
      character(len=4), intent(in)  :: type1, type3 
      character(len=2), intent(in)  :: type2
      real*8, intent(in)            :: v1(*), v2(*)        
      real*8, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8)
      integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop3_12, iop3_34
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1, fact2        

      ! --- 自动映射字符变量 (Auto Index Mapper) ---
      character(len=4) :: t1_map, t3_map
      character(len=2) :: t2_map
      character(len=1) :: map_from(5), map_to(5)
      integer :: n_mapped, i, j
      logical :: is_contracted

      ! --- 合法性校验辅助变量 ---
      integer :: pop1(4,8), pop2(2,8)
      integer :: pos1_c, pos2_c
      logical :: pair1_intact, pair2_intact, pair3_intact
      logical :: out1_intact, out2_intact, iop_valid

      ! =============================================================
      ! 【安全墙 1：全局对称性守恒校验】
      ! =============================================================
      if (dirprd(sym1, sym2) /= sym3) then 
         write(6,*) "========================================================="
         write(6,*) " FATAL ERROR in vv2v_dpd_424 (Dispatcher):"
         write(6,*) " 张量收缩全局对称性守恒遭到破坏！"
         write(6,*) " Sym1(", sym1, ") * Sym2(", sym2, ") != Sym3(", sym3, ")"
         write(6,*) "========================================================="
         stop
      end if

      ! =============================================================
      ! 【组件：自由度自动映射翻译器 (Automatic Index Mapper)】
      ! =============================================================
      
      ! 1. 强制将结果矩阵 v3 (type3) 的四个外部指标标定为 '1', '2', '3', '4'
      map_from(1) = type3(1:1); map_to(1) = '1'
      map_from(2) = type3(2:2); map_to(2) = '2'
      map_from(3) = type3(3:3); map_to(3) = '3'
      map_from(4) = type3(4:4); map_to(4) = '4'

      ! 2. 遍历 v1 (type1)，找出唯一的一个“收缩指标” (存在于 v1 但不在 v3 中)
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
              n_mapped = 5
              map_from(5) = type1(i:i)
              map_to(5) = '5'
          end if
      end do

      ! 防御：如果未能精确找到 1 个收缩指标，立即拦截
      if (n_mapped /= 5) then
         write(6,*) " FATAL ERROR: vv2v_dpd_424 requires EXACTLY 1 contracted index."
         write(6,*) " Check your inputs: ", type1, " * ", type2, " -> ", type3
         stop
      end if

      ! 3. 将推导出的映射规则应用到局部的 t1_map, t2_map, t3_map
      t1_map = type1(1:4); t2_map = type2(1:2); t3_map = type3(1:4)
      
      do i = 1, 4
          do j = 1, 5
              if (type1(i:i) == map_from(j)) then
                  t1_map(i:i) = map_to(j); exit
              end if
          end do
          do j = 1, 5
              if (type3(i:i) == map_from(j)) then
                  t3_map(i:i) = map_to(j); exit
              end if
          end do
      end do
      do i = 1, 2
          do j = 1, 5
              if (type2(i:i) == map_from(j)) then
                  t2_map(i:i) = map_to(j); exit
              end if
          end do
      end do

      ! =============================================================
      ! 【安全墙 2：同源轨道维数一致性校验】
      ! =============================================================
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:)

      pos1_c = 0; pos2_c = 0
      do i = 1, 4; if (t1_map(i:i) == '5') pos1_c = i; end do
      do i = 1, 2; if (t2_map(i:i) == '5') pos2_c = i; end do

      if (pos1_c /= 0 .and. pos2_c /= 0) then
          if (any(pop1(pos1_c, :) /= pop2(pos2_c, :))) then
              write(6,*) "========================================================="
              write(6,*) " FATAL ERROR in vv2v_dpd_424 (Dispatcher):"
              write(6,*) " 收缩轨道维度验证失败！"
              write(6,*) " 张量 V1 和 V2 在参与收缩的指标 '", map_from(5), "' 上轨道数不一致。"
              write(6,*) " 请检查调用组合: ", type1, " * ", type2, " -> ", type3
              write(6,*) "========================================================="
              stop
          end if
      else
          write(6,*) " FATAL ERROR: Contracted index '", map_from(5), "' missing in V2."
          stop
      end if

      ! =============================================================
      ! 【安全墙 3：智能 IOP 探测墙 (Smart IOP Barrier)】
      ! 基于 DPD 复合维度特性：
      ! 只要一个 Pair 在收缩映射前后保持绝对完整 (顺序不变且未拆分)，
      ! 则其可被安全视为一个连续的大维度，放行 IOP=0/2。
      ! 若涉及拆分、重组或转置 (例如 pqrm*ms -> prqs 中的 pq 被拆开)，则强校验拦截。
      ! =============================================================
      
      ! 验证输入端 Pair 是否能完整流向输出端
      pair1_intact = (t1_map(1:2) == t3_map(1:2) .and. iop1_12 == iop3_12) .or. &
                     (t1_map(1:2) == t3_map(3:4) .and. iop1_12 == iop3_34)
      pair2_intact = (t1_map(3:4) == t3_map(1:2) .and. iop1_34 == iop3_12) .or. &
                     (t1_map(3:4) == t3_map(3:4) .and. iop1_34 == iop3_34)
      pair3_intact = (t2_map(1:2) == t3_map(1:2) .and. iop2_12 == iop3_12) .or. &
                     (t2_map(1:2) == t3_map(3:4) .and. iop2_12 == iop3_34)
                     
      ! 验证输出端 Pair 是否直接来源于输入端的完整块
      out1_intact = (t3_map(1:2) == t1_map(1:2) .and. iop3_12 == iop1_12) .or. &
                    (t3_map(1:2) == t1_map(3:4) .and. iop3_12 == iop1_34) .or. &
                    (t3_map(1:2) == t2_map(1:2) .and. iop3_12 == iop2_12)
      out2_intact = (t3_map(3:4) == t1_map(1:2) .and. iop3_34 == iop1_12) .or. &
                    (t3_map(3:4) == t1_map(3:4) .and. iop3_34 == iop1_34) .or. &
                    (t3_map(3:4) == t2_map(1:2) .and. iop3_34 == iop2_12)

      iop_valid = .true.
      if (iop1_12 /= 1 .and. .not. pair1_intact) iop_valid = .false.
      if (iop1_34 /= 1 .and. .not. pair2_intact) iop_valid = .false.
      if (iop2_12 /= 1 .and. .not. pair3_intact) iop_valid = .false.
      if (iop3_12 /= 1 .and. .not. out1_intact)  iop_valid = .false.
      if (iop3_34 /= 1 .and. .not. out2_intact)  iop_valid = .false.

      if (.not. iop_valid) then
          write(6,*) "========================================================="
          write(6,*) " FATAL ERROR in vv2v_dpd_424 (Dispatcher):"
          write(6,*) " 智能 IOP 空间校验拦截！"
          write(6,*) " 只有在收缩前后保持完整的 Pair (例如 pq -> pq) 允许 IOP=0/2。"
          write(6,*) " 涉及收缩、拆分或内部转置的维度 (例如 rm, ms, pr) 必须解压为 IOP=1。"
          write(6,*) " 请检查您的运算张量，或使用 sstgenso1 将其转化为 Full Space。"
          write(6,*) " 逻辑维流向: ", t1_map, " * ", t2_map, " -> ", t3_map
          write(6,*) "========================================================="
          stop
      end if

      ! =============================================================
      ! 核心路由分发器 (基于全局系统参数 prec 和 igpu)
      ! =============================================================
      if (prec == 8 .and. igpu == 0) then 
         ! CPU 双精度核心
         call vv2v_dpd_424_d(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 4 .and. igpu == 0) then 
         ! CPU 单精度核心 (统一传入 real*8，由底层截断)
         call vv2v_dpd_424_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 8 .and. igpu == 1) then 
         ! GPU 双精度核心
         call vv2v_dpd_424_dg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
                             
      elseif (prec == 4 .and. igpu == 1) then 
         ! GPU 单精度核心 (统一传入 real*8，由底层截断)
         call vv2v_dpd_424_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, t3_map, fact1, fact2)
      endif

      return
      end subroutine vv2v_dpd_424
