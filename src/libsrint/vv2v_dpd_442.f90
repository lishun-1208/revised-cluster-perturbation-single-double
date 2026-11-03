! =======================================================================
! 终极工业级张量引擎 - 统一调度入口 (Dispatcher) [442 模式]
!
! 【模块功能】
!   作为 4-index * 4-index -> 2-index 收缩运算的唯一前端入口。
!   负责隔离外层算法与底层硬件，根据全局环境变量 (prec, igpu) 自动分发任务：
!     - vv2v_dpd_442_d  (CPU, 双精度)
!     - vv2v_dpd_442_s  (CPU, 单精度)
!     - vv2v_dpd_442_dg (GPU, 双精度)
!     - vv2v_dpd_442_sg (GPU, 单精度)
!
! 【核心特性 1：自由度自动映射翻译器 (Auto Index Mapper)】
!   允许外部传入任意别名组合（如 pqrm * pqrn -> mn）。
!   引擎会自动将其提取，将结果张量映射为 '12'，将 3 个收缩指标分别映射为 '3', '4', '5'，
!   并向下透传给底层数字敏感的执行核心。
!
! 【核心特性 2：智能 IOP 探测墙 (Smart IOP Barrier)】
!   原生支持 Macro-Index 宏维收缩拓扑！
!   若检测到未拆分的宏维 Pair (例如 pqrm * pqrn 中的 pq)，则放行 IOP=0 或 2 的压缩态。
!   对于独立散乱的收缩和生成维度，强制拦截并要求 IOP=1。
! =======================================================================
      subroutine vv2v_dpd_442(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, sym3, type3, fact1, fact2)
      
      use mod_size  ! 导入全局硬件与精度配置变量 (prec, igpu)
      use mod_symm  ! 对称性校验
      use mod_iop
      implicit none

      ! --- 接口参数声明 (对外维持标准的 Double Precision 签名) ---
      character(len=4), intent(in)  :: type1, type2 
      character(len=2), intent(in)  :: type3
      real*8, intent(in)            :: v1(*), v2(*)        
      real*8, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12
      integer, intent(in)           :: sym1, sym2, sym3     
      real*8, intent(in)            :: fact1, fact2        

      ! --- 自动映射字符变量 (Auto Index Mapper) ---
      character(len=4) :: t1_map, t2_map
      character(len=2) :: t3_map
      character(len=1) :: map_from(5), map_to(5)
      integer :: n_mapped, i, j
      logical :: is_contracted

      ! --- 合法性校验辅助变量 ---
      integer :: pop1(4,8), pop2(4,8)
      integer :: pos1_c1, pos1_c2, pos1_c3
      integer :: pos2_c1, pos2_c2, pos2_c3
      logical :: is_macro_LL, is_macro_RR, is_macro_RL ! [修改] 新增 is_macro_RL
      character(len=1) :: c_r, c_m, c_n

      ! =============================================================
      ! 【安全墙 1：全局对称性守恒校验】
      ! =============================================================
      if (dirprd(sym1, sym2) /= sym3) then 
         write(6,*) "========================================================="
         write(6,*) " FATAL ERROR in vv2v_dpd_442 (Dispatcher):"
         write(6,*) " 张量收缩全局对称性守恒遭到破坏！"
         write(6,*) " Sym1(", sym1, ") * Sym2(", sym2, ") != Sym3(", sym3, ")"
         write(6,*) "========================================================="
         stop
      end if

      ! =============================================================
      ! 【组件：自由度自动映射翻译器 (Automatic Index Mapper)】
      ! =============================================================
      
      ! 1. 强制将结果矩阵 v3 (type3) 的两个外部指标标定为 '1', '2'
      map_from(1) = type3(1:1); map_to(1) = '1'
      map_from(2) = type3(2:2); map_to(2) = '2'

      ! 2. 遍历 v1 (type1)，找出 3 个“收缩指标” (存在于 v1 但不在 v3 中)
      n_mapped = 2
      do i = 1, 4
          is_contracted = .true.
          do j = 1, 2
              if (type1(i:i) == type3(j:j)) then
                  is_contracted = .false.
                  exit
              end if
          end do
          if (is_contracted) then
              n_mapped = n_mapped + 1
              map_from(n_mapped) = type1(i:i)
              ! 动态分配为 '3', '4', '5'
              if (n_mapped == 3) map_to(3) = '3'
              if (n_mapped == 4) map_to(4) = '4'
              if (n_mapped == 5) map_to(5) = '5'
          end if
      end do

      ! 防御：442 模式严格要求有且仅有 3 个收缩指标
      if (n_mapped /= 5) then
         write(6,*) " FATAL ERROR: vv2v_dpd_442 requires EXACTLY 3 contracted indices."
         write(6,*) " Check your inputs: ", type1, " * ", type2, " -> ", type3
         stop
      end if

      ! 3. 将推导出的映射规则应用到局部的 t1_map, t2_map, t3_map
      t1_map = type1(1:4); t2_map = type2(1:4); t3_map = type3(1:2)
      
      do i = 1, 4
          do j = 1, 5
              if (type1(i:i) == map_from(j)) then
                  t1_map(i:i) = map_to(j); exit
              end if
          end do
          do j = 1, 5
              if (type2(i:i) == map_from(j)) then
                  t2_map(i:i) = map_to(j); exit
              end if
          end do
      end do
      do i = 1, 2
          do j = 1, 5
              if (type3(i:i) == map_from(j)) then
                  t3_map(i:i) = map_to(j); exit
              end if
          end do
      end do

      ! =============================================================
      ! 【安全墙 2：同源轨道维数一致性校验】
      ! =============================================================
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:); pop2(3,:) = p2_3(:); pop2(4,:) = p2_4(:)

      pos1_c1 = 0; pos1_c2 = 0; pos1_c3 = 0
      pos2_c1 = 0; pos2_c2 = 0; pos2_c3 = 0

      do i = 1, 4
         if (t1_map(i:i) == '3') pos1_c1 = i
         if (t1_map(i:i) == '4') pos1_c2 = i
         if (t1_map(i:i) == '5') pos1_c3 = i
         if (t2_map(i:i) == '3') pos2_c1 = i
         if (t2_map(i:i) == '4') pos2_c2 = i
         if (t2_map(i:i) == '5') pos2_c3 = i
      end do

      ! 验证收缩指标 '3'
      if (pos1_c1 /= 0 .and. pos2_c1 /= 0) then
          if (any(pop1(pos1_c1, :) /= pop2(pos2_c1, :))) stop " FATAL ERROR: Dimension mismatch on contracted index 1"
      else
          stop " FATAL ERROR: Contracted index 1 missing in V2."
      end if
      ! 验证收缩指标 '4'
      if (pos1_c2 /= 0 .and. pos2_c2 /= 0) then
          if (any(pop1(pos1_c2, :) /= pop2(pos2_c2, :))) stop " FATAL ERROR: Dimension mismatch on contracted index 2"
      else
          stop " FATAL ERROR: Contracted index 2 missing in V2."
      end if
      ! 验证收缩指标 '5'
      if (pos1_c3 /= 0 .and. pos2_c3 /= 0) then
          if (any(pop1(pos1_c3, :) /= pop2(pos2_c3, :))) stop " FATAL ERROR: Dimension mismatch on contracted index 3"
      else
          stop " FATAL ERROR: Contracted index 3 missing in V2."
      end if

      ! =============================================================
      ! 【安全墙 3：智能 IOP 探测墙 (Smart IOP Barrier)】
      ! 原生探测 442 特有的宏维收缩拓扑。
      ! =============================================================
      is_macro_LL = .false.; is_macro_RR = .false.; is_macro_RL = .false. ! [修改]
      c_r = ' '

      ! 探测左侧宏维 (Macro-LL)
      if (t1_map(1:2) == t2_map(1:2) .and. iop1_12 == iop2_12) then
          if (t1_map(3:3) == t2_map(3:3)) then
              c_r=t1_map(3:3); c_m=t1_map(4:4); c_n=t2_map(4:4)
          else if (t1_map(3:3) == t2_map(4:4)) then
              c_r=t1_map(3:3); c_m=t1_map(4:4); c_n=t2_map(3:3)
          else if (t1_map(4:4) == t2_map(3:3)) then
              c_r=t1_map(4:4); c_m=t1_map(3:3); c_n=t2_map(4:4)
          else if (t1_map(4:4) == t2_map(4:4)) then
              c_r=t1_map(4:4); c_m=t1_map(3:3); c_n=t2_map(3:3)
          end if
          if (c_r /= ' ' .and. ((t3_map(1:1)==c_m .and. t3_map(2:2)==c_n) .or. (t3_map(1:1)==c_n .and. t3_map(2:2)==c_m))) is_macro_LL = .true.
      ! 探测右侧宏维 (Macro-RR)
      else if (t1_map(3:4) == t2_map(3:4) .and. iop1_34 == iop2_34) then
          if (t1_map(1:1) == t2_map(1:1)) then
              c_r=t1_map(1:1); c_m=t1_map(2:2); c_n=t2_map(2:2)
          else if (t1_map(1:1) == t2_map(2:2)) then
              c_r=t1_map(1:1); c_m=t1_map(2:2); c_n=t2_map(1:1)
          else if (t1_map(2:2) == t2_map(1:1)) then
              c_r=t1_map(2:2); c_m=t1_map(1:1); c_n=t2_map(2:2)
          else if (t1_map(2:2) == t2_map(2:2)) then
              c_r=t1_map(2:2); c_m=t1_map(1:1); c_n=t2_map(1:1)
          end if
          if (c_r /= ' ' .and. ((t3_map(1:1)==c_m .and. t3_map(2:2)==c_n) .or. (t3_map(1:1)==c_n .and. t3_map(2:2)==c_m))) is_macro_RR = .true.
      ! [新增] 探测右-左交叉宏维 (Macro-RL)
      else if (t1_map(3:4) == t2_map(1:2) .and. iop1_34 == iop2_12) then
          if (t1_map(1:1) == t2_map(3:3)) then
              c_r=t1_map(1:1); c_m=t1_map(2:2); c_n=t2_map(4:4)
          else if (t1_map(1:1) == t2_map(4:4)) then
              c_r=t1_map(1:1); c_m=t1_map(2:2); c_n=t2_map(3:3)
          else if (t1_map(2:2) == t2_map(3:3)) then
              c_r=t1_map(2:2); c_m=t1_map(1:1); c_n=t2_map(4:4)
          else if (t1_map(2:2) == t2_map(4:4)) then
              c_r=t1_map(2:2); c_m=t1_map(1:1); c_n=t2_map(3:3)
          end if
          if (c_r /= ' ' .and. ((t3_map(1:1)==c_m .and. t3_map(2:2)==c_n) .or. (t3_map(1:1)==c_n .and. t3_map(2:2)==c_m))) is_macro_RL = .true.
      end if

      ! 拦截越权 IOP 配置
      if (is_macro_LL) then
          if (iop1_34 /= 1 .or. iop2_34 /= 1 .or. iop3_12 /= 1) then
              write(6,*) " FATAL ERROR: 442 Macro-LL must have IOP=1 for split indices."
              stop
          end if
      else if (is_macro_RR) then
          if (iop1_12 /= 1 .or. iop2_12 /= 1 .or. iop3_12 /= 1) then
              write(6,*) " FATAL ERROR: 442 Macro-RR must have IOP=1 for split indices."
              stop
          end if
      ! [新增] 校验 Macro-RL
      else if (is_macro_RL) then
          if (iop1_12 /= 1 .or. iop2_34 /= 1 .or. iop3_12 /= 1) then
              write(6,*) " FATAL ERROR: 442 Macro-RL must have IOP=1 for split indices."
              stop
          end if
      else
          if (iop1_12/=1 .or. iop1_34/=1 .or. iop2_12/=1 .or. iop2_34/=1 .or. iop3_12/=1) then
              write(6,*) "========================================================="
              write(6,*) " FATAL ERROR in vv2v_dpd_442 (Dispatcher):"
              write(6,*) " 智能 IOP 空间校验拦截！"
              write(6,*) " 纯 442 模式严格要求所有维度均为解压态 (IOP=1)。"
              write(6,*) " 逻辑维流向: ", t1_map, " * ", t2_map, " -> ", t3_map
              write(6,*) "========================================================="
              stop
          end if
      end if

      ! =============================================================
      ! 核心路由分发器 (基于全局系统参数 prec 和 igpu)
      ! =============================================================
      if (prec == 8 .and. igpu == 0) then 
         ! CPU 双精度核心
         call vv2v_dpd_442_d(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 4 .and. igpu == 0) then 
         ! CPU 单精度核心
         call vv2v_dpd_442_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                             v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                             v3, p3_1, p3_2, iop3_12, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 8 .and. igpu == 1) then 
         ! GPU 双精度核心
         call vv2v_dpd_442_dg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, sym3, t3_map, fact1, fact2)
                           
      elseif (prec == 4 .and. igpu == 1) then 
         ! GPU 单精度核心
         call vv2v_dpd_442_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, t1_map, &
                              v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, t2_map, &
                              v3, p3_1, p3_2, iop3_12, sym3, t3_map, fact1, fact2)
      endif

      return
      end subroutine vv2v_dpd_442
