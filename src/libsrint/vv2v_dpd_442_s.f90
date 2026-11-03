! =======================================================================
! 终极工业级张量引擎 - CPU 单精度核心 (Single Precision, CPU)
! [核心能力] 442 模式: 4-index * 4-index -> 2-index 
! [适用场景] 3指标收缩运算 (例如: pqrm * pqrn -> mn)
! [架构设计] 1. 原生支持 Macro-Index (宏维) 探测，完美支持未拆分 Pair 的 IOP=0/2 压缩态。
!            2. 彻底消灭竞态条件！在任务分配上按 V3 块严格正交隔离，内层累加 R。
!            3. 混合并行线控 (MKL/OpenBLAS)，防止多线程嵌套爆炸。
!            4. 专门针对单精度 real*4 进行了内存占用和 Block Size 的重调优。
! =======================================================================
      subroutine vv2v_dpd_442_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                                v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
                                v3, p3_1, p3_2, iop3_12, sym3, type3, fact1, fact2)
      use mod_symm
      use omp_lib
#if defined(USE_OPENBLAS)
      use iso_c_binding
#endif
      implicit none

      ! =============================================================
      ! 底层线程控制宏定义
      ! =============================================================
#ifdef USE_MKL
      interface
         subroutine mkl_set_num_threads(num)
            integer, intent(in) :: num
         end subroutine mkl_set_num_threads
         integer function mkl_get_max_threads()
         end function mkl_get_max_threads
      end interface
      integer :: old_threads
#elif defined(USE_OPENBLAS)
      interface
         subroutine openblas_set_num_threads(num) bind(c, name="openblas_set_num_threads")
            import :: c_int
            integer(c_int), value :: num
         end subroutine openblas_set_num_threads
         function openblas_get_num_threads() bind(c, name="openblas_get_num_threads")
            import :: c_int
            integer(c_int) :: openblas_get_num_threads
         end function openblas_get_num_threads
      end interface
      integer(c_int) :: old_threads
#else
      integer :: old_threads
#endif

      ! --- 接口参数声明 (注意: 数据数组为 real*4, 但 fact 为 real*8 保持签名兼容) ---
      character(len=4), intent(in)  :: type1, type2 
      character(len=2), intent(in)  :: type3        
      real*4, intent(in)            :: v1(*), v2(*)        
      real*4, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1, fact2        

      ! 外部接口 (替换为单精度 xgemms)
      external :: xgemms
      integer, external :: irpdso  

      ! --- 路由解析变量 ---
      character(len=1) :: c_con(3), c_r, c_m, c_n
      integer :: pos1_m, pos2_n, pos1_k(3), pos2_k(3), v3_map_m, v3_map_n
      integer :: pos1_r, pos2_r
      integer :: i, j, idx, max_threads
      logical :: is_macro_LL, is_macro_RR, is_macro_RL, found_r
      integer :: pop1(4,8), pop2(4,8), pop3(2,8)
      
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposlft1_fb(33,33), iposrht1(33,33)
      integer :: inumlft2(33), inumrht2(33), iabsoff2(33), iposlft2_fb(33,33), iposrht2(33,33)
      integer :: iabsoff3(33)
      
      integer :: num_tasks, task_list(10000)
      integer*8 :: max_mem_bytes_per_thread, req_mem_bytes
      integer*8 :: global_max_m, global_max_n, global_max_k

      ! --- OpenMP 线程私有变量 (局部缓存设为 real*4) ---
      integer :: t_idx, irr_M, irr_N, irr_R, irr_PQ, irr_l1, irr_r1, irr_l2, irr_r2, irr_v3_2, irr_k1, irr_k2, irr_k3
      integer :: cur_M, cur_N, cur_K, BS, dim_PQ, dim_M, dim_N, dim_R, dim_k(3), irr_logi1(4), irr_logi2(4)
      integer :: m_start, n_start, k_start, cur_bs_m, cur_bs_n, cur_bs_k
      integer :: m_l, n_l, k_l, m_g, n_g, k_g, i_pq, i_r, k_v(3), rem, idx_v3
      integer :: base_v1, base_v2, base_v3, b_v1_a, b_v2_a, b_v3_a, s1(4), s2(4), s3_m, s3_n, s_phys1(4), s_phys2(4)
      integer :: irr_pos_m, irr_pos_n
      integer*8 :: req_mem_bytes_local
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)

      ! =============================================================
      ! 【阶段 1：宏维 (Macro-Index) 探测逻辑】
      ! =============================================================
      is_macro_LL = .false.; is_macro_RR = .false.; is_macro_RL = .false.
      c_r = ' '

      if (type1(1:2) == type2(1:2) .and. iop1_12 == iop2_12) then
          if (type1(3:3) == type2(3:3)) then
              c_r=type1(3:3); c_m=type1(4:4); c_n=type2(4:4); pos1_r=3; pos1_m=4; pos2_r=3; pos2_n=4
          else if (type1(3:3) == type2(4:4)) then
              c_r=type1(3:3); c_m=type1(4:4); c_n=type2(3:3); pos1_r=3; pos1_m=4; pos2_r=4; pos2_n=3
          else if (type1(4:4) == type2(3:3)) then
              c_r=type1(4:4); c_m=type1(3:3); c_n=type2(4:4); pos1_r=4; pos1_m=3; pos2_r=3; pos2_n=4
          else if (type1(4:4) == type2(4:4)) then
              c_r=type1(4:4); c_m=type1(3:3); c_n=type2(3:3); pos1_r=4; pos1_m=3; pos2_r=4; pos2_n=3
          end if
          if (c_r /= ' ' .and. ((type3(1:1)==c_m .and. type3(2:2)==c_n) .or. (type3(1:1)==c_n .and. type3(2:2)==c_m))) is_macro_LL = .true.
      else if (type1(3:4) == type2(3:4) .and. iop1_34 == iop2_34) then
          if (type1(1:1) == type2(1:1)) then
              c_r=type1(1:1); c_m=type1(2:2); c_n=type2(2:2); pos1_r=1; pos1_m=2; pos2_r=1; pos2_n=2
          else if (type1(1:1) == type2(2:2)) then
              c_r=type1(1:1); c_m=type1(2:2); c_n=type2(1:1); pos1_r=1; pos1_m=2; pos2_r=2; pos2_n=1
          else if (type1(2:2) == type2(1:1)) then
              c_r=type1(2:2); c_m=type1(1:1); c_n=type2(2:2); pos1_r=2; pos1_m=1; pos2_r=1; pos2_n=2
          else if (type1(2:2) == type2(2:2)) then
              c_r=type1(2:2); c_m=type1(1:1); c_n=type2(1:1); pos1_r=2; pos1_m=1; pos2_r=2; pos2_n=1
          end if
          if (c_r /= ' ' .and. ((type3(1:1)==c_m .and. type3(2:2)==c_n) .or. (type3(1:1)==c_n .and. type3(2:2)==c_m))) is_macro_RR = .true.
      else if (type1(3:4) == type2(1:2) .and. iop1_34 == iop2_12) then
          if (type1(1:1) == type2(3:3)) then
              c_r=type1(1:1); c_m=type1(2:2); c_n=type2(4:4); pos1_r=1; pos1_m=2; pos2_r=3; pos2_n=4
          else if (type1(1:1) == type2(4:4)) then
              c_r=type1(1:1); c_m=type1(2:2); c_n=type2(3:3); pos1_r=1; pos1_m=2; pos2_r=4; pos2_n=3
          else if (type1(2:2) == type2(3:3)) then
              c_r=type1(2:2); c_m=type1(1:1); c_n=type2(4:4); pos1_r=2; pos1_m=1; pos2_r=3; pos2_n=4
          else if (type1(2:2) == type2(4:4)) then
              c_r=type1(2:2); c_m=type1(1:1); c_n=type2(3:3); pos1_r=2; pos1_m=1; pos2_r=4; pos2_n=3
          end if
          if (c_r /= ' ' .and. ((type3(1:1)==c_m .and. type3(2:2)==c_n) .or. (type3(1:1)==c_n .and. type3(2:2)==c_m))) is_macro_RL = .true.
      end if

      if (is_macro_LL) then
          if (iop1_34 /= 1 .or. iop2_34 /= 1 .or. iop3_12 /= 1) stop " FATAL: 442 Macro-LL must have IOP=1 for split indices."
      else if (is_macro_RR) then
          if (iop1_12 /= 1 .or. iop2_12 /= 1 .or. iop3_12 /= 1) stop " FATAL: 442 Macro-RR must have IOP=1 for split indices."
      else if (is_macro_RL) then
          if (iop1_12 /= 1 .or. iop2_34 /= 1 .or. iop3_12 /= 1) stop " FATAL: 442 Macro-RL must have IOP=1 for split indices."
      else
          if (iop1_12/=1 .or. iop1_34/=1 .or. iop2_12/=1 .or. iop2_34/=1 .or. iop3_12/=1) &
             stop " FATAL: Pure 442 strictly requires all IOP=1 as indices are scattered."
      end if

      ! =============================================================
      ! 【阶段 2：构建 DPD 物理表】
      ! =============================================================
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:); pop2(3,:) = p2_3(:); pop2(4,:) = p2_4(:)
      pop3(1,:) = p3_1(:); pop3(2,:) = p3_2(:)

      do i = 1, nirrep
         inumlft1(i)=irpdso(i,p1_1,p1_2,iop1_12); inumrht1(i)=irpdso(i,p1_3,p1_4,iop1_34)
         inumlft2(i)=irpdso(i,p2_1,p2_2,iop2_12); inumrht2(i)=irpdso(i,p2_3,p2_4,iop2_34)
      end do
      
      iabsoff1(1)=1; iabsoff2(1)=1; iabsoff3(1)=1
      do i = 1, nirrep
         iabsoff1(i+1)=iabsoff1(i)+inumlft1(dirprd(i,sym1))*inumrht1(i)
         iabsoff2(i+1)=iabsoff2(i)+inumlft2(dirprd(i,sym2))*inumrht2(i)
         iposlft1_fb(1,i)=0; iposrht1(1,i)=0; iposlft2_fb(1,i)=0; iposrht2(1,i)=0
         do j = 1, nirrep-1
            iposlft1_fb(j+1,i)=iposlft1_fb(j,i)+pop1(1,dirprd(j,i))*pop1(2,j)
            iposrht1(j+1,i)=iposrht1(j,i)+pop1(3,dirprd(j,i))*pop1(4,j)
            iposlft2_fb(j+1,i)=iposlft2_fb(j,i)+pop2(1,dirprd(j,i))*pop2(2,j)
            iposrht2(j+1,i)=iposrht2(j,i)+pop2(3,dirprd(j,i))*pop2(4,j)
         end do
      end do
      do i = 1, nirrep-1; iabsoff3(i+1)=iabsoff3(i)+pop3(1,dirprd(i,sym3))*pop3(2,i); end do

      ! =============================================================
      ! 【阶段 3：线程管理与任务分发探测】
      ! =============================================================
#ifdef USE_MKL
      old_threads = mkl_get_max_threads()
      call mkl_set_num_threads(1)
#elif defined(USE_OPENBLAS)
      old_threads = openblas_get_num_threads()
      call openblas_set_num_threads(1_c_int)
#endif
      max_threads = omp_get_max_threads()
      max_mem_bytes_per_thread = 3500000000_8 / max_threads

      num_tasks = 0; global_max_m=0; global_max_n=0; global_max_k=0
      
      if (is_macro_LL .or. is_macro_RR .or. is_macro_RL) then
          do i = 1, nirrep
             j = dirprd(sym3, i)
             if (pop1(pos1_m, i)==0 .or. pop2(pos2_n, j)==0) cycle
             
             found_r = .false.
             do idx = 1, nirrep
                if (pop1(pos1_r, idx)==0) cycle
                if (is_macro_LL .and. inumlft1(dirprd(sym1, dirprd(idx, i))) == 0) cycle
                if (is_macro_RR .and. inumrht1(dirprd(sym1, dirprd(idx, i))) == 0) cycle
                if (is_macro_RL .and. inumrht1(dirprd(sym1, dirprd(idx, i))) == 0) cycle
                
                found_r = .true.
                if (is_macro_LL) global_max_k = max(global_max_k, int(inumlft1(dirprd(sym1,dirprd(idx,i))),8)*pop1(pos1_r,idx))
                if (is_macro_RR) global_max_k = max(global_max_k, int(inumrht1(dirprd(sym1,dirprd(idx,i))),8)*pop1(pos1_r,idx))
                if (is_macro_RL) global_max_k = max(global_max_k, int(inumrht1(dirprd(sym1,dirprd(idx,i))),8)*pop1(pos1_r,idx))
             end do
             
             if (found_r) then
                num_tasks = num_tasks + 1; task_list(num_tasks) = i
                global_max_m = max(global_max_m, int(pop1(pos1_m, i), 8))
                global_max_n = max(global_max_n, int(pop2(pos2_n, j), 8))
             end if
          end do
      else
          do i = 1, 2
             if (index(type1, type3(i:i)) > 0) then; pos1_m = index(type1, type3(i:i)); v3_map_m = i; end if
             if (index(type2, type3(i:i)) > 0) then; pos2_n = index(type2, type3(i:i)); v3_map_n = i; end if
          end do
          idx=0; do i=1,4; if(i/=pos1_m) then; idx=idx+1; c_con(idx)=type1(i:i); pos1_k(idx)=i; end if; end do
          do i=1,3; pos2_k(i)=index(type2, c_con(i)); end do
          
          do i = 1, nirrep
             if (pop3(2, i) == 0 .or. pop3(1, dirprd(i, sym3)) == 0) cycle
             num_tasks = num_tasks + 1; task_list(num_tasks) = i
             
             irr_m = i; irr_n = dirprd(i, sym3)
             if (v3_map_m == 1) then; irr_m = dirprd(i, sym3); irr_n = i; end if
             
             global_max_m = max(global_max_m, int(pop1(pos1_m, irr_m), 8))
             global_max_n = max(global_max_n, int(pop2(pos2_n, irr_n), 8))
             
             cur_K = 0
             do irr_k1=1,nirrep; do irr_k2=1,nirrep; irr_k3=dirprd(sym1,dirprd(irr_m,dirprd(irr_k1,irr_k2)))
                irr_logi2(pos2_n)=irr_n; irr_logi2(pos2_k(1))=irr_k1; irr_logi2(pos2_k(2))=irr_k2; irr_logi2(pos2_k(3))=irr_k3
                if(dirprd(dirprd(irr_logi2(1),irr_logi2(2)),dirprd(irr_logi2(3),irr_logi2(4)))/=sym2) cycle
                if(pop1(pos1_k(1),irr_k1)==0 .or. pop1(pos1_k(2),irr_k2)==0 .or. pop1(pos1_k(3),irr_k3)==0) cycle
                if(pop2(pos2_k(1),irr_k1)==0 .or. pop2(pos2_k(2),irr_k2)==0 .or. pop2(pos2_k(3),irr_k3)==0) cycle
                cur_K = cur_K + pop1(pos1_k(1),irr_k1)*pop1(pos1_k(2),irr_k2)*pop1(pos1_k(3),irr_k3)
             end do; end do
             global_max_k = max(global_max_k, int(cur_K, 8))
          end do
      end if
      if (num_tasks == 0) goto 999

      ! =============================================================
      ! 【阶段 4：执行引擎】
      ! =============================================================
!$omp parallel default(shared) &
!$omp private(t_idx, irr_M, irr_N, irr_R, irr_PQ, irr_l1, irr_r1, irr_l2, irr_r2, irr_v3_2, irr_k1, irr_k2, irr_k3, &
!$omp         cur_M, cur_N, cur_K, BS, dim_PQ, dim_M, dim_N, dim_R, dim_k, irr_logi1, irr_logi2, &
!$omp         m_start, n_start, k_start, cur_bs_m, cur_bs_n, cur_bs_k, &
!$omp         m_l, n_l, k_l, m_g, n_g, k_g, i_pq, i_r, k_v, rem, idx_v3, &
!$omp         base_v1, base_v2, base_v3, b_v1_a, b_v2_a, b_v3_a, s1, s2, s3_m, s3_n, s_phys1, s_phys2, &
!$omp         irr_pos_m, irr_pos_n, req_mem_bytes_local, buf_A, buf_B, buf_C)

         ! 针对单精度 real*4 重调优内存体积与瓦片 BS (除以 12.0d0 代替 24.0d0)
         req_mem_bytes_local = (global_max_m*global_max_k + global_max_k*global_max_n + global_max_m*global_max_n) * 4_8
         if (req_mem_bytes_local > max_mem_bytes_per_thread) then
            BS = int(sqrt(max_mem_bytes_per_thread / 12.0d0))
         else
            BS = max(int(global_max_m), int(global_max_n), int(global_max_k))
         end if
         if (BS > 12000) BS = 12000; if (BS < 100) BS = 100
         allocate(buf_A(BS, BS), buf_B(BS, BS), buf_C(BS, BS))

!$omp do schedule(dynamic, 1)
         do t_idx = 1, num_tasks
            ! --- 分支 A/B/D: Macro-Index 宏维极速收缩 ---
            if (is_macro_LL .or. is_macro_RR .or. is_macro_RL) then
                irr_M = task_list(t_idx)
                if (type3(1:1)==c_m) then; irr_N = dirprd(sym3, irr_M); base_v3=iabsoff3(irr_N)
                else; irr_N = irr_M; irr_M = dirprd(sym3, irr_N); base_v3=iabsoff3(irr_M); end if
                
                dim_M = pop1(pos1_m, irr_M); dim_N = pop2(pos2_n, irr_N)
                if (type3(1:1)==c_m) then; s3_m=1; s3_n=dim_M; else; s3_n=1; s3_m=dim_N; end if
                b_v3_a = base_v3 - s3_m - s3_n
                
                cur_M = dim_M; cur_N = dim_N
                
                do m_start = 1, cur_M, BS; cur_bs_m = min(BS, cur_M - m_start + 1)
                do n_start = 1, cur_N, BS; cur_bs_n = min(BS, cur_N - n_start + 1)
                   
                   buf_C(1:cur_bs_m, 1:cur_bs_n) = 0.0e0 
                   
                   do irr_R = 1, nirrep
                      dim_R = pop1(pos1_r, irr_R)
                      if (dim_R == 0) cycle
                      
                      if (is_macro_LL) then
                         irr_PQ = dirprd(sym1, dirprd(irr_R, irr_M)); dim_PQ = inumlft1(irr_PQ)
                         if (dim_PQ == 0) cycle
                         irr_r1=dirprd(irr_R,irr_M); irr_r2=dirprd(irr_R,irr_N)
                         irr_pos_m = irr_M; if(pos1_r==4) irr_pos_m=irr_R
                         base_v1 = iabsoff1(irr_r1) + dim_PQ * iposrht1(irr_pos_m, irr_r1); s1(1)=1
                         if(pos1_r==3) then; s1(2)=dim_PQ; s1(3)=dim_PQ*dim_R; else; s1(3)=dim_PQ; s1(2)=dim_PQ*dim_M; end if
                         irr_pos_n = irr_N; if(pos2_r==4) irr_pos_n=irr_R
                         base_v2 = iabsoff2(irr_r2) + dim_PQ * iposrht2(irr_pos_n, irr_r2); s2(1)=1
                         if(pos2_r==3) then; s2(2)=dim_PQ; s2(3)=dim_PQ*dim_R; else; s2(3)=dim_PQ; s2(2)=dim_PQ*dim_N; end if
                      
                      else if (is_macro_RR) then
                         irr_PQ = dirprd(sym1, dirprd(irr_R, irr_M)); dim_PQ = inumrht1(irr_PQ)
                         if (dim_PQ == 0) cycle
                         irr_l1=dirprd(irr_R,irr_M); irr_l2=dirprd(irr_R,irr_N)
                         irr_pos_m = irr_M; if(pos1_r==2) irr_pos_m=irr_R
                         base_v1 = iabsoff1(irr_PQ) + iposlft1_fb(irr_pos_m, irr_l1); s1(1)=inumlft1(irr_l1)
                         if(pos1_r==1) then; s1(2)=1; s1(3)=dim_R; else; s1(3)=1; s1(2)=dim_M; end if
                         irr_pos_n = irr_N; if(pos2_r==2) irr_pos_n=irr_R
                         base_v2 = iabsoff2(irr_PQ) + iposlft2_fb(irr_pos_n, irr_l2); s2(1)=inumlft2(irr_l2)
                         if(pos2_r==1) then; s2(2)=1; s2(3)=dim_R; else; s2(3)=1; s2(2)=dim_N; end if
                      
                      else if (is_macro_RL) then
                         irr_PQ = dirprd(sym1, dirprd(irr_R, irr_M))
                         dim_PQ = inumrht1(irr_PQ)
                         if (dim_PQ == 0) cycle
                         
                         ! --- V1 采取右侧宏维逻辑 ---
                         irr_l1 = dirprd(irr_R, irr_M)
                         irr_pos_m = irr_M
                         if (pos1_r == 2) irr_pos_m = irr_R
                         base_v1 = iabsoff1(irr_PQ) + iposlft1_fb(irr_pos_m, irr_l1)
                         s1(1) = inumlft1(irr_l1)
                         if (pos1_r == 1) then; s1(2) = 1; s1(3) = dim_R
                         else; s1(3) = 1; s1(2) = dim_M; end if

                         ! --- V2 采取左侧宏维逻辑 ---
                         irr_r2 = dirprd(irr_R, irr_N)
                         irr_pos_n = irr_N
                         if (pos2_r == 4) irr_pos_n = irr_R
                         base_v2 = iabsoff2(irr_r2) + dim_PQ * iposrht2(irr_pos_n, irr_r2)
                         s2(1) = 1
                         if (pos2_r == 3) then; s2(2) = dim_PQ; s2(3) = dim_PQ * dim_R
                         else; s2(3) = dim_PQ; s2(2) = dim_PQ * dim_N; end if
                      end if
                      
                      b_v1_a = base_v1-s1(1)-s1(2)-s1(3); b_v2_a = base_v2-s2(1)-s2(2)-s2(3)
                      cur_K = dim_PQ * dim_R

                      do k_start = 1, cur_K, BS; cur_bs_k = min(BS, cur_K - k_start + 1)
                         do m_l = 1, cur_bs_m; m_g = m_start + m_l - 1
                         do k_l = 1, cur_bs_k; k_g = k_start + k_l - 1; i_pq = mod(k_g-1, dim_PQ)+1; i_r = (k_g-1)/dim_PQ+1
                            buf_A(k_l, m_l) = v1(b_v1_a + i_pq*s1(1) + i_r*s1(2) + m_g*s1(3)); end do; end do
                         
                         do n_l = 1, cur_bs_n; n_g = n_start + n_l - 1
                         do k_l = 1, cur_bs_k; k_g = k_start + k_l - 1; i_pq = mod(k_g-1, dim_PQ)+1; i_r = (k_g-1)/dim_PQ+1
                            buf_B(k_l, n_l) = v2(b_v2_a + i_pq*s2(1) + i_r*s2(2) + n_g*s2(3)); end do; end do
                            
                         ! 单精度计算: xgemms
                         call xgemms('T', 'N', cur_bs_m, cur_bs_n, cur_bs_k, 1.0e0, buf_A, BS, buf_B, BS, 1.0e0, buf_C, BS)
                      end do
                   end do 

                   do n_l = 1, cur_bs_n; n_g = n_start + n_l - 1
                   do m_l = 1, cur_bs_m; m_g = m_start + m_l - 1; idx_v3 = b_v3_a + m_g*s3_m + n_g*s3_n
                      ! 注意：安全强制类型转换 real(..., 4) 防止实数精度越界截断出错
                      if (fact2 == 0.0d0) then
                          v3(idx_v3) = real(fact1, 4) * buf_C(m_l, n_l)
                      else
                          v3(idx_v3) = real(fact2, 4) * v3(idx_v3) + real(fact1, 4) * buf_C(m_l, n_l)
                      end if
                   end do; end do
                   
                end do; end do

            ! --- 分支 C: 标准 442 完全独立收缩 (K=k1*k2*k3) ---
            else
                irr_v3_2 = task_list(t_idx)
                if (v3_map_m == 2) then; irr_m = irr_v3_2; irr_n = dirprd(irr_v3_2, sym3)
                else; irr_n = irr_v3_2; irr_m = dirprd(irr_v3_2, sym3); end if
                
                cur_M = pop1(pos1_m, irr_m); cur_N = pop2(pos2_n, irr_n)
                base_v3 = iabsoff3(irr_v3_2)
                if (v3_map_m == 1) then; s3_m = 1; s3_n = cur_M; else; s3_m = cur_N; s3_n = 1; end if
                b_v3_a = base_v3 - s3_m - s3_n

                do m_start = 1, cur_M, BS; cur_bs_m = min(BS, cur_M - m_start + 1)
                do n_start = 1, cur_N, BS; cur_bs_n = min(BS, cur_N - n_start + 1)
                   buf_C(1:cur_bs_m, 1:cur_bs_n) = 0.0e0
                   do irr_k1=1,nirrep; do irr_k2=1,nirrep; irr_k3=dirprd(sym1,dirprd(irr_m,dirprd(irr_k1,irr_k2)))
                      irr_logi1(pos1_m)=irr_m; irr_logi1(pos1_k(1))=irr_k1; irr_logi1(pos1_k(2))=irr_k2; irr_logi1(pos1_k(3))=irr_k3
                      irr_logi2(pos2_n)=irr_n; irr_logi2(pos2_k(1))=irr_k1; irr_logi2(pos2_k(2))=irr_k2; irr_logi2(pos2_k(3))=irr_k3
                      if(dirprd(dirprd(irr_logi2(1),irr_logi2(2)),dirprd(irr_logi2(3),irr_logi2(4)))/=sym2) cycle
                      
                      dim_k(1)=pop1(pos1_k(1),irr_k1); dim_k(2)=pop1(pos1_k(2),irr_k2); dim_k(3)=pop1(pos1_k(3),irr_k3)
                      if(dim_k(1)==0 .or. dim_k(2)==0 .or. dim_k(3)==0) cycle
                      if(pop2(pos2_k(1),irr_k1)==0 .or. pop2(pos2_k(2),irr_k2)==0 .or. pop2(pos2_k(3),irr_k3)==0) cycle
                      cur_K = dim_k(1)*dim_k(2)*dim_k(3)
                     
                      irr_r1=dirprd(irr_logi1(3),irr_logi1(4)); irr_l1=dirprd(irr_logi1(1),irr_logi1(2))
                      irr_r2=dirprd(irr_logi2(3),irr_logi2(4)); irr_l2=dirprd(irr_logi2(1),irr_logi2(2))
                      
                      base_v1 = iabsoff1(irr_r1) + inumlft1(irr_l1)*iposrht1(irr_logi1(4),irr_r1) + iposlft1_fb(irr_logi1(2),irr_l1)
                      s_phys1(1)=1; s_phys1(2)=pop1(1,irr_logi1(1)); s_phys1(3)=inumlft1(irr_l1); s_phys1(4)=inumlft1(irr_l1)*pop1(3,irr_logi1(3))
                      
                      base_v2 = iabsoff2(irr_r2) + inumlft2(irr_l2)*iposrht2(irr_logi2(4),irr_r2) + iposlft2_fb(irr_logi2(2),irr_l2)
                      s_phys2(1)=1; s_phys2(2)=pop2(1,irr_logi2(1)); s_phys2(3)=inumlft2(irr_l2); s_phys2(4)=inumlft2(irr_l2)*pop2(3,irr_logi2(3))
                      
                      b_v1_a = base_v1-s_phys1(pos1_m)-s_phys1(pos1_k(1))-s_phys1(pos1_k(2))-s_phys1(pos1_k(3))
                      b_v2_a = base_v2-s_phys2(pos2_n)-s_phys2(pos2_k(1))-s_phys2(pos2_k(2))-s_phys2(pos2_k(3))
                      
                      do k_start = 1, cur_K, BS; cur_bs_k = min(BS, cur_K - k_start + 1)
                         do m_l = 1, cur_bs_m; m_g = m_start + m_l - 1
                         do k_l = 1, cur_bs_k; k_g = k_start + k_l - 1
                            k_v(1)=mod(k_g-1, dim_k(1))+1; rem=(k_g-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                            buf_A(k_l, m_l) = v1(b_v1_a + m_g*s_phys1(pos1_m) + k_v(1)*s_phys1(pos1_k(1)) + k_v(2)*s_phys1(pos1_k(2)) + k_v(3)*s_phys1(pos1_k(3)))
                         end do; end do
                         
                         do n_l = 1, cur_bs_n; n_g = n_start + n_l - 1
                         do k_l = 1, cur_bs_k; k_g = k_start + k_l - 1
                            k_v(1)=mod(k_g-1, dim_k(1))+1; rem=(k_g-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                            buf_B(k_l, n_l) = v2(b_v2_a + n_g*s_phys2(pos2_n) + k_v(1)*s_phys2(pos2_k(1)) + k_v(2)*s_phys2(pos2_k(2)) + k_v(3)*s_phys2(pos2_k(3)))
                         end do; end do
                         
                         call xgemms('T', 'N', cur_bs_m, cur_bs_n, cur_bs_k, 1.0e0, buf_A, BS, buf_B, BS, 1.0e0, buf_C, BS)
                      end do
             
                   end do; end do
                   
                   do n_l = 1, cur_bs_n; n_g = n_start + n_l - 1
                   do m_l = 1, cur_bs_m; m_g = m_start + m_l - 1; idx_v3 = b_v3_a + m_g*s3_m + n_g*s3_n
                      if (fact2 == 0.0d0) then
                          v3(idx_v3) = real(fact1, 4) * buf_C(m_l, n_l)
                      else
                          v3(idx_v3) = real(fact2, 4) * v3(idx_v3) + real(fact1, 4) * buf_C(m_l, n_l)
                      end if
                   end do; end do
                end do; end do
            end if
         end do
!$omp end do
         deallocate(buf_A, buf_B, buf_C)
!$omp end parallel

999   continue
#ifdef USE_MKL
      call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
      call openblas_set_num_threads(old_threads)
#endif
      return
      end subroutine vv2v_dpd_442_s
