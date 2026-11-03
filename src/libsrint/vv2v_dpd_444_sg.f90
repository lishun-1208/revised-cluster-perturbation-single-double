! =======================================================================
! 终极工业级张量引擎 - GPU 单精度核心 (Single Precision, GPU + OpenACC)
! 架构：动态双擎路由 + 极限降维映射 + 双精度防截断装甲
! 特性：
!   1. 降维打击：彻底消灭内层 if！通过 Flat-to-ND 索引映射，榨干访存速度。
!   2. 精度满血：保留 KxM 连续布局与 transA='T'，统一 cuBLAS 底层加法树。
!   3. 防截断装甲：引入 buf_C_fp64，彻底吸收 Tiling 跨块累加时的单精度截断误差！
!   4. 严谨编译：全面补全所有临时调度变量声明，完美通过严格强类型检查。
!   5. 动态内存自适应：完美支持任意字符串别名路由带来的非连续索引维度提取。
! =======================================================================
      subroutine vv2v_dpd_444_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                               v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
                               v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1in, fact2in)
 
      use mod_symm
      implicit none

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type1, type2, type3 
      real*4, intent(in)            :: v1(*), v2(*)        
      real*4, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12, iop3_34
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1in, fact2in        

      ! 外部接口
      external :: xgemmsg
      integer, external :: irpdso  

      ! --- 静态架构配置 ---
      integer, parameter :: TILE_DIM = 4000 

      ! --- Fast Path 变量 ---
      logical :: is_fast_path
      character(len=1) :: trans1, trans2
      integer :: m_size, n_size, k_size, lda1, lda2
      integer :: inumpq(33), inumrs(33), inummn(33)
      integer :: iabsoff1_fp(33), iabsoff2_fp(33), iabsoff3_fp(33)
      integer :: irr_pq, irr_rs, irr_mn, irr_L1_fp, irr_L2_fp, irr_L3_fp
      integer :: iop_mn, ioff1, ioff2, ioff3
      integer*8 :: m_size8, n_size8, k_size8, lda1_8, lda2_8, lda3_8
      real*8 :: alpha8_gemm, beta8_gemm, beta8_accum
      real*8 :: d_fact1, d_fact2  

      ! --- Fallback 引擎映射与调度变量 ---
      character(len=1) :: c_ext(4), c_dum(2)
      integer :: pos1(6), pos2(6), pos3(6) 
      integer :: i, j, idx, cnt_m, cnt_n
      logical :: found, skip_block
      integer :: idx_m_in_v3(2), idx_n_in_v3(2)
      integer :: pop1(4,8), pop2(4,8), pop3(4,8)
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposlft1(33,33), iposrht1(33,33)
      integer :: inumlft2(33), inumrht2(33), iabsoff2(33), iposlft2(33,33), iposrht2(33,33)
      integer :: inumlft3(33), inumrht3(33), iabsoff3(33), iposlft3(33,33), iposrht3(33,33)
      integer :: num_tasks, task_list(3, 10000)

      ! --- Fallback GPU 线程私有计算变量 ---
      integer :: task_idx, irr_rht3, irr_lft3, irr_L4, irr_L3, irr_L2, irr_L1
      integer :: dim_L1, dim_L2, dim_L3, dim_L4, dim_X_arr(4)
      integer :: cur_M, cur_N, cur_K, cur_irr_D1
      integer :: irr_phys1(4), irr_phys2(4), irr_phys3(4), irr_logi(6)
      integer :: base_v3, s_phys(4), s_v3(6), base_v3_adj
      integer :: d_D1(8), d_D2(8), irr_D1, irr_D2, irr_lft1, irr_rht1, irr_lft2, irr_rht2
      integer :: base_v1, base_v2, s_v1_arr(6), s_v2_arr(6)
      integer :: m_dim1, m_dim2, n_dim1, n_dim2
      
      ! --- 双擎动态路由与防截断装甲变量 ---
      integer :: max_m, max_n, max_k, k_offset
      integer*8 :: global_max_m, global_max_n, global_max_k
      integer*8 :: req_mem_mb
      integer :: b_v1, b_v2, s_m1, s_m2, s_n1, s_n2, s_k1, s_k2, k_dim1, k_dim2
      integer :: s_v3_m1, s_v3_m2, s_v3_n1, s_v3_n2
      integer :: m_2, m_1, n_2, n_1, k_2, k_1, v1_idx, v2_idx, idx_v3
      integer :: m_start, n_start, k_start, m_chunk, n_chunk, k_chunk
      integer :: m_loc, n_loc, k_loc, m_flat, n_flat, k_flat
      real*4  :: fact1, fact2
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)
      real*8, allocatable :: buf_C_fp64(:,:)  

      fact1 = sngl(fact1in); fact2 = sngl(fact2in)

      ! =============================================================
      ! 【阶段 1：智能路由探测器】
      ! =============================================================
      is_fast_path = .false.
      if (type1 == '3412' .and. type2 == '3456' .and. type3 == '1256') then
          if (iop1_12 == iop2_12 .and. iop1_34 == iop3_12 .and. iop2_34 == iop3_34) then
              is_fast_path = .true.; trans1 = 'T'; trans2 = 'N'; iop_mn = iop1_12
          end if
      else if (type1 == '1234' .and. type2 == '3456' .and. type3 == '1256') then
          if (iop1_34 == iop2_12 .and. iop1_12 == iop3_12 .and. iop2_34 == iop3_34) then
              is_fast_path = .true.; trans1 = 'N'; trans2 = 'N'; iop_mn = iop1_34
          end if
      else if (type1 == '1234' .and. type2 == '5634' .and. type3 == '1256') then
          if (iop1_34 == iop2_34 .and. iop1_12 == iop3_12 .and. iop2_12 == iop3_34) then
              is_fast_path = .true.; trans1 = 'N'; trans2 = 'T'; iop_mn = iop1_34
          end if
      else if (type1 == '3412' .and. type2 == '5634' .and. type3 == '1256') then
          if (iop1_12 == iop2_34 .and. iop1_34 == iop3_12 .and. iop2_12 == iop3_34) then
              is_fast_path = .true.; trans1 = 'T'; trans2 = 'T'; iop_mn = iop1_12
          end if
      end if

      ! >>> Fast Path 极速计算分支 <<<
      if (is_fast_path) then
          inumpq = 0; inumrs = 0; inummn = 0
          do i = 1, nirrep
              inumpq(i) = irpdso(i, p3_1, p3_2, iop3_12)
              inumrs(i) = irpdso(i, p3_3, p3_4, iop3_34)
              if (trans1 == 'T') then; inummn(i) = irpdso(i, p1_1, p1_2, iop_mn)
              else; inummn(i) = irpdso(i, p1_3, p1_4, iop_mn)
              end if
          end do

          iabsoff1_fp(1) = 1; iabsoff2_fp(1) = 1; iabsoff3_fp(1) = 1
          do i = 1, nirrep - 1
              irr_L1_fp = dirprd(i, sym1)
              if (trans1 == 'T') then; iabsoff1_fp(i+1) = iabsoff1_fp(i) + inummn(irr_L1_fp) * inumpq(i)
              else; iabsoff1_fp(i+1) = iabsoff1_fp(i) + inumpq(irr_L1_fp) * inummn(i)
              end if
              irr_L2_fp = dirprd(i, sym2)
              if (trans2 == 'N') then; iabsoff2_fp(i+1) = iabsoff2_fp(i) + inummn(irr_L2_fp) * inumrs(i)
              else; iabsoff2_fp(i+1) = iabsoff2_fp(i) + inumrs(irr_L2_fp) * inummn(i)
              end if
              irr_L3_fp = dirprd(i, sym3)
              iabsoff3_fp(i+1) = iabsoff3_fp(i) + inumpq(irr_L3_fp) * inumrs(i)
          end do

          ! 严密防护：使用原始的 fact1in, fact2in 传递给底层，防止早期单精度截断
          alpha8_gemm = fact1in; beta8_gemm = fact2in
          
          do irr_rs = 1, nirrep
              irr_pq = dirprd(irr_rs, sym3); irr_mn = dirprd(irr_pq, sym1)
              if (trans1 == 'T') then; ioff1 = iabsoff1_fp(irr_pq)
              else; ioff1 = iabsoff1_fp(irr_mn); end if
              if (trans2 == 'N') then; ioff2 = iabsoff2_fp(irr_rs)
              else; ioff2 = iabsoff2_fp(irr_mn); end if
              ioff3 = iabsoff3_fp(irr_rs)
              
              m_size = inumpq(irr_pq); n_size = inumrs(irr_rs); k_size = inummn(irr_mn)
              
              ! 【核心修复：Fast Path】分离 M/N=0 与 K=0 判定，防止事实上的跳过清零
              if (m_size == 0 .or. n_size == 0) cycle
              
              if (k_size == 0) then
                  if (fact2in == 0.0d0) then
!$acc parallel loop deviceptr(v3)
                      do i = 0, m_size*n_size - 1
                          v3(ioff3 + i) = 0.0e0
                      end do
                  else if (fact2in /= 1.0d0) then
!$acc parallel loop deviceptr(v3)
                      do i = 0, m_size*n_size - 1
                          v3(ioff3 + i) = sngl(fact2in) * v3(ioff3 + i)
                      end do
                  end if
                  cycle
              end if
              
              lda1 = k_size; if (trans1 == 'N') lda1 = m_size
              lda2 = k_size; if (trans2 == 'T') lda2 = n_size
              
              m_size8 = m_size; n_size8 = n_size; k_size8 = k_size
              lda1_8  = lda1; lda2_8  = lda2
              
              call xgemmsg(trans1, trans2, m_size8, n_size8, k_size8, alpha8_gemm, &
                         v1(ioff1), lda1_8, v2(ioff2), lda2_8, beta8_gemm, v3(ioff3), m_size8)
          end do
          return
      end if

      ! =============================================================
      ! 【阶段 2：Fallback - 预调度与 DPD 表构建】
      ! =============================================================
      if (iop1_12 /= 1 .or. iop1_34 /= 1 .or. iop2_12 /= 1 .or. &
          iop2_34 /= 1 .or. iop3_12 /= 1 .or. iop3_34 /= 1) then
          write(6,*) " FATAL ERROR in vv2v_general_sg: IOP=0/2 not allowed in slow path."
          stop
      end if

      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:); pop2(3,:) = p2_3(:); pop2(4,:) = p2_4(:)
      pop3(1,:) = p3_1(:); pop3(2,:) = p3_2(:); pop3(3,:) = p3_3(:); pop3(4,:) = p3_4(:)

      c_ext(1) = type3(1:1); c_ext(2) = type3(2:2); c_ext(3) = type3(3:3); c_ext(4) = type3(4:4)
      idx = 0
      do i = 1, 4
         found = .false.
         do j = 1, 4; if (type1(i:i) == c_ext(j)) found = .true.; end do
         if (.not. found) then; idx = idx + 1; if(idx <= 2) c_dum(idx) = type1(i:i); end if
      end do

      pos1 = 0; pos2 = 0; pos3 = 0
      do i = 1, 4
         do j = 1, 4
            if (type1(j:j) == c_ext(i)) pos1(i) = j
            if (type2(j:j) == c_ext(i)) pos2(i) = j
            if (type3(j:j) == c_ext(i)) pos3(i) = j
         end do
      end do
 
      do i = 1, 2; do j = 1, 4
         if (type1(j:j) == c_dum(i)) pos1(4+i) = j
         if (type2(j:j) == c_dum(i)) pos2(4+i) = j
      end do; end do

      cnt_m = 0; cnt_n = 0; idx_m_in_v3 = 0; idx_n_in_v3 = 0
      do i = 1, 4
         if (pos1(i) /= 0) then; cnt_m = cnt_m + 1; idx_m_in_v3(cnt_m) = i; end if
         if (pos2(i) /= 0) then; cnt_n = cnt_n + 1; idx_n_in_v3(cnt_n) = i; end if
      end do

      call build_dpd_tables(pop1, sym1, inumlft1, inumrht1, iabsoff1, iposlft1, iposrht1)
      call build_dpd_tables(pop2, sym2, inumlft2, inumrht2, iabsoff2, iposlft2, iposrht2)
      call build_dpd_tables(pop3, sym3, inumlft3, inumrht3, iabsoff3, iposlft3, iposrht3)

      num_tasks = 0
      do irr_rht3 = 1, nirrep
         irr_lft3 = dirprd(irr_rht3, sym3)      
         do irr_L4 = 1, nirrep
            irr_L3 = dirprd(irr_rht3, irr_L4)     
            do irr_L2 = 1, nirrep
               irr_L1 = dirprd(irr_lft3, irr_L2)    
               dim_L1 = pop3(pos3(1), irr_L1); dim_L2 = pop3(pos3(2), irr_L2)
               dim_L3 = pop3(pos3(3), irr_L3); dim_L4 = pop3(pos3(4), irr_L4)
               if (dim_L1 == 0 .or. dim_L2 == 0 .or. dim_L3 == 0 .or. dim_L4 == 0) cycle
               num_tasks = num_tasks + 1
               task_list(1, num_tasks) = irr_rht3
               task_list(2, num_tasks) = irr_L4
               task_list(3, num_tasks) = irr_L2
            end do
         end do
      end do
      if (num_tasks == 0) return

      ! =============================================================
      ! 【阶段 3：全局极速预扫描，探测显存需求】
      ! =============================================================
      global_max_m = 0; global_max_n = 0; global_max_k = 0
      do task_idx = 1, num_tasks
         irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
         irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)
         
         dim_X_arr(1) = pop3(pos3(1), irr_L1)
         dim_X_arr(2) = pop3(pos3(2), irr_L2)
         dim_X_arr(3) = pop3(pos3(3), irr_L3)
         dim_X_arr(4) = pop3(pos3(4), irr_L4)
         
         cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2))
         cur_N = dim_X_arr(idx_n_in_v3(1)) * dim_X_arr(idx_n_in_v3(2))
         
         if (cur_M > global_max_m) global_max_m = cur_M
         if (cur_N > global_max_n) global_max_n = cur_N
         max_k = 0
         
         do irr_D2 = 1, nirrep
            irr_D1 = dirprd(sym1, irr_D2)
            do i = 1, 4; if (pos1(i)/=0) then         
               if (i==1) irr_D1 = dirprd(irr_D1, irr_L1); if (i==2) irr_D1 = dirprd(irr_D1, irr_L2)
               if (i==3) irr_D1 = dirprd(irr_D1, irr_L3); if (i==4) irr_D1 = dirprd(irr_D1, irr_L4)
            end if; end do
            irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /); irr_logi(5:6) = (/ irr_D1, irr_D2 /)
            skip_block = .false.
            do i = 1, 6
               if (pos1(i)/=0 .and. pop1(pos1(i), irr_logi(i))==0) skip_block = .true.
               if (pos2(i)/=0 .and. pop2(pos2(i), irr_logi(i))==0) skip_block = .true.
            end do
            if (skip_block) cycle
            max_k = max_k + (pop1(pos1(5), irr_D1) * pop1(pos1(6), irr_D2))
         end do
         if (max_k > global_max_k) global_max_k = max_k
      end do

      if (global_max_m == 0 .or. global_max_n == 0 .or. global_max_k == 0) return

      req_mem_mb = (global_max_m * global_max_k + &
                    global_max_k * global_max_n + &
                    global_max_m * global_max_n) * 4_8 / (1024_8 * 1024_8)

      d_fact1 = fact1in
      d_fact2 = fact2in

      ! =============================================================
      ! 【阶段 4：动态双擎路由执行】
      ! =============================================================
      if (req_mem_mb <= 500_8) then
          ! *************************************************************
          ! [引擎 A：狂飙模式] 显存充足，单次全量 GEMM
          ! *************************************************************
          allocate(buf_A(global_max_k, global_max_m), buf_B(global_max_k, global_max_n), buf_C(global_max_m, global_max_n))
!$acc data create(buf_A, buf_B, buf_C)

          do task_idx = 1, num_tasks
             irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
             irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)    
             irr_phys3(pos3(1):pos3(4)) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /)
             dim_X_arr(1:4) = (/ pop3(pos3(1), irr_L1), pop3(pos3(2), irr_L2), pop3(pos3(3), irr_L3), pop3(pos3(4), irr_L4) /)
             cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2))
             cur_N = dim_X_arr(idx_n_in_v3(1)) * dim_X_arr(idx_n_in_v3(2))
           
             base_v3 = iabsoff3(irr_rht3) + iposrht3(irr_phys3(4), irr_rht3) + iposlft3(irr_phys3(2), irr_lft3)
             s_phys(1) = 1; s_phys(2) = pop3(1, irr_phys3(1)); s_phys(3) = inumlft3(irr_lft3); s_phys(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))
             do i = 1, 4; s_v3(i) = s_phys(pos3(i)); end do
             base_v3_adj = base_v3 - sum(s_v3(1:4)) + 1
             m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2))
             n_dim1 = dim_X_arr(idx_n_in_v3(1)); n_dim2 = dim_X_arr(idx_n_in_v3(2))
             max_k = 0; k_offset = 0

             do irr_D2 = 1, nirrep
                irr_D1 = dirprd(sym1, irr_D2)
                do i = 1, 4; if(pos1(i)/=0) then
                   if(i==1) irr_D1=dirprd(irr_D1,irr_L1); if(i==2) irr_D1=dirprd(irr_D1,irr_L2)
                   if(i==3) irr_D1=dirprd(irr_D1,irr_L3); if(i==4) irr_D1=dirprd(irr_D1,irr_L4)
                end if; end do
                irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /); irr_logi(5:6) = (/ irr_D1, irr_D2 /)
                skip_block = .false.
                do i = 1, 6; if(pos1(i)/=0.and.pop1(pos1(i),irr_logi(i))==0) skip_block=.true.
                             if(pos2(i)/=0.and.pop2(pos2(i),irr_logi(i))==0) skip_block=.true.; end do
                if (skip_block) cycle
                
                k_dim1 = pop1(pos1(5), irr_D1); k_dim2 = pop1(pos1(6), irr_D2); cur_K = k_dim1 * k_dim2
                if (cur_K == 0) cycle

                do i = 1, 6; if(pos1(i)/=0) irr_phys1(pos1(i))=irr_logi(i); if(pos2(i)/=0) irr_phys2(pos2(i))=irr_logi(i); end do
                irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2)); irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                irr_lft2 = dirprd(irr_phys2(1), irr_phys2(2)); irr_rht2 = dirprd(irr_phys2(3), irr_phys2(4))
                base_v1 = iabsoff1(irr_rht1) + iposrht1(irr_phys1(4), irr_rht1) + iposlft1(irr_phys1(2), irr_lft1)
                base_v2 = iabsoff2(irr_rht2) + iposrht2(irr_phys2(4), irr_rht2) + iposlft2(irr_phys2(2), irr_lft2)
                s_phys(1)=1; s_phys(2)=pop1(1,irr_phys1(1)); s_phys(3)=inumlft1(irr_lft1); s_phys(4)=inumlft1(irr_lft1)*pop1(3,irr_phys1(3))
                s_v1_arr(:)=0; do i=1,6; if(pos1(i)/=0) s_v1_arr(i)=s_phys(pos1(i)); end do
                s_phys(1)=1; s_phys(2)=pop2(1,irr_phys2(1)); s_phys(3)=inumlft2(irr_lft2); s_phys(4)=inumlft2(irr_lft2)*pop2(3,irr_phys2(3))
                s_v2_arr(:)=0; do i=1,6; if(pos2(i)/=0) s_v2_arr(i)=s_phys(pos2(i)); end do
                b_v1 = base_v1 - sum(s_v1_arr) + 1; b_v2 = base_v2 - sum(s_v2_arr) + 1
                s_m1 = s_v1_arr(idx_m_in_v3(1)); s_m2 = s_v1_arr(idx_m_in_v3(2))
                s_n1 = s_v2_arr(idx_n_in_v3(1)); s_n2 = s_v2_arr(idx_n_in_v3(2))
                s_k1 = s_v1_arr(5); s_k2 = s_v1_arr(6)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, m_1, m_2, k_1, k_2, v1_idx)
                do k_loc = 1, cur_K
                   do m_loc = 1, cur_M
                      m_flat = m_loc
                      k_flat = k_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      m_2 = (m_flat - 1) / m_dim1 + 1
                      k_1 = mod(k_flat - 1, k_dim1) + 1
                      k_2 = (k_flat - 1) / k_dim1 + 1
                      v1_idx = b_v1 + m_2*s_m2 + m_1*s_m1 + k_1*s_k1 + k_2*s_k2
                      buf_A(k_offset + k_loc, m_loc) = v1(v1_idx)
                   end do
                end do

                s_k1 = s_v2_arr(5); s_k2 = s_v2_arr(6)

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, n_1, n_2, k_1, k_2, v2_idx)
                do n_loc = 1, cur_N
                   do k_loc = 1, cur_K
                      n_flat = n_loc
                      k_flat = k_loc
                      n_1 = mod(n_flat - 1, n_dim1) + 1
                      n_2 = (n_flat - 1) / n_dim1 + 1
                      k_1 = mod(k_flat - 1, k_dim1) + 1
                      k_2 = (k_flat - 1) / k_dim1 + 1
                      v2_idx = b_v2 + n_2*s_n2 + n_1*s_n1 + k_1*s_k1 + k_2*s_k2
                      buf_B(k_offset + k_loc, n_loc) = v2(v2_idx)
                   end do
                end do

                k_offset = k_offset + cur_K
             end do
             max_k = k_offset
             
             ! 【核心修复：Fallback Engine A】 当 K=0 (max_k=0) 且外维不为0时，单独为 v3 清零或缩放！
             if (max_k == 0) then
                 s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2))
                 s_v3_n1 = s_v3(idx_n_in_v3(1)); s_v3_n2 = s_v3(idx_n_in_v3(2))
                 if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          m_flat = m_loc; n_flat = n_loc
                          m_1 = mod(m_flat - 1, m_dim1) + 1; m_2 = (m_flat - 1) / m_dim1 + 1
                          n_1 = mod(n_flat - 1, n_dim1) + 1; n_2 = (n_flat - 1) / n_dim1 + 1
                          idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                          v3(idx_v3) = 0.0e0
                       end do
                    end do
                 else if (fact2in /= 1.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          m_flat = m_loc; n_flat = n_loc
                          m_1 = mod(m_flat - 1, m_dim1) + 1; m_2 = (m_flat - 1) / m_dim1 + 1
                          n_1 = mod(n_flat - 1, n_dim1) + 1; n_2 = (n_flat - 1) / n_dim1 + 1
                          idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                          v3(idx_v3) = sngl(fact2in) * v3(idx_v3)
                       end do
                    end do
                 end if
                 cycle
             end if
             
!$acc wait

             m_size8 = cur_M; n_size8 = cur_N; k_size8 = max_k
             lda1_8 = global_max_k; lda2_8 = global_max_k; lda3_8 = global_max_m
             alpha8_gemm = 1.0d0; beta8_accum = 0.0d0 

!$acc host_data use_device(buf_A, buf_B, buf_C)
             call xgemmsg('T', 'N', m_size8, n_size8, k_size8, alpha8_gemm, &
                          buf_A, lda1_8, buf_B, lda2_8, beta8_accum, buf_C, lda3_8)
!$acc end host_data

             s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2))
             s_v3_n1 = s_v3(idx_n_in_v3(1)); s_v3_n2 = s_v3(idx_n_in_v3(2))

             ! 严格判定传入 fact2in 为 0.0d0
             if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      m_flat = m_loc
                      n_flat = n_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      m_2 = (m_flat - 1) / m_dim1 + 1
                      n_1 = mod(n_flat - 1, n_dim1) + 1
                      n_2 = (n_flat - 1) / n_dim1 + 1
                      idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                      v3(idx_v3) = sngl(d_fact1 * dble(buf_C(m_loc, n_loc)))
                   end do
                end do
             else
!$acc parallel loop collapse(2) present(buf_C) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      m_flat = m_loc
                      n_flat = n_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      m_2 = (m_flat - 1) / m_dim1 + 1
                      n_1 = mod(n_flat - 1, n_dim1) + 1
                      n_2 = (n_flat - 1) / n_dim1 + 1
                      idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                      v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * dble(buf_C(m_loc, n_loc)))
                   end do
                end do
             end if
          end do
!$acc end data
          deallocate(buf_A, buf_B, buf_C)

      else
          ! *************************************************************
          ! [引擎 B：生存模式] (超大内存需求 C1 体系)
          ! *************************************************************
          allocate(buf_A(TILE_DIM, TILE_DIM), buf_B(TILE_DIM, TILE_DIM), buf_C(TILE_DIM, TILE_DIM), buf_C_fp64(TILE_DIM, TILE_DIM))
!$acc data create(buf_A, buf_B, buf_C, buf_C_fp64)

          do task_idx = 1, num_tasks
             irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
             irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)    
             irr_phys3(pos3(1):pos3(4)) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /)
             dim_X_arr(1:4) = (/ pop3(pos3(1), irr_L1), pop3(pos3(2), irr_L2), pop3(pos3(3), irr_L3), pop3(pos3(4), irr_L4) /)
             cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2))
             cur_N = dim_X_arr(idx_n_in_v3(1)) * dim_X_arr(idx_n_in_v3(2))
           
             base_v3 = iabsoff3(irr_rht3) + iposrht3(irr_phys3(4), irr_rht3) + iposlft3(irr_phys3(2), irr_lft3)
             s_phys(1) = 1; s_phys(2) = pop3(1, irr_phys3(1)); s_phys(3) = inumlft3(irr_lft3); s_phys(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))
             do i = 1, 4; s_v3(i) = s_phys(pos3(i)); end do
             base_v3_adj = base_v3 - sum(s_v3(1:4)) + 1
             m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2))
             n_dim1 = dim_X_arr(idx_n_in_v3(1)); n_dim2 = dim_X_arr(idx_n_in_v3(2))

             do m_start = 1, cur_M, TILE_DIM
                m_chunk = min(TILE_DIM, cur_M - m_start + 1)
                do n_start = 1, cur_N, TILE_DIM
                   n_chunk = min(TILE_DIM, cur_N - n_start + 1)

!$acc kernels present(buf_C_fp64)
                   buf_C_fp64 = 0.0d0
!$acc end kernels

                   do irr_D2 = 1, nirrep
                      irr_D1 = dirprd(sym1, irr_D2)
                      do i = 1, 4; if(pos1(i)/=0) then
                         if(i==1) irr_D1=dirprd(irr_D1,irr_L1); if(i==2) irr_D1=dirprd(irr_D1,irr_L2)
                         if(i==3) irr_D1=dirprd(irr_D1,irr_L3); if(i==4) irr_D1=dirprd(irr_D1,irr_L4)
                      end if; end do
                      irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /); irr_logi(5:6) = (/ irr_D1, irr_D2 /)
                      skip_block = .false.
                      do i = 1, 6; if(pos1(i)/=0.and.pop1(pos1(i),irr_logi(i))==0) skip_block=.true.
                                   if(pos2(i)/=0.and.pop2(pos2(i),irr_logi(i))==0) skip_block=.true.; end do
                      if (skip_block) cycle
                      
                      k_dim1 = pop1(pos1(5), irr_D1); k_dim2 = pop1(pos1(6), irr_D2); cur_K = k_dim1 * k_dim2
                      if (cur_K == 0) cycle

                      do i = 1, 6; if(pos1(i)/=0) irr_phys1(pos1(i))=irr_logi(i); if(pos2(i)/=0) irr_phys2(pos2(i))=irr_logi(i); end do
                      irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2)); irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                      irr_lft2 = dirprd(irr_phys2(1), irr_phys2(2)); irr_rht2 = dirprd(irr_phys2(3), irr_phys2(4))
                      base_v1 = iabsoff1(irr_rht1) + iposrht1(irr_phys1(4), irr_rht1) + iposlft1(irr_phys1(2), irr_lft1)
                      base_v2 = iabsoff2(irr_rht2) + iposrht2(irr_phys2(4), irr_rht2) + iposlft2(irr_phys2(2), irr_lft2)
                      s_phys(1)=1; s_phys(2)=pop1(1,irr_phys1(1)); s_phys(3)=inumlft1(irr_lft1); s_phys(4)=inumlft1(irr_lft1)*pop1(3,irr_phys1(3))
                      s_v1_arr(:)=0; do i=1,6; if(pos1(i)/=0) s_v1_arr(i)=s_phys(pos1(i)); end do
                      s_phys(1)=1; s_phys(2)=pop2(1,irr_phys2(1)); s_phys(3)=inumlft2(irr_lft2); s_phys(4)=inumlft2(irr_lft2)*pop2(3,irr_phys2(3))
                      s_v2_arr(:)=0; do i=1,6; if(pos2(i)/=0) s_v2_arr(i)=s_phys(pos2(i)); end do
                      b_v1 = base_v1 - sum(s_v1_arr) + 1; b_v2 = base_v2 - sum(s_v2_arr) + 1
                      s_m1 = s_v1_arr(idx_m_in_v3(1)); s_m2 = s_v1_arr(idx_m_in_v3(2))
                      s_n1 = s_v2_arr(idx_n_in_v3(1)); s_n2 = s_v2_arr(idx_n_in_v3(2))

                      do k_start = 1, cur_K, TILE_DIM
                         k_chunk = min(TILE_DIM, cur_K - k_start + 1)
                         s_k1 = s_v1_arr(5); s_k2 = s_v1_arr(6)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, m_1, m_2, k_1, k_2, v1_idx)
                         do k_loc = 1, k_chunk
                            do m_loc = 1, m_chunk
                               m_flat = m_start + m_loc - 1
                               k_flat = k_start + k_loc - 1
                               m_1 = mod(m_flat - 1, m_dim1) + 1
                               m_2 = (m_flat - 1) / m_dim1 + 1
                               k_1 = mod(k_flat - 1, k_dim1) + 1
                               k_2 = (k_flat - 1) / k_dim1 + 1
                               v1_idx = b_v1 + m_2*s_m2 + m_1*s_m1 + k_1*s_k1 + k_2*s_k2
                               buf_A(k_loc, m_loc) = v1(v1_idx)
                            end do
                         end do

                         s_k1 = s_v2_arr(5); s_k2 = s_v2_arr(6)
!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, n_1, n_2, k_1, k_2, v2_idx)
                         do n_loc = 1, n_chunk
                            do k_loc = 1, k_chunk
                               n_flat = n_start + n_loc - 1
                               k_flat = k_start + k_loc - 1
                               n_1 = mod(n_flat - 1, n_dim1) + 1
                               n_2 = (n_flat - 1) / n_dim1 + 1
                               k_1 = mod(k_flat - 1, k_dim1) + 1
                               k_2 = (k_flat - 1) / k_dim1 + 1
                               v2_idx = b_v2 + n_2*s_n2 + n_1*s_n1 + k_1*s_k1 + k_2*s_k2
                               buf_B(k_loc, n_loc) = v2(v2_idx)
                            end do
                         end do
!$acc wait

                         m_size8 = m_chunk; n_size8 = n_chunk; k_size8 = k_chunk
                         lda1_8 = TILE_DIM; lda2_8 = TILE_DIM; lda3_8 = TILE_DIM
                         alpha8_gemm = 1.0d0; beta8_accum = 0.0d0 
!$acc host_data use_device(buf_A, buf_B, buf_C)
                         call xgemmsg('T', 'N', m_size8, n_size8, k_size8, alpha8_gemm, &
                                      buf_A, lda1_8, buf_B, lda2_8, beta8_accum, buf_C, lda3_8)
!$acc end host_data

!$acc parallel loop collapse(2) present(buf_C, buf_C_fp64)
                         do n_loc = 1, n_chunk
                            do m_loc = 1, m_chunk
                               buf_C_fp64(m_loc, n_loc) = buf_C_fp64(m_loc, n_loc) + dble(buf_C(m_loc, n_loc))
                            end do
                         end do

                      end do ! end k_start
                   end do ! end irr_D2

                   s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2))
                   s_v3_n1 = s_v3(idx_n_in_v3(1)); s_v3_n2 = s_v3(idx_n_in_v3(2))

                   if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1
                            n_flat = n_start + n_loc - 1
                            m_1 = mod(m_flat - 1, m_dim1) + 1
                            m_2 = (m_flat - 1) / m_dim1 + 1
                            n_1 = mod(n_flat - 1, n_dim1) + 1
                            n_2 = (n_flat - 1) / n_dim1 + 1
                            idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                            v3(idx_v3) = sngl(d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   else
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, n_1, n_2, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1
                            n_flat = n_start + n_loc - 1
                            m_1 = mod(m_flat - 1, m_dim1) + 1
                            m_2 = (m_flat - 1) / m_dim1 + 1
                            n_1 = mod(n_flat - 1, n_dim1) + 1
                            n_2 = (n_flat - 1) / n_dim1 + 1
                            idx_v3 = base_v3_adj + n_1*s_v3_n1 + n_2*s_v3_n2 + m_1*s_v3_m1 + m_2*s_v3_m2
                            v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   end if

                end do ! end n_start
             end do ! end m_start
          end do ! end task_idx

!$acc end data
          deallocate(buf_A, buf_B, buf_C, buf_C_fp64)
      end if

      return
      end subroutine vv2v_dpd_444_sg
