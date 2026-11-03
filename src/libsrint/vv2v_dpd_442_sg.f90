! =======================================================================
! 终极工业级张量引擎 - GPU 单精度核心 (Single Precision, GPU + OpenACC)
! [核心能力] 442 模式: 4-index * 4-index -> 2-index 
! [适用场景] 3指标收缩运算 (例如: pqrm * pqrn -> mn)
! [架构设计] 1. 彻底消灭竞态！按 V3 分块正交隔离任务，防止 irr_R 子块结果互相覆盖。
!            2. K 维级联拼接：在 Engine A 模式下，将所有 irr_R 在显存中横向拼接，
!               实现 1 次超大 GEMM，彻底榨干 GPU 算力。
!            3. FP64装甲：Engine B 的 Tiling 模式下使用 real*8 瓦片无损累加。
!            4. 动态双擎：单精度按 4Bytes 评估内存，小数据走狂飙，大数据走 Tiling。
! =======================================================================
      subroutine vv2v_dpd_442_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                                 v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
                                 v3, p3_1, p3_2, iop3_12, sym3, type3, fact1in, fact2in)
      
      use mod_symm
      implicit none

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type1, type2 
      character(len=2), intent(in)  :: type3        
      real*4, intent(in)            :: v1(*), v2(*)        
      real*4, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1in, fact2in        

      ! 外部接口
      external :: xgemmsg
      integer, external :: irpdso  

      ! --- 静态架构配置 ---
      integer, parameter :: TILE_DIM = 4000 

      ! --- 路由解析变量 ---
      character(len=1) :: c_con(3), c_r, c_m, c_n
      integer :: pos1_m, pos2_n, pos1_k(3), pos2_k(3), v3_map_m, v3_map_n
      integer :: pos1_r, pos2_r
      integer :: i, j, idx
      logical :: is_macro_LL, is_macro_RR, found_r
      integer :: pop1(4,8), pop2(4,8), pop3(2,8)
      
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposlft1_fb(33,33), iposrht1(33,33)
      integer :: inumlft2(33), inumrht2(33), iabsoff2(33), iposlft2_fb(33,33), iposrht2(33,33)
      integer :: iabsoff3(33)
      
      integer :: num_tasks, task_list(10000)

      ! --- GPU 调度私有变量 ---
      integer :: t_idx, irr_M, irr_N, irr_R, irr_PQ, irr_l1, irr_r1, irr_l2, irr_r2, irr_v3_2, irr_k1, irr_k2, irr_k3
      integer :: cur_M, cur_N, cur_K, cur_K_block, dim_PQ, dim_M, dim_N, dim_R, dim_k(3), irr_logi1(4), irr_logi2(4)
      integer :: m_start, n_start, k_start, m_chunk, n_chunk, k_chunk
      integer :: m_loc, n_loc, k_loc, m_flat, n_flat, k_flat, i_pq, i_r, k_v(3), rem
      integer :: base_v1, base_v2, base_v3, b_v1_a, b_v2_a, b_v3_a, s1(4), s2(4), s3_m, s3_n, s_phys1(4), s_phys2(4)
      integer :: irr_pos_m, irr_pos_n
      
      ! --- 双擎动态路由变量 ---
      integer :: k_offset
      integer*8 :: global_max_m, global_max_n, global_max_k, cur_K_sum
      integer*8 :: req_mem_mb, v1_idx, v2_idx, idx_v3
      integer*8 :: m_size8, n_size8, k_size8, lda1_8, lda2_8, lda3_8
      real*8 :: alpha8_gemm, beta8_accum
      real*8 :: d_fact1, d_fact2
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)
      real*8, allocatable :: buf_C_fp64(:,:)

      d_fact1 = fact1in
      d_fact2 = fact2in

      ! =============================================================
      ! 【阶段 1：宏维 (Macro-Index) 探测逻辑】
      ! =============================================================
      is_macro_LL = .false.; is_macro_RR = .false.
      c_r = ' '

      ! 检查左侧宏维: pqrm * pqrn = mn
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
      
      ! 检查右侧宏维: rmpq * rnpq = mn
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
      end if

      ! --- IOP 拦截墙 ---
      if (is_macro_LL) then
          if (iop1_34 /= 1 .or. iop2_34 /= 1 .or. iop3_12 /= 1) stop " FATAL: 442 Macro-LL must have IOP=1 for split indices."
      else if (is_macro_RR) then
          if (iop1_12 /= 1 .or. iop2_12 /= 1 .or. iop3_12 /= 1) stop " FATAL: 442 Macro-RR must have IOP=1 for split indices."
      else
          if (iop1_12/=1 .or. iop1_34/=1 .or. iop2_12/=1 .or. iop2_34/=1 .or. iop3_12/=1) stop " FATAL: Pure 442 strictly requires all IOP=1."
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
      ! 【阶段 3：全局极速预扫描，探测显存需求】
      ! =============================================================
      num_tasks = 0; global_max_m = 0; global_max_n = 0; global_max_k = 0
      
      if (is_macro_LL .or. is_macro_RR) then
          do i = 1, nirrep ! irr_M
             j = dirprd(sym3, i) ! irr_N
             if (pop1(pos1_m, i)==0 .or. pop2(pos2_n, j)==0) cycle
             
             found_r = .false.
             cur_K_sum = 0_8
             do idx = 1, nirrep ! irr_R
                if (pop1(pos1_r, idx)==0) cycle
                if (is_macro_LL .and. inumlft1(dirprd(sym1, dirprd(idx, i))) == 0) cycle
                if (is_macro_RR .and. inumrht1(dirprd(sym1, dirprd(idx, i))) == 0) cycle
                found_r = .true.
                if (is_macro_LL) cur_K_sum = cur_K_sum + int(inumlft1(dirprd(sym1,dirprd(idx,i))),8)*pop1(pos1_r,idx)
                if (is_macro_RR) cur_K_sum = cur_K_sum + int(inumrht1(dirprd(sym1,dirprd(idx,i))),8)*pop1(pos1_r,idx)
             end do
             
             if (found_r) then
                num_tasks = num_tasks + 1; task_list(num_tasks) = i ! 仅记录 irr_M，保证写入互斥
                global_max_m = max(global_max_m, int(pop1(pos1_m, i), 8))
                global_max_n = max(global_max_n, int(pop2(pos2_n, j), 8))
                global_max_k = max(global_max_k, cur_K_sum) ! 获取拼接后的全量 K
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
             
             cur_K_sum = 0_8
             do irr_k1=1,nirrep; do irr_k2=1,nirrep; irr_k3=dirprd(sym1,dirprd(irr_m,dirprd(irr_k1,irr_k2)))
                irr_logi2(pos2_n)=irr_n; irr_logi2(pos2_k(1))=irr_k1; irr_logi2(pos2_k(2))=irr_k2; irr_logi2(pos2_k(3))=irr_k3
                if(dirprd(dirprd(irr_logi2(1),irr_logi2(2)),dirprd(irr_logi2(3),irr_logi2(4)))/=sym2) cycle
                if(pop1(pos1_k(1),irr_k1)==0 .or. pop1(pos1_k(2),irr_k2)==0 .or. pop1(pos1_k(3),irr_k3)==0) cycle
                if(pop2(pos2_k(1),irr_k1)==0 .or. pop2(pos2_k(2),irr_k2)==0 .or. pop2(pos2_k(3),irr_k3)==0) cycle
                cur_K_sum = cur_K_sum + int(pop1(pos1_k(1),irr_k1),8)*pop1(pos1_k(2),irr_k2)*pop1(pos1_k(3),irr_k3)
             end do; end do
             global_max_k = max(global_max_k, cur_K_sum)
          end do
      end if
      if (num_tasks == 0 .or. global_max_m == 0 .or. global_max_n == 0) return

      ! 单精度体系按照 4 Bytes 计算
      req_mem_mb = (global_max_m * global_max_k + global_max_k * global_max_n + global_max_m * global_max_n) * 4_8 / (1024_8 * 1024_8)

      ! =============================================================
      ! 【阶段 4：动态双擎路由执行】
      ! =============================================================
      if (req_mem_mb <= 500_8) then
          ! *************************************************************
          ! [引擎 A：狂飙模式] 显存充足，支持 K 维级联的单次超大单精度 GEMM
          ! *************************************************************
          allocate(buf_A(global_max_k, global_max_m), buf_B(global_max_k, global_max_n), buf_C(global_max_m, global_max_n))
!$acc data create(buf_A, buf_B, buf_C)

          do t_idx = 1, num_tasks
             if (is_macro_LL .or. is_macro_RR) then
                irr_M = task_list(t_idx)
                if (type3(1:1)==c_m) then; irr_N = dirprd(sym3, irr_M); base_v3=iabsoff3(irr_N)
                else; irr_N = irr_M; irr_M = dirprd(sym3, irr_N); base_v3=iabsoff3(irr_M); end if
                
                dim_M = pop1(pos1_m, irr_M); dim_N = pop2(pos2_n, irr_N)
                if (type3(1:1)==c_m) then; s3_m=1; s3_n=dim_M; else; s3_n=1; s3_m=dim_N; end if
                b_v3_a = base_v3 - s3_m - s3_n
                cur_M = dim_M; cur_N = dim_N
                
                k_offset = 0
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
                   else
                      irr_PQ = dirprd(sym1, dirprd(irr_R, irr_M)); dim_PQ = inumrht1(irr_PQ)
                      if (dim_PQ == 0) cycle
                      irr_l1=dirprd(irr_R,irr_M); irr_l2=dirprd(irr_R,irr_N)
                      irr_pos_m = irr_M; if(pos1_r==2) irr_pos_m=irr_R
                      base_v1 = iabsoff1(irr_PQ) + iposlft1_fb(irr_pos_m, irr_l1); s1(1)=inumlft1(irr_l1)
                      if(pos1_r==1) then; s1(2)=1; s1(3)=dim_R; else; s1(3)=1; s1(2)=dim_M; end if
                      irr_pos_n = irr_N; if(pos2_r==2) irr_pos_n=irr_R
                      base_v2 = iabsoff2(irr_PQ) + iposlft2_fb(irr_pos_n, irr_l2); s2(1)=inumlft2(irr_l2)
                      if(pos2_r==1) then; s2(2)=1; s2(3)=dim_R; else; s2(3)=1; s2(2)=dim_N; end if
                   end if
                   
                   b_v1_a = base_v1-s1(1)-s1(2)-s1(3); b_v2_a = base_v2-s2(1)-s2(2)-s2(3)
                   cur_K_block = dim_PQ * dim_R

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(k_loc, m_loc, i_pq, i_r, v1_idx)
                   do k_loc = 1, cur_K_block
                      do m_loc = 1, cur_M
                         i_pq = mod(k_loc - 1, dim_PQ) + 1; i_r = (k_loc - 1) / dim_PQ + 1
                         v1_idx = int(b_v1_a, 8) + int(i_pq, 8)*s1(1) + int(i_r, 8)*s1(2) + int(m_loc, 8)*s1(3)
                         buf_A(k_offset + k_loc, m_loc) = v1(v1_idx)
                      end do
                   end do

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(k_loc, n_loc, i_pq, i_r, v2_idx)
                   do n_loc = 1, cur_N
                      do k_loc = 1, cur_K_block
                         i_pq = mod(k_loc - 1, dim_PQ) + 1; i_r = (k_loc - 1) / dim_PQ + 1
                         v2_idx = int(b_v2_a, 8) + int(i_pq, 8)*s2(1) + int(i_r, 8)*s2(2) + int(n_loc, 8)*s2(3)
                         buf_B(k_offset + k_loc, n_loc) = v2(v2_idx)
                      end do
                   end do
                   k_offset = k_offset + cur_K_block
                end do
                cur_K = k_offset
             else
                ! --- 标准分支 C ---
                irr_v3_2 = task_list(t_idx)
                if (v3_map_m == 2) then; irr_m = irr_v3_2; irr_n = dirprd(irr_v3_2, sym3)
                else; irr_n = irr_v3_2; irr_m = dirprd(irr_v3_2, sym3); end if
                
                cur_M = pop1(pos1_m, irr_m); cur_N = pop2(pos2_n, irr_n)
                base_v3 = iabsoff3(irr_v3_2)
                if (v3_map_m == 1) then; s3_m = 1; s3_n = cur_M; else; s3_m = cur_N; s3_n = 1; end if
                b_v3_a = base_v3 - s3_m - s3_n
                
                k_offset = 0
                do irr_k1=1,nirrep; do irr_k2=1,nirrep; irr_k3=dirprd(sym1,dirprd(irr_m,dirprd(irr_k1,irr_k2)))
                   irr_logi1(pos1_m)=irr_m; irr_logi1(pos1_k(1))=irr_k1; irr_logi1(pos1_k(2))=irr_k2; irr_logi1(pos1_k(3))=irr_k3
                   irr_logi2(pos2_n)=irr_n; irr_logi2(pos2_k(1))=irr_k1; irr_logi2(pos2_k(2))=irr_k2; irr_logi2(pos2_k(3))=irr_k3
                   if(dirprd(dirprd(irr_logi2(1),irr_logi2(2)),dirprd(irr_logi2(3),irr_logi2(4)))/=sym2) cycle
                   
                   dim_k(1)=pop1(pos1_k(1),irr_k1); dim_k(2)=pop1(pos1_k(2),irr_k2); dim_k(3)=pop1(pos1_k(3),irr_k3)
                   if(dim_k(1)==0 .or. dim_k(2)==0 .or. dim_k(3)==0) cycle
                   if(pop2(pos2_k(1),irr_k1)==0 .or. pop2(pos2_k(2),irr_k2)==0 .or. pop2(pos2_k(3),irr_k3)==0) cycle
                   cur_K_block = dim_k(1)*dim_k(2)*dim_k(3)
                   
                   irr_r1=dirprd(irr_logi1(3),irr_logi1(4)); irr_l1=dirprd(irr_logi1(1),irr_logi1(2))
                   irr_r2=dirprd(irr_logi2(3),irr_logi2(4)); irr_l2=dirprd(irr_logi2(1),irr_logi2(2))
                   
                   base_v1 = iabsoff1(irr_r1) + inumlft1(irr_l1)*iposrht1(irr_logi1(4),irr_r1) + iposlft1_fb(irr_logi1(2),irr_l1)
                   s_phys1(1)=1; s_phys1(2)=pop1(1,irr_logi1(1)); s_phys1(3)=inumlft1(irr_l1); s_phys1(4)=inumlft1(irr_l1)*pop1(3,irr_logi1(3))
                   
                   base_v2 = iabsoff2(irr_r2) + inumlft2(irr_l2)*iposrht2(irr_logi2(4),irr_r2) + iposlft2_fb(irr_logi2(2),irr_l2)
                   s_phys2(1)=1; s_phys2(2)=pop2(1,irr_logi2(1)); s_phys2(3)=inumlft2(irr_l2); s_phys2(4)=inumlft2(irr_l2)*pop2(3,irr_logi2(3))
                   
                   b_v1_a = base_v1-s_phys1(pos1_m)-s_phys1(pos1_k(1))-s_phys1(pos1_k(2))-s_phys1(pos1_k(3))
                   b_v2_a = base_v2-s_phys2(pos2_n)-s_phys2(pos2_k(1))-s_phys2(pos2_k(2))-s_phys2(pos2_k(3))

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(k_loc, m_loc, k_v, rem, v1_idx)
                   do k_loc = 1, cur_K_block
                      do m_loc = 1, cur_M
                         k_v(1)=mod(k_loc-1, dim_k(1))+1; rem=(k_loc-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                         v1_idx = int(b_v1_a, 8) + int(m_loc, 8)*s_phys1(pos1_m) + int(k_v(1), 8)*s_phys1(pos1_k(1)) + int(k_v(2), 8)*s_phys1(pos1_k(2)) + int(k_v(3), 8)*s_phys1(pos1_k(3))
                         buf_A(k_offset + k_loc, m_loc) = v1(v1_idx)
                      end do
                   end do

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(k_loc, n_loc, k_v, rem, v2_idx)
                   do n_loc = 1, cur_N
                      do k_loc = 1, cur_K_block
                         k_v(1)=mod(k_loc-1, dim_k(1))+1; rem=(k_loc-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                         v2_idx = int(b_v2_a, 8) + int(n_loc, 8)*s_phys2(pos2_n) + int(k_v(1), 8)*s_phys2(pos2_k(1)) + int(k_v(2), 8)*s_phys2(pos2_k(2)) + int(k_v(3), 8)*s_phys2(pos2_k(3))
                         buf_B(k_offset + k_loc, n_loc) = v2(v2_idx)
                      end do
                   end do
                   k_offset = k_offset + cur_K_block
                end do; end do
                cur_K = k_offset
             end if
             
             if (cur_K == 0) then
                 if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_loc, n_loc, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          idx_v3 = int(b_v3_a, 8) + int(m_loc, 8)*s3_m + int(n_loc, 8)*s3_n
                          v3(idx_v3) = 0.0e0
                       end do
                    end do
                 else if (fact2in /= 1.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_loc, n_loc, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          idx_v3 = int(b_v3_a, 8) + int(m_loc, 8)*s3_m + int(n_loc, 8)*s3_n
                          v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)))
                       end do
                    end do
                 end if
                 cycle
             end if
             
!$acc wait

             m_size8 = cur_M; n_size8 = cur_N; k_size8 = cur_K
             lda1_8 = global_max_k; lda2_8 = global_max_k; lda3_8 = global_max_m
             alpha8_gemm = 1.0d0; beta8_accum = 0.0d0 

!$acc host_data use_device(buf_A, buf_B, buf_C)
             call xgemmsg('T', 'N', m_size8, n_size8, k_size8, alpha8_gemm, &
                          buf_A, lda1_8, buf_B, lda2_8, beta8_accum, buf_C, lda3_8)
!$acc end host_data

             if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C) deviceptr(v3) private(m_loc, n_loc, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      idx_v3 = int(b_v3_a, 8) + int(m_loc, 8)*s3_m + int(n_loc, 8)*s3_n
                      v3(idx_v3) = sngl(d_fact1 * dble(buf_C(m_loc, n_loc)))
                   end do
                end do
             else
!$acc parallel loop collapse(2) present(buf_C) deviceptr(v3) private(m_loc, n_loc, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      idx_v3 = int(b_v3_a, 8) + int(m_loc, 8)*s3_m + int(n_loc, 8)*s3_n
                      v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * dble(buf_C(m_loc, n_loc)))
                   end do
                end do
             end if
          end do
!$acc end data
          deallocate(buf_A, buf_B, buf_C)

      else
          ! *************************************************************
          ! [引擎 B：生存模式] (超大内存需求体系的 Tiling 双层 FP64 防截断版本)
          ! *************************************************************
          allocate(buf_A(TILE_DIM, TILE_DIM), buf_B(TILE_DIM, TILE_DIM), buf_C(TILE_DIM, TILE_DIM), buf_C_fp64(TILE_DIM, TILE_DIM))
!$acc data create(buf_A, buf_B, buf_C, buf_C_fp64)

          do t_idx = 1, num_tasks
             if (is_macro_LL .or. is_macro_RR) then
                irr_M = task_list(t_idx)
                if (type3(1:1)==c_m) then; irr_N = dirprd(sym3, irr_M); base_v3=iabsoff3(irr_N)
                else; irr_N = irr_M; irr_M = dirprd(sym3, irr_N); base_v3=iabsoff3(irr_M); end if
                
                dim_M = pop1(pos1_m, irr_M); dim_N = pop2(pos2_n, irr_N)
                if (type3(1:1)==c_m) then; s3_m=1; s3_n=dim_M; else; s3_n=1; s3_m=dim_N; end if
                b_v3_a = base_v3 - s3_m - s3_n
                cur_M = dim_M; cur_N = dim_N

                do m_start = 1, cur_M, TILE_DIM; m_chunk = min(TILE_DIM, cur_M - m_start + 1)
                do n_start = 1, cur_N, TILE_DIM; n_chunk = min(TILE_DIM, cur_N - n_start + 1)

!$acc kernels present(buf_C_fp64)
                   buf_C_fp64 = 0.0d0
!$acc end kernels

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
                      else
                         irr_PQ = dirprd(sym1, dirprd(irr_R, irr_M)); dim_PQ = inumrht1(irr_PQ)
                         if (dim_PQ == 0) cycle
                         irr_l1=dirprd(irr_R,irr_M); irr_l2=dirprd(irr_R,irr_N)
                         irr_pos_m = irr_M; if(pos1_r==2) irr_pos_m=irr_R
                         base_v1 = iabsoff1(irr_PQ) + iposlft1_fb(irr_pos_m, irr_l1); s1(1)=inumlft1(irr_l1)
                         if(pos1_r==1) then; s1(2)=1; s1(3)=dim_R; else; s1(3)=1; s1(2)=dim_M; end if
                         irr_pos_n = irr_N; if(pos2_r==2) irr_pos_n=irr_R
                         base_v2 = iabsoff2(irr_PQ) + iposlft2_fb(irr_pos_n, irr_l2); s2(1)=inumlft2(irr_l2)
                         if(pos2_r==1) then; s2(2)=1; s2(3)=dim_R; else; s2(3)=1; s2(2)=dim_N; end if
                      end if
                      
                      b_v1_a = base_v1-s1(1)-s1(2)-s1(3); b_v2_a = base_v2-s2(1)-s2(2)-s2(3)
                      cur_K = dim_PQ * dim_R

                      do k_start = 1, cur_K, TILE_DIM
                         k_chunk = min(TILE_DIM, cur_K - k_start + 1)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, i_pq, i_r, v1_idx)
                         do k_loc = 1, k_chunk
                            do m_loc = 1, m_chunk
                               m_flat = m_start + m_loc - 1; k_flat = k_start + k_loc - 1
                               i_pq = mod(k_flat - 1, dim_PQ) + 1; i_r = (k_flat - 1) / dim_PQ + 1
                               v1_idx = int(b_v1_a, 8) + int(i_pq, 8)*s1(1) + int(i_r, 8)*s1(2) + int(m_flat, 8)*s1(3)
                               buf_A(k_loc, m_loc) = v1(v1_idx)
                            end do
                         end do

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, i_pq, i_r, v2_idx)
                         do n_loc = 1, n_chunk
                            do k_loc = 1, k_chunk
                               n_flat = n_start + n_loc - 1; k_flat = k_start + k_loc - 1
                               i_pq = mod(k_flat - 1, dim_PQ) + 1; i_r = (k_flat - 1) / dim_PQ + 1
                               v2_idx = int(b_v2_a, 8) + int(i_pq, 8)*s2(1) + int(i_r, 8)*s2(2) + int(n_flat, 8)*s2(3)
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
                   end do ! end irr_R

                   if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            idx_v3 = int(b_v3_a, 8) + int(m_flat, 8)*s3_m + int(n_flat, 8)*s3_n
                            v3(idx_v3) = sngl(d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   else
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            idx_v3 = int(b_v3_a, 8) + int(m_flat, 8)*s3_m + int(n_flat, 8)*s3_n
                            v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   end if
                end do; end do
             else
                ! --- 标准分支 C Tiling ---
                irr_v3_2 = task_list(t_idx)
                if (v3_map_m == 2) then; irr_m = irr_v3_2; irr_n = dirprd(irr_v3_2, sym3)
                else; irr_n = irr_v3_2; irr_m = dirprd(irr_v3_2, sym3); end if
                
                cur_M = pop1(pos1_m, irr_m); cur_N = pop2(pos2_n, irr_n)
                base_v3 = iabsoff3(irr_v3_2)
                if (v3_map_m == 1) then; s3_m = 1; s3_n = cur_M; else; s3_m = cur_N; s3_n = 1; end if
                b_v3_a = base_v3 - s3_m - s3_n

                do m_start = 1, cur_M, TILE_DIM; m_chunk = min(TILE_DIM, cur_M - m_start + 1)
                do n_start = 1, cur_N, TILE_DIM; n_chunk = min(TILE_DIM, cur_N - n_start + 1)

!$acc kernels present(buf_C_fp64)
                   buf_C_fp64 = 0.0d0
!$acc end kernels

                   do irr_k1=1,nirrep; do irr_k2=1,nirrep; irr_k3=dirprd(sym1,dirprd(irr_m,dirprd(irr_k1,irr_k2)))
                      irr_logi1(pos1_m)=irr_m; irr_logi1(pos1_k(1))=irr_k1; irr_logi1(pos1_k(2))=irr_k2; irr_logi1(pos1_k(3))=irr_k3
                      irr_logi2(pos2_n)=irr_n; irr_logi2(pos2_k(1))=irr_k1; irr_logi2(pos2_k(2))=irr_k2; irr_logi2(pos2_k(3))=irr_k3
                      if(dirprd(dirprd(irr_logi2(1),irr_logi2(2)),dirprd(irr_logi2(3),irr_logi2(4)))/=sym2) cycle
                      
                      dim_k(1)=pop1(pos1_k(1),irr_k1); dim_k(2)=pop1(pos1_k(2),irr_k2); dim_k(3)=pop1(pos1_k(3),irr_k3)
                      if(dim_k(1)==0 .or. dim_k(2)==0 .or. dim_k(3)==0) cycle
                      if(pop2(pos2_k(1),irr_k1)==0 .or. pop2(pos2_k(2),irr_k2)==0 .or. pop2(pos2_k(3),irr_k3)==0) cycle
                      cur_K_block = dim_k(1)*dim_k(2)*dim_k(3)
                      
                      irr_r1=dirprd(irr_logi1(3),irr_logi1(4)); irr_l1=dirprd(irr_logi1(1),irr_logi1(2))
                      irr_r2=dirprd(irr_logi2(3),irr_logi2(4)); irr_l2=dirprd(irr_logi2(1),irr_logi2(2))
                      
                      base_v1 = iabsoff1(irr_r1) + inumlft1(irr_l1)*iposrht1(irr_logi1(4),irr_r1) + iposlft1_fb(irr_logi1(2),irr_l1)
                      s_phys1(1)=1; s_phys1(2)=pop1(1,irr_logi1(1)); s_phys1(3)=inumlft1(irr_l1); s_phys1(4)=inumlft1(irr_l1)*pop1(3,irr_logi1(3))
                      
                      base_v2 = iabsoff2(irr_r2) + inumlft2(irr_l2)*iposrht2(irr_logi2(4),irr_r2) + iposlft2_fb(irr_logi2(2),irr_l2)
                      s_phys2(1)=1; s_phys2(2)=pop2(1,irr_logi2(1)); s_phys2(3)=inumlft2(irr_l2); s_phys2(4)=inumlft2(irr_l2)*pop2(3,irr_logi2(3))
                      
                      b_v1_a = base_v1-s_phys1(pos1_m)-s_phys1(pos1_k(1))-s_phys1(pos1_k(2))-s_phys1(pos1_k(3))
                      b_v2_a = base_v2-s_phys2(pos2_n)-s_phys2(pos2_k(1))-s_phys2(pos2_k(2))-s_phys2(pos2_k(3))
                      
                      do k_start = 1, cur_K_block, TILE_DIM
                         k_chunk = min(TILE_DIM, cur_K_block - k_start + 1)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, k_v, rem, v1_idx)
                         do k_loc = 1, k_chunk
                            do m_loc = 1, m_chunk
                               m_flat = m_start + m_loc - 1; k_flat = k_start + k_loc - 1
                               k_v(1)=mod(k_flat-1, dim_k(1))+1; rem=(k_flat-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                               v1_idx = int(b_v1_a, 8) + int(m_flat, 8)*s_phys1(pos1_m) + int(k_v(1), 8)*s_phys1(pos1_k(1)) + int(k_v(2), 8)*s_phys1(pos1_k(2)) + int(k_v(3), 8)*s_phys1(pos1_k(3))
                               buf_A(k_loc, m_loc) = v1(v1_idx)
                            end do
                         end do

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, k_v, rem, v2_idx)
                         do n_loc = 1, n_chunk
                            do k_loc = 1, k_chunk
                               n_flat = n_start + n_loc - 1; k_flat = k_start + k_loc - 1
                               k_v(1)=mod(k_flat-1, dim_k(1))+1; rem=(k_flat-1)/dim_k(1); k_v(2)=mod(rem, dim_k(2))+1; k_v(3)=rem/dim_k(2)+1
                               v2_idx = int(b_v2_a, 8) + int(n_flat, 8)*s_phys2(pos2_n) + int(k_v(1), 8)*s_phys2(pos2_k(1)) + int(k_v(2), 8)*s_phys2(pos2_k(2)) + int(k_v(3), 8)*s_phys2(pos2_k(3))
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
                   end do; end do ! end irr_D1
                   
                   if (fact2in == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            idx_v3 = int(b_v3_a, 8) + int(m_flat, 8)*s3_m + int(n_flat, 8)*s3_n
                            v3(idx_v3) = sngl(d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   else
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            idx_v3 = int(b_v3_a, 8) + int(m_flat, 8)*s3_m + int(n_flat, 8)*s3_n
                            v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   end if
                end do; end do
             end if
          end do
!$acc end data
          deallocate(buf_A, buf_B, buf_C, buf_C_fp64)
      end if

      return
      end subroutine vv2v_dpd_442_sg
