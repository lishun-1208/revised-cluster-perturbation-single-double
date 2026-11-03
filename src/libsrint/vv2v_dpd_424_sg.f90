! =======================================================================
! 终极工业级张量引擎 - GPU 单精度核心 (Single Precision, GPU + OpenACC)
! [核心能力] 424 模式: 4-index * 2-index -> 4-index 
! [架构特性]
!   1. 极速通道回归：完美支持 IOP=0/2 的 Intact Pair 压缩态快速收缩。
!   2. 防截断装甲：引入 buf_C_fp64，彻底吸收 Tiling 跨块累加时的单精度截断误差。
!   3. 降维打击：彻底消灭内层 if！通过 Flat-to-ND 索引映射，榨干访存速度。
!   4. 算力解放：连续化内存后，强制调用 cuBLAS (xgemmsg)，满血发挥 GPU FP32 单元。
! [稳定性修复] 在保留原版 800 行双擎框架与所有 OpenACC 指令的前提下，
!              完美修复了极速通道 (Fast Path) 中导致虚部/非全对称块
!              越界写入的物理偏移量 (iposlft) 和跨度步长 (inumlft) Bug！
! =======================================================================
      subroutine vv2v_dpd_424_sg(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                                 v2, p2_1, p2_2, iop2_12, sym2, type2, &
                                 v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1in, fact2in)
      
      use mod_symm
      implicit none

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type1, type3 
      character(len=2), intent(in)  :: type2        
      real*4, intent(in)            :: v1(*), v2(*)        
      real*4, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8) 
      integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop3_12, iop3_34
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1in, fact2in        

      ! 外部接口
      external :: xgemmsg
      integer, external :: irpdso  

      ! --- 内部计算系数 ---
      real*8 :: d_fact1, d_fact2
      real*4 :: alpha4_gemm, beta4_gemm

      ! --- Fast Path (IOP=0/2 压缩通道) 变量 ---
      logical :: is_fast_path, is_scatter
      character(len=1) :: transA, transB
      integer :: irr_lft, irr_fast
      integer :: irr_12, irr_34, irr_1, irr_2, irr_3, irr_4, irr_5
      integer :: irr_lft1, irr_rht1, irr_lft3, irr_rht3
      integer :: dim_12, dim_34, dim_1, dim_2, dim_3, dim_4, dim_5
      integer :: off_v1, off_v2, off_v3
      integer*8 :: M8, N8, K8, lda8, ldb8, ldc8
     
      integer :: M, N, K, lda, ldb
      integer :: idx_A, i1, i3, i4, i5, i12, i34
      integer*8 :: idx_C, total_size8, i12_8, i34_8
      real*4, allocatable :: fp_buf_A(:), fp_buf_C(:)

      ! --- Fallback 引擎映射与调度变量 ---
      integer, parameter :: TILE_DIM = 4000 
      character(len=1) :: c_ext(4), c_dum(1)
      integer :: pos1(5), pos2(5), pos3(4) 
      integer :: i, j, idx, cnt_m, cnt_n
      logical :: found, skip_block
      integer :: idx_m_in_v3(3), idx_n_in_v3(1)
      integer :: pop1(4,8), pop2(2,8), pop3(4,8)
      
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposlft1_fb(33,33), iposrht1(33,33)
      integer :: iabsoff2(33)
      integer :: inumlft3(33), inumrht3(33), iabsoff3(33), iposlft3_fb(33,33), iposrht3(33,33)
      
      integer :: num_tasks, task_list(3, 10000)
      integer :: task_idx, irr_L4, irr_L3, irr_L2, irr_L1
      integer :: dim_L1, dim_L2, dim_L3, dim_L4, dim_X_arr(4)
      integer :: cur_M, cur_N, cur_K
      integer :: irr_phys1(4), irr_phys2(2), irr_phys3(4), irr_logi(5)
      integer :: base_v3, s_phys3(4), s_v3(4), base_v3_adj
      integer :: irr_D1, irr_ext_v1, irr_ext_v2, irr_L1_tmp
    
      integer :: base_v1, base_v2, s_phys1(4), s_phys2(2)
      integer :: s_v1_arr(5), s_v2_arr(5)
      integer :: m_dim1, m_dim2, m_dim3, n_dim1
      
      ! --- 双擎动态路由变量 ---
      integer :: max_k, k_offset
      integer*8 :: global_max_m, global_max_n, global_max_k
      integer*8 :: req_mem_mb
      integer :: b_v1, b_v2, s_m1, s_m2, s_m3, s_n1, s_k1, s_k2
      integer :: s_v3_m1, s_v3_m2, s_v3_m3, s_v3_n1 
      integer :: m_1, m_2, m_3, n_1, k_1, rem
      integer*8 :: v1_idx, v2_idx, idx_v3
      integer :: m_start, n_start, k_start, m_chunk, n_chunk, k_chunk
      integer :: m_loc, n_loc, k_loc, m_flat, n_flat, k_flat
      
      integer*8 :: m_size8, n_size8, k_size8, lda1_8, lda2_8, lda3_8
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)
      real*8, allocatable :: buf_C_fp64(:,:)  ! 单精度防截断装甲

      ! 浮点系数转换
      d_fact1 = fact1in
      d_fact2 = fact2in
      alpha4_gemm = sngl(fact1in)
      beta4_gemm  = sngl(fact2in)

      ! =============================================================
      ! 【阶段 0：映射重组与物理表构建】
      ! =============================================================
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:)
      pop3(1,:) = p3_1(:); pop3(2,:) = p3_2(:); pop3(3,:) = p3_3(:); pop3(4,:) = p3_4(:)

      inumlft3 = 0; inumrht3 = 0
      do i = 1, nirrep
         inumlft3(i) = irpdso(i, p3_1, p3_2, iop3_12)
         inumrht3(i) = irpdso(i, p3_3, p3_4, iop3_34)
      end do
      iabsoff3(1) = 1
      do i = 1, nirrep
         irr_lft = dirprd(i, sym3)
         iabsoff3(i+1) = iabsoff3(i) + inumlft3(irr_lft) * inumrht3(i)
        
         iposlft3_fb(1, i) = 0
         iposrht3(1, i) = 0
         do j = 1, nirrep - 1
            irr_fast = dirprd(j, i)
            iposlft3_fb(j+1, i) = iposlft3_fb(j, i) + pop3(1, irr_fast) * pop3(2, j)
            iposrht3(j+1, i) = iposrht3(j, i) + pop3(3, irr_fast) * pop3(4, j)
         end do
      end do
      
      inumlft1 = 0; inumrht1 = 0
      do i = 1, nirrep
         inumlft1(i) = irpdso(i, p1_1, p1_2, iop1_12)
         inumrht1(i) = irpdso(i, p1_3, p1_4, iop1_34)
      end do
      iabsoff1(1) = 1
      do i = 1, nirrep
         irr_lft = dirprd(i, sym1)
         iabsoff1(i+1) = iabsoff1(i) + inumlft1(irr_lft) * inumrht1(i)
        
         iposlft1_fb(1, i) = 0
         iposrht1(1, i) = 0
         do j = 1, nirrep - 1
            irr_fast = dirprd(j, i)
            iposlft1_fb(j+1, i) = iposlft1_fb(j, i) + pop1(1, irr_fast) * pop1(2, j)
            iposrht1(j+1, i) = iposrht1(j, i) + pop1(3, irr_fast) * pop1(4, j)
         end do
      end do
      
      iabsoff2(1) = 1
      do i = 1, nirrep - 1
         irr_L1_tmp = dirprd(i, sym2)
         iabsoff2(i+1) = iabsoff2(i) + pop2(1, irr_L1_tmp) * pop2(2, i)
      end do

      ! =============================================================
      ! 【阶段 1：智能路由探测器 (完美解决 IOP=0/2)】
      ! =============================================================
      is_fast_path = .false.
      is_scatter   = .false.
      if (type3 == '1234') then
          if (type1(1:2) == '12') then
              if (type1 == '1235' .and. type2 == '54') then; is_fast_path = .true.; end if
              if (type1 == '1235' .and. type2 == '45') then; is_fast_path = .true.; end if
              if (type1 == '1253' .and. type2 == '54') then; is_fast_path = .true.; end if
              if (type1 == '1253' .and. type2 == '45') then; is_fast_path = .true.; end if
          else if (type1(3:4) == '12') then
              if (type1 == '3512' .and. type2 == '54') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3512' .and. type2 == '45') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '5312' .and. type2 == '54') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '5312' .and. type2 == '45') then; is_fast_path = .true.; is_scatter = .true.; end if
          else if (type1(3:4) == '34') then
              if (type1 == '1534' .and. type2 == '52') then; is_fast_path = .true.; end if
              if (type1 == '1534' .and. type2 == '25') then; is_fast_path = .true.; end if
              if (type1 == '5134' .and. type2 == '52') then; is_fast_path = .true.; end if
              if (type1 == '5134' .and. type2 == '25') then; is_fast_path = .true.; end if
          else if (type1(1:2) == '34') then
              if (type1 == '3415' .and. type2 == '52') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3415' .and. type2 == '25') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3451' .and. type2 == '52') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3451' .and. type2 == '25') then; is_fast_path = .true.; is_scatter = .true.; end if
          end if
      end if

      ! >>> 极速通道 (Fast Path) GPU 版本 <<<
      if (is_fast_path) then
          if (is_scatter) then
             allocate(fp_buf_A(64000000), fp_buf_C(64000000))
!$acc enter data create(fp_buf_A, fp_buf_C)
          end if
          
          if (type1(1:2) == '12' .or. type1(3:4) == '12') then
              do irr_4 = 1, nirrep
                 irr_5 = dirprd(irr_4, sym2)
                 do irr_3 = 1, nirrep
                    irr_12 = dirprd(sym3, dirprd(irr_3, irr_4))
                    dim_12 = inumlft3(irr_12)
                    dim_3  = pop3(3, irr_3)
                    dim_4  = pop3(4, irr_4)
                    
                    dim_5  = 0
                    if (type2 == '54') dim_5 = pop2(1, irr_5)
                    if (type2 == '45') dim_5 = pop2(2, irr_5)
                    
                    irr_rht3 = dirprd(irr_3, irr_4)
                    off_v3 = iabsoff3(irr_rht3) + dim_12 * iposrht3(irr_4, irr_rht3)

                    if (dim_12 == 0 .or. dim_3 == 0 .or. dim_4 == 0) cycle
                    if (dim_5 == 0) then
                        total_size8 = int(dim_12, 8) * int(dim_3, 8) * int(dim_4, 8)
                        if (d_fact2 == 0.0d0) then
!$acc parallel loop deviceptr(v3)
                           do i12_8 = 0, total_size8 - 1
                              v3(off_v3 + i12_8) = 0.0e0
                           end do
                        else if (d_fact2 /= 1.0d0) then
!$acc parallel loop deviceptr(v3)
                           do i12_8 = 0, total_size8 - 1
                              v3(off_v3 + i12_8) = sngl(d_fact2 * dble(v3(off_v3 + i12_8)))
                           end do
                        end if
                        cycle
                    end if
                    
                    if (type2 == '54') off_v2 = iabsoff2(irr_4)
                    if (type2 == '45') off_v2 = iabsoff2(irr_5)

                    if (type1 == '1235' .or. type1 == '1253') then
                       irr_rht1 = dirprd(irr_3, irr_5)
                       off_v1 = iabsoff1(irr_rht1)
                       if (type1 == '1235') off_v1 = off_v1 + dim_12 * iposrht1(irr_5, irr_rht1)
                       if (type1 == '1253') off_v1 = off_v1 + dim_12 * iposrht1(irr_3, irr_rht1)
                       
                       if (type1 == '1235') then
                          M = dim_12 * dim_3; N = dim_4; K = dim_5; lda = M
                          if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                          M8 = M; N8 = N; K8 = K; lda8 = lda; ldb8 = ldb
                          call xgemmsg('N', transB, M8, N8, K8, alpha4_gemm, v1(off_v1), lda8, v2(off_v2), ldb8, beta4_gemm, v3(off_v3), M8)
                       else
                          M = dim_12; N = dim_4; K = dim_5; lda = M
                          if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                          M8 = M; N8 = N; K8 = K; lda8 = lda; ldb8 = ldb; ldc8 = M8 * int(dim_3, 8)
                          do i3 = 0, dim_3 - 1
                             call xgemmsg('N', transB, M8, N8, K8, alpha4_gemm, v1(off_v1 + i3*M*K), lda8, v2(off_v2), ldb8, beta4_gemm, v3(off_v3 + i3*M), ldc8)
                          end do
                       end if
                       
                    else if (type1 == '3512' .or. type1 == '5312') then
                       ! 【BUGFIX 1】: 补充内部对称性偏移量
                       irr_lft1 = dirprd(irr_3, irr_5)
                       off_v1 = iabsoff1(irr_12)
                       if (type1 == '3512') off_v1 = off_v1 + iposlft1_fb(irr_5, irr_lft1)
                       if (type1 == '5312') off_v1 = off_v1 + iposlft1_fb(irr_3, irr_lft1)

                       if (type1 == '3512') then; M = dim_3; K = dim_5; lda = dim_3; transA = 'N'; else; M = dim_3; K = dim_5; lda = dim_5; transA = 'T'; end if
                       N = dim_4
                       if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                       M8 = M; N8 = N; K8 = K; lda8 = lda; ldb8 = ldb
                       
                       do i12 = 0, dim_12 - 1
                          ! 【BUGFIX 2】: 修正 v1 的步长使用全局宏维度跨度 inumlft1
                          call xgemmsg(transA, transB, M8, N8, K8, 1.0e0, v1(off_v1 + i12*inumlft1(irr_lft1)), lda8, v2(off_v2), ldb8, 0.0e0, fp_buf_C, M8)
                          
!$acc parallel loop collapse(2) present(fp_buf_C) deviceptr(v3) private(idx_C)
                          do i4 = 0, dim_4 - 1
                             do i3 = 0, dim_3 - 1
                                idx_C = int(off_v3, 8) + int(i12, 8) + int(i3, 8)*dim_12 + int(i4, 8)*dim_12*dim_3
                                if (d_fact2 == 0.0d0) then
                                   v3(idx_C) = sngl(d_fact1 * dble(fp_buf_C(1 + i3 + i4*dim_3)))
                                else
                                   v3(idx_C) = sngl(d_fact2 * dble(v3(idx_C)) + d_fact1 * dble(fp_buf_C(1 + i3 + i4*dim_3)))
                                end if
                             end do
                          end do
                       end do
                    end if
                 end do
              end do
              
          else
              ! 34xx 系列
              do irr_2 = 1, nirrep
                 irr_5 = dirprd(irr_2, sym2)
                 do irr_1 = 1, nirrep
                    irr_12 = dirprd(irr_1, irr_2)
                    irr_34 = dirprd(sym3, irr_12)
                    
                    dim_34 = inumrht3(irr_34)
                    dim_1  = pop3(1, irr_1)
                    dim_2  = pop3(2, irr_2)
                    
                    dim_5  = 0
                    if (type2 == '52') dim_5 = pop2(1, irr_5)
                    if (type2 == '25') dim_5 = pop2(2, irr_5)
                    
                    ! 【BUGFIX 3】: 修正 v3 内部块物理偏移
                    off_v3 = iabsoff3(irr_34) + iposlft3_fb(irr_2, irr_12)
                    
                    if (dim_34 == 0 .or. dim_1 == 0 .or. dim_2 == 0) cycle
                    if (dim_5 == 0) then
                        total_size8 = int(dim_34, 8) * int(dim_1, 8) * int(dim_2, 8)
                        if (d_fact2 == 0.0d0) then
!$acc parallel loop deviceptr(v3)
                           do i34_8 = 0, total_size8 - 1
                              v3(off_v3 + i34_8) = 0.0e0
                           end do
                        else if (d_fact2 /= 1.0d0) then
!$acc parallel loop deviceptr(v3)
                           do i34_8 = 0, total_size8 - 1
                              v3(off_v3 + i34_8) = sngl(d_fact2 * dble(v3(off_v3 + i34_8)))
                           end do
                        end if
                        cycle
                    end if
                    
                    if (type2 == '52') off_v2 = iabsoff2(irr_2)
                    if (type2 == '25') off_v2 = iabsoff2(irr_5)
                    
                    if (type1 == '1534' .or. type1 == '5134') then
                       irr_lft1 = dirprd(irr_1, irr_5)
                       off_v1 = iabsoff1(irr_34)
                       if (type1 == '1534') then
                           off_v1 = off_v1 + iposlft1_fb(irr_5, irr_lft1)
                           M = dim_1; K = dim_5; N = dim_2; lda = M; transA = 'N'
                       else
                           off_v1 = off_v1 + iposlft1_fb(irr_1, irr_lft1)
                           M = dim_1; K = dim_5; N = dim_2; lda = K; transA = 'T'
                       end if
                       if (type2 == '52') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                       M8 = M; N8 = N; K8 = K; lda8 = lda; ldb8 = ldb
                       do i34 = 0, dim_34 - 1
                          ! 【BUGFIX 4】: 修正 34xx 系列调用 xgemm 时的跨度步长 (inumlft)
                          call xgemmsg(transA, transB, M8, N8, K8, alpha4_gemm, v1(off_v1 + i34*inumlft1(irr_lft1)), lda8, v2(off_v2), ldb8, beta4_gemm, v3(off_v3 + i34*inumlft3(irr_12)), M8)
                       end do
                       
                    else if (type1 == '3415' .or. type1 == '3451') then
                       irr_rht1 = dirprd(irr_1, irr_5)
                       off_v1 = iabsoff1(irr_rht1)
                       if (type1 == '3415') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_5, irr_rht1)
                       if (type1 == '3451') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_1, irr_rht1)
                       
                       if (type2 == '52') then; ldb = dim_5; transB = 'N'; else; ldb = dim_2; transB = 'T'; end if
                       ldb8 = ldb
                       
                       do i34 = 0, dim_34 - 1
                          if (type1 == '3415') then
!$acc parallel loop collapse(2) present(fp_buf_A) deviceptr(v1)
                             do i5 = 0, dim_5 - 1
                                do i1 = 0, dim_1 - 1
                                   fp_buf_A(1 + i1 + i5*dim_1) = v1(off_v1 + i34 + i1*inumlft1(irr_34) + i5*inumlft1(irr_34)*dim_1)
                                end do
                             end do
                             M = dim_1; K = dim_5; lda = dim_1; transA = 'N'
                          else
!$acc parallel loop collapse(2) present(fp_buf_A) deviceptr(v1)
                             do i1 = 0, dim_1 - 1
                                do i5 = 0, dim_5 - 1
                                   fp_buf_A(1 + i5 + i1*dim_5) = v1(off_v1 + i34 + i5*inumlft1(irr_34) + i1*inumlft1(irr_34)*dim_5)
                                end do
                             end do
                             M = dim_1; K = dim_5; lda = dim_5; transA = 'T'
                          end if
                          M8 = M; N8 = dim_2; K8 = K; lda8 = lda
                          call xgemmsg(transA, transB, M8, N8, K8, alpha4_gemm, fp_buf_A, lda8, v2(off_v2), ldb8, beta4_gemm, v3(off_v3 + i34*inumlft3(irr_12)), M8)
                       end do
                    end if
                 end do
              end do
          end if
          
          if (is_scatter) then
!$acc exit data delete(fp_buf_A, fp_buf_C)
             deallocate(fp_buf_A, fp_buf_C)
          end if
          return
      end if

      ! =============================================================
      ! 【阶段 2：Fallback 阻断控制墙】
      ! 所有合法 IOP=0/2 已被极速通道拦截。走到此步的若含 IOP=0/2 必为物理拆分异常项。
      ! =============================================================
      if (iop1_12 /= 1 .or. iop1_34 /= 1 .or. iop2_12 /= 1 .or. &
          iop3_12 /= 1 .or. iop3_34 /= 1) then
          write(6,*) " FATAL ERROR in vv2v_dpd_424_sg: Invalid physical scramble for IOP=0/2 compressed space."
          stop
      end if

      ! =============================================================
      ! 【阶段 3：预调度探测器】
      ! =============================================================
      do i = 1, 4
         c_ext(i) = type3(i:i)
      end do
      
      idx = 0
      do i = 1, 4
         found = .false.
         do j = 1, 4
            if (type1(i:i) == c_ext(j)) found = .true.
         end do
         if (.not. found) c_dum(1) = type1(i:i)
      end do

      pos1 = 0; pos2 = 0; pos3 = 0
      do i = 1, 4
         do j = 1, 4
            if (type1(j:j) == c_ext(i)) pos1(i) = j
            if (type3(j:j) == c_ext(i)) pos3(i) = j
         end do
         do j = 1, 2
            if (type2(j:j) == c_ext(i)) pos2(i) = j
         end do
      end do
      
      do j = 1, 4
         if (type1(j:j) == c_dum(1)) pos1(5) = j
      end do
      do j = 1, 2
         if (type2(j:j) == c_dum(1)) pos2(5) = j
      end do

      cnt_m = 0; cnt_n = 0; idx_m_in_v3 = 0; idx_n_in_v3 = 0
      do i = 1, 4
         if (pos1(i) /= 0) then
            cnt_m = cnt_m + 1
            idx_m_in_v3(cnt_m) = i
         end if
         if (pos2(i) /= 0) then
            cnt_n = cnt_n + 1
            idx_n_in_v3(cnt_n) = i
         end if
      end do

      num_tasks = 0
      do irr_rht3 = 1, nirrep
         irr_lft3 = dirprd(irr_rht3, sym3)      
         do irr_L4 = 1, nirrep
            irr_L3 = dirprd(irr_rht3, irr_L4)    
            do irr_L2 = 1, nirrep
               irr_L1 = dirprd(irr_lft3, irr_L2)    
               dim_L1 = pop3(pos3(1), irr_L1)
               dim_L2 = pop3(pos3(2), irr_L2)
               dim_L3 = pop3(pos3(3), irr_L3)
               dim_L4 = pop3(pos3(4), irr_L4)
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
      ! 【阶段 4：显存极速预扫描】
      ! =============================================================
      global_max_m = 0; global_max_n = 0; global_max_k = 0
      do task_idx = 1, num_tasks
         irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
         irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)
         
         dim_X_arr(1) = pop3(pos3(1), irr_L1)
         dim_X_arr(2) = pop3(pos3(2), irr_L2)
         dim_X_arr(3) = pop3(pos3(3), irr_L3)
         dim_X_arr(4) = pop3(pos3(4), irr_L4)
         
         cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2)) * dim_X_arr(idx_m_in_v3(3))
         cur_N = dim_X_arr(idx_n_in_v3(1))
         
         if (cur_M > global_max_m) global_max_m = cur_M
         if (cur_N > global_max_n) global_max_n = cur_N
         max_k = 0
         
         do irr_D1 = 1, nirrep
            irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /)
            irr_logi(5)   = irr_D1
            
            irr_ext_v2 = irr_logi(idx_n_in_v3(1))
            if (dirprd(irr_ext_v2, irr_D1) /= sym2) cycle

            irr_ext_v1 = dirprd(irr_logi(idx_m_in_v3(1)), irr_logi(idx_m_in_v3(2)))
            irr_ext_v1 = dirprd(irr_ext_v1, irr_logi(idx_m_in_v3(3)))
            if (dirprd(irr_ext_v1, irr_D1) /= sym1) cycle

            skip_block = .false.
            do i = 1, 5
               if (pos1(i)/=0 .and. pop1(pos1(i), irr_logi(i))==0) skip_block = .true.
               if (pos2(i)/=0 .and. pop2(pos2(i), irr_logi(i))==0) skip_block = .true.
            end do
            if (skip_block) cycle
            max_k = max_k + pop1(pos1(5), irr_D1)
         end do
         if (max_k > global_max_k) global_max_k = max_k
      end do

      if (global_max_m == 0 .or. global_max_n == 0) return
      req_mem_mb = (global_max_m * global_max_k + global_max_k * global_max_n + global_max_m * global_max_n) * 8_8 / (1024_8 * 1024_8)

      ! =============================================================
      ! 【阶段 5：Fallback 动态双擎路由执行 (搭载 fp64 防截断装甲)】
      ! =============================================================
      if (req_mem_mb <= 500_8) then
          ! *************************************************************
          ! [引擎 A：狂飙模式]
          ! *************************************************************
          allocate(buf_A(global_max_k, global_max_m), buf_B(global_max_k, global_max_n), buf_C(global_max_m, global_max_n))
          allocate(buf_C_fp64(global_max_m, global_max_n))
!$acc data create(buf_A, buf_B, buf_C, buf_C_fp64)

          do task_idx = 1, num_tasks
             irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
             irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)    
             
             irr_phys3(pos3(1)) = irr_L1
             irr_phys3(pos3(2)) = irr_L2
             irr_phys3(pos3(3)) = irr_L3
             irr_phys3(pos3(4)) = irr_L4
             
             dim_X_arr(1) = pop3(pos3(1), irr_L1)
             dim_X_arr(2) = pop3(pos3(2), irr_L2)
             dim_X_arr(3) = pop3(pos3(3), irr_L3)
             dim_X_arr(4) = pop3(pos3(4), irr_L4)
             
             cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2)) * dim_X_arr(idx_m_in_v3(3))
             cur_N = dim_X_arr(idx_n_in_v3(1))
           
             base_v3 = iabsoff3(irr_rht3) + inumlft3(irr_lft3) * iposrht3(irr_phys3(4), irr_rht3) + iposlft3_fb(irr_phys3(2), irr_lft3)
             s_phys3(1) = 1; s_phys3(2) = pop3(1, irr_phys3(1)); s_phys3(3) = inumlft3(irr_lft3); s_phys3(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))
             do i = 1, 4; s_v3(i) = s_phys3(pos3(i)); end do
             base_v3_adj = base_v3 - sum(s_v3(1:4))
             
             m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2)); m_dim3 = dim_X_arr(idx_m_in_v3(3))
             n_dim1 = dim_X_arr(idx_n_in_v3(1))
             max_k = 0; k_offset = 0

             ! [安全归零] fp64 防截断装甲块初始化
!$acc parallel loop collapse(2) present(buf_C_fp64)
             do n_loc = 1, cur_N
                do m_loc = 1, cur_M
                   buf_C_fp64(m_loc, n_loc) = 0.0d0
                end do
             end do

             do irr_D1 = 1, nirrep
                irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /)
                irr_logi(5)   = irr_D1
                
                irr_ext_v2 = irr_logi(idx_n_in_v3(1))
                if (dirprd(irr_ext_v2, irr_D1) /= sym2) cycle

                irr_ext_v1 = dirprd(irr_logi(idx_m_in_v3(1)), irr_logi(idx_m_in_v3(2)))
                irr_ext_v1 = dirprd(irr_ext_v1, irr_logi(idx_m_in_v3(3)))
                if (dirprd(irr_ext_v1, irr_D1) /= sym1) cycle

                skip_block = .false.
                do i = 1, 5
                   if (pos1(i)/=0 .and. pop1(pos1(i), irr_logi(i))==0) skip_block = .true.
                   if (pos2(i)/=0 .and. pop2(pos2(i), irr_logi(i))==0) skip_block = .true.
                end do
                if (skip_block) cycle
                
                cur_K = pop1(pos1(5), irr_D1)
                if (cur_K == 0) cycle

                do i = 1, 5
                   if(pos1(i)/=0) irr_phys1(pos1(i))=irr_logi(i)
                   if(pos2(i)/=0) irr_phys2(pos2(i))=irr_logi(i)
                end do
                
                irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2))
                irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                
                base_v1 = iabsoff1(irr_rht1) + inumlft1(irr_lft1) * iposrht1(irr_phys1(4), irr_rht1) + iposlft1_fb(irr_phys1(2), irr_lft1)
                base_v2 = iabsoff2(irr_phys2(2))
                
                s_phys1(1)=1; s_phys1(2)=pop1(1,irr_phys1(1)); s_phys1(3)=inumlft1(irr_lft1); s_phys1(4)=inumlft1(irr_lft1)*pop1(3,irr_phys1(3))
                s_v1_arr(:)=0; do i=1,5; if(pos1(i)/=0) s_v1_arr(i)=s_phys1(pos1(i)); end do
                
                s_phys2(1)=1; s_phys2(2)=pop2(1,irr_phys2(1))
                s_v2_arr(:)=0; do i=1,5; if(pos2(i)/=0) s_v2_arr(i)=s_phys2(pos2(i)); end do
                
                b_v1 = base_v1 - sum(s_v1_arr)
                b_v2 = base_v2 - sum(s_v2_arr)
                
                s_m1 = s_v1_arr(idx_m_in_v3(1)); s_m2 = s_v1_arr(idx_m_in_v3(2)); s_m3 = s_v1_arr(idx_m_in_v3(3))
                s_n1 = s_v2_arr(idx_n_in_v3(1))
                s_k1 = s_v1_arr(5)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, m_1, m_2, m_3, rem, k_1, v1_idx)
                do k_loc = 1, cur_K
                   do m_loc = 1, cur_M
                      m_flat = m_loc; k_flat = k_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      rem = (m_flat - 1) / m_dim1
                      m_2 = mod(rem, m_dim2) + 1
                      m_3 = rem / m_dim2 + 1
                      k_1 = k_flat
                      v1_idx = int(b_v1, 8) + int(m_3, 8)*s_m3 + int(m_2, 8)*s_m2 + int(m_1, 8)*s_m1 + int(k_1, 8)*s_k1
                      buf_A(k_offset + k_loc, m_loc) = v1(v1_idx)
                   end do
                end do

                s_k2 = s_v2_arr(5)

!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, n_1, k_1, v2_idx)
                do n_loc = 1, cur_N
                   do k_loc = 1, cur_K
                      n_flat = n_loc; k_flat = k_loc
                      n_1 = n_flat; k_1 = k_flat
                      v2_idx = int(b_v2, 8) + int(n_1, 8)*s_n1 + int(k_1, 8)*s_k2
                      buf_B(k_offset + k_loc, n_loc) = v2(v2_idx)
                   end do
                end do

                k_offset = k_offset + cur_K
             end do
             max_k = k_offset
             
             if (max_k == 0) then
                 s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2)); s_v3_m3 = s_v3(idx_m_in_v3(3))
                 s_v3_n1 = s_v3(idx_n_in_v3(1))
                 if (d_fact2 == 0.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          m_flat = m_loc; n_flat = n_loc
                          m_1 = mod(m_flat - 1, m_dim1) + 1
                          rem = (m_flat - 1) / m_dim1
                          m_2 = mod(rem, m_dim2) + 1
                          m_3 = rem / m_dim2 + 1
                          n_1 = n_flat
                          idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
                          v3(idx_v3) = 0.0e0
                       end do
                    end do
                 else if (d_fact2 /= 1.0d0) then
!$acc parallel loop collapse(2) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                    do n_loc = 1, cur_N
                       do m_loc = 1, cur_M
                          m_flat = m_loc; n_flat = n_loc
                          m_1 = mod(m_flat - 1, m_dim1) + 1
                          rem = (m_flat - 1) / m_dim1
                          m_2 = mod(rem, m_dim2) + 1
                          m_3 = rem / m_dim2 + 1
                          n_1 = n_flat
                          idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
                          v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)))
                       end do
                    end do
                 end if
                 cycle
             end if
             
!$acc wait

             m_size8 = cur_M; n_size8 = cur_N; k_size8 = max_k
             lda1_8 = global_max_k; lda2_8 = global_max_k; lda3_8 = global_max_m

!$acc host_data use_device(buf_A, buf_B, buf_C)
             call xgemmsg('T', 'N', m_size8, n_size8, k_size8, 1.0e0, &
                          buf_A, lda1_8, buf_B, lda2_8, 0.0e0, buf_C, lda3_8)
!$acc end host_data

             ! [截断吸收] 双精度高阶累加
!$acc parallel loop collapse(2) present(buf_C, buf_C_fp64)
             do n_loc = 1, cur_N
                do m_loc = 1, cur_M
                   buf_C_fp64(m_loc, n_loc) = buf_C_fp64(m_loc, n_loc) + dble(buf_C(m_loc, n_loc))
                end do
             end do

             s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2)); s_v3_m3 = s_v3(idx_m_in_v3(3))
             s_v3_n1 = s_v3(idx_n_in_v3(1))

             if (d_fact2 == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      m_flat = m_loc; n_flat = n_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      rem = (m_flat - 1) / m_dim1
                      m_2 = mod(rem, m_dim2) + 1
                      m_3 = rem / m_dim2 + 1
                      n_1 = n_flat
                      idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
                      v3(idx_v3) = sngl(d_fact1 * buf_C_fp64(m_loc, n_loc))
                   end do
                end do
             else
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                do n_loc = 1, cur_N
                   do m_loc = 1, cur_M
                      m_flat = m_loc; n_flat = n_loc
                      m_1 = mod(m_flat - 1, m_dim1) + 1
                      rem = (m_flat - 1) / m_dim1
                      m_2 = mod(rem, m_dim2) + 1
                      m_3 = rem / m_dim2 + 1
                      n_1 = n_flat
                      idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
                      v3(idx_v3) = sngl(d_fact2 * dble(v3(idx_v3)) + d_fact1 * buf_C_fp64(m_loc, n_loc))
                   end do
                end do
             end if
          end do
!$acc end data
          deallocate(buf_A, buf_B, buf_C, buf_C_fp64)

      else
          ! *************************************************************
          ! [引擎 B：生存模式 (Tiling 分块)] 
          ! *************************************************************
          allocate(buf_A(TILE_DIM, TILE_DIM), buf_B(TILE_DIM, TILE_DIM), buf_C(TILE_DIM, TILE_DIM))
          allocate(buf_C_fp64(TILE_DIM, TILE_DIM))
!$acc data create(buf_A, buf_B, buf_C, buf_C_fp64)

          do task_idx = 1, num_tasks
             irr_rht3 = task_list(1, task_idx); irr_L4 = task_list(2, task_idx); irr_L2 = task_list(3, task_idx)
             irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)    
             
             irr_phys3(pos3(1)) = irr_L1
             irr_phys3(pos3(2)) = irr_L2
             irr_phys3(pos3(3)) = irr_L3
             irr_phys3(pos3(4)) = irr_L4
             
             dim_X_arr(1) = pop3(pos3(1), irr_L1)
             dim_X_arr(2) = pop3(pos3(2), irr_L2)
             dim_X_arr(3) = pop3(pos3(3), irr_L3)
             dim_X_arr(4) = pop3(pos3(4), irr_L4)
             
             cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2)) * dim_X_arr(idx_m_in_v3(3))
             cur_N = dim_X_arr(idx_n_in_v3(1))
           
             base_v3 = iabsoff3(irr_rht3) + inumlft3(irr_lft3) * iposrht3(irr_phys3(4), irr_rht3) + iposlft3_fb(irr_phys3(2), irr_lft3)
             s_phys3(1) = 1; s_phys3(2) = pop3(1, irr_phys3(1)); s_phys3(3) = inumlft3(irr_lft3); s_phys3(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))
             do i = 1, 4; s_v3(i) = s_phys3(pos3(i)); end do
             base_v3_adj = base_v3 - sum(s_v3(1:4))
             
             m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2)); m_dim3 = dim_X_arr(idx_m_in_v3(3))
             n_dim1 = dim_X_arr(idx_n_in_v3(1))

             do m_start = 1, cur_M, TILE_DIM
                m_chunk = min(TILE_DIM, cur_M - m_start + 1)
                do n_start = 1, cur_N, TILE_DIM
                  n_chunk = min(TILE_DIM, cur_N - n_start + 1)

!$acc kernels present(buf_C_fp64)
                   buf_C_fp64 = 0.0d0
!$acc end kernels

                   do irr_D1 = 1, nirrep
                      irr_logi(1:4) = (/ irr_L1, irr_L2, irr_L3, irr_L4 /)
                      irr_logi(5)   = irr_D1
                      
                      irr_ext_v2 = irr_logi(idx_n_in_v3(1))
                      if (dirprd(irr_ext_v2, irr_D1) /= sym2) cycle

                      irr_ext_v1 = dirprd(irr_logi(idx_m_in_v3(1)), irr_logi(idx_m_in_v3(2)))
                      irr_ext_v1 = dirprd(irr_ext_v1, irr_logi(idx_m_in_v3(3)))
                      if (dirprd(irr_ext_v1, irr_D1) /= sym1) cycle

                      skip_block = .false.
                      do i = 1, 5
                         if (pos1(i)/=0 .and. pop1(pos1(i), irr_logi(i))==0) skip_block = .true.
                         if (pos2(i)/=0 .and. pop2(pos2(i), irr_logi(i))==0) skip_block = .true.
                      end do
                      if (skip_block) cycle
                      
                      cur_K = pop1(pos1(5), irr_D1)
                      if (cur_K == 0) cycle

                      do i = 1, 5
                         if(pos1(i)/=0) irr_phys1(pos1(i))=irr_logi(i)
                         if(pos2(i)/=0) irr_phys2(pos2(i))=irr_logi(i)
                      end do
                      
                      irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2))
                      irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                     
                      base_v1 = iabsoff1(irr_rht1) + inumlft1(irr_lft1) * iposrht1(irr_phys1(4), irr_rht1) + iposlft1_fb(irr_phys1(2), irr_lft1)
                      base_v2 = iabsoff2(irr_phys2(2))
                      
                      s_phys1(1)=1; s_phys1(2)=pop1(1,irr_phys1(1)); s_phys1(3)=inumlft1(irr_lft1); s_phys1(4)=inumlft1(irr_lft1)*pop1(3,irr_phys1(3))
                      s_v1_arr(:)=0; do i=1,5; if(pos1(i)/=0) s_v1_arr(i)=s_phys1(pos1(i)); end do
                      
                      s_phys2(1)=1; s_phys2(2)=pop2(1,irr_phys2(1))
                      s_v2_arr(:)=0; do i=1,5; if(pos2(i)/=0) s_v2_arr(i)=s_phys2(pos2(i)); end do
                      
                      b_v1 = base_v1 - sum(s_v1_arr)
                      b_v2 = base_v2 - sum(s_v2_arr)
                     
                      s_m1 = s_v1_arr(idx_m_in_v3(1)); s_m2 = s_v1_arr(idx_m_in_v3(2)); s_m3 = s_v1_arr(idx_m_in_v3(3))
                      s_n1 = s_v2_arr(idx_n_in_v3(1))

                      do k_start = 1, cur_K, TILE_DIM
                         k_chunk = min(TILE_DIM, cur_K - k_start + 1)
                         s_k1 = s_v1_arr(5)

!$acc parallel loop collapse(2) present(buf_A) deviceptr(v1) async(1) private(m_flat, k_flat, m_1, m_2, m_3, rem, k_1, v1_idx)
                         do k_loc = 1, k_chunk
                            do m_loc = 1, m_chunk
                               m_flat = m_start + m_loc - 1; k_flat = k_start + k_loc - 1
                               m_1 = mod(m_flat - 1, m_dim1) + 1
                               rem = (m_flat - 1) / m_dim1
                               m_2 = mod(rem, m_dim2) + 1
                               m_3 = rem / m_dim2 + 1
                               k_1 = k_flat
                               v1_idx = int(b_v1, 8) + int(m_3, 8)*s_m3 + int(m_2, 8)*s_m2 + int(m_1, 8)*s_m1 + int(k_1, 8)*s_k1
                               buf_A(k_loc, m_loc) = v1(v1_idx)
                            end do
                         end do

                         s_k2 = s_v2_arr(5)
!$acc parallel loop collapse(2) present(buf_B) deviceptr(v2) async(2) private(n_flat, k_flat, n_1, k_1, v2_idx)
                         do n_loc = 1, n_chunk
                           do k_loc = 1, k_chunk
                               n_flat = n_start + n_loc - 1; k_flat = k_start + k_loc - 1
                               n_1 = n_flat; k_1 = k_flat
                               v2_idx = int(b_v2, 8) + int(n_1, 8)*s_n1 + int(k_1, 8)*s_k2
                               buf_B(k_loc, n_loc) = v2(v2_idx)
                            end do
                         end do
!$acc wait

                         m_size8 = m_chunk; n_size8 = n_chunk; k_size8 = k_chunk
                         lda1_8 = TILE_DIM; lda2_8 = TILE_DIM; lda3_8 = TILE_DIM

!$acc host_data use_device(buf_A, buf_B, buf_C)
                         call xgemmsg('T', 'N', m_size8, n_size8, k_size8, 1.0e0, &
                                      buf_A, lda1_8, buf_B, lda2_8, 0.0e0, buf_C, lda3_8)
!$acc end host_data

!$acc parallel loop collapse(2) present(buf_C, buf_C_fp64)
                         do n_loc = 1, n_chunk
                            do m_loc = 1, m_chunk
                               buf_C_fp64(m_loc, n_loc) = buf_C_fp64(m_loc, n_loc) + dble(buf_C(m_loc, n_loc))
                            end do
                         end do

                      end do ! end k_start
                   end do ! end irr_D1

                   s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2)); s_v3_m3 = s_v3(idx_m_in_v3(3))
                   s_v3_n1 = s_v3(idx_n_in_v3(1))

                   if (d_fact2 == 0.0d0) then
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                      do n_loc = 1, n_chunk
                          do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            m_1 = mod(m_flat - 1, m_dim1) + 1
                            rem = (m_flat - 1) / m_dim1
                            m_2 = mod(rem, m_dim2) + 1
                            m_3 = rem / m_dim2 + 1
                            n_1 = n_flat
                            idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
                            v3(idx_v3) = sngl(d_fact1 * buf_C_fp64(m_loc, n_loc))
                         end do
                      end do
                   else
!$acc parallel loop collapse(2) present(buf_C_fp64) deviceptr(v3) private(m_flat, n_flat, m_1, m_2, m_3, rem, n_1, idx_v3)
                      do n_loc = 1, n_chunk
                         do m_loc = 1, m_chunk
                            m_flat = m_start + m_loc - 1; n_flat = n_start + n_loc - 1
                            m_1 = mod(m_flat - 1, m_dim1) + 1
                            rem = (m_flat - 1) / m_dim1
                            m_2 = mod(rem, m_dim2) + 1
                            m_3 = rem / m_dim2 + 1
                            n_1 = n_flat
                            idx_v3 = int(base_v3_adj, 8) + int(n_1, 8)*s_v3_n1 + int(m_1, 8)*s_v3_m1 + int(m_2, 8)*s_v3_m2 + int(m_3, 8)*s_v3_m3
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
      end subroutine vv2v_dpd_424_sg
