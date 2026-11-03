! =======================================================================
! 终极精简版张量引擎 - CPU 单精度核心 (Single Precision, CPU)
! [核心能力] 424 模式: 4-index * 2-index -> 4-index 
! [优化特性] 拥有 26 条 Intact Pair 极速通道，原生支持 IOP=0/2
! [重大升级] 1. 补充了 ebij/beij 和 ijeb/ijbe 的极速路由，防 Fallback 崩溃。
!            2. 修复了 1253 (ijmb) 通道非连续步长 (LDC) 导致的交叉覆盖 Bug，
!               引入 fp_buf_C 保证内存块映射绝对安全。
!            3. 新增 1254/53,35 通道，支持 abmi*mj/jm -> abji full。
! =======================================================================
      subroutine vv2v_dpd_424_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                                v2, p2_1, p2_2, iop2_12, sym2, type2, &
                                v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1in, fact2in)
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
      integer :: old_threads, new_threads
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
      integer(c_int) :: old_threads, new_threads
#else
      integer :: old_threads, new_threads
#endif

      ! --- 接口参数声明 (单精度 real*4；外部 fact 保持 real*8 兼容前端) ---
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
      real*4                         :: fact1, fact2        

      external :: xgemmso
      integer, external :: irpdso  

      ! --- Fast Path 相关变量 ---
      logical :: is_fast_path, is_scatter
      logical :: is_intact12_path, is_intact34_to12_path
      character(len=1) :: transA, transB
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposrht1(33,33)
      integer :: iabsoff2(33)
      integer :: inumlft3(33), inumrht3(33), iabsoff3(33), iposrht3(33,33)
      integer :: irr_lft, irr_fast
      integer :: irr_12, irr_34, irr_1, irr_2, irr_3, irr_4, irr_5
      integer :: irr_lft1, irr_rht1, irr_lft3, irr_rht3
      integer :: dim_12, dim_34, dim_1, dim_2, dim_3, dim_4, dim_5
      integer :: off_v1, off_v2, off_v3, M, N, K, lda, ldb, ldc
      integer :: idx_A, i1, i2, i3, i4, i5, i12, i34
      integer :: cv1_3, cv1_4, cv2_1, cv2_2
      integer :: s1_3, s1_4, s2_1, s2_2
      integer(8) :: idx_C
      real*4, allocatable :: fp_buf_A(:), fp_buf_C(:)

      ! --- Fallback 引擎映射与调度变量 ---
      character(len=1) :: c_ext(4), c_dum(1)
      integer :: pos1(5), pos2(5), pos3(4) 
      integer :: i, j, idx, cnt_m, cnt_n
      logical :: found, skip_block
      integer :: idx_m_in_v3(3), idx_n_in_v3(1)
      integer :: pop1(4,8), pop2(2,8), pop3(4,8)
      integer :: iposlft1_fb(33,33), iposlft3_fb(33,33)
      
      integer :: num_tasks, task_list(3, 10000)
      integer :: num_threads, BS
      logical :: run_parallel
      integer(8) :: max_mem_bytes
      integer :: task_idx, irr_L4, irr_L3, irr_L2, irr_L1
      integer :: dim_L1, dim_L2, dim_L3, dim_L4, dim_X_arr(4)
      integer :: cur_M, cur_N, cur_K, cur_bs_m, cur_bs_n, cur_bs_k
      integer :: irr_phys1(4), irr_phys2(2), irr_phys3(4), irr_logi(5)
      integer :: irr_ext_v1, irr_ext_v2, irr_L1_tmp
      integer :: base_v3, s_phys1(4), s_phys2(2), s_phys3(4), s_v3(4), base_v3_adj
      integer :: d_D1(8), d_D2(8), irr_D1
      integer :: base_v1, base_v2, b_v1_adj, b_v2_adj
      integer :: s_v1_arr(5), s_v2_arr(5)
      integer :: m_dim1, m_dim2, m_dim3, n_dim1
      integer :: s_m1, s_m2, s_m3, s_n1, s_k_v1, s_k_v2
      integer :: s_v3_m1, s_v3_m2, s_v3_m3, s_v3_n1
      integer :: m1, m2, m3, rem, base_k, base_n, idx_v3
      integer :: m_start, n_start, k_start, m_global, n_global, k_global
      integer :: m_loc, n_loc, k_loc
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)

      fact1 = sngl(fact1in)
      fact2 = sngl(fact2in)

#ifdef USE_MKL
      old_threads = mkl_get_max_threads()
#elif defined(USE_OPENBLAS)
      old_threads = openblas_get_num_threads()
#else
      old_threads = 1
#endif

      ! 映射整合
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:)
      pop3(1,:) = p3_1(:); pop3(2,:) = p3_2(:); pop3(3,:) = p3_3(:); pop3(4,:) = p3_4(:)

      ! =============================================================
      ! 【原生构建 DPD 偏移表】
      ! =============================================================
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
      ! 【阶段 1：智能路由 (Fast Path Router)】
      ! =============================================================
      is_fast_path = .false.
      is_scatter   = .false.
      if (type3 == '1234') then
          if (type1(1:2) == '12') then
              if (type1 == '1235' .and. type2 == '54') then; is_fast_path = .true.; end if
              if (type1 == '1235' .and. type2 == '45') then; is_fast_path = .true.; end if
              ! 【修复】为 1253 配置 is_scatter = .true. 确保为其分配安全缓存 fp_buf_C
              if (type1 == '1253' .and. type2 == '54') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '1253' .and. type2 == '45') then; is_fast_path = .true.; is_scatter = .true.; end if
              ! 【新增】1254: abmi * mj/jm -> abji(full)，输出右 pair 反序
              if (type1 == '1254' .and. type2 == '53') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '1254' .and. type2 == '35') then; is_fast_path = .true.; is_scatter = .true.; end if
          else if (type1(3:4) == '12') then
              if (type1 == '3512' .and. type2 == '54') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3512' .and. type2 == '45') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '5312' .and. type2 == '54') then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '5312' .and. type2 == '45') then; is_fast_path = .true.; is_scatter = .true.; end if
          else if (type1(3:4) == '34') then
              ! 原有的 aeij 极速通道
              if (type1 == '1534' .and. (type2 == '52' .or. type2 == '25')) then; is_fast_path = .true.; end if
              if (type1 == '5134' .and. (type2 == '52' .or. type2 == '25')) then; is_fast_path = .true.; end if
              ! [新增] ebij / beij (i<j) * ae -> abij 的极速路由
              if (type1 == '2534' .and. (type2 == '15' .or. type2 == '51')) then; is_fast_path = .true.; end if
              if (type1 == '5234' .and. (type2 == '15' .or. type2 == '51')) then; is_fast_path = .true.; end if
          else if (type1(1:2) == '34') then
              ! 原有的 ijeb -> 极速散列通道
              if (type1 == '3415' .and. (type2 == '52' .or. type2 == '25')) then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3451' .and. (type2 == '52' .or. type2 == '25')) then; is_fast_path = .true.; is_scatter = .true.; end if
              ! [新增] ijbe / ijeb * ae -> abij 的极速路由
              if (type1 == '3425' .and. (type2 == '15' .or. type2 == '51')) then; is_fast_path = .true.; is_scatter = .true.; end if
              if (type1 == '3452' .and. (type2 == '15' .or. type2 == '51')) then; is_fast_path = .true.; is_scatter = .true.; end if
          end if
      end if

      ! >>> 激活极速通道 <<<
      if (is_fast_path) then
          num_threads = omp_get_max_threads()
          if (num_threads < 1) num_threads = 1
#ifdef USE_MKL
          call mkl_set_num_threads(num_threads)
#elif defined(USE_OPENBLAS)
          new_threads = int(num_threads, c_int)
          call openblas_set_num_threads(new_threads)
#endif
          ! 只要 is_scatter 为 true (包含 1253/ijmb)，就会分配缓存
          if (is_scatter) allocate(fp_buf_A(64000000), fp_buf_C(64000000))

          if (type1(1:2) == '12' .or. type1(3:4) == '12') then
              do irr_4 = 1, nirrep
                 do irr_3 = 1, nirrep
                    ! type2=54/45: contracted index 5 pairs with output index 4
                    ! type2=53/35: contracted index 5 pairs with output index 3
                    if (type2 == '54' .or. type2 == '45') then
                       irr_5 = dirprd(irr_4, sym2)
                    else if (type2 == '53' .or. type2 == '35') then
                       irr_5 = dirprd(irr_3, sym2)
                    else
                       cycle
                    end if

                    irr_12 = dirprd(sym3, dirprd(irr_3, irr_4))
                    dim_12 = inumlft3(irr_12)
                    dim_3  = p3_3(irr_3)
                    dim_4  = p3_4(irr_4)
                    
                    dim_5  = 0
                    if (type2 == '54') dim_5 = p2_1(irr_5)
                    if (type2 == '45') dim_5 = p2_2(irr_5)
                    if (type2 == '53') dim_5 = p2_1(irr_5)
                    if (type2 == '35') dim_5 = p2_2(irr_5)
                    if (dim_12 == 0 .or. dim_3 == 0 .or. dim_4 == 0 .or. dim_5 == 0) cycle
                    
                    irr_rht3 = dirprd(irr_3, irr_4)
                    off_v3 = iabsoff3(irr_rht3) + dim_12 * iposrht3(irr_4, irr_rht3)
                    
                    if (type2 == '54') off_v2 = iabsoff2(irr_4)
                    if (type2 == '45') off_v2 = iabsoff2(irr_5)
                    if (type2 == '53') off_v2 = iabsoff2(irr_3)
                    if (type2 == '35') off_v2 = iabsoff2(irr_5)

                    if (type1 == '1235' .or. type1 == '1253' .or. type1 == '1254') then
                       if (type1 == '1254') then
                          irr_rht1 = dirprd(irr_5, irr_4)
                       else
                          irr_rht1 = dirprd(irr_3, irr_5)
                       end if

                       off_v1 = iabsoff1(irr_rht1)
                       if (type1 == '1235') off_v1 = off_v1 + dim_12 * iposrht1(irr_5, irr_rht1)
                       if (type1 == '1253') off_v1 = off_v1 + dim_12 * iposrht1(irr_3, irr_rht1)
                       if (type1 == '1254') off_v1 = off_v1 + dim_12 * iposrht1(irr_4, irr_rht1)

                       if (type1 == '1235') then
                          ! 1235 模式: (12,3,5) * (5,4) -> (12,3,4)，可直接 xgemmso
                          M = dim_12 * dim_3; N = dim_4; K = dim_5; lda = M
                          if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                          call xgemmso('N', transB, M, N, K, fact1in, v1(off_v1), lda, v2(off_v2), ldb, fact2in, v3(off_v3), M)

                       else if (type1 == '1253') then
                          ! 1253 模式: (12,5,3) * (5,4) -> (12,3,4)
                          ! 右 pair 中 3 是慢变块，需要逐 i3 做 GEMM 后 scatter 回 full (3,4)。
                          M = dim_12; N = dim_4; K = dim_5; lda = M
                          if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if

                          do i3 = 0, dim_3 - 1
                             call xgemmso('N', transB, M, N, K, 1.0d0, v1(off_v1 + i3*M*K), lda, v2(off_v2), ldb, 0.0d0, fp_buf_C, M)

                             do i4 = 0, dim_4 - 1
                                do i12 = 0, dim_12 - 1
                                   idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                                   if (fact2 == 0.0d0) then
                                      v3(idx_C) = fact1 * fp_buf_C(i12 + i4*M + 1)
                                   else
                                      v3(idx_C) = fact2 * v3(idx_C) + fact1 * fp_buf_C(i12 + i4*M + 1)
                                   end if
                                end do
                             end do
                          end do

                       else if (type1 == '1254') then
                          ! 1254 模式: (12,5,4) * (5,3) -> (12,3,4)
                          ! 对应原始指标 abmi * mj -> abji(full)。
                          ! 这里固定输出第 4 指标 i，GEMM 生成第 3 指标 j，再 scatter 回 v3(12,3,4)。
                          M = dim_12; N = dim_3; K = dim_5; lda = M
                          if (type2 == '53') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if

                          do i4 = 0, dim_4 - 1
                             call xgemmso('N', transB, M, N, K, 1.0d0, v1(off_v1 + i4*M*K), lda, v2(off_v2), ldb, 0.0d0, fp_buf_C, M)

                             do i3 = 0, dim_3 - 1
                                do i12 = 0, dim_12 - 1
                                   idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                                   if (fact2 == 0.0d0) then
                                      v3(idx_C) = fact1 * fp_buf_C(i12 + i3*M + 1)
                                   else
                                      v3(idx_C) = fact2 * v3(idx_C) + fact1 * fp_buf_C(i12 + i3*M + 1)
                                   end if
                                end do
                             end do
                          end do
                       end if

                    else if (type1 == '3512' .or. type1 == '5312') then
                       irr_lft1 = dirprd(irr_3, irr_5)
                       off_v1 = iabsoff1(irr_12)
                       
                       if (type1 == '3512') off_v1 = off_v1 + iposlft1_fb(irr_5, irr_lft1)
                       if (type1 == '5312') off_v1 = off_v1 + iposlft1_fb(irr_3, irr_lft1)
                       
                       if (type1 == '3512') then; M = dim_3; K = dim_5; lda = dim_3; transA = 'N'; else; M = dim_3; K = dim_5; lda = dim_5; transA = 'T'; end if
                       N = dim_4
                       if (type2 == '54') then; ldb = K; transB = 'N'; else; ldb = N; transB = 'T'; end if
                       
                       do i12 = 0, dim_12 - 1
                          call xgemmso(transA, transB, M, N, K, 1.0d0, v1(off_v1 + i12*inumlft1(irr_lft1)), lda, v2(off_v2), ldb, 0.0d0, fp_buf_C, M)
                          do i4 = 0, dim_4 - 1
                             do i3 = 0, dim_3 - 1
                                idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                                if (fact2 == 0.0d0) then
                                   v3(idx_C) = fact1 * fp_buf_C(i3 + i4*dim_3 + 1)
                                else
                                   v3(idx_C) = fact2 * v3(idx_C) + fact1 * fp_buf_C(i3 + i4*dim_3 + 1)
                                end if
                             end do
                          end do
                       end do
                    end if
                 end do
              end do
              
          else
              ! [重构] 为 ebij / beij 增加专属动态路由环
              do irr_1 = 1, nirrep
                 do irr_2 = 1, nirrep
                    if (type2 == '52' .or. type2 == '25') then
                        irr_5 = dirprd(irr_2, sym2)
                    else if (type2 == '15' .or. type2 == '51') then
                        irr_5 = dirprd(irr_1, sym2)
                    end if
                    
                    irr_12 = dirprd(irr_1, irr_2)
                    irr_34 = dirprd(sym3, irr_12)
                    
                    dim_34 = inumrht3(irr_34)
                    dim_1  = p3_1(irr_1)
                    dim_2  = p3_2(irr_2)
                    
                    dim_5 = 0
                    if (type2 == '52' .or. type2 == '51') dim_5 = p2_1(irr_5)
                    if (type2 == '25' .or. type2 == '15') dim_5 = p2_2(irr_5)
                    if (dim_34 == 0 .or. dim_1 == 0 .or. dim_2 == 0 .or. dim_5 == 0) cycle
                    
                    off_v3 = iabsoff3(irr_34) + iposlft3_fb(irr_2, irr_12)
                    
                    if (type2 == '52') off_v2 = iabsoff2(irr_2)
                    if (type2 == '25') off_v2 = iabsoff2(irr_5)
                    if (type2 == '15') off_v2 = iabsoff2(irr_5)
                    if (type2 == '51') off_v2 = iabsoff2(irr_1)
                    
                    if (type1(3:4) == '34') then
                       if (type1 == '1534' .or. type1 == '5134') then
                           irr_lft1 = dirprd(irr_1, irr_5)
                           off_v1 = iabsoff1(irr_34)
                           
                           if (type1 == '1534') then
                               off_v1 = off_v1 + iposlft1_fb(irr_5, irr_lft1)
                               M = dim_1; K = dim_5; N = dim_2; lda = dim_1; transA = 'N'
                           else
                               off_v1 = off_v1 + iposlft1_fb(irr_1, irr_lft1)
                               M = dim_1; K = dim_5; N = dim_2; lda = dim_5; transA = 'T'
                           end if
                           if (type2 == '52') then; ldb = dim_5; transB = 'N'; else; ldb = dim_2; transB = 'T'; end if
                           
                           do i34 = 0, dim_34 - 1
                              call xgemmso(transA, transB, M, N, K, fact1in, v1(off_v1 + i34*inumlft1(irr_lft1)), lda, v2(off_v2), ldb, fact2in, v3(off_v3 + i34*inumlft3(irr_12)), dim_1)
                           end do
                           
                       ! [新增] 处理 ebij / beij
                       else if (type1 == '2534' .or. type1 == '5234') then
                           irr_lft1 = dirprd(irr_2, irr_5)
                           off_v1 = iabsoff1(irr_34)
                           
                           if (type1 == '2534') then
                               off_v1 = off_v1 + iposlft1_fb(irr_5, irr_lft1)
                           else
                               off_v1 = off_v1 + iposlft1_fb(irr_2, irr_lft1)
                           end if
                           
                           M = dim_1; N = dim_2; K = dim_5
                           if (type2 == '15') then; lda = dim_1; transA = 'N'; else; lda = dim_5; transA = 'T'; end if
                           if (type1 == '2534') then; ldb = dim_2; transB = 'T'; else; ldb = dim_5; transB = 'N'; end if
                           
                           do i34 = 0, dim_34 - 1
                              call xgemmso(transA, transB, M, N, K, fact1in, v2(off_v2), lda, v1(off_v1 + i34*inumlft1(irr_lft1)), ldb, fact2in, v3(off_v3 + i34*inumlft3(irr_12)), M)
                           end do
                       end if
                       
                    else if (type1(1:2) == '34') then
                       if (type1 == '3415' .or. type1 == '3451') then
                           irr_rht1 = dirprd(irr_1, irr_5)
                           off_v1 = iabsoff1(irr_rht1)
                           if (type1 == '3415') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_5, irr_rht1)
                           if (type1 == '3451') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_1, irr_rht1)
                           if (type2 == '52') then; ldb = dim_5; transB = 'N'; else; ldb = dim_2; transB = 'T'; end if
                           
                           do i34 = 0, dim_34 - 1
                              idx_A = 1
                              if (type1 == '3415') then
                                 do i5 = 0, dim_5 - 1; do i1 = 0, dim_1 - 1
                                    fp_buf_A(idx_A) = v1(off_v1 + i34 + i1*inumlft1(irr_34) + i5*inumlft1(irr_34)*dim_1); idx_A = idx_A + 1
                                 end do; end do
                                 M = dim_1; K = dim_5; lda = dim_1; transA = 'N'
                              else
                                 do i1 = 0, dim_1 - 1; do i5 = 0, dim_5 - 1
                                    fp_buf_A(idx_A) = v1(off_v1 + i34 + i5*inumlft1(irr_34) + i1*inumlft1(irr_34)*dim_5); idx_A = idx_A + 1
                                 end do; end do
                                 M = dim_1; K = dim_5; lda = dim_5; transA = 'T'
                              end if
                              call xgemmso(transA, transB, M, dim_2, K, fact1in, fp_buf_A, lda, v2(off_v2), ldb, fact2in, v3(off_v3 + i34*inumlft3(irr_12)), dim_1)
                           end do
                           
                       ! [新增] 处理 ijeb / ijbe
                       else if (type1 == '3425' .or. type1 == '3452') then
                           irr_rht1 = dirprd(irr_2, irr_5)
                           off_v1 = iabsoff1(irr_rht1)
                           if (type1 == '3425') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_5, irr_rht1)
                           if (type1 == '3452') off_v1 = off_v1 + inumlft1(irr_34) * iposrht1(irr_2, irr_rht1)
                           
                           M = dim_1; N = dim_2; K = dim_5
                           if (type2 == '15') then; lda = dim_1; transA = 'N'; else; lda = dim_5; transA = 'T'; end if
                           if (type1 == '3425') then; transB = 'T'; ldb = dim_2; else; transB = 'N'; ldb = dim_5; end if
                           
                           do i34 = 0, dim_34 - 1
                              idx_A = 1
                              if (type1 == '3425') then
                                 do i5 = 0, dim_5 - 1; do i2 = 0, dim_2 - 1
                                    fp_buf_A(idx_A) = v1(off_v1 + i34 + i2*inumlft1(irr_34) + i5*inumlft1(irr_34)*dim_2); idx_A = idx_A + 1
                                 end do; end do
                              else
                                 do i2 = 0, dim_2 - 1; do i5 = 0, dim_5 - 1
                                    fp_buf_A(idx_A) = v1(off_v1 + i34 + i5*inumlft1(irr_34) + i2*inumlft1(irr_34)*dim_5); idx_A = idx_A + 1
                                 end do; end do
                              end if
                              call xgemmso(transA, transB, M, N, K, fact1in, v2(off_v2), lda, fp_buf_A, ldb, fact2in, v3(off_v3 + i34*inumlft3(irr_12)), dim_1)
                           end do
                       end if
                    end if
                 end do
              end do
          end if
          
          if (is_scatter) deallocate(fp_buf_A, fp_buf_C)
#ifdef USE_MKL
          call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
          call openblas_set_num_threads(old_threads)
#endif
          return
      end if

      ! =============================================================
      ! 【阶段 2：Fallback 全动态瓦片装配引擎】
      ! =============================================================
      do i = 1, 4; c_ext(i) = type3(i:i); end do
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
      
      do j = 1, 4; if (type1(j:j) == c_dum(1)) pos1(5) = j; end do
      do j = 1, 2; if (type2(j:j) == c_dum(1)) pos2(5) = j; end do

      cnt_m = 0; cnt_n = 0; idx_m_in_v3 = 0; idx_n_in_v3 = 0
      do i = 1, 4
         if (pos1(i) /= 0) then; cnt_m = cnt_m + 1; idx_m_in_v3(cnt_m) = i; end if
         if (pos2(i) /= 0) then; cnt_n = cnt_n + 1; idx_n_in_v3(cnt_n) = i; end if
      end do

      if (cnt_m /= 3 .or. cnt_n /= 1) stop "FATAL ERROR: 424 kernel expects precisely 1 contracted index."

      ! =============================================================
      ! 【阶段 2A：Intact-left-pair 专用正确性通道】
      !
      ! 目标修复：type1='abim', type2='mj', type3='abji' 这类
      ! a<b/a<=b 压缩 pair 完整保留在左 pair，而右 pair 发生 ij/ji
      ! 重排的 424 情况。
      !
      ! 旧 fallback 会把 packed left pair 拆成 full a,b 双循环，导致
      ! iop=0/2 时写入位置错误。这里把 left pair 当作一个 opaque
      ! packed 维度 i12，只在右 pair 和 dummy 指标上显式循环。
      ! =============================================================
      is_intact12_path = .false.
      if (type1(1:2) == type3(1:2)) then
         if (iop1_12 == iop3_12 .and. iop1_34 == 1 .and. iop2_12 == 1 .and. iop3_34 == 1) then
            if (pos1(1) == 1 .and. pos1(2) == 2 .and. pos1(5) /= 0 .and. pos2(5) /= 0) then
               is_intact12_path = .true.
            end if
         end if
      end if

      if (is_intact12_path) then
#ifdef USE_MKL
         call mkl_set_num_threads(1)
#elif defined(USE_OPENBLAS)
         new_threads = 1_c_int
         call openblas_set_num_threads(new_threads)
#endif

         do irr_rht3 = 1, nirrep
            irr_lft3 = dirprd(irr_rht3, sym3)
            do irr_4 = 1, nirrep
               irr_3 = dirprd(irr_rht3, irr_4)

               dim_12 = inumlft3(irr_lft3)
               dim_3  = p3_3(irr_3)
               dim_4  = p3_4(irr_4)
               if (dim_12 == 0 .or. dim_3 == 0 .or. dim_4 == 0) cycle

               off_v3 = iabsoff3(irr_rht3) + dim_12 * iposrht3(irr_4, irr_rht3)

               do i4 = 0, dim_4 - 1
                  do i3 = 0, dim_3 - 1
                     do i12 = 0, dim_12 - 1
                        idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                        if (fact2 == 0.0d0) then
                           v3(idx_C) = 0.0d0
                        else
                           v3(idx_C) = fact2 * v3(idx_C)
                        end if
                     end do
                  end do
               end do

               do irr_5 = 1, nirrep
                  ! Build physical irreps for V2 from type2 and output type3.
                  do i = 1, 2
                     if (type2(i:i) == c_dum(1)) then
                        irr_phys2(i) = irr_5
                     else if (type2(i:i) == type3(3:3)) then
                        irr_phys2(i) = irr_3
                     else if (type2(i:i) == type3(4:4)) then
                        irr_phys2(i) = irr_4
                     else
                        irr_phys2(i) = 0
                     end if
                  end do
                  if (irr_phys2(1) == 0 .or. irr_phys2(2) == 0) cycle
                  if (dirprd(irr_phys2(1), irr_phys2(2)) /= sym2) cycle

                  ! Build physical irreps for the right pair of V1.
                  do i = 3, 4
                     if (type1(i:i) == c_dum(1)) then
                        irr_phys1(i) = irr_5
                     else if (type1(i:i) == type3(3:3)) then
                        irr_phys1(i) = irr_3
                     else if (type1(i:i) == type3(4:4)) then
                        irr_phys1(i) = irr_4
                     else
                        irr_phys1(i) = 0
                     end if
                  end do
                  if (irr_phys1(3) == 0 .or. irr_phys1(4) == 0) cycle

                  irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                  irr_lft1 = dirprd(irr_rht1, sym1)
                  if (irr_lft1 /= irr_lft3) cycle

                  dim_5 = pop1(pos1(5), irr_5)
                  if (dim_5 == 0) cycle
                  if (pop2(pos2(5), irr_5) /= dim_5) cycle

                  off_v1 = iabsoff1(irr_rht1) + dim_12 * iposrht1(irr_phys1(4), irr_rht1)
                  off_v2 = iabsoff2(irr_phys2(2))

                  s1_3 = dim_12
                  s1_4 = dim_12 * pop1(3, irr_phys1(3))
                  s2_1 = 1
                  s2_2 = pop2(1, irr_phys2(1))

                  do i5 = 0, dim_5 - 1
                     do i4 = 0, dim_4 - 1
                        do i3 = 0, dim_3 - 1
                           ! Coordinates for V1 right pair positions 3 and 4.
                           if (type1(3:3) == c_dum(1)) then
                              cv1_3 = i5
                           else if (type1(3:3) == type3(3:3)) then
                              cv1_3 = i3
                           else
                              cv1_3 = i4
                           end if
                           if (type1(4:4) == c_dum(1)) then
                              cv1_4 = i5
                           else if (type1(4:4) == type3(3:3)) then
                              cv1_4 = i3
                           else
                              cv1_4 = i4
                           end if

                           ! Coordinates for V2 positions 1 and 2.
                           if (type2(1:1) == c_dum(1)) then
                              cv2_1 = i5
                           else if (type2(1:1) == type3(3:3)) then
                              cv2_1 = i3
                           else
                              cv2_1 = i4
                           end if
                           if (type2(2:2) == c_dum(1)) then
                              cv2_2 = i5
                           else if (type2(2:2) == type3(3:3)) then
                              cv2_2 = i3
                           else
                              cv2_2 = i4
                           end if

                           do i12 = 0, dim_12 - 1
                              idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                              v3(idx_C) = v3(idx_C) + fact1 * &
                                   v1(off_v1 + i12 + cv1_3*s1_3 + cv1_4*s1_4) * &
                                   v2(off_v2 + cv2_1*s2_1 + cv2_2*s2_2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

#ifdef USE_MKL
         call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
         call openblas_set_num_threads(old_threads)
#endif
         return
      end if


      ! =============================================================
      ! 【阶段 2B：Intact-right-pair -> left-pair 专用正确性通道】
      !
      ! 目标修复：type1='miab', type2='jm', type3='abji' 这类
      ! a<b/a<=b 压缩 pair 在输入 V1 的右 pair，但在输出 V3
      ! 的左 pair。归一化后典型为：
      !      V1(5,4,1,2) * V2(3,5) -> V3(1,2,3,4)
      ! 即 P3-04: miab * jm -> abji。
      !
      ! 旧 fallback 会把输出左侧 packed pair 当作 full a,b 双循环。
      ! 这里把 V1 右 pair / V3 左 pair 当作同一个 opaque packed
      ! 维度 i12，只显式循环 V3 的右 pair 和 dummy 收缩指标。
      ! =============================================================
      is_intact34_to12_path = .false.
      if (type1(3:4) == type3(1:2)) then
         if (iop1_34 == iop3_12 .and. iop1_12 == 1 .and. iop2_12 == 1 .and. iop3_34 == 1) then
            if (pos1(1) == 3 .and. pos1(2) == 4 .and. pos1(5) /= 0 .and. pos2(5) /= 0) then
               is_intact34_to12_path = .true.
            end if
         end if
      end if

      if (is_intact34_to12_path) then
#ifdef USE_MKL
         call mkl_set_num_threads(1)
#elif defined(USE_OPENBLAS)
         new_threads = 1_c_int
         call openblas_set_num_threads(new_threads)
#endif

         do irr_rht3 = 1, nirrep
            irr_lft3 = dirprd(irr_rht3, sym3)
            do irr_4 = 1, nirrep
               irr_3 = dirprd(irr_rht3, irr_4)

               dim_12 = inumlft3(irr_lft3)
               dim_3  = p3_3(irr_3)
               dim_4  = p3_4(irr_4)
               if (dim_12 == 0 .or. dim_3 == 0 .or. dim_4 == 0) cycle

               off_v3 = iabsoff3(irr_rht3) + dim_12 * iposrht3(irr_4, irr_rht3)

               do i4 = 0, dim_4 - 1
                  do i3 = 0, dim_3 - 1
                     do i12 = 0, dim_12 - 1
                        idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                        if (fact2 == 0.0d0) then
                           v3(idx_C) = 0.0d0
                        else
                           v3(idx_C) = fact2 * v3(idx_C)
                        end if
                     end do
                  end do
               end do

               do irr_5 = 1, nirrep
                  ! Build physical irreps for V2 from type2 and output type3.
                  do i = 1, 2
                     if (type2(i:i) == c_dum(1)) then
                        irr_phys2(i) = irr_5
                     else if (type2(i:i) == type3(3:3)) then
                        irr_phys2(i) = irr_3
                     else if (type2(i:i) == type3(4:4)) then
                        irr_phys2(i) = irr_4
                     else
                        irr_phys2(i) = 0
                     end if
                  end do
                  if (irr_phys2(1) == 0 .or. irr_phys2(2) == 0) cycle
                  if (dirprd(irr_phys2(1), irr_phys2(2)) /= sym2) cycle

                  ! Build physical irreps for V1 left pair; V1 right pair is the intact packed 12 pair.
                  do i = 1, 2
                     if (type1(i:i) == c_dum(1)) then
                        irr_phys1(i) = irr_5
                     else if (type1(i:i) == type3(3:3)) then
                        irr_phys1(i) = irr_3
                     else if (type1(i:i) == type3(4:4)) then
                        irr_phys1(i) = irr_4
                     else
                        irr_phys1(i) = 0
                     end if
                  end do
                  if (irr_phys1(1) == 0 .or. irr_phys1(2) == 0) cycle

                  irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2))
                  irr_rht1 = irr_lft3
                  if (dirprd(irr_lft1, irr_rht1) /= sym1) cycle

                  dim_5 = pop1(pos1(5), irr_5)
                  if (dim_5 == 0) cycle
                  if (pop2(pos2(5), irr_5) /= dim_5) cycle

                  off_v1 = iabsoff1(irr_rht1) + iposlft1_fb(irr_phys1(2), irr_lft1)
                  off_v2 = iabsoff2(irr_phys2(2))

                  s1_3 = 1
                  s1_4 = pop1(1, irr_phys1(1))
                  s2_1 = 1
                  s2_2 = pop2(1, irr_phys2(1))

                  do i5 = 0, dim_5 - 1
                     do i4 = 0, dim_4 - 1
                        do i3 = 0, dim_3 - 1
                           ! Coordinates for V1 left pair positions 1 and 2.
                           if (type1(1:1) == c_dum(1)) then
                              cv1_3 = i5
                           else if (type1(1:1) == type3(3:3)) then
                              cv1_3 = i3
                           else
                              cv1_3 = i4
                           end if
                           if (type1(2:2) == c_dum(1)) then
                              cv1_4 = i5
                           else if (type1(2:2) == type3(3:3)) then
                              cv1_4 = i3
                           else
                              cv1_4 = i4
                           end if

                           ! Coordinates for V2 positions 1 and 2.
                           if (type2(1:1) == c_dum(1)) then
                              cv2_1 = i5
                           else if (type2(1:1) == type3(3:3)) then
                              cv2_1 = i3
                           else
                              cv2_1 = i4
                           end if
                           if (type2(2:2) == c_dum(1)) then
                              cv2_2 = i5
                           else if (type2(2:2) == type3(3:3)) then
                              cv2_2 = i3
                           else
                              cv2_2 = i4
                           end if

                           do i12 = 0, dim_12 - 1
                              idx_C = off_v3 + i12 + i3*dim_12 + i4*dim_12*dim_3
                              v3(idx_C) = v3(idx_C) + fact1 * &
                                   v1(off_v1 + cv1_3*s1_3 + cv1_4*s1_4 + i12*inumlft1(irr_lft1)) * &
                                   v2(off_v2 + cv2_1*s2_1 + cv2_2*s2_2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

#ifdef USE_MKL
         call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
         call openblas_set_num_threads(old_threads)
#endif
         return
      end if

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

      num_threads = omp_get_max_threads()
      if (num_threads < 1) num_threads = 1
      if (num_tasks >= num_threads * 2) then
         run_parallel = .true.
         max_mem_bytes = 3500000000_8 / num_threads 
         BS = int(sqrt(real(max_mem_bytes, 8) / 24.0d0))
#ifdef USE_MKL
         call mkl_set_num_threads(1)
#elif defined(USE_OPENBLAS)
         new_threads = 1_c_int
         call openblas_set_num_threads(new_threads)
#endif
      else
         run_parallel = .false.
         max_mem_bytes = 3500000000_8 
         BS = int(sqrt(real(max_mem_bytes, 8) / 24.0d0))
#ifdef USE_MKL
         call mkl_set_num_threads(num_threads)
#elif defined(USE_OPENBLAS)
         new_threads = int(num_threads, c_int)
         call openblas_set_num_threads(new_threads)
#endif
      end if
      
      if (BS > 12000) BS = 12000
      if (BS < 256) BS = 256

!$omp parallel if(run_parallel) private(task_idx, irr_rht3, irr_lft3, irr_L4, irr_L3, irr_L2, irr_L1, &
!$omp dim_L1, dim_L2, dim_L3, dim_L4, dim_X_arr, cur_M, cur_N, cur_K, &
!$omp irr_phys1, irr_phys2, irr_phys3, base_v3, s_phys1, s_phys2, s_phys3, s_v3, base_v3_adj, &
!$omp irr_D1, irr_ext_v2, irr_ext_v1, irr_logi, skip_block, i, j, &
!$omp irr_lft1, irr_rht1, base_v1, base_v2, b_v1_adj, b_v2_adj, &
!$omp s_v1_arr, s_v2_arr, buf_A, buf_B, buf_C, m_dim1, m_dim2, m_dim3, n_dim1, &
!$omp s_m1, s_m2, s_m3, s_n1, s_k_v1, s_k_v2, s_v3_m1, s_v3_m2, s_v3_m3, s_v3_n1, &
!$omp m1, m2, m3, rem, base_k, base_n, idx_v3, &
!$omp m_start, n_start, k_start, cur_bs_m, cur_bs_n, cur_bs_k, &
!$omp m_global, n_global, k_global, m_loc, n_loc, k_loc)

      allocate(buf_A(BS, BS), buf_B(BS, BS), buf_C(BS, BS))

!$omp do schedule(dynamic, 1)
      do task_idx = 1, num_tasks
         irr_rht3 = task_list(1, task_idx)
         irr_L4   = task_list(2, task_idx)
         irr_L2   = task_list(3, task_idx)

         irr_lft3 = dirprd(irr_rht3, sym3)
         irr_L3   = dirprd(irr_rht3, irr_L4)
         irr_L1   = dirprd(irr_lft3, irr_L2)    

         irr_phys3(pos3(1)) = irr_L1; irr_phys3(pos3(2)) = irr_L2
         irr_phys3(pos3(3)) = irr_L3; irr_phys3(pos3(4)) = irr_L4

         dim_X_arr(1) = pop3(pos3(1), irr_L1); dim_X_arr(2) = pop3(pos3(2), irr_L2)
         dim_X_arr(3) = pop3(pos3(3), irr_L3); dim_X_arr(4) = pop3(pos3(4), irr_L4)
         
         cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2)) * dim_X_arr(idx_m_in_v3(3))
         cur_N = dim_X_arr(idx_n_in_v3(1))

         base_v3 = iabsoff3(irr_rht3) + inumlft3(irr_lft3) * iposrht3(irr_phys3(4), irr_rht3) + iposlft3_fb(irr_phys3(2), irr_lft3)
         
         s_phys3(1) = 1; s_phys3(2) = pop3(1, irr_phys3(1))
         s_phys3(3) = inumlft3(irr_lft3); s_phys3(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))

         do i = 1, 4; s_v3(i) = s_phys3(pos3(i)); end do
         
         b_v1_adj = 0 
         b_v2_adj = 0 
         
         base_v3_adj = base_v3 - sum(s_v3(1:4)) 

         m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2)); m_dim3 = dim_X_arr(idx_m_in_v3(3))
         n_dim1 = dim_X_arr(idx_n_in_v3(1))

         do n_start = 1, cur_N, BS
            cur_bs_n = min(n_start + BS - 1, cur_N) - n_start + 1

            do m_start = 1, cur_M, BS
               cur_bs_m = min(m_start + BS - 1, cur_M) - m_start + 1

               do j = 1, cur_bs_n; do i = 1, cur_bs_m; buf_C(i, j) = 0.0d0; end do; end do

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
                     if (pos1(i) /= 0) irr_phys1(pos1(i)) = irr_logi(i)
                     if (pos2(i) /= 0) irr_phys2(pos2(i)) = irr_logi(i)
                  end do

                  irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2))
                  irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                  
                  base_v1 = iabsoff1(irr_rht1) + inumlft1(irr_lft1) * iposrht1(irr_phys1(4), irr_rht1) + iposlft1_fb(irr_phys1(2), irr_lft1)
                  base_v2 = iabsoff2(irr_phys2(2)) 
                  
                  s_phys1(1) = 1; s_phys1(2) = pop1(1, irr_phys1(1))
                  s_phys1(3) = inumlft1(irr_lft1); s_phys1(4) = inumlft1(irr_lft1) * pop1(3, irr_phys1(3))
                  s_v1_arr(:) = 0
                  do i = 1, 5; if(pos1(i)/=0) s_v1_arr(i) = s_phys1(pos1(i)); end do

                  s_phys2(1) = 1; s_phys2(2) = pop2(1, irr_phys2(1))
                  s_v2_arr(:) = 0
                  do i = 1, 5; if(pos2(i)/=0) s_v2_arr(i) = s_phys2(pos2(i)); end do

                  b_v1_adj = base_v1 - sum(s_v1_arr) 
                  b_v2_adj = base_v2 - sum(s_v2_arr) 

                  s_m1 = s_v1_arr(idx_m_in_v3(1)); s_m2 = s_v1_arr(idx_m_in_v3(2)); s_m3 = s_v1_arr(idx_m_in_v3(3))
                  s_n1 = s_v2_arr(idx_n_in_v3(1))
                  s_k_v1 = s_v1_arr(5); s_k_v2 = s_v2_arr(5)

                  do k_start = 1, cur_K, BS
                     cur_bs_k = min(k_start + BS - 1, cur_K) - k_start + 1

                     do k_loc = 1, cur_bs_k
                        k_global = k_start + k_loc - 1
                        do m_loc = 1, cur_bs_m
                           m_global = m_start + m_loc - 1
                           m3 = (m_global - 1) / (m_dim1 * m_dim2) + 1
                           rem = mod(m_global - 1, m_dim1 * m_dim2)
                           m2 = rem / m_dim1 + 1; m1 = mod(rem, m_dim1) + 1
                           buf_A(m_loc, k_loc) = v1(b_v1_adj + m1*s_m1 + m2*s_m2 + m3*s_m3 + k_global*s_k_v1)
                        end do
                     end do

                     do n_loc = 1, cur_bs_n
                        n_global = n_start + n_loc - 1
                        do k_loc = 1, cur_bs_k
                           k_global = k_start + k_loc - 1
                           buf_B(k_loc, n_loc) = v2(b_v2_adj + n_global*s_n1 + k_global*s_k_v2)
                        end do
                     end do

                     call xgemmso('N', 'N', cur_bs_m, cur_bs_n, cur_bs_k, 1.0d0, buf_A, BS, buf_B, BS, 1.0d0, buf_C, BS)
                  end do 
               end do 

               s_v3_m1 = s_v3(idx_m_in_v3(1)); s_v3_m2 = s_v3(idx_m_in_v3(2)); s_v3_m3 = s_v3(idx_m_in_v3(3))
               s_v3_n1 = s_v3(idx_n_in_v3(1))

               do n_loc = 1, cur_bs_n
                  n_global = n_start + n_loc - 1
                  do m_loc = 1, cur_bs_m
                     m_global = m_start + m_loc - 1
                     m3 = (m_global - 1) / (m_dim1 * m_dim2) + 1
                     rem = mod(m_global - 1, m_dim1 * m_dim2)
                     m2 = rem / m_dim1 + 1; m1 = mod(rem, m_dim1) + 1
                     
                     idx_v3 = base_v3_adj + n_global*s_v3_n1 + m1*s_v3_m1 + m2*s_v3_m2 + m3*s_v3_m3
                     
                     if (fact2 == 0.0d0) then
                        v3(idx_v3) = fact1 * buf_C(m_loc, n_loc)
                     else
                        v3(idx_v3) = fact2 * v3(idx_v3) + fact1 * buf_C(m_loc, n_loc)
                     end if
                  end do
               end do

            end do 
         end do 
      end do 
!$omp end do
      deallocate(buf_A, buf_B, buf_C)
!$omp end parallel

#ifdef USE_MKL
      call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
      call openblas_set_num_threads(old_threads)
#endif

      return
      end subroutine vv2v_dpd_424_s



