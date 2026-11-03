! =======================================================================
! 终极精简版张量引擎 - CPU 单精度核心 (Single Precision, CPU)
! =======================================================================
      subroutine vv2v_dpd_444_s(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                               v2, p2_1, p2_2, iop2_12, p2_3, p2_4, iop2_34, sym2, type2, &
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

      ! --- 接口参数声明 (单精度 real*4) ---
      character(len=4), intent(in)  :: type1, type2, type3 
      real*4, intent(in)            :: v1(*), v2(*)        
      real*4, intent(inout)         :: v3(*)               
      integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
      integer, intent(in)           :: p2_1(8), p2_2(8), p2_3(8), p2_4(8)
      integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
      integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop2_34, iop3_12, iop3_34
      integer, intent(in)           :: sym1, sym2, sym3    
      real*8, intent(in)            :: fact1in, fact2in        

      integer, external :: irpdso  

      ! --- Fast Path (智能路由) 变量 ---
      logical :: is_fast_path
      character(len=1) :: trans1, trans2
      integer :: m_size, n_size, k_size, lda1, lda2
      integer :: inumpq(33), inumrs(33), inummn(33)
      integer :: iabsoff1_fp(33), iabsoff2_fp(33), iabsoff3_fp(33)
      integer :: irr_pq, irr_rs, irr_mn, irr_L1_fp, irr_L2_fp, irr_L3_fp
      integer :: iop_mn
      integer :: ioff1, ioff2, ioff3

      ! --- Fallback 引擎映射与调度变量 ---
      character(len=1) :: c_ext(4), c_dum(2)
      integer :: pos1(6), pos2(6), pos3(6) 
      integer :: i, j, idx, cnt_m, cnt_n
      logical :: found
      integer :: idx_m_in_v3(2), idx_n_in_v3(2)
      integer :: pop1(4,8), pop2(4,8), pop3(4,8)
      integer :: inumlft1(33), inumrht1(33), iabsoff1(33), iposlft1(33,33), iposrht1(33,33)
      integer :: inumlft2(33), inumrht2(33), iabsoff2(33), iposlft2(33,33), iposrht2(33,33)
      integer :: inumlft3(33), inumrht3(33), iabsoff3(33), iposlft3(33,33), iposrht3(33,33)
      
      integer :: num_tasks
      integer :: task_list(3, 10000)
      integer :: num_threads, BS
      logical :: run_parallel
      integer(8) :: max_mem_bytes

      ! --- Fallback 线程私有计算变量 ---
      integer :: task_idx, irr_rht3, irr_lft3, irr_L4, irr_L3, irr_L2, irr_L1
      integer :: dim_L1, dim_L2, dim_L3, dim_L4, dim_X_arr(4)
      integer :: cur_M, cur_N, cur_K, cur_irr_D1
      integer :: irr_phys1(4), irr_phys2(4), irr_phys3(4), irr_logi(6)
      integer :: base_v3, s_phys(4), s_v3(6), base_v3_adj
      logical :: skip_block
      integer :: d_D1(8), d_D2(8), irr_D1, irr_D2, irr_lft1, irr_rht1, irr_lft2, irr_rht2
      integer :: base_v1, base_v2, b_v1_adj(8), b_v2_adj(8)
      integer :: s_v1_arr(6,8), s_v2_arr(6,8)
      integer :: m_dim1, m_dim2, n_dim1, n_dim2, k_dim1, k_dim2
      integer :: s_m1, s_m2, s_n1, s_n2, s_k1_v1, s_k2_v1, s_k1_v2, s_k2_v2
      integer :: k1, k2, m1, m2, n1, n2, base_k, base_n, base_n_v3, idx_v3
      integer :: m_start, n_start, k_start, k_end
      integer :: cur_bs_m, cur_bs_n, cur_bs_k
      integer :: m_global, n_global, k_global
      integer :: m_loc, n_loc, k_loc
      real*4  :: fact1, fact2
       
      ! 缓冲数组 (单精度 real*4)
      real*4, allocatable :: buf_A(:,:), buf_B(:,:), buf_C(:,:)

      ! =============================================================
      ! 初始系统环境备份
      ! =============================================================
#ifdef USE_MKL
      old_threads = mkl_get_max_threads()
#elif defined(USE_OPENBLAS)
      old_threads = openblas_get_num_threads()
#else
      old_threads = 1
#endif
      
      ! 安全截断：将外部传入的双精度转为单精度，用于局部的纯四字节计算
      fact1 = sngl(fact1in)
      fact2 = sngl(fact2in) 

      ! =============================================================
      ! 【阶段 1：智能路由探测器 (Smart Pattern Router)】
      ! =============================================================
      is_fast_path = .false.

      if (type1 == '3412' .and. type2 == '3456' .and. type3 == '1256') then
          if (iop1_12 == iop2_12 .and. iop1_34 == iop3_12 .and. iop2_34 == iop3_34) then
              is_fast_path = .true.; trans1 = 'T'; trans2 = 'N'
              iop_mn = iop1_12
          end if
      else if (type1 == '1234' .and. type2 == '3456' .and. type3 == '1256') then
          if (iop1_34 == iop2_12 .and. iop1_12 == iop3_12 .and. iop2_34 == iop3_34) then
              is_fast_path = .true.; trans1 = 'N'; trans2 = 'N'
              iop_mn = iop1_34
          end if
      else if (type1 == '1234' .and. type2 == '5634' .and. type3 == '1256') then
          if (iop1_34 == iop2_34 .and. iop1_12 == iop3_12 .and. iop2_12 == iop3_34) then
              is_fast_path = .true.; trans1 = 'N'; trans2 = 'T'
              iop_mn = iop1_34
          end if
      else if (type1 == '3412' .and. type2 == '5634' .and. type3 == '1256') then
          if (iop1_12 == iop2_34 .and. iop1_34 == iop3_12 .and. iop2_12 == iop3_34) then
              is_fast_path = .true.; trans1 = 'T'; trans2 = 'T'
              iop_mn = iop1_12
          end if
      end if

      ! =============================================================
      ! >>> Fast Path 极速计算分支 <<<
      ! =============================================================
      if (is_fast_path) then
          num_threads = omp_get_max_threads()
          if (num_threads < 1) num_threads = 1
#ifdef USE_MKL
          call mkl_set_num_threads(num_threads)
#elif defined(USE_OPENBLAS)
          new_threads = int(num_threads, c_int)
          call openblas_set_num_threads(new_threads)
#endif
          
          inumpq = 0; inumrs = 0; inummn = 0
          do i = 1, nirrep
              inumpq(i) = irpdso(i, p3_1, p3_2, iop3_12)
              inumrs(i) = irpdso(i, p3_3, p3_4, iop3_34)
              if (trans1 == 'T') then
                  inummn(i) = irpdso(i, p1_1, p1_2, iop_mn)
              else
                  inummn(i) = irpdso(i, p1_3, p1_4, iop_mn)
              end if
          end do

          iabsoff1_fp(1) = 1; iabsoff2_fp(1) = 1; iabsoff3_fp(1) = 1
          do i = 1, nirrep - 1
              irr_L1_fp = dirprd(i, sym1)
              if (trans1 == 'T') then
                  iabsoff1_fp(i+1) = iabsoff1_fp(i) + inummn(irr_L1_fp) * inumpq(i)
              else
                  iabsoff1_fp(i+1) = iabsoff1_fp(i) + inumpq(irr_L1_fp) * inummn(i)
              end if
              
              irr_L2_fp = dirprd(i, sym2)
              if (trans2 == 'N') then
                  iabsoff2_fp(i+1) = iabsoff2_fp(i) + inummn(irr_L2_fp) * inumrs(i)
              else
                  iabsoff2_fp(i+1) = iabsoff2_fp(i) + inumrs(irr_L2_fp) * inummn(i)
              end if
              
              irr_L3_fp = dirprd(i, sym3)
              iabsoff3_fp(i+1) = iabsoff3_fp(i) + inumpq(irr_L3_fp) * inumrs(i)
          end do

          do irr_rs = 1, nirrep
              irr_pq = dirprd(irr_rs, sym3)
              irr_mn = dirprd(irr_pq, sym1)
              
              if (trans1 == 'T') then; ioff1 = iabsoff1_fp(irr_pq)
              else; ioff1 = iabsoff1_fp(irr_mn)
              end if
              
              if (trans2 == 'N') then; ioff2 = iabsoff2_fp(irr_rs)
              else; ioff2 = iabsoff2_fp(irr_mn)
              end if
              
              ioff3 = iabsoff3_fp(irr_rs)
              
              m_size = inumpq(irr_pq)
              n_size = inumrs(irr_rs)
              k_size = inummn(irr_mn)
              
              ! 分离 M/N=0 与 K=0 判定，使用 DO 循环安全清理
              if (m_size == 0 .or. n_size == 0) cycle
              
              if (k_size == 0) then
                  if (fact2 == 0.0e0) then
                      do i = 0, m_size*n_size - 1
                          v3(ioff3 + i) = 0.0e0
                      end do
                  else if (fact2 /= 1.0e0) then
                      do i = 0, m_size*n_size - 1
                          v3(ioff3 + i) = fact2 * v3(ioff3 + i)
                      end do
                  end if
                  cycle
              end if
              
              lda1 = k_size; if (trans1 == 'N') lda1 = m_size
              lda2 = k_size; if (trans2 == 'T') lda2 = n_size
              
              ! 保证传递给底层 xgemms 的标量为双精度 REAL*8，以防底层强制取 8 字节边界引发越界
              call xgemmso(trans1, trans2, m_size, n_size, k_size, fact1in, &
                         v1(ioff1), lda1, v2(ioff2), lda2, fact2in, v3(ioff3), m_size)
          end do

#ifdef USE_MKL
          call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
          call openblas_set_num_threads(old_threads)
#endif
          return
      end if

      ! =============================================================
      ! 【阶段 2：严格 IOP 控制墙 (Strict Barrier)】
      ! =============================================================
      if (iop1_12 /= 1 .or. iop1_34 /= 1 .or. iop2_12 /= 1 .or. &
          iop2_34 /= 1 .or. iop3_12 /= 1 .or. iop3_34 /= 1) then
          write(6,*) "========================================================="
          write(6,*) " FATAL ERROR in vv2v_general_s:"
          write(6,*) " You are attempting to contract tensors with crossed or"
          write(6,*) " interleaved dimensions (e.g. 3142) while using compressed"
          write(6,*) " IOP settings (IOP=0 or 2). This is mathematically ambiguous"
          write(6,*) " and blocked by the engine."
          write(6,*) ""
          write(6,*) " ACTION REQUIRED:"
          write(6,*) " Please use 'sstgenso1' to unpack your tensors into full"
          write(6,*) " space (IOP=1) before passing them to this Slow Path."
          write(6,*) "========================================================="
          stop
      end if

      ! =============================================================
      ! 【阶段 3：Fallback 瓦片装配引擎】
      ! =============================================================
      pop1(1,:) = p1_1(:); pop1(2,:) = p1_2(:); pop1(3,:) = p1_3(:); pop1(4,:) = p1_4(:)
      pop2(1,:) = p2_1(:); pop2(2,:) = p2_2(:); pop2(3,:) = p2_3(:); pop2(4,:) = p2_4(:)
      pop3(1,:) = p3_1(:); pop3(2,:) = p3_2(:); pop3(3,:) = p3_3(:); pop3(4,:) = p3_4(:)

      c_ext(1) = type3(1:1); c_ext(2) = type3(2:2)
      c_ext(3) = type3(3:3); c_ext(4) = type3(4:4)
      idx = 0
      do i = 1, 4
         found = .false.
         do j = 1, 4; if (type1(i:i) == c_ext(j)) found = .true.; end do
         if (.not. found) then
            idx = idx + 1; if(idx <= 2) c_dum(idx) = type1(i:i)
         end if
      end do

      pos1 = 0; pos2 = 0; pos3 = 0
      do i = 1, 4
         do j = 1, 4
            if (type1(j:j) == c_ext(i)) pos1(i) = j
            if (type2(j:j) == c_ext(i)) pos2(i) = j
            if (type3(j:j) == c_ext(i)) pos3(i) = j
         end do
      end do
 
      do i = 1, 2
         do j = 1, 4
            if (type1(j:j) == c_dum(i)) pos1(4+i) = j
            if (type2(j:j) == c_dum(i)) pos2(4+i) = j
         end do
      end do

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

      num_threads = omp_get_max_threads()
      if (num_threads < 1) num_threads = 1

      ! 显存评估调整。单精度只有 3*4 = 12 字节开销，从而让 BS 更大，吞吐更强
      if (num_tasks >= num_threads * 2) then
         run_parallel = .true.
         max_mem_bytes = 3500000000_8 / num_threads 
         BS = int(sqrt(real(max_mem_bytes, 8) / 12.0d0))
#ifdef USE_MKL
         call mkl_set_num_threads(1)
#elif defined(USE_OPENBLAS)
         new_threads = 1_c_int
         call openblas_set_num_threads(new_threads)
#endif
      else
         run_parallel = .false.
         max_mem_bytes = 3500000000_8 
         BS = int(sqrt(real(max_mem_bytes, 8) / 12.0d0))
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
!$omp irr_phys1, irr_phys2, irr_phys3, base_v3, s_phys, s_v3, base_v3_adj, &
!$omp irr_D1, irr_D2, irr_logi, skip_block, d_D1, d_D2, i, j, &
!$omp irr_lft1, irr_rht1, irr_lft2, irr_rht2, base_v1, base_v2, b_v1_adj, b_v2_adj, &
!$omp s_v1_arr, s_v2_arr, buf_A, buf_B, buf_C, m_dim1, m_dim2, n_dim1, n_dim2, &
!$omp k_dim1, k_dim2, s_m1, s_m2, s_n1, s_n2, s_k1_v1, s_k2_v1, s_k1_v2, s_k2_v2, &
!$omp k1, k2, m1, m2, n1, n2, base_k, base_n, base_n_v3, idx_v3, &
!$omp m_start, n_start, k_start, k_end, cur_bs_m, cur_bs_n, cur_bs_k, &
!$omp m_global, n_global, k_global, m_loc, n_loc, k_loc)

      allocate(buf_A(BS, BS), buf_B(BS, BS), buf_C(BS, BS))

!$omp do schedule(dynamic, 1)
      do task_idx = 1, num_tasks
         irr_rht3 = task_list(1, task_idx)
         irr_L4   = task_list(2, task_idx)
         irr_L2   = task_list(3, task_idx)

         irr_lft3 = dirprd(irr_rht3, sym3); irr_L3 = dirprd(irr_rht3, irr_L4); irr_L1 = dirprd(irr_lft3, irr_L2)    

         irr_phys3(pos3(1)) = irr_L1; irr_phys3(pos3(2)) = irr_L2; irr_phys3(pos3(3)) = irr_L3; irr_phys3(pos3(4)) = irr_L4

         dim_L1 = pop3(pos3(1), irr_L1); dim_L2 = pop3(pos3(2), irr_L2)
         dim_L3 = pop3(pos3(3), irr_L3); dim_L4 = pop3(pos3(4), irr_L4)
         
         dim_X_arr(1) = dim_L1; dim_X_arr(2) = dim_L2; dim_X_arr(3) = dim_L3; dim_X_arr(4) = dim_L4
         cur_M = dim_X_arr(idx_m_in_v3(1)) * dim_X_arr(idx_m_in_v3(2))
         cur_N = dim_X_arr(idx_n_in_v3(1)) * dim_X_arr(idx_n_in_v3(2))

         base_v3 = iabsoff3(irr_rht3) + iposrht3(irr_phys3(4), irr_rht3) + iposlft3(irr_phys3(2), irr_lft3)
         s_phys(1) = 1; s_phys(2) = pop3(1, irr_phys3(1))
         s_phys(3) = inumlft3(irr_lft3); s_phys(4) = inumlft3(irr_lft3) * pop3(3, irr_phys3(3))

         do i = 1, 4; s_v3(i) = s_phys(pos3(i)); end do
         base_v3_adj = base_v3 - sum(s_v3(1:4)) + 1

         m_dim1 = dim_X_arr(idx_m_in_v3(1)); m_dim2 = dim_X_arr(idx_m_in_v3(2))
         n_dim1 = dim_X_arr(idx_n_in_v3(1)); n_dim2 = dim_X_arr(idx_n_in_v3(2))

         do n_start = 1, cur_N, BS
            cur_bs_n = min(n_start + BS - 1, cur_N) - n_start + 1

            do m_start = 1, cur_M, BS
               cur_bs_m = min(m_start + BS - 1, cur_M) - m_start + 1

               do j = 1, cur_bs_n
                  do i = 1, cur_bs_m
                     buf_C(i, j) = 0.0e0
                  end do
               end do

               do irr_D1 = 1, nirrep
                  irr_D2 = dirprd(sym1, irr_D1)
                  do i = 1, 4
                     if (pos1(i)/=0) then         
                        if (i==1) irr_D2 = dirprd(irr_D2, irr_L1)
                        if (i==2) irr_D2 = dirprd(irr_D2, irr_L2)
                        if (i==3) irr_D2 = dirprd(irr_D2, irr_L3)
                        if (i==4) irr_D2 = dirprd(irr_D2, irr_L4)
                     end if
                  end do

                  irr_logi(1) = irr_L1; irr_logi(2) = irr_L2; irr_logi(3) = irr_L3
                  irr_logi(4) = irr_L4; irr_logi(5) = irr_D1; irr_logi(6) = irr_D2

                  skip_block = .false.
                  do i = 1, 6
                     if (pos1(i)/=0 .and. pop1(pos1(i), irr_logi(i))==0) skip_block = .true.
                     if (pos2(i)/=0 .and. pop2(pos2(i), irr_logi(i))==0) skip_block = .true.
                  end do
                  if (skip_block) cycle

                  d_D1(irr_D1) = pop1(pos1(5), irr_D1)
                  d_D2(irr_D1) = pop1(pos1(6), irr_D2)
                  cur_K = d_D1(irr_D1) * d_D2(irr_D1)
                  if (cur_K == 0) cycle

                  do i = 1, 6
                     if (pos1(i) /= 0) irr_phys1(pos1(i)) = irr_logi(i)
                     if (pos2(i) /= 0) irr_phys2(pos2(i)) = irr_logi(i)
                  end do

                  irr_lft1 = dirprd(irr_phys1(1), irr_phys1(2)); irr_rht1 = dirprd(irr_phys1(3), irr_phys1(4))
                  irr_lft2 = dirprd(irr_phys2(1), irr_phys2(2)); irr_rht2 = dirprd(irr_phys2(3), irr_phys2(4))
                  
                  base_v1 = iabsoff1(irr_rht1) + iposrht1(irr_phys1(4), irr_rht1) + iposlft1(irr_phys1(2), irr_lft1)
                  base_v2 = iabsoff2(irr_rht2) + iposrht2(irr_phys2(4), irr_rht2) + iposlft2(irr_phys2(2), irr_lft2)

                  s_phys(1) = 1; s_phys(2) = pop1(1, irr_phys1(1)); s_phys(3) = inumlft1(irr_lft1); s_phys(4) = inumlft1(irr_lft1) * pop1(3, irr_phys1(3))
                  s_v1_arr(:, irr_D1) = 0
                  do i = 1, 6; if(pos1(i)/=0) s_v1_arr(i, irr_D1) = s_phys(pos1(i)); end do

                  s_phys(1) = 1; s_phys(2) = pop2(1, irr_phys2(1)); s_phys(3) = inumlft2(irr_lft2); s_phys(4) = inumlft2(irr_lft2) * pop2(3, irr_phys2(3))
                  s_v2_arr(:, irr_D1) = 0
                  do i = 1, 6; if(pos2(i)/=0) s_v2_arr(i, irr_D1) = s_phys(pos2(i)); end do

                  b_v1_adj(irr_D1) = base_v1 - sum(s_v1_arr(:, irr_D1)) + 1
                  b_v2_adj(irr_D1) = base_v2 - sum(s_v2_arr(:, irr_D1)) + 1

                  s_m1 = s_v1_arr(idx_m_in_v3(1), irr_D1); s_m2 = s_v1_arr(idx_m_in_v3(2), irr_D1)
                  s_n1 = s_v2_arr(idx_n_in_v3(1), irr_D1); s_n2 = s_v2_arr(idx_n_in_v3(2), irr_D1)
                  k_dim1 = d_D1(irr_D1)
                  s_k1_v1 = s_v1_arr(5, irr_D1); s_k2_v1 = s_v1_arr(6, irr_D1)
                  s_k1_v2 = s_v2_arr(5, irr_D1); s_k2_v2 = s_v2_arr(6, irr_D1)

                  do k_start = 1, cur_K, BS
                     cur_bs_k = min(k_start + BS - 1, cur_K) - k_start + 1

                     k_global = k_start
                     k2 = (k_global - 1) / k_dim1 + 1; k1 = mod(k_global - 1, k_dim1) + 1
                     do k_loc = 1, cur_bs_k
                        base_k = b_v1_adj(irr_D1) + k1*s_k1_v1 + k2*s_k2_v1
                        m_global = m_start
                        m2 = (m_global - 1) / m_dim1 + 1; m1 = mod(m_global - 1, m_dim1) + 1
                        do m_loc = 1, cur_bs_m
                           buf_A(m_loc, k_loc) = v1(base_k + m2*s_m2 + m1*s_m1)
                           m1 = m1 + 1; if (m1 > m_dim1) then; m1 = 1; m2 = m2 + 1; end if
                        end do
                        k1 = k1 + 1; if (k1 > k_dim1) then; k1 = 1; k2 = k2 + 1; end if
                     end do

                     n_global = n_start
                     n2 = (n_global - 1) / n_dim1 + 1; n1 = mod(n_global - 1, n_dim1) + 1
                     do n_loc = 1, cur_bs_n
                        base_n = b_v2_adj(irr_D1) + n1*s_n1 + n2*s_n2
                        k_global = k_start
                        k2 = (k_global - 1) / k_dim1 + 1; k1 = mod(k_global - 1, k_dim1) + 1
                        do k_loc = 1, cur_bs_k
                           buf_B(k_loc, n_loc) = v2(base_n + k1*s_k1_v2 + k2*s_k2_v2)
                           k1 = k1 + 1; if (k1 > k_dim1) then; k1 = 1; k2 = k2 + 1; end if
                        end do
                        n1 = n1 + 1; if (n1 > n_dim1) then; n1 = 1; n2 = n2 + 1; end if
                     end do

                     ! 【终极隐患排除】：强制使用 REAL*8 的 1.0d0 避免底层越位读取导致单精度核心计算乱码
                     call xgemmso('N', 'N', cur_bs_m, cur_bs_n, cur_bs_k, 1.0d0, &
                                buf_A, BS, buf_B, BS, 1.0d0, buf_C, BS)

                   end do 
               end do 

               n_global = n_start
               n2 = (n_global - 1) / n_dim1 + 1; n1 = mod(n_global - 1, n_dim1) + 1
               do n_loc = 1, cur_bs_n
                  base_n_v3 = base_v3_adj + n1*s_v3(idx_n_in_v3(1)) + n2*s_v3(idx_n_in_v3(2))
                  m_global = m_start
                  m2 = (m_global - 1) / m_dim1 + 1; m1 = mod(m_global - 1, m_dim1) + 1
                  
                  do m_loc = 1, cur_bs_m
                     idx_v3 = base_n_v3 + m2*s_v3(idx_m_in_v3(2)) + m1*s_v3(idx_m_in_v3(1))
                     
                     if (fact2 == 0.0e0) then
                        v3(idx_v3) = fact1 * buf_C(m_loc, n_loc)
                     else
                        v3(idx_v3) = fact2 * v3(idx_v3) + fact1 * buf_C(m_loc, n_loc)
                     end if
                     
                     m1 = m1 + 1; if (m1 > m_dim1) then; m1 = 1; m2 = m2 + 1; end if
                  end do
                  n1 = n1 + 1; if (n1 > n_dim1) then; n1 = 1; n2 = n2 + 1; end if
               end do

            end do 
         end do 
      end do 
!$omp end do
      deallocate(buf_A, buf_B, buf_C)
!$omp end parallel

      ! =============================================================
      ! 恢复状态
      ! =============================================================
#ifdef USE_MKL
      call mkl_set_num_threads(old_threads)
#elif defined(USE_OPENBLAS)
      call openblas_set_num_threads(old_threads)
#endif

      return
      end subroutine vv2v_dpd_444_s
