
      subroutine sstgenso(win, wout, nsize, pop1, pop2, pop3, pop4, ipos, irrepx, type)
      use mod_iop
      use mod_symm                 ! 引入对称性直积运算模块
      use omp_lib                  ! 引入 openmp 多线程库
      implicit none                ! 强制显式声明所有变量

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type     ! 转置模式，例如 '3142'
      real*4, intent(in)  :: win(*)   ! 输入的一维连续张量 w(p,q,r,s)
      real*4, intent(out) :: wout(*)  ! 输出的一维连续张量 w(新排布)
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8) ! 四个指标的各irrep基函数数
      integer, intent(in)           :: irrepx   ! 张量总空间对称性
      integer, intent(in) :: nsize 
      real*4, intent(in) :: ipos(*)
      call sstgenso1(win, wout, pop1, pop2, pop3, pop4, irrepx, type)

      return
      end


      subroutine sstgenso1(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      use mod_iop
      use mod_symm                 ! 引入对称性直积运算模块
      use omp_lib                  ! 引入 openmp 多线程库
      implicit none                ! 强制显式声明所有变量

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type     ! 转置模式，例如 '3142'
      real*4, intent(in)  :: win(*)   ! 输入的一维连续张量 w(p,q,r,s)
      real*4, intent(out) :: wout(*)  ! 输出的一维连续张量 w(新排布)
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8) ! 四个指标的各irrep基函数数
      integer, intent(in)           :: irrepx   ! 张量总空间对称性

      if(prec==8.and.igpu==0) then 
         call sstgend_opt(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      elseif(prec==4.and.igpu==0) then 
         call sstgens_opt(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      elseif(prec==8.and.igpu==1) then 
         call sstgend_opt_g(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      elseif(prec==4.and.igpu==1) then 
         call sstgens_opt_g(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      endif
      return
      end

subroutine sstgend_opt(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      use mod_symm                 ! 引入对称性直积运算模块
      use omp_lib                  ! 引入 openmp 多线程库
      implicit none                ! 强制显式声明所有变量

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type     ! 转置模式，例如 '3142'
      double precision, intent(in)  :: win(*)   ! 输入的一维连续张量 w(p,q,r,s)
      double precision, intent(out) :: wout(*)  ! 输出的一维连续张量 w(新排布)
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8) ! 四个指标的各irrep基函数数
      integer, intent(in)           :: irrepx   ! 张量总空间对称性

      ! --- 局部数组：对称性元数据 ---
      integer*8 :: inumlftn(9), inumrhtn(9), inumlfto(9), inumrhto(9) ! 新旧张量的左右半区总行/列数
      integer*8 :: iabsoffn(9), iabsoffo(9)                           ! 新旧张量按右半区irrep划分的大块基址
      integer*8 :: iposrhtn(9,9), iposlftn(9,9), iposrhto(9,9), iposlfto(9,9) ! 新旧张量在块内的相对偏移前缀和

      ! --- 局部变量：通用映射与步长 ---
      integer   :: pop_in(4, 8), pop_out(4, 8) ! 二维化pop数组：(指标序号1~4, 对称性1~8)
      integer   :: p(4), n_p(4), irr_p(4)      ! p:位置映射数组; n_p:当前维度长; irr_p:当前维度irrep
      integer   :: i, irrep, irr_r, irr_l, irr_l2 ! 循环与对称性直积相关的临时变量
      integer   :: irrrht, irrlft, irr4, irr3, irr2, irr1 ! dpd四层嵌套的各指标irrep
      integer   :: irr_lft_in, irr_rht_in, irr_lft_out, irr_rht_out ! 新旧张量左右半区的总体irrep
      integer   :: i1, i2, i3, i4, n1, n2, n3, n4 ! 实际数据存取的循环索引与边界
      integer*8 :: ilfto_val, ilftn_val        ! 当前irrep下，新旧张量左半区(列主序的行)的总容量
      integer*8 :: base_old, base_new, base_old_adj, base_new_adj ! 基础物理地址与消除+1偏移后的调整地址
      integer*8 :: stride_in(4), stride_out(4) ! 新旧张量在四个维度上的物理内存跨度(步长)
      integer*8 :: ioffold, ioffnew            ! 最内层展开后的一维物理内存指针

      ! -------------------------------------------------------------------
      ! 1. 解析 type 字符串，建立维度映射关系
      ! -------------------------------------------------------------------
      do i = 1, 4
         p(i) = ichar(type(i:i)) - ichar('0') ! 将字符转化为映射索引，如 '3' -> 3
         if (p(i) < 1 .or. p(i) > 4) then     ! 越界检查，防止接口传参错位导致崩溃
            print *, "fatal error in sstgend_opt: invalid type received: ", type
            stop
         end if
      end do

      pop_in(1, :) = pop1(:) ! 将离散的旧指标1的pop数组打包成二维矩阵第1行
      pop_in(2, :) = pop2(:) ! 将离散的旧指标2的pop数组打包成二维矩阵第2行
      pop_in(3, :) = pop3(:) ! 将离散的旧指标3的pop数组打包成二维矩阵第3行
      pop_in(4, :) = pop4(:) ! 将离散的旧指标4的pop数组打包成二维矩阵第4行

      do i = 1, 4
         pop_out(i, :) = pop_in(p(i), :) ! 核心：根据p数组，动态生成新张量四个指标的pop分布
      end do

      ! -------------------------------------------------------------------
      ! 2. 计算给定对称性下的左半侧(指标1,2)和右半侧(指标3,4)的基函数总数
      ! -------------------------------------------------------------------
      inumlftn = 0_8; inumrhtn = 0_8 ! 初始化新张量(n)的左右块容量
      inumlfto = 0_8; inumrhto = 0_8 ! 初始化旧张量(o)的左右块容量

      do irrep = 1, nirrep ! 遍历半区的总体对称性 (左半侧或右半侧的irrep)
         do irr_r = 1, nirrep ! 遍历该半区内“第二个指标”的irrep (如指标2或指标4)
            irr_l = dirprd(irr_r, irrep) ! 反推该半区内“第一个指标”的irrep，确保直积等于irrep
            inumlfto(irrep) = inumlfto(irrep) + int(pop_in(1, irr_l), 8) * pop_in(2, irr_r) ! 累加旧张量左半区(1,2)容量
            inumrhto(irrep) = inumrhto(irrep) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) ! 累加旧张量右半区(3,4)容量
            inumlftn(irrep) = inumlftn(irrep) + int(pop_out(1, irr_l), 8) * pop_out(2, irr_r) ! 累加新张量左半区容量
            inumrhtn(irrep) = inumrhtn(irrep) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) ! 累加新张量右半区容量
         end do
      end do

      ! -------------------------------------------------------------------
      ! 3. 计算 dpd 三级跳寻址的前置偏移量表 (大块起点与块内前缀和)
      ! -------------------------------------------------------------------
      iabsoffn(1) = 0_8 ! 新张量第一个超级大块(右半侧irrep=1)的绝对起点永远是0
      iabsoffo(1) = 0_8 ! 旧张量第一个超级大块的绝对起点永远是0

      do irrrht = 1, nirrep ! 遍历右半区总体对称性(控制超级大块)
         irrlft = dirprd(irrrht, irrepx)  ! 反推左半侧要求，必须满足：左侧irrep * 右侧irrep = 张量总irrep
         
         iabsoffo(irrrht+1) = iabsoffo(irrrht) + inumlfto(irrlft) * inumrhto(irrrht) ! 旧张量下一个大块的起点
         iabsoffn(irrrht+1) = iabsoffn(irrrht) + inumlftn(irrlft) * inumrhtn(irrrht) ! 新张量下一个大块的起点
         
         iposrhto(1, irrrht) = 0_8; iposlfto(1, irrlft) = 0_8 ! 旧张量内部偏移前缀和的起点初始化
         iposrhtn(1, irrrht) = 0_8; iposlftn(1, irrlft) = 0_8 ! 新张量内部偏移前缀和的起点初始化

         do irr_r = 1, nirrep ! 遍历组成当前半区的“第二个指标”的对称性
            irr_l = dirprd(irr_r, irrrht) ! 反推组成右半区的“第一个指标”(如指标3)的对称性
            iposrhto(irr_r+1, irrrht) = iposrhto(irr_r, irrrht) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) * inumlfto(irrlft) ! 旧张量右半区内部偏移累加
            iposrhtn(irr_r+1, irrrht) = iposrhtn(irr_r, irrrht) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) * inumlftn(irrlft) ! 新张量右半区内部偏移累加

            irr_l2 = dirprd(irr_r, irrlft) ! 反推组成左半区的“第一个指标”(如指标1)的对称性
            iposlfto(irr_r+1, irrlft) = iposlfto(irr_r, irrlft) + int(pop_in(1, irr_l2), 8) * pop_in(2, irr_r) ! 旧张量左半区内部偏移累加
            iposlftn(irr_r+1, irrlft) = iposlftn(irr_r, irrlft) + int(pop_out(1, irr_l2), 8) * pop_out(2, irr_r) ! 新张量左半区内部偏移累加
         end do
      end do

      ! -------------------------------------------------------------------
      ! 4. 张量转置核心循环：遍历所有合法的 irrep 组合并搬运数据
      ! -------------------------------------------------------------------
      do irrrht = 1, nirrep ! 第1层：固定右半区(旧指标3,4)总体对称性
         irrlft = dirprd(irrrht, irrepx) ! 反推左半区(旧指标1,2)总体对称性
         do irr4 = 1, nirrep ! 第2层：固定旧指标4对称性
            irr3 = dirprd(irrrht, irr4) ! 反推旧指标3对称性
            do irr2 = 1, nirrep ! 第3层：固定旧指标2对称性
               irr1 = dirprd(irrlft, irr2) ! 反推旧指标1对称性

               n1 = pop_in(1, irr1) ! 获取当前irrep组合下，旧指标1的维度大小
               n2 = pop_in(2, irr2) ! 获取当前irrep组合下，旧指标2的维度大小
               n3 = pop_in(3, irr3) ! 获取当前irrep组合下，旧指标3的维度大小
               n4 = pop_in(4, irr4) ! 获取当前irrep组合下，旧指标4的维度大小

               if (n1 == 0 .or. n2 == 0 .or. n3 == 0 .or. n4 == 0) cycle ! 如果任何一维为0，该dpd块为空，直接跳过

               n_p(1) = n1;   n_p(2) = n2;   n_p(3) = n3;   n_p(4) = n4 ! 将旧维度打包为数组，供通用映射p读取
               irr_p(1) = irr1; irr_p(2) = irr2; irr_p(3) = irr3; irr_p(4) = irr4 ! 将旧irrep打包为数组，供通用映射p读取
               
               irr_lft_in = irrlft ! 明确旧张量的左半区总体对称性
               irr_rht_in = irrrht ! 明确旧张量的右半区总体对称性
               irr_lft_out = dirprd(irr_p(p(1)), irr_p(p(2))) ! 根据p映射，动态算出新张量左半侧(新指标1,2)的总体对称性
               irr_rht_out = dirprd(irr_p(p(3)), irr_p(p(4))) ! 根据p映射，动态算出新张量右半侧(新指标3,4)的总体对称性

               ilfto_val = inumlfto(irr_lft_in)  ! 获取旧张量在当前左侧irrep下的一整行(左半区)的数据量大小
               ilftn_val = inumlftn(irr_lft_out) ! 获取新张量在当前左侧irrep下的一整行(左半区)的数据量大小

               ! 查表获取块级起始物理地址 = 大块基址 + 右侧子块局部偏移 + 左侧子块局部偏移
               base_old = iabsoffo(irr_rht_in)  + iposrhto(irr4, irr_rht_in)      + iposlfto(irr2, irr_lft_in) ! 旧张量(顺序存取)的块级基址
               base_new = iabsoffn(irr_rht_out) + iposrhtn(irr_p(p(4)), irr_rht_out) + iposlftn(irr_p(p(2)), irr_lft_out) ! 新张量(按p映射提取对称性后查表)的块级基址

               ! 预计算多维数组退化为一维数组时的物理步长(stride)
               stride_in(4) = int(n3, 8) * ilfto_val ! 旧指标4每增加1，需跨过整个指标3和整个左半区
               stride_in(3) = ilfto_val              ! 旧指标3每增加1，需跨过一个完整的左半区
               stride_in(2) = int(n1, 8)             ! 旧指标2每增加1，需跨过n1个元素
               stride_in(1) = 1_8                    ! 旧指标1在最内层，内存严格连续，跨度为1

               ! 将新张量映射回旧维度循环体的步长
               stride_out(p(4)) = int(n_p(p(3)), 8) * ilftn_val ! 映射过去的新指标4带来的内存跨度
               stride_out(p(3)) = ilftn_val                     ! 映射过去的新指标3带来的内存跨度
               stride_out(p(2)) = int(n_p(p(1)), 8)             ! 映射过去的新指标2带来的内存跨度
               stride_out(p(1)) = 1_8                           ! 映射过去的新指标1带来的内存跨度

               ! 彻底消除 i 从 1 开始引起的代数偏移：减去所有步长并补偿1，保证 i=1 时公式完美闭合到 base
               base_old_adj = base_old - stride_in(4) - stride_in(3) - stride_in(2) - stride_in(1) + 1_8
               base_new_adj = base_new - stride_out(4) - stride_out(3) - stride_out(2) - stride_out(1) + 1_8

!$omp parallel do collapse(3) private(i3, i2, i1, ioffold, ioffnew) ! 多线程并行最外层i4(临时私有化内部变量，严防冲突)
               do i4 = 1, n4 ! 遍历旧指标4
                  do i3 = 1, n3 ! 遍历旧指标3
                     do i2 = 1, n2 ! 遍历旧指标2
                        ioffold = base_old_adj + i4*stride_in(4) + i3*stride_in(3) + i2*stride_in(2) ! 根据外侧三层步长，累加出旧一维指针
                        ioffnew = base_new_adj + i4*stride_out(4) + i3*stride_out(3) + i2*stride_out(2) ! 同理，累加出目标新一维指针
                        
                        do i1 = 1, n1 ! 最内层循环遍历旧指标1 (交由编译器自动向量化simd)
                           wout(ioffnew + i1*stride_out(1)) = win(ioffold + i1) ! 旧指标1连续读取，新指标按映射步长离散/连续写入
                        end do
                     end do
                  end do
               end do
!$omp end parallel do

            end do
         end do
      end do

      return
      end subroutine sstgend_opt

      subroutine sstgens_opt(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      use mod_symm                 ! 引入对称性直积运算模块
      use omp_lib                  ! 引入 openmp 多线程库
      implicit none                ! 强制显式声明所有变量

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type     ! 转置模式，例如 '3142'
      real*4, intent(in)  :: win(*)   ! 输入的一维连续张量 w(p,q,r,s)
      real*4, intent(out) :: wout(*)  ! 输出的一维连续张量 w(新排布)
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8) ! 四个指标的各irrep基函数数
      integer, intent(in)           :: irrepx   ! 张量总空间对称性

      ! --- 局部数组：对称性元数据 ---
      integer*8 :: inumlftn(9), inumrhtn(9), inumlfto(9), inumrhto(9) ! 新旧张量的左右半区总行/列数
      integer*8 :: iabsoffn(9), iabsoffo(9)                           ! 新旧张量按右半区irrep划分的大块基址
      integer*8 :: iposrhtn(9,9), iposlftn(9,9), iposrhto(9,9), iposlfto(9,9) ! 新旧张量在块内的相对偏移前缀和

      ! --- 局部变量：通用映射与步长 ---
      integer   :: pop_in(4, 8), pop_out(4, 8) ! 二维化pop数组：(指标序号1~4, 对称性1~8)
      integer   :: p(4), n_p(4), irr_p(4)      ! p:位置映射数组; n_p:当前维度长; irr_p:当前维度irrep
      integer   :: i, irrep, irr_r, irr_l, irr_l2 ! 循环与对称性直积相关的临时变量
      integer   :: irrrht, irrlft, irr4, irr3, irr2, irr1 ! dpd四层嵌套的各指标irrep
      integer   :: irr_lft_in, irr_rht_in, irr_lft_out, irr_rht_out ! 新旧张量左右半区的总体irrep
      integer   :: i1, i2, i3, i4, n1, n2, n3, n4 ! 实际数据存取的循环索引与边界
      integer*8 :: ilfto_val, ilftn_val        ! 当前irrep下，新旧张量左半区(列主序的行)的总容量
      integer*8 :: base_old, base_new, base_old_adj, base_new_adj ! 基础物理地址与消除+1偏移后的调整地址
      integer*8 :: stride_in(4), stride_out(4) ! 新旧张量在四个维度上的物理内存跨度(步长)
      integer*8 :: ioffold, ioffnew            ! 最内层展开后的一维物理内存指针

      ! -------------------------------------------------------------------
      ! 1. 解析 type 字符串，建立维度映射关系
      ! -------------------------------------------------------------------
      do i = 1, 4
         p(i) = ichar(type(i:i)) - ichar('0') ! 将字符转化为映射索引，如 '3' -> 3
         if (p(i) < 1 .or. p(i) > 4) then     ! 越界检查，防止接口传参错位导致崩溃
            print *, "fatal error in sstgend_opt: invalid type received: ", type
            stop
         end if
      end do

      pop_in(1, :) = pop1(:) ! 将离散的旧指标1的pop数组打包成二维矩阵第1行
      pop_in(2, :) = pop2(:) ! 将离散的旧指标2的pop数组打包成二维矩阵第2行
      pop_in(3, :) = pop3(:) ! 将离散的旧指标3的pop数组打包成二维矩阵第3行
      pop_in(4, :) = pop4(:) ! 将离散的旧指标4的pop数组打包成二维矩阵第4行

      do i = 1, 4
         pop_out(i, :) = pop_in(p(i), :) ! 核心：根据p数组，动态生成新张量四个指标的pop分布
      end do

      ! -------------------------------------------------------------------
      ! 2. 计算给定对称性下的左半侧(指标1,2)和右半侧(指标3,4)的基函数总数
      ! -------------------------------------------------------------------
      inumlftn = 0_8; inumrhtn = 0_8 ! 初始化新张量(n)的左右块容量
      inumlfto = 0_8; inumrhto = 0_8 ! 初始化旧张量(o)的左右块容量

      do irrep = 1, nirrep ! 遍历半区的总体对称性 (左半侧或右半侧的irrep)
         do irr_r = 1, nirrep ! 遍历该半区内“第二个指标”的irrep (如指标2或指标4)
            irr_l = dirprd(irr_r, irrep) ! 反推该半区内“第一个指标”的irrep，确保直积等于irrep
            inumlfto(irrep) = inumlfto(irrep) + int(pop_in(1, irr_l), 8) * pop_in(2, irr_r) ! 累加旧张量左半区(1,2)容量
            inumrhto(irrep) = inumrhto(irrep) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) ! 累加旧张量右半区(3,4)容量
            inumlftn(irrep) = inumlftn(irrep) + int(pop_out(1, irr_l), 8) * pop_out(2, irr_r) ! 累加新张量左半区容量
            inumrhtn(irrep) = inumrhtn(irrep) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) ! 累加新张量右半区容量
         end do
      end do

      ! -------------------------------------------------------------------
      ! 3. 计算 dpd 三级跳寻址的前置偏移量表 (大块起点与块内前缀和)
      ! -------------------------------------------------------------------
      iabsoffn(1) = 0_8 ! 新张量第一个超级大块(右半侧irrep=1)的绝对起点永远是0
      iabsoffo(1) = 0_8 ! 旧张量第一个超级大块的绝对起点永远是0

      do irrrht = 1, nirrep ! 遍历右半区总体对称性(控制超级大块)
         irrlft = dirprd(irrrht, irrepx)  ! 反推左半侧要求，必须满足：左侧irrep * 右侧irrep = 张量总irrep
         
         iabsoffo(irrrht+1) = iabsoffo(irrrht) + inumlfto(irrlft) * inumrhto(irrrht) ! 旧张量下一个大块的起点
         iabsoffn(irrrht+1) = iabsoffn(irrrht) + inumlftn(irrlft) * inumrhtn(irrrht) ! 新张量下一个大块的起点
         
         iposrhto(1, irrrht) = 0_8; iposlfto(1, irrlft) = 0_8 ! 旧张量内部偏移前缀和的起点初始化
         iposrhtn(1, irrrht) = 0_8; iposlftn(1, irrlft) = 0_8 ! 新张量内部偏移前缀和的起点初始化

         do irr_r = 1, nirrep ! 遍历组成当前半区的“第二个指标”的对称性
            irr_l = dirprd(irr_r, irrrht) ! 反推组成右半区的“第一个指标”(如指标3)的对称性
            iposrhto(irr_r+1, irrrht) = iposrhto(irr_r, irrrht) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) * inumlfto(irrlft) ! 旧张量右半区内部偏移累加
            iposrhtn(irr_r+1, irrrht) = iposrhtn(irr_r, irrrht) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) * inumlftn(irrlft) ! 新张量右半区内部偏移累加

            irr_l2 = dirprd(irr_r, irrlft) ! 反推组成左半区的“第一个指标”(如指标1)的对称性
            iposlfto(irr_r+1, irrlft) = iposlfto(irr_r, irrlft) + int(pop_in(1, irr_l2), 8) * pop_in(2, irr_r) ! 旧张量左半区内部偏移累加
            iposlftn(irr_r+1, irrlft) = iposlftn(irr_r, irrlft) + int(pop_out(1, irr_l2), 8) * pop_out(2, irr_r) ! 新张量左半区内部偏移累加
         end do
      end do

      ! -------------------------------------------------------------------
      ! 4. 张量转置核心循环：遍历所有合法的 irrep 组合并搬运数据
      ! -------------------------------------------------------------------
      do irrrht = 1, nirrep ! 第1层：固定右半区(旧指标3,4)总体对称性
         irrlft = dirprd(irrrht, irrepx) ! 反推左半区(旧指标1,2)总体对称性
         do irr4 = 1, nirrep ! 第2层：固定旧指标4对称性
            irr3 = dirprd(irrrht, irr4) ! 反推旧指标3对称性
            do irr2 = 1, nirrep ! 第3层：固定旧指标2对称性
               irr1 = dirprd(irrlft, irr2) ! 反推旧指标1对称性

               n1 = pop_in(1, irr1) ! 获取当前irrep组合下，旧指标1的维度大小
               n2 = pop_in(2, irr2) ! 获取当前irrep组合下，旧指标2的维度大小
               n3 = pop_in(3, irr3) ! 获取当前irrep组合下，旧指标3的维度大小
               n4 = pop_in(4, irr4) ! 获取当前irrep组合下，旧指标4的维度大小

               if (n1 == 0 .or. n2 == 0 .or. n3 == 0 .or. n4 == 0) cycle ! 如果任何一维为0，该dpd块为空，直接跳过

               n_p(1) = n1;   n_p(2) = n2;   n_p(3) = n3;   n_p(4) = n4 ! 将旧维度打包为数组，供通用映射p读取
               irr_p(1) = irr1; irr_p(2) = irr2; irr_p(3) = irr3; irr_p(4) = irr4 ! 将旧irrep打包为数组，供通用映射p读取
               
               irr_lft_in = irrlft ! 明确旧张量的左半区总体对称性
               irr_rht_in = irrrht ! 明确旧张量的右半区总体对称性
               irr_lft_out = dirprd(irr_p(p(1)), irr_p(p(2))) ! 根据p映射，动态算出新张量左半侧(新指标1,2)的总体对称性
               irr_rht_out = dirprd(irr_p(p(3)), irr_p(p(4))) ! 根据p映射，动态算出新张量右半侧(新指标3,4)的总体对称性

               ilfto_val = inumlfto(irr_lft_in)  ! 获取旧张量在当前左侧irrep下的一整行(左半区)的数据量大小
               ilftn_val = inumlftn(irr_lft_out) ! 获取新张量在当前左侧irrep下的一整行(左半区)的数据量大小

               ! 查表获取块级起始物理地址 = 大块基址 + 右侧子块局部偏移 + 左侧子块局部偏移
               base_old = iabsoffo(irr_rht_in)  + iposrhto(irr4, irr_rht_in)      + iposlfto(irr2, irr_lft_in) ! 旧张量(顺序存取)的块级基址
               base_new = iabsoffn(irr_rht_out) + iposrhtn(irr_p(p(4)), irr_rht_out) + iposlftn(irr_p(p(2)), irr_lft_out) ! 新张量(按p映射提取对称性后查表)的块级基址

               ! 预计算多维数组退化为一维数组时的物理步长(stride)
               stride_in(4) = int(n3, 8) * ilfto_val ! 旧指标4每增加1，需跨过整个指标3和整个左半区
               stride_in(3) = ilfto_val              ! 旧指标3每增加1，需跨过一个完整的左半区
               stride_in(2) = int(n1, 8)             ! 旧指标2每增加1，需跨过n1个元素
               stride_in(1) = 1_8                    ! 旧指标1在最内层，内存严格连续，跨度为1

               ! 将新张量映射回旧维度循环体的步长
               stride_out(p(4)) = int(n_p(p(3)), 8) * ilftn_val ! 映射过去的新指标4带来的内存跨度
               stride_out(p(3)) = ilftn_val                     ! 映射过去的新指标3带来的内存跨度
               stride_out(p(2)) = int(n_p(p(1)), 8)             ! 映射过去的新指标2带来的内存跨度
               stride_out(p(1)) = 1_8                           ! 映射过去的新指标1带来的内存跨度

               ! 彻底消除 i 从 1 开始引起的代数偏移：减去所有步长并补偿1，保证 i=1 时公式完美闭合到 base
               base_old_adj = base_old - stride_in(4) - stride_in(3) - stride_in(2) - stride_in(1) + 1_8
               base_new_adj = base_new - stride_out(4) - stride_out(3) - stride_out(2) - stride_out(1) + 1_8

!$omp parallel do collapse(3) private(i3, i2, i1, ioffold, ioffnew) ! 多线程并行最外层i4(临时私有化内部变量，严防冲突)
               do i4 = 1, n4 ! 遍历旧指标4
                  do i3 = 1, n3 ! 遍历旧指标3
                     do i2 = 1, n2 ! 遍历旧指标2
                        ioffold = base_old_adj + i4*stride_in(4) + i3*stride_in(3) + i2*stride_in(2) ! 根据外侧三层步长，累加出旧一维指针
                        ioffnew = base_new_adj + i4*stride_out(4) + i3*stride_out(3) + i2*stride_out(2) ! 同理，累加出目标新一维指针
                        
                        do i1 = 1, n1 ! 最内层循环遍历旧指标1 (交由编译器自动向量化simd)
                           wout(ioffnew + i1*stride_out(1)) = win(ioffold + i1) ! 旧指标1连续读取，新指标按映射步长离散/连续写入
                        end do
                     end do
                  end do
               end do
!$omp end parallel do

            end do
         end do
      end do

      return
      end subroutine sstgens_opt

      subroutine sstgend_opt_g(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      use mod_symm                 ! 引入对称性直积运算模块
      implicit none                ! 强制显式声明

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type     ! 转置模式，例如 '3142'
      double precision, intent(in)  :: win(*)   ! 纯设备端(GPU)指针：输入张量
      double precision, intent(out) :: wout(*)  ! 纯设备端(GPU)指针：输出张量
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8) ! CPU端的各irrep基函数数
      integer, intent(in)           :: irrepx   ! 张量总空间对称性

      ! --- 局部变量 ---
      integer*8 :: inumlftn(9), inumrhtn(9), inumlfto(9), inumrhto(9) ! 新旧张量左右半区容量
      integer*8 :: iabsoffn(9), iabsoffo(9)                           ! 大块基址
      integer*8 :: iposrhtn(9,9), iposlftn(9,9), iposrhto(9,9), iposlfto(9,9) ! 块内偏移
      integer   :: pop_in(4, 8), pop_out(4, 8) ! 二维化pop数组
      integer   :: p(4), n_p(4), irr_p(4)      ! 通用映射与维度
      integer   :: i, irrep, irr_r, irr_l, irr_l2
      integer   :: irrrht, irrlft, irr4, irr3, irr2, irr1
      integer   :: irr_lft_in, irr_rht_in, irr_lft_out, irr_rht_out
      integer   :: i1, i2, i3, i4, n1, n2, n3, n4
      integer*8 :: ilfto_val, ilftn_val
      integer*8 :: base_old, base_new, base_old_adj, base_new_adj
      integer*8 :: stride_in(4), stride_out(4) ! 新旧张量物理内存步长
      integer*8 :: ioffold, ioffnew            ! 一维物理指针

      ! --- 1. 解析映射与预计算元数据 (在CPU端执行，计算量极小) ---
      do i = 1, 4
         p(i) = ichar(type(i:i)) - ichar('0')
         if (p(i) < 1 .or. p(i) > 4) stop ! 防御性检查
      end do

      ! 打包pop数组
      pop_in(1, :) = pop1(:); pop_in(2, :) = pop2(:)
      pop_in(3, :) = pop3(:); pop_in(4, :) = pop4(:)

      ! 根据映射生成输出的pop分布
      do i = 1, 4
         pop_out(i, :) = pop_in(p(i), :)
      end do

      ! 计算左右半区容量
      inumlftn = 0_8; inumrhtn = 0_8; inumlfto = 0_8; inumrhto = 0_8
      do irrep = 1, nirrep
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrep)
            inumlfto(irrep) = inumlfto(irrep) + int(pop_in(1, irr_l), 8) * pop_in(2, irr_r)
            inumrhto(irrep) = inumrhto(irrep) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r)
            inumlftn(irrep) = inumlftn(irrep) + int(pop_out(1, irr_l), 8) * pop_out(2, irr_r)
            inumrhtn(irrep) = inumrhtn(irrep) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r)
         end do
      end do

      ! 计算基址与偏移表
      iabsoffn(1) = 0_8; iabsoffo(1) = 0_8
      do irrrht = 1, nirrep
         irrlft = dirprd(irrrht, irrepx)
         iabsoffo(irrrht+1) = iabsoffo(irrrht) + inumlfto(irrlft) * inumrhto(irrrht)
         iabsoffn(irrrht+1) = iabsoffn(irrrht) + inumlftn(irrlft) * inumrhtn(irrrht)
         iposrhto(1, irrrht) = 0_8; iposlfto(1, irrlft) = 0_8
         iposrhtn(1, irrrht) = 0_8; iposlftn(1, irrlft) = 0_8
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrrht)
            iposrhto(irr_r+1, irrrht) = iposrhto(irr_r, irrrht) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) * inumlfto(irrlft)
            iposrhtn(irr_r+1, irrrht) = iposrhtn(irr_r, irrrht) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) * inumlftn(irrlft)
            irr_l2 = dirprd(irr_r, irrlft)
            iposlfto(irr_r+1, irrlft) = iposlfto(irr_r, irrlft) + int(pop_in(1, irr_l2), 8) * pop_in(2, irr_r)
            iposlftn(irr_r+1, irrlft) = iposlftn(irr_r, irrlft) + int(pop_out(1, irr_l2), 8) * pop_out(2, irr_r)
         end do
      end do

      ! --- 2. 核心搬运循环 (将计算下发至GPU) ---
      do irrrht = 1, nirrep
         irrlft = dirprd(irrrht, irrepx)
         do irr4 = 1, nirrep
            irr3 = dirprd(irrrht, irr4)
            do irr2 = 1, nirrep
               irr1 = dirprd(irrlft, irr2)

               ! 获取当前块四个维度的长度
               n1 = pop_in(1, irr1); n2 = pop_in(2, irr2); n3 = pop_in(3, irr3); n4 = pop_in(4, irr4)
               if (n1 == 0 .or. n2 == 0 .or. n3 == 0 .or. n4 == 0) cycle ! 空块跳过

               n_p(1) = n1;   n_p(2) = n2;   n_p(3) = n3;   n_p(4) = n4
               irr_p(1) = irr1; irr_p(2) = irr2; irr_p(3) = irr3; irr_p(4) = irr4

               irr_lft_in = irrlft; irr_rht_in = irrrht
               irr_lft_out = dirprd(irr_p(p(1)), irr_p(p(2)))
               irr_rht_out = dirprd(irr_p(p(3)), irr_p(p(4)))

               ilfto_val = inumlfto(irr_lft_in)
               ilftn_val = inumlftn(irr_lft_out)

               ! 计算当前张量块在新旧内存中的基址
               base_old = iabsoffo(irr_rht_in)  + iposrhto(irr4, irr_rht_in)      + iposlfto(irr2, irr_lft_in)
               base_new = iabsoffn(irr_rht_out) + iposrhtn(irr_p(p(4)), irr_rht_out) + iposlftn(irr_p(p(2)), irr_lft_out)

               ! 计算物理内存寻址步长
               stride_in(4) = int(n3, 8) * ilfto_val; stride_in(3) = ilfto_val
               stride_in(2) = int(n1, 8);             stride_in(1) = 1_8
               stride_out(p(4)) = int(n_p(p(3)), 8) * ilftn_val; stride_out(p(3)) = ilftn_val
               stride_out(p(2)) = int(n_p(p(1)), 8);             stride_out(p(1)) = 1_8

               ! 调整基址以抵消 Fortran 中 i 从 1 开始带来的偏移
               base_old_adj = base_old - stride_in(4) - stride_in(3) - stride_in(2) - stride_in(1) + 1_8
               base_new_adj = base_new - stride_out(4) - stride_out(3) - stride_out(2) - stride_out(1) + 1_8

! ---------------------------------------------------------------------------------
! OpenACC 核心区域：
! 1. deviceptr(win, wout): 明确告诉编译器这两个是裸设备指针，直接用！
! 2. copyin(stride_in, stride_out): 将 CPU 刚算好的步长数据拷入 GPU 以供读取。
! 3. collapse(3): 将外三层循环合并，分配给 GPU 的 Grid/Block 级并行。
! ---------------------------------------------------------------------------------
!$acc parallel loop collapse(3) deviceptr(win, wout) &
!$acc copyin(stride_in, stride_out) private(i1, ioffold, ioffnew)
               do i4 = 1, n4
                  do i3 = 1, n3
                     do i2 = 1, n2
                        ! 在 GPU 线程内部独立计算内存指针
                        ioffold = base_old_adj + i4*stride_in(4) + i3*stride_in(3) + i2*stride_in(2)
                        ioffnew = base_new_adj + i4*stride_out(4) + i3*stride_out(3) + i2*stride_out(2)

                        ! 提示编译器内层循环使用 GPU 的 SIMT 向量化核心执行
                        !$acc loop vector
                        do i1 = 1, n1
                           wout(ioffnew + i1*stride_out(1)) = win(ioffold + i1)
                        end do
                     end do
                  end do
               end do
!$acc end parallel loop

            end do
         end do
      end do

      return
      end subroutine sstgend_opt_g

      ! =======================================================================
      ! 单精度版本：逻辑与双精度完全一致，仅改了接口变量 win 和 wout 的类型
      ! =======================================================================
      subroutine sstgens_opt_g(win, wout, pop1, pop2, pop3, pop4, irrepx, type)
      use mod_symm
      implicit none

      ! --- 接口参数声明 ---
      character(len=4), intent(in)  :: type
      real*4, intent(in)            :: win(*)   ! 纯设备端(GPU)单精度指针
      real*4, intent(out)           :: wout(*)  ! 纯设备端(GPU)单精度指针
      integer, intent(in)           :: pop1(8), pop2(8), pop3(8), pop4(8)
      integer, intent(in)           :: irrepx

      ! --- 局部变量 ---
      integer*8 :: inumlftn(9), inumrhtn(9), inumlfto(9), inumrhto(9)
      integer*8 :: iabsoffn(9), iabsoffo(9)
      integer*8 :: iposrhtn(9,9), iposlftn(9,9), iposrhto(9,9), iposlfto(9,9)
      integer   :: pop_in(4, 8), pop_out(4, 8)
      integer   :: p(4), n_p(4), irr_p(4)
      integer   :: i, irrep, irr_r, irr_l, irr_l2
      integer   :: irrrht, irrlft, irr4, irr3, irr2, irr1
      integer   :: irr_lft_in, irr_rht_in, irr_lft_out, irr_rht_out
      integer   :: i1, i2, i3, i4, n1, n2, n3, n4
      integer*8 :: ilfto_val, ilftn_val
      integer*8 :: base_old, base_new, base_old_adj, base_new_adj
      integer*8 :: stride_in(4), stride_out(4)
      integer*8 :: ioffold, ioffnew

      ! 前置解析映射逻辑与基址计算 (CPU端运行)
      do i = 1, 4
         p(i) = ichar(type(i:i)) - ichar('0')
         if (p(i) < 1 .or. p(i) > 4) stop
      end do

      pop_in(1, :) = pop1(:); pop_in(2, :) = pop2(:)
      pop_in(3, :) = pop3(:); pop_in(4, :) = pop4(:)
      do i = 1, 4
         pop_out(i, :) = pop_in(p(i), :)
      end do

      inumlftn = 0_8; inumrhtn = 0_8; inumlfto = 0_8; inumrhto = 0_8
      do irrep = 1, nirrep
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrep)
            inumlfto(irrep) = inumlfto(irrep) + int(pop_in(1, irr_l), 8) * pop_in(2, irr_r)
            inumrhto(irrep) = inumrhto(irrep) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r)
            inumlftn(irrep) = inumlftn(irrep) + int(pop_out(1, irr_l), 8) * pop_out(2, irr_r)
            inumrhtn(irrep) = inumrhtn(irrep) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r)
         end do
      end do

      iabsoffn(1) = 0_8; iabsoffo(1) = 0_8
      do irrrht = 1, nirrep
         irrlft = dirprd(irrrht, irrepx)
         iabsoffo(irrrht+1) = iabsoffo(irrrht) + inumlfto(irrlft) * inumrhto(irrrht)
         iabsoffn(irrrht+1) = iabsoffn(irrrht) + inumlftn(irrlft) * inumrhtn(irrrht)
         iposrhto(1, irrrht) = 0_8; iposlfto(1, irrlft) = 0_8
         iposrhtn(1, irrrht) = 0_8; iposlftn(1, irrlft) = 0_8
         do irr_r = 1, nirrep
            irr_l = dirprd(irr_r, irrrht)
            iposrhto(irr_r+1, irrrht) = iposrhto(irr_r, irrrht) + int(pop_in(3, irr_l), 8) * pop_in(4, irr_r) * inumlfto(irrlft)
            iposrhtn(irr_r+1, irrrht) = iposrhtn(irr_r, irrrht) + int(pop_out(3, irr_l), 8) * pop_out(4, irr_r) * inumlftn(irrlft)
            irr_l2 = dirprd(irr_r, irrlft)
            iposlfto(irr_r+1, irrlft) = iposlfto(irr_r, irrlft) + int(pop_in(1, irr_l2), 8) * pop_in(2, irr_r)
            iposlftn(irr_r+1, irrlft) = iposlftn(irr_r, irrlft) + int(pop_out(1, irr_l2), 8) * pop_out(2, irr_r)
         end do
      end do

      ! 核心 GPU 搬运循环
      do irrrht = 1, nirrep
         irrlft = dirprd(irrrht, irrepx)
         do irr4 = 1, nirrep
            irr3 = dirprd(irrrht, irr4)
            do irr2 = 1, nirrep
               irr1 = dirprd(irrlft, irr2)

               n1 = pop_in(1, irr1); n2 = pop_in(2, irr2); n3 = pop_in(3, irr3); n4 = pop_in(4, irr4)
               if (n1 == 0 .or. n2 == 0 .or. n3 == 0 .or. n4 == 0) cycle

               n_p(1) = n1;   n_p(2) = n2;   n_p(3) = n3;   n_p(4) = n4
               irr_p(1) = irr1; irr_p(2) = irr2; irr_p(3) = irr3; irr_p(4) = irr4

               irr_lft_in = irrlft; irr_rht_in = irrrht
               irr_lft_out = dirprd(irr_p(p(1)), irr_p(p(2)))
               irr_rht_out = dirprd(irr_p(p(3)), irr_p(p(4)))

               ilfto_val = inumlfto(irr_lft_in); ilftn_val = inumlftn(irr_lft_out)

               base_old = iabsoffo(irr_rht_in)  + iposrhto(irr4, irr_rht_in)      + iposlfto(irr2, irr_lft_in)
               base_new = iabsoffn(irr_rht_out) + iposrhtn(irr_p(p(4)), irr_rht_out) + iposlftn(irr_p(p(2)), irr_lft_out)

               stride_in(4) = int(n3, 8) * ilfto_val; stride_in(3) = ilfto_val
               stride_in(2) = int(n1, 8);             stride_in(1) = 1_8
               stride_out(p(4)) = int(n_p(p(3)), 8) * ilftn_val; stride_out(p(3)) = ilftn_val
               stride_out(p(2)) = int(n_p(p(1)), 8);             stride_out(p(1)) = 1_8

               base_old_adj = base_old - stride_in(4) - stride_in(3) - stride_in(2) - stride_in(1) + 1_8
               base_new_adj = base_new - stride_out(4) - stride_out(3) - stride_out(2) - stride_out(1) + 1_8

! 使用 deviceptr 声明输入输出为纯显存指针
!$acc parallel loop collapse(3) deviceptr(win, wout) &
!$acc copyin(stride_in, stride_out) private(i1, ioffold, ioffnew)
               do i4 = 1, n4
                  do i3 = 1, n3
                     do i2 = 1, n2
                        ioffold = base_old_adj + i4*stride_in(4) + i3*stride_in(3) + i2*stride_in(2)
                        ioffnew = base_new_adj + i4*stride_out(4) + i3*stride_out(3) + i2*stride_out(2)

                        ! 提示编译器内层循环向量化
                        !$acc loop vector
                        do i1 = 1, n1
                           wout(ioffnew + i1*stride_out(1)) = win(ioffold + i1)
                        end do
                     end do
                  end do
               end do
!$acc end parallel loop

            end do
         end do
      end do

      return
      end subroutine sstgens_opt_g
