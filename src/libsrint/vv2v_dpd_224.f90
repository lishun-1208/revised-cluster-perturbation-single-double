! =======================================================================
! 终极工业级张量引擎 - 封装层 (Wrapper) [224 模式]
! 
! 【功能说明】
! 将 2-index * 2-index -> 4-index (例如 Vpq * Vrs -> Vpqrs，即外积) 
! 映射并转换为 4-index * 4-index -> 4-index (vv2v_dpd_444) 统一处理。
!
! 【零开销外积机制 (Outer Product)】
! 自动提取未使用的虚拟字符 (如 'x', 'y')。
! 构造映射: V1(pq) -> V1(pqxy), V2(rs) -> V2(xyrs), V3(pqrs) -> V3(pqrs)
! 
! 【IOP 完美继承特性】
! V1(pq) 占据 Pair 1，V2(rs) 占据 Pair 2。物理维度未被拆分，
! 因此该模式下原生支持 ioppq = 0, 1, 2 及 ioprs = 0, 1, 2 的任意压缩态！
! =======================================================================
      subroutine vv2v_dpd_224(v1, p1_1, p1_2, iop1_12, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1, fact2)
        implicit none

        ! --- 接口参数声明 ---
        character(len=2), intent(in)  :: type1, type2
        character(len=4), intent(in)  :: type3
        real*8, intent(in)            :: v1(*), v2(*)
        real*8, intent(inout)         :: v3(*)
        integer, intent(in)           :: p1_1(8), p1_2(8)
        integer, intent(in)           :: p2_1(8), p2_2(8)
        integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8), p3_4(8)
        integer, intent(in)           :: iop1_12, iop2_12, iop3_12, iop3_34
        integer, intent(in)           :: sym1, sym2, sym3
        real*8, intent(in)            :: fact1, fact2

        ! --- 内部扩展变量 ---
        integer                       :: pdum(8)
        character(len=4)              :: type1_444, type2_444
        character(len=1)              :: c1, c2
        integer                       :: i
        character(len=26), parameter  :: alphabet = 'abcdefghijklmnopqrstuvwxyz'
        logical                       :: found_c1, found_c2

        ! --- 外部例程声明 ---
        external :: vv2v_dpd_444

        ! =================================================================
        ! 1. 构造虚拟维度 (大小恒为 1，仅存在于 Irrep 1，全展开态 IOP=1)
        ! =================================================================
        pdum = 0
        pdum(1) = 1

        ! =================================================================
        ! 2. 自动寻找 2 个未被使用的字符，作为虚拟维度的名称标识
        !    - c1 (例如 x): 虚拟收缩维度 1
        !    - c2 (例如 y): 虚拟收缩维度 2
        ! =================================================================
        c1 = ' '; c2 = ' '
        found_c1 = .false.; found_c2 = .false.
        
        do i = 1, 26
            if (index(type1, alphabet(i:i)) == 0 .and. &
                index(type2, alphabet(i:i)) == 0 .and. &
                index(type3, alphabet(i:i)) == 0) then
                if (.not. found_c1) then
                    c1 = alphabet(i:i)
                    found_c1 = .true.
                else if (.not. found_c2) then
                    c2 = alphabet(i:i)
                    found_c2 = .true.
                    exit
                end if
            end if
        end do

        if (.not. found_c1 .or. .not. found_c2) then
            write(6,*) " FATAL ERROR in vv2v_dpd_224: Cannot find enough available dummy characters."
            stop
        end if

        ! =================================================================
        ! 3. 极简的 Pair 组合映射：
        !    将虚拟维度分别塞入 V1 的后半区 (Pair 2) 和 V2 的前半区 (Pair 1)。
        !    这样物理的 type1(pq) 和 type2(rs) 都不会被截断！
        ! =================================================================
        type1_444 = type1(1:2) // c1 // c2  ! 例如: pqxy
        type2_444 = c1 // c2 // type2(1:2)  ! 例如: xyrs

        ! =================================================================
        ! 4. 路由转发：调用 444 统一入口引擎
        !    V1 物理态传给 p1_1, p1_2 (IOP 原封不动)；虚拟态传给 p1_3, p1_4 (IOP=1)
        !    V2 虚拟态传给 p2_1, p2_2 (IOP=1)；物理态传给 p2_3, p2_4 (IOP 原封不动)
        ! =================================================================
        call vv2v_dpd_444(v1, p1_1, p1_2, iop1_12, pdum, pdum, 1, sym1, type1_444, &
                          v2, pdum, pdum, 1, p2_1, p2_2, iop2_12, sym2, type2_444, &
                          v3, p3_1, p3_2, iop3_12, p3_3, p3_4, iop3_34, sym3, type3, fact1, fact2)

      end subroutine vv2v_dpd_224
