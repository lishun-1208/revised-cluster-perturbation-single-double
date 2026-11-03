! =======================================================================
! 终极工业级张量引擎 - 封装层 (Wrapper) [323 模式]
! 
! 【功能说明】
! 将 3-index * 2-index -> 3-index (例如 Vpqm * Vmr -> Vpqr) 
! 映射并转换为 4-index * 4-index -> 4-index (vv2v_dpd_444) 统一处理。
!
! 【架构约定】
! 323 模式下的 3 指标张量 (V1, V3)，默认前两个指标组成 Pair，第三个指标为 Independent。
!
! 【高速通道交错封装策略 (Interleaved Dummy Strategy)】
! V1(pqm) -> V1(pq, m*c1*)
! V2(mr)  -> V2(m*c1*, r*c2*)
! V3(pqr) -> V3(pq, r*c2*)
! 这种交错策略完美凑成了 444 引擎中最标准的连续收缩形式 (1234 * 3456 -> 1256)。
! 
! 【智能 IOP 探测墙 (Smart IOP Barrier)】
! 如果收缩前后 pair 保持绝对完整 (例如 pqm * mr -> pqr)，则放行 IOP=0/2。
! 如果 pair 发生颠倒、拆分或移位 (例如 pqm * mr -> rpq)，则强制要求 IOP=1。
! =======================================================================
      subroutine vv2v_dpd_323(v1, p1_1, p1_2, iop1_12, p1_3, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, p3_3, sym3, type3, fact1, fact2)
        implicit none

        ! --- 接口参数声明 ---
        character(len=3), intent(in)  :: type1, type3
        character(len=2), intent(in)  :: type2
        real*8, intent(in)            :: v1(*), v2(*)
        real*8, intent(inout)         :: v3(*)
        integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8)
        integer, intent(in)           :: p2_1(8), p2_2(8)
        integer, intent(in)           :: p3_1(8), p3_2(8), p3_3(8)
        integer, intent(in)           :: iop1_12, iop2_12, iop3_12
        integer, intent(in)           :: sym1, sym2, sym3
        real*8, intent(in)            :: fact1, fact2

        ! --- 内部扩展变量 ---
        integer                       :: pdum(8)
        character(len=4)              :: type1_444, type2_444, type3_444
        character(len=1)              :: c1, c2, char_c
        integer                       :: i
        character(len=26), parameter  :: alphabet = 'abcdefghijklmnopqrstuvwxyz'
        logical                       :: found_c1, found_c2, is_intact

        ! --- 外部例程声明 ---
        external :: vv2v_dpd_444

        ! =================================================================
        ! 1. 提取特征与物理完整性判定
        ! =================================================================
        ! 获取 V1 的收缩指标 (即第 3 个独立指标)
        char_c = type1(3:3)
        
        ! 确保该指标在 V3 中不存在 (否则说明收缩了 pair 中的指标)
        if (index(type3, char_c) /= 0) then
            write(6,*) "========================================================="
            write(6,*) " FATAL ERROR in vv2v_dpd_323 (Wrapper):"
            write(6,*) " 非法收缩维度！"
            write(6,*) " 323 模式强制要求 V1 的第 3 个指标 (", char_c, ") 为收缩基。"
            write(6,*) "========================================================="
            stop
        end if
        
        ! 判断 Pair (前两个字符) 是否原封不动地流向了 V3 的前两个字符
        is_intact = (type1(1:2) == type3(1:2))
        
        ! =================================================================
        ! 2. 智能 IOP 拦截墙
        ! =================================================================
        if (iop2_12 /= 1) then
            write(6,*) " FATAL ERROR in vv2v_dpd_323: 2-index tensor V2 must be uncompressed (IOP=1)."
            stop
        end if
        
        if (.not. is_intact) then
            if (iop1_12 /= 1 .or. iop3_12 /= 1) then
                write(6,*) "========================================================="
                write(6,*) " FATAL ERROR in vv2v_dpd_323 (Wrapper):"
                write(6,*) " 智能 IOP 空间校验拦截！"
                write(6,*) " 收缩前后 Pair 发生破坏或移位 (", type1(1:2), " -> ", type3(1:2), ")。"
                write(6,*) " 涉及重组乱序的张量强制要求输入/输出为解压态 (IOP=1)。"
                write(6,*) "========================================================="
                stop
            end if
        end if

        ! =================================================================
        ! 3. 构造虚拟维度 (大小恒为 1，仅存在于 Irrep 1)
        ! =================================================================
        pdum = 0
        pdum(1) = 1

        ! 自动寻找 2 个未被使用的字母，作为虚拟维度的名称标识
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
            write(6,*) " FATAL ERROR in vv2v_dpd_323: Cannot find enough available dummy characters."
            stop
        end if

        ! =================================================================
        ! 4. 交错式构造传递给 444 引擎的扩展 type 字符串
        ! =================================================================
        ! V1(pqm) -> V1(pq, m c1)
        type1_444 = type1(1:2) // type1(3:3) // c1
        
        ! V3(pqr) -> V3(pq, r c2)   (即使 V3 为 rpq, 此映射依然能使 444 引擎完美执行)
        type3_444 = type3(1:2) // type3(3:3) // c2
        
        ! 智能拼装 V2，确保它包含来自 V1 和 V3 的相关虚/实伴随指标
        ! 保证收缩块 (m c1) 始终成对出现
        if (type2(1:1) == char_c) then
            type2_444 = type2(1:1) // c1 // type2(2:2) // c2
        else if (type2(2:2) == char_c) then
            type2_444 = type2(1:1) // c2 // type2(2:2) // c1
        else
            write(6,*) " FATAL ERROR in vv2v_dpd_323: Contracted index '", char_c, "' not found in V2."
            stop
        end if

        ! =================================================================
        ! 5. 路由转发：调用 444 统一入口引擎
        ! 原生 Pair 指标正常传递，虚拟指标强制 IOP=1。
        ! =================================================================
        call vv2v_dpd_444(v1, p1_1, p1_2, iop1_12, p1_3, pdum, 1, sym1, type1_444, &
                          v2, p2_1, pdum, 1,       p2_2, pdum, 1, sym2, type2_444, &
                          v3, p3_1, p3_2, iop3_12, p3_3, pdum, 1, sym3, type3_444, fact1, fact2)

      end subroutine vv2v_dpd_323
