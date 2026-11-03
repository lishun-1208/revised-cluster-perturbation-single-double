! =======================================================================
! 终极工业级张量引擎 - 封装层 (Wrapper) [222 模式]
! 
! 【功能说明】
! 将 2-index * 2-index -> 2-index (例如 Vpm * Vmr -> Vpr) 
! 映射并转换为 4-index * 4-index -> 4-index (vv2v_dpd_444) 统一处理。
!
! 【高速通道交错封装策略 (Interleaved Dummy Strategy)】
! 用户指定策略：V1(pm) -> V1(px, my)
! V1: (p, m) -> (p, x, m, y)
! V2: (m, r) -> (m, y, r, z)
! V3: (p, r) -> (p, x, r, z)
! 这种交错策略极其精妙，它恰好凑成了 444 引擎中最标准的连续收缩形式 (1234 * 3456 -> 1256)。
! 
! 【严格的 IOP 约束】
! 222 模式下，由于指标被拆分重组到不同的 Pair 中 (p 和 m 分离)，
! 原张量必须是完全展开态 (IOP=1)。如果传入 0 或 2，将被强制拦截并报错。
! =======================================================================
      subroutine vv2v_dpd_222(v1, p1_1, p1_2, iop1_12, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, sym3, type3, fact1, fact2)
        implicit none

        ! --- 接口参数声明 ---
        character(len=2), intent(in)  :: type1, type2, type3
        real*8, intent(in)            :: v1(*), v2(*)
        real*8, intent(inout)         :: v3(*)
        integer, intent(in)           :: p1_1(8), p1_2(8)
        integer, intent(in)           :: p2_1(8), p2_2(8)
        integer, intent(in)           :: p3_1(8), p3_2(8)
        integer, intent(in)           :: iop1_12, iop2_12, iop3_12
        integer, intent(in)           :: sym1, sym2, sym3
        real*8, intent(in)            :: fact1, fact2

        ! --- 内部扩展变量 ---
        integer                       :: pdum(8)
        character(len=4)              :: type1_444, type2_444, type3_444
        character(len=1)              :: c1, c2, c3
        integer                       :: i
        character(len=26), parameter  :: alphabet = 'abcdefghijklmnopqrstuvwxyz'
        logical                       :: found_c1, found_c2, found_c3

        ! --- 外部例程声明 ---
        external :: vv2v_dpd_444

        ! =================================================================
        ! 0. 强制检查 IOP 状态 (222模式必须为展开态)
        ! =================================================================
        if (iop1_12 /= 1 .or. iop2_12 /= 1 .or. iop3_12 /= 1) then
            write(6,*) " FATAL ERROR in vv2v_dpd_222: iop=0/2 not allowed. 2-index tensors must be uncompressed (IOP=1)."
            stop
        end if

        ! =================================================================
        ! 1. 构造虚拟维度 (大小恒为 1，仅存在于 Irrep 1)
        ! =================================================================
        pdum = 0
        pdum(1) = 1

        ! =================================================================
        ! 2. 自动寻找 3 个未被使用的字符，作为虚拟维度的名称标识
        !    - c1 (对应 x): V1 和 V3 的前置保留伴随维度
        !    - c2 (对应 y): V1 和 V2 的收缩伴随维度
        !    - c3 (对应 z): V2 和 V3 的后置保留伴随维度
        ! =================================================================
        c1 = ' '; c2 = ' '; c3 = ' '
        found_c1 = .false.; found_c2 = .false.; found_c3 = .false.
        
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
                else if (.not. found_c3) then
                    c3 = alphabet(i:i)
                    found_c3 = .true.
                    exit
                end if
            end if
        end do

        if (.not. found_c1 .or. .not. found_c2 .or. .not. found_c3) then
            write(6,*) " FATAL ERROR in vv2v_dpd_222: Cannot find enough available dummy characters."
            stop
        end if

        ! =================================================================
        ! 3. 交错式构造传递给 444 引擎的扩展 type 字符串
        !    例如：V1(pm) -> V1(px, my)   即 type1(1:1)//c1//type1(2:2)//c2
        !          V2(mr) -> V2(my, rz)   即 type2(1:1)//c2//type2(2:2)//c3
        !          V3(pr) -> V3(px, rz)   即 type3(1:1)//c1//type3(2:2)//c3
        !    完美凑成 444 的标准收缩形态！
        ! =================================================================
        type1_444 = type1(1:1) // c1 // type1(2:2) // c2
        type2_444 = type2(1:1) // c2 // type2(2:2) // c3
        type3_444 = type3(1:1) // c1 // type3(2:2) // c3

        ! =================================================================
        ! 4. 路由转发：调用 444 统一入口引擎
        !    原生指标使用真实的 pX_X，虚拟指标使用 pdum
        !    所有构成 pair 的内部 IOP 固定设为 1
        ! =================================================================
        call vv2v_dpd_444(v1, p1_1, pdum, 1, p1_2, pdum, 1, sym1, type1_444, &
                          v2, p2_1, pdum, 1, p2_2, pdum, 1, sym2, type2_444, &
                          v3, p3_1, pdum, 1, p3_2, pdum, 1, sym3, type3_444, fact1, fact2)

      end subroutine vv2v_dpd_222
