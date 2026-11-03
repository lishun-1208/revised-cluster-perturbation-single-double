! =======================================================================
! 终极工业级张量引擎 - 封装层 (Wrapper) [422 模式]
! 
! 【功能说明】
! 将 4-index * 2-index -> 2-index (例如 Vpqmn * Vmn -> Vpq) 
! 映射并转换为 4-index * 4-index -> 4-index (vv2v_dpd_444) 统一处理。
!
! 【零开销降维机制】
! 动态寻找不冲突的字符（如 'x','y'）补齐 type 字符串。
! 传入 pdum(1)=1 构造虚拟的 1x1 虚维度，对称性始终保持全对称 (Irrep 1)。
! 利用乘数为 1 的 Flat-mapping 属性，原封不动复用物理内存，完美无开销调用 444 引擎。
!
! 【BLAS 标准缩并系数 (Scaling Factors)】
!   fact1 (alpha): 乘法缩放系数。控制 V1 与 V2 收缩结果的权重。
!   fact2 (beta):  累加缩放系数。控制结果张量 V3 原有数据的保留权重。
!   物理方程: V3 = fact1 * (V1 * V2) + fact2 * V3
!   - 常用场景 1 (直接赋值) : fact1 = 1.0, fact2 = 0.0
!   - 常用场景 2 (叠加增量) : fact1 = 1.0, fact2 = 1.0
!   - 常用场景 3 (带系数叠加): fact1 = fact, fact2 = 1.0
! =======================================================================
      subroutine vv2v_dpd_422(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                              v2, p2_1, p2_2, iop2_12, sym2, type2, &
                              v3, p3_1, p3_2, iop3_12, sym3, type3, fact1, fact2)
        implicit none

        ! --- 接口参数声明 ---
        character(len=4), intent(in)  :: type1
        character(len=2), intent(in)  :: type2, type3
        real*8, intent(in)            :: v1(*), v2(*)
        real*8, intent(inout)         :: v3(*)
        integer, intent(in)           :: p1_1(8), p1_2(8), p1_3(8), p1_4(8)
        integer, intent(in)           :: p2_1(8), p2_2(8)
        integer, intent(in)           :: p3_1(8), p3_2(8)
        integer, intent(in)           :: iop1_12, iop1_34, iop2_12, iop3_12
        integer, intent(in)           :: sym1, sym2, sym3
        real*8, intent(in)            :: fact1, fact2

        ! --- 内部扩展变量 ---
        integer                       :: pdum(8)
        character(len=4)              :: type2_444, type3_444
        character(len=1)              :: c1, c2
        integer                       :: i
        character(len=26), parameter  :: alphabet = 'abcdefghijklmnopqrstuvwxyz'
        logical                       :: found_c1, found_c2

        ! --- 外部例程声明 ---
        external :: vv2v_dpd_444

        ! =================================================================
        ! 1. 构造虚拟维度 (大小为 1，仅存在于 Irrep 1)
        ! =================================================================
        pdum = 0
        pdum(1) = 1

        ! =================================================================
        ! 2. 自动寻找两个未被使用的字符，作为虚拟维度的名称标识
        !    保证调度器的自由度提取器不会产生自相交。
        ! =================================================================
        c1 = ' '
        found_c1 = .false.
        do i = 1, 26
            if (index(type1, alphabet(i:i)) == 0 .and. &
                index(type2, alphabet(i:i)) == 0 .and. &
                index(type3, alphabet(i:i)) == 0) then
                c1 = alphabet(i:i)
                found_c1 = .true.
                exit
            end if
        end do

        c2 = ' '
        found_c2 = .false.
        do i = 1, 26
            if (index(type1, alphabet(i:i)) == 0 .and. &
                index(type2, alphabet(i:i)) == 0 .and. &
                index(type3, alphabet(i:i)) == 0 .and. &
                alphabet(i:i) /= c1) then
                c2 = alphabet(i:i)
                found_c2 = .true.
                exit
            end if
        end do

        if (.not. found_c1 .or. .not. found_c2) then
            write(6,*) " FATAL ERROR in vv2v_dpd_422: Cannot find available dummy characters."
            stop
        end if

        ! =================================================================
        ! 3. 构造传递给 444 引擎的扩展 type 字符串 (追加 dummy 维度)
        ! =================================================================
        type2_444 = type2(1:2) // c1 // c2
        type3_444 = type3(1:2) // c1 // c2

        ! =================================================================
        ! 4. 路由转发：完美伪装并调用 444 统一入口引擎
        !    将 V2 扩展为 Vmnxy, V3 扩展为 Vpqxy
        !    pdum 作为虚拟的 p2_3, p2_4, p3_3, p3_4
        !    虚拟维度的 IOP 固定设为 1 (Full space)
        ! =================================================================
        call vv2v_dpd_444(v1, p1_1, p1_2, iop1_12, p1_3, p1_4, iop1_34, sym1, type1, &
                          v2, p2_1, p2_2, iop2_12, pdum, pdum, 1, sym2, type2_444, &
                          v3, p3_1, p3_2, iop3_12, pdum, pdum, 1, sym3, type3_444, fact1, fact2)

      end subroutine vv2v_dpd_422
