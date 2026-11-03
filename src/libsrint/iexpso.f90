!================================================================================
! 模块说明: 高性能对称性积分展开模块 (iexpso_opt)
! 
! 功能描述:
!   该模块用于将压缩存储的积分 (如 p<q 或 r<=s) 展开为全矩阵格式 (p,q 或 r,s)。
!   支持单/双精度切换、CPU(OpenMP)/GPU(OpenACC) 硬件加速。
!
! 子程序命名规则:
!   - iexpso: 主驱动入口，根据 ioppq/ioprs 自动分发任务。
!   - exp_pq  系列: 处理 p < q  (对应 ioppq=0, 原 symexp2so)。
!   - exp_pq2 系列: 处理 p <= q (对应 ioppq=2, 原 exp4symso)。
!   - exp_rs  系列: 处理 r < s  (对应 ioprs=0, 原 symexpso)。
!   - exp_rs2 系列: 处理 r <= s (对应 ioprs=2, 原 exp3symso)。
!
! 核心安全特性:
!   1. 强制维度检查: 在展开某侧时，严格检查两侧维度是否一致 (popp == popq 或 popr == pops)。
!   2. FPE 崩溃保护: 针对对角线和非对角线进行解耦，添加 >=2 的边界锁，杜绝空迭代除零错误。
!================================================================================

!================================================================================
! 主驱动入口: iexpso
! 根据 ioppq 和 ioprs 自动分发到对应的单侧或双侧展开引擎
!================================================================================
subroutine iexpso(iin, iout, popp, popq, ioppq, iexppq, &
                  popr, pops, ioprs, iexprs, sym)
    use mod_symm
    use mod_iop
    implicit none
    integer, external :: irpdso

    integer, intent(in) :: popp(8), popq(8), popr(8), pops(8)
    integer, intent(in) :: ioppq, iexppq, ioprs, iexprs, sym
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)

    integer :: irrepr, irrepl, numl, numr, numl2, numr2
    integer :: nl_c, nr_c, nl_f, nr_f
    integer :: ioff, ioff1

    ! --- 控制字合法性检查 ---
    if (ioppq /= 0 .and. ioppq /= 1 .and. ioppq /= 2) stop 'fatal error: ioppq'
    if (ioprs /= 0 .and. ioprs /= 1 .and. ioprs /= 2) stop 'fatal error: ioprs'

    ! --- 维度一致性安全检查 (物理对角线防线) ---
    if (iexppq == 1 .and. ioppq /= 1 .and. any(popp /= popq)) stop 'fatal error: popp /= popq'
    if (iexprs == 1 .and. ioprs /= 1 .and. any(popr /= pops)) stop 'fatal error: popr /= pops'

    ! ========================================================================
    ! 路由分发 (Routing Dispatcher)
    ! ========================================================================
    
    ! 1. PQ 侧单边展开 (p<q -> p,q)
    if (iexppq == 1 .and. ioppq == 0 .and. iexprs == 0) then
        ioff = 1; ioff1 = 1
        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, sym)
            numl = irpdso(irrepl, popp, popq, ioppq); numr = irpdso(irrepr, popr, pops, ioprs)
            numl2 = irpdso(irrepl, popp, popq, iexppq)
            call exp_pq(irrepl, popp, numl2, numl, numr, iout(ioff1), iin(ioff))
            ioff = ioff + numl*numr*isd; ioff1 = ioff1 + numl2*numr*isd
        end do

    ! 2. RS 侧单边展开 (r<s -> r,s)
    else if (ioprs == 0 .and. iexprs == 1 .and. iexppq == 0) then
        ioff = 1; ioff1 = 1
        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, sym)
            numl = irpdso(irrepl, popp, popq, ioppq); numr = irpdso(irrepr, popr, pops, ioprs)
            numr2 = irpdso(irrepr, popr, pops, iexprs)
            call exp_rs(irrepr, popr, numl, iout(ioff1), iin(ioff))
            ioff = ioff + numl*numr*isd; ioff1 = ioff1 + numl*numr2*isd
        end do

    ! 3. PQ 侧单边展开 (p<=q -> p,q)
    else if (iexppq == 1 .and. ioppq == 2 .and. iexprs == 0) then
        ioff = 1; ioff1 = 1
        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, sym)
            numr = irpdso(irrepr, popr, pops, ioprs)
            numl = irpdso(irrepl, popp, popq, ioppq); numl2 = irpdso(irrepl, popp, popq, iexppq) 
            call exp_pq2(irrepl, popp, numl, numl2, numr, iin(ioff), iout(ioff1))
            ioff = ioff + numl*numr*isd; ioff1 = ioff1 + numl2*numr*isd
        end do

    ! 4. RS 侧单边展开 (r<=s -> r,s)
    else if (ioprs == 2 .and. iexprs == 1 .and. iexppq == 0) then
        ioff = 1; ioff1 = 1
        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, sym)
            numl = irpdso(irrepl, popp, popq, ioppq); numr = irpdso(irrepr, popr, pops, ioprs)
            numr2 = irpdso(irrepr, popr, pops, iexprs)
            call exp_rs2(irrepr, popr, numl, iin(ioff), iout(ioff1))
            ioff = ioff + numl*numr*isd; ioff1 = ioff1 + numl*numr2*isd
        end do

    ! ========================================================================
    ! 5. 双侧同步解压 (PQ & RS simultaneously: 00, 22, 02, 20 全兼容)
    ! ========================================================================
    else if (iexppq == 1 .and. ioppq /= 1 .and. iexprs == 1 .and. ioprs /= 1) then
        ioff = 1; ioff1 = 1
        do irrepr = 1, nirrep
            irrepl = dirprd(irrepr, sym)
            nl_f = irpdso(irrepl, popp, popq, 1);     nr_f = irpdso(irrepr, popr, pops, 1)
            nl_c = irpdso(irrepl, popp, popq, ioppq); nr_c = irpdso(irrepr, popr, pops, ioprs)
            
            ! 防御空块
            if (nl_f == 0 .or. nr_f == 0) cycle
            
            ! 统一分发至万能双侧引擎
            call exp_pqrs(ioppq, ioprs, irrepl, irrepr, popp, popr, nl_f, nr_f, nl_c, nr_c, iin(ioff), iout(ioff1))
            
            ! 地址步长缩放
            ioff = ioff + nl_c*nr_c*isd; ioff1 = ioff1 + nl_f*nr_f*isd
        end do

    else
        write(6,*) 'fatal error: iexpso 必须且只能展开其中一侧，或两侧控制字同时匹配双侧展开。'
        stop
    end if
end subroutine iexpso
!================================================================================
! 子程序系列: EXP_PQ (包含 _s, _d, _sg, _dg 四种精度与设备版本)
! 功能描述: 将压缩存储的严格反对称积分 I(p < q, L) 展开为全矩阵格式 I(p, q, L)。
! 数学映射关系:
!   - 下三角映射: IOUT(p, q) =  IIN(p < q)
!   - 上三角映射: IOUT(q, p) = -IIN(p < q)  (反对称展开，带负号)
!   - 对角线置零: IOUT(q, q) =  0.0
!================================================================================
subroutine exp_pq(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
    use mod_iop
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)

    if (prec == 8) then
        if (igpu == 0) then
            call exp_pq_d(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
        else if (igpu == 1) then
            call exp_pq_dg(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
        end if
    else if (prec == 4) then
        if (igpu == 0) then
            call exp_pq_s(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
        else if (igpu == 1) then
            call exp_pq_sg(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
        end if
    end if
end subroutine exp_pq

! -------------------------------------------------------------------------------
! 单精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq_s(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 0)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            ! FPE 安全锁：大于等于 2 才能执行三角展开
            if (numq >= 2) then
                !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
                do rs = 1, numrs
                    do q = 2, numq
                        base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * (q - 2) / 2
                        base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                        base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                        do p = 1, q - 1
                            iout(base_out + p) = iin(base_in + p)
                            iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                        end do
                    end do
                end do
            end if
            
            ! 对角线必定置零
            !$omp parallel do private(q) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    iout((rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump + q) = 0.0
                end do
            end do
        else
            !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq_s

! -------------------------------------------------------------------------------
! 双精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq_d(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ, NUMRS
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 0)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            if (numq >= 2) then
                !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
                do rs = 1, numrs
                    do q = 2, numq
                        base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * (q - 2) / 2
                        base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                        base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                        do p = 1, q - 1
                            iout(base_out + p) = iin(base_in + p)
                            iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                        end do
                    end do
                end do
            end if
            
            !$omp parallel do private(q) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    iout((rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump + q) = 0.0d0
                end do
            end do
        else
            !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq_d

! -------------------------------------------------------------------------------
! 单精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq_sg(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
    use mod_symm
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 0)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            if (numq >= 2) then
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
                do rs = 1, numrs
                    do q = 2, numq
                        base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * (q - 2) / 2
                        base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                        base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                        !$acc loop independent
                        do p = 1, q - 1
                            iout(base_out + p) = iin(base_in + p)
                            iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                        end do
                    end do
                end do
            end if
            
            !$acc parallel loop collapse(2) deviceptr(iin, iout)
            do rs = 1, numrs
                do q = 1, numq
                    iout((rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump + q) = 0.0
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq_sg

! -------------------------------------------------------------------------------
! 双精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq_dg(IRREPPQ, POPP, NUMPQ2, NUMPQ, NUMRS, IOUT, IIN)
    use mod_symm
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ2, NUMPQ, NUMRS
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 0)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            if (numq >= 2) then
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
                do rs = 1, numrs
                    do q = 2, numq
                        base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * (q - 2) / 2
                        base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                        base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                        !$acc loop independent
                        do p = 1, q - 1
                            iout(base_out + p) = iin(base_in + p)
                            iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                        end do
                    end do
                end do
            end if
            
            !$acc parallel loop collapse(2) deviceptr(iin, iout)
            do rs = 1, numrs
                do q = 1, numq
                    iout((rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump + q) = 0.0d0
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = -iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq_dg


!================================================================================
! 子程序系列: EXP_PQ2 (包含 _s, _d, _sg, _dg 四种精度与设备版本)
! 功能描述: 将对称积分 I(p <= q, L) 展开为全矩阵格式 I(p, q, L)。(上下三角均正数)
!================================================================================
subroutine exp_pq2(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
    use mod_iop
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ, NUMPQ2, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)

    if (prec == 8) then
        if (igpu == 0) then
            call exp_pq2_d(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
        else if (igpu == 1) then
            call exp_pq2_dg(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
        end if
    else if (prec == 4) then
        if (igpu == 0) then
            call exp_pq2_s(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
        else if (igpu == 1) then
            call exp_pq2_sg(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
        end if
    end if
end subroutine exp_pq2

! -------------------------------------------------------------------------------
! 单精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq2_s(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ, NUMPQ2, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 2)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            !$omp parallel do private(q, p, base_in, base_out) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in  = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    do p = 1, q
                        iout(base_out + p) = iin(base_in + p)
                    end do
                end do
            end do

            !$omp parallel do private(q, p, base_in, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, q
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        else
            !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq2_s

! -------------------------------------------------------------------------------
! 双精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq2_d(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ, NUMPQ2, NUMRS
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 2)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            !$omp parallel do private(q, p, base_in, base_out) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in  = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    do p = 1, q
                        iout(base_out + p) = iin(base_in + p)
                    end do
                end do
            end do

            !$omp parallel do private(q, p, base_in, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, q
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        else
            !$omp parallel do private(q, p, base_in, base_out, base_out2) collapse(2) schedule(static)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq2_d

! -------------------------------------------------------------------------------
! 单精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq2_sg(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
    use mod_symm
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ, NUMPQ2, NUMRS
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 2)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out)
            do rs = 1, numrs
                do q = 1, numq
                    base_in  = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    !$acc loop independent
                    do p = 1, q
                        iout(base_out + p) = iin(base_in + p)
                    end do
                end do
            end do

            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, q
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq2_sg

! -------------------------------------------------------------------------------
! 双精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_pq2_dg(IRREPPQ, POPP, NUMPQ, NUMPQ2, NUMRS, IIN, IOUT)
    use mod_symm
    implicit none
    integer, intent(in) :: popp(8), IRREPPQ, NUMPQ, NUMPQ2, NUMRS
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irrepp, irrepq, nump, numq, p, q, rs
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irrepq = 1, nirrep
        irrepp = dirprd(irrepq, irreppq)
        if (irrepp > irrepq) cycle
        numq = popp(irrepq)
        nump = popp(irrepp)
        if (numq == 0 .or. nump == 0 .or. numrs == 0) cycle

        vstarto  = isymoffso(irrepq, irreppq, popp, popp, 2)
        vstartn  = isymoffso(irrepq, irreppq, popp, popp, 1)
        vstartn2 = isymoffso(irrepp, irreppq, popp, popp, 1)

        if (irrepp == irrepq) then
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out)
            do rs = 1, numrs
                do q = 1, numq
                    base_in  = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    !$acc loop independent
                    do p = 1, q
                        iout(base_out + p) = iin(base_in + p)
                    end do
                end do
            end do

            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * q / 2
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, q
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do rs = 1, numrs
                do q = 1, numq
                    base_in   = (rs - 1) * numpq  + vstarto - 1 + (q - 1) * nump
                    base_out  = (rs - 1) * numpq2 + vstartn - 1 + (q - 1) * nump
                    base_out2 = (rs - 1) * numpq2 + vstartn2 - 1 + q
                    !$acc loop independent
                    do p = 1, nump
                        iout(base_out + p) = iin(base_in + p)
                        iout(base_out2 + (p - 1) * numq) = iin(base_in + p)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_pq2_dg


!================================================================================
! 子程序系列: EXP_RS (包含 _s, _d, _sg, _dg 四种精度与设备版本)
! 功能描述: 将压缩存储的严格反对称积分 I(L, r < s) 展开为全矩阵格式 I(L, r, s)。
! 数学映射关系:
!   - 下三角映射: IOUT(r, s) =  IIN(r < s)
!   - 上三角映射: IOUT(s, r) = -IIN(r < s)  (反对称展开，带负号)
!   - 对角线置零: IOUT(s, s) =  0.0
!================================================================================
subroutine exp_rs(sym, POPR, NUMPQ, IOUT, IIN)
    use mod_iop
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)

    if (prec == 8) then
        if (igpu == 0) then
            call exp_rs_d(sym, POPR, NUMPQ, IOUT, IIN)
        else if (igpu == 1) then
            call exp_rs_dg(sym, POPR, NUMPQ, IOUT, IIN)
        end if
    else if (prec == 4) then
        if (igpu == 0) then
            call exp_rs_s(sym, POPR, NUMPQ, IOUT, IIN)
        else if (igpu == 1) then
            call exp_rs_sg(sym, POPR, NUMPQ, IOUT, IIN)
        end if
    end if
end subroutine exp_rs

! -------------------------------------------------------------------------------
! 单精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs_s(sym, POPR, NUMPQ, IOUT, IIN)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 0)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            ! FPE 安全锁：大于等于 2 才能执行三角展开
            if (nums >= 2) then
                !$omp parallel do private(r, pq, base_in, base_out, base_out2) schedule(static)
                do s = 2, nums
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * (s - 2) / 2 + r - 1) * numpq
                        base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        do pq = 1, numpq
                            iout(base_out + pq) = iin(base_in + pq)
                            iout(base_out2 + pq) = -iin(base_in + pq)
                        end do
                    end do
                end do
            end if
            
            ! 对角线置零
            !$omp parallel do private(pq) collapse(2) schedule(static)
            do s = 1, nums
                do pq = 1, numpq
                    iout((vstartn - 1 + (s - 1) * numr + s - 1) * numpq + pq) = 0.0
                end do
            end do
        else
            !$omp parallel do private(r, pq, base_in, base_out, base_out2) collapse(2) schedule(static)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = -iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs_s

! -------------------------------------------------------------------------------
! 双精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs_d(sym, POPR, NUMPQ, IOUT, IIN)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 0)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            if (nums >= 2) then
                !$omp parallel do private(r, pq, base_in, base_out, base_out2) schedule(static)
                do s = 2, nums
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * (s - 2) / 2 + r - 1) * numpq
                        base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        do pq = 1, numpq
                            iout(base_out + pq) = iin(base_in + pq)
                            iout(base_out2 + pq) = -iin(base_in + pq)
                        end do
                    end do
                end do
            end if
            
            !$omp parallel do private(pq) collapse(2) schedule(static)
            do s = 1, nums
                do pq = 1, numpq
                    iout((vstartn - 1 + (s - 1) * numr + s - 1) * numpq + pq) = 0.0d0
                end do
            end do
        else
            !$omp parallel do private(r, pq, base_in, base_out, base_out2) collapse(2) schedule(static)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = -iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs_d

! -------------------------------------------------------------------------------
! 单精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs_sg(sym, POPR, NUMPQ, IOUT, IIN)
    use mod_symm
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 0)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            if (nums >= 2) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 2, nums
                    !$acc loop private(base_in, base_out, base_out2)
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * (s - 2) / 2 + r - 1) * numpq
                        base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        !$acc loop independent
                        do pq = 1, numpq
                            iout(base_out + pq) = iin(base_in + pq)
                            iout(base_out2 + pq) = -iin(base_in + pq)
                        end do
                    end do
                end do
            end if
            
            !$acc parallel loop collapse(2) deviceptr(iin, iout)
            do s = 1, nums
                do pq = 1, numpq
                    iout((vstartn - 1 + (s - 1) * numr + s - 1) * numpq + pq) = 0.0
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = -iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs_sg

! -------------------------------------------------------------------------------
! 双精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs_dg(sym, POPR, NUMPQ, IOUT, IIN)
    use mod_symm
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 0)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            if (nums >= 2) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 2, nums
                    !$acc loop private(base_in, base_out, base_out2)
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * (s - 2) / 2 + r - 1) * numpq
                        base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        !$acc loop independent
                        do pq = 1, numpq
                            iout(base_out + pq) = iin(base_in + pq)
                            iout(base_out2 + pq) = -iin(base_in + pq)
                        end do
                    end do
                end do
            end if
            
            !$acc parallel loop collapse(2) deviceptr(iin, iout)
            do s = 1, nums
                do pq = 1, numpq
                    iout((vstartn - 1 + (s - 1) * numr + s - 1) * numpq + pq) = 0.0d0
                end do
            end do
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = -iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs_dg


!================================================================================
! 子程序系列: EXP_RS2 (包含 _s, _d, _sg, _dg 四种精度与设备版本)
! 功能描述: 将包含对角线的下三角对称张量 I(L, r <= s) 展开为全矩阵 I(L, r, s)。
!================================================================================
subroutine exp_rs2(sym, POPR, NUMPQ, IIN, IOUT)
    use mod_iop
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)

    if (prec == 8) then
        if (igpu == 0) then
            call exp_rs2_d(sym, POPR, NUMPQ, IIN, IOUT)
        else if (igpu == 1) then
            call exp_rs2_dg(sym, POPR, NUMPQ, IIN, IOUT)
        end if
    else if (prec == 4) then
        if (igpu == 0) then
            call exp_rs2_s(sym, POPR, NUMPQ, IIN, IOUT)
        else if (igpu == 1) then
            call exp_rs2_sg(sym, POPR, NUMPQ, IIN, IOUT)
        end if
    end if
end subroutine exp_rs2

! -------------------------------------------------------------------------------
! 单精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs2_s(sym, POPR, NUMPQ, IIN, IOUT)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 2)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            !$omp parallel do private(r, pq, base_in, base_out) schedule(static)
            do s = 1, nums
                do r = 1, s
                    base_in  = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                    base_out = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                    end do
                end do
            end do
            
            if (nums >= 2) then
                !$omp parallel do private(r, pq, base_in, base_out2) schedule(static)
                do s = 2, nums
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        do pq = 1, numpq
                            iout(base_out2 + pq) = iin(base_in + pq)
                        end do
                    end do
                end do
            end if
        else
            !$omp parallel do private(r, pq, base_in, base_out, base_out2) collapse(2) schedule(static)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs2_s

! -------------------------------------------------------------------------------
! 双精度 CPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs2_d(sym, POPR, NUMPQ, IIN, IOUT)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 2)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            !$omp parallel do private(r, pq, base_in, base_out) schedule(static)
            do s = 1, nums
                do r = 1, s
                    base_in  = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                    base_out = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                    end do
                end do
            end do
            
            if (nums >= 2) then
                !$omp parallel do private(r, pq, base_in, base_out2) schedule(static)
                do s = 2, nums
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        do pq = 1, numpq
                            iout(base_out2 + pq) = iin(base_in + pq)
                        end do
                    end do
                end do
            end if
        else
            !$omp parallel do private(r, pq, base_in, base_out, base_out2) collapse(2) schedule(static)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs2_d

! -------------------------------------------------------------------------------
! 单精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs2_sg(sym, POPR, NUMPQ, IIN, IOUT)
    use mod_symm
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 2)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            !$acc parallel loop deviceptr(iin, iout)
            do s = 1, nums
                !$acc loop private(base_in, base_out)
                do r = 1, s
                    base_in  = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                    base_out = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                    end do
                end do
            end do
            
            if (nums >= 2) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 2, nums
                    !$acc loop private(base_in, base_out2)
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        !$acc loop independent
                        do pq = 1, numpq
                            iout(base_out2 + pq) = iin(base_in + pq)
                        end do
                    end do
                end do
            end if
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs2_sg

! -------------------------------------------------------------------------------
! 双精度 GPU 版本
! -------------------------------------------------------------------------------
subroutine exp_rs2_dg(sym, POPR, NUMPQ, IIN, IOUT)
    use mod_symm
    implicit none
    integer, intent(in) :: popr(8), sym, numpq
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, numr, nums, r, s, pq
    integer :: vstarto, vstartn, vstartn2
    integer :: base_in, base_out, base_out2
    integer, external :: isymoffso

    do irreps = 1, nirrep
        irrepr = dirprd(sym, irreps)
        if (irrepr > irreps) cycle
        numr = popr(irrepr)
        nums = popr(irreps)
        if (numr == 0 .or. nums == 0 .or. numpq == 0) cycle

        vstarto  = isymoffso(irreps, sym, popr, popr, 2)
        vstartn  = isymoffso(irreps, sym, popr, popr, 1)
        vstartn2 = isymoffso(irrepr, sym, popr, popr, 1)

        if (irrepr == irreps) then
            !$acc parallel loop deviceptr(iin, iout)
            do s = 1, nums
                !$acc loop private(base_in, base_out)
                do r = 1, s
                    base_in  = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                    base_out = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                    end do
                end do
            end do
            
            if (nums >= 2) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 2, nums
                    !$acc loop private(base_in, base_out2)
                    do r = 1, s - 1
                        base_in   = (vstarto - 1 + (s - 1) * s / 2 + r - 1) * numpq
                        base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                        !$acc loop independent
                        do pq = 1, numpq
                            iout(base_out2 + pq) = iin(base_in + pq)
                        end do
                    end do
                end do
            end if
        else
            !$acc parallel loop collapse(2) deviceptr(iin, iout) private(base_in, base_out, base_out2)
            do s = 1, nums
                do r = 1, numr
                    base_in   = (vstarto - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out  = (vstartn - 1 + (s - 1) * numr + r - 1) * numpq
                    base_out2 = (vstartn2 - 1 + (r - 1) * nums + s - 1) * numpq
                    !$acc loop independent
                    do pq = 1, numpq
                        iout(base_out + pq) = iin(base_in + pq)
                        iout(base_out2 + pq) = iin(base_in + pq)
                    end do
                end do
            end do
        end if
    end do
end subroutine exp_rs2_dg
! ==============================================================================
! 统一双侧展开包装器 (Unified Wrapper for Dual-Expansion)
! ==============================================================================
subroutine exp_pqrs(ioppq, ioprs, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout)
    use mod_iop
    implicit none
    integer, intent(in) :: ioppq, ioprs, sympq, symrs, numpq_f, numrs_f, numpq_c, numrs_c, popp(8), popr(8)
    real*4, intent(in)  :: iin(*); real*4, intent(out) :: iout(*)
    
    if (prec == 8) then
        if (igpu == 0) then; call exp_pqrs_d(ioppq, ioprs, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout); endif
        if (igpu == 1) then; call exp_pqrs_dg(ioppq, ioprs, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout); endif
    else
        if (igpu == 0) then; call exp_pqrs_s(ioppq, ioprs, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout); endif
        if (igpu == 1) then; call exp_pqrs_sg(ioppq, ioprs, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout); endif
    endif
end subroutine exp_pqrs


! ==============================================================================
! 万能全量解压底层核心 (Master Unified Engine for PQRS)
! ==============================================================================

! -------------------------------------------------------------------------------
! 单精度 CPU 核心 (_s)
! -------------------------------------------------------------------------------
subroutine exp_pqrs_s(IOPPQ, IOPRS, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: IOPPQ, IOPRS, sympq, symrs, numpq_f, numrs_f, numpq_c, numrs_c
    integer, intent(in) :: popp(8), popr(8)
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    
    integer :: irreps, irrepr, irrepq, irrepp, nums, numr, numq, nump
    integer :: rs_off_c, rs_n1, rs_n2, pq_off_c, pq_n1, pq_n2
    integer :: s, r, q, p, i, base_in, b11, b21, b12, b22
    real*4  :: spq, srs
    integer, external :: isymoffso
    
    spq = 1.0; if (IOPPQ == 0) spq = -1.0
    srs = 1.0; if (IOPRS == 0) srs = -1.0

    ! 【对角线安全锁】: 提前将当前块全面清零，配合 cycle 完美消除幽灵垃圾
    !$omp parallel do schedule(static)
    do i = 1, numpq_f * numrs_f; iout(i) = 0.0; end do

    do irreps = 1, nirrep
        irrepr = dirprd(irreps, symrs)
        if (irrepr > irreps) cycle
        nums = popr(irreps); numr = popr(irrepr)
        if (nums == 0 .or. numr == 0) cycle

        rs_off_c = isymoffso(irreps, symrs, popr, popr, IOPRS)
        rs_n1 = isymoffso(irreps, symrs, popr, popr, 1); rs_n2 = isymoffso(irrepr, symrs, popr, popr, 1)

        do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, sympq)
            if (irrepp > irrepq) cycle
            numq = popp(irrepq); nump = popp(irrepp)
            if (numq == 0 .or. nump == 0) cycle

            pq_off_c = isymoffso(irrepq, sympq, popp, popp, IOPPQ)
            pq_n1 = isymoffso(irrepq, sympq, popp, popp, 1); pq_n2 = isymoffso(irrepp, sympq, popp, popp, 1)

            if (irrepr == irreps .and. irrepp == irrepq) then
                !$omp parallel do private(r, q, p, base_in, b11, b21, b12, b22) schedule(static)
                do s = 1, nums
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr == irreps .and. irrepp < irrepq) then
                !$omp parallel do private(r, q, p, base_in, b11, b21, b12, b22) schedule(static)
                do s = 1, nums
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, nump
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                base_in = base_in + pq_off_c - 1 + (q-1)*nump + p
                                
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr < irreps .and. irrepp == irrepq) then
                !$omp parallel do private(q, p, base_in, b11, b21, b12, b22) collapse(2) schedule(static)
                do s = 1, nums
                    do r = 1, numr
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else
                !$omp parallel do private(q, p, base_in, b11, b21, b12, b22) collapse(2) schedule(static)
                do s = 1, nums
                    do r = 1, numr
                        do q = 1, numq
                            do p = 1, nump
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c + pq_off_c - 1 + (q-1)*nump + p
                                
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            end if
        end do
    end do
end subroutine exp_pqrs_s

! -------------------------------------------------------------------------------
! 双精度 CPU 核心 (_d)
! -------------------------------------------------------------------------------
subroutine exp_pqrs_d(IOPPQ, IOPRS, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout)
    use mod_symm
    use omp_lib
    implicit none
    integer, intent(in) :: IOPPQ, IOPRS, sympq, symrs, numpq_f, numrs_f, numpq_c, numrs_c
    integer, intent(in) :: popp(8), popr(8)
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, irrepq, irrepp, nums, numr, numq, nump
    integer :: rs_off_c, rs_n1, rs_n2, pq_off_c, pq_n1, pq_n2
    integer :: s, r, q, p, i, base_in, b11, b21, b12, b22
    real*8  :: spq, srs
    integer, external :: isymoffso
    
    spq = 1.0d0; if (IOPPQ == 0) spq = -1.0d0
    srs = 1.0d0; if (IOPRS == 0) srs = -1.0d0

    !$omp parallel do schedule(static)
    do i = 1, numpq_f * numrs_f; iout(i) = 0.0d0; end do

    do irreps = 1, nirrep
        irrepr = dirprd(irreps, symrs)
        if (irrepr > irreps) cycle
        nums = popr(irreps); numr = popr(irrepr)
        if (nums == 0 .or. numr == 0) cycle

        rs_off_c = isymoffso(irreps, symrs, popr, popr, IOPRS)
        rs_n1 = isymoffso(irreps, symrs, popr, popr, 1); rs_n2 = isymoffso(irrepr, symrs, popr, popr, 1)

        do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, sympq)
            if (irrepp > irrepq) cycle
            numq = popp(irrepq); nump = popp(irrepp)
            if (numq == 0 .or. nump == 0) cycle

            pq_off_c = isymoffso(irrepq, sympq, popp, popp, IOPPQ)
            pq_n1 = isymoffso(irrepq, sympq, popp, popp, 1); pq_n2 = isymoffso(irrepp, sympq, popp, popp, 1)

            if (irrepr == irreps .and. irrepp == irrepq) then
                !$omp parallel do private(r, q, p, base_in, b11, b21, b12, b22) schedule(static)
                do s = 1, nums
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr == irreps .and. irrepp < irrepq) then
                !$omp parallel do private(r, q, p, base_in, b11, b21, b12, b22) schedule(static)
                do s = 1, nums
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, nump
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                base_in = base_in + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr < irreps .and. irrepp == irrepq) then
                !$omp parallel do private(q, p, base_in, b11, b21, b12, b22) collapse(2) schedule(static)
                do s = 1, nums
                    do r = 1, numr
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else
                !$omp parallel do private(q, p, base_in, b11, b21, b12, b22) collapse(2) schedule(static)
                do s = 1, nums
                    do r = 1, numr
                        do q = 1, numq
                            do p = 1, nump
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            end if
        end do
    end do
end subroutine exp_pqrs_d

! -------------------------------------------------------------------------------
! 单精度 GPU 核心 (_sg)
! -------------------------------------------------------------------------------
subroutine exp_pqrs_sg(IOPPQ, IOPRS, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout)
    use mod_symm
    implicit none
    integer, intent(in) :: IOPPQ, IOPRS, sympq, symrs, numpq_f, numrs_f, numpq_c, numrs_c
    integer, intent(in) :: popp(8), popr(8)
    real*4, intent(in)  :: iin(*)
    real*4, intent(out) :: iout(*)
    integer :: irreps, irrepr, irrepq, irrepp, nums, numr, numq, nump
    integer :: rs_off_c, rs_n1, rs_n2, pq_off_c, pq_n1, pq_n2
    integer :: s, r, q, p, i, base_in, b11, b21, b12, b22
    real*4  :: spq, srs
    integer, external :: isymoffso
    
    spq = 1.0; if (IOPPQ == 0) spq = -1.0
    srs = 1.0; if (IOPRS == 0) srs = -1.0

    !$acc parallel loop deviceptr(iout)
    do i = 1, numpq_f * numrs_f; iout(i) = 0.0; end do

    do irreps = 1, nirrep
        irrepr = dirprd(irreps, symrs)
        if (irrepr > irreps) cycle
        nums = popr(irreps); numr = popr(irrepr)
        if (nums == 0 .or. numr == 0) cycle
        rs_off_c = isymoffso(irreps, symrs, popr, popr, IOPRS)
        rs_n1 = isymoffso(irreps, symrs, popr, popr, 1); rs_n2 = isymoffso(irrepr, symrs, popr, popr, 1)

        do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, sympq)
            if (irrepp > irrepq) cycle
            numq = popp(irrepq); nump = popp(irrepp)
            if (numq == 0 .or. nump == 0) cycle
            pq_off_c = isymoffso(irrepq, sympq, popp, popp, IOPPQ)
            pq_n1 = isymoffso(irrepq, sympq, popp, popp, 1); pq_n2 = isymoffso(irrepp, sympq, popp, popp, 1)

            if (irrepr == irreps .and. irrepp == irrepq) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 1, nums
                    !$acc loop private(q, p, base_in, b11, b21, b12, b22)
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr == irreps .and. irrepp < irrepq) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 1, nums
                    !$acc loop private(q, p, base_in, b11, b21, b12, b22)
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, nump
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                base_in = base_in + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr < irreps .and. irrepp == irrepq) then
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(q, p, base_in, b11, b21, b12, b22)
                do s = 1, nums
                    do r = 1, numr
                        !$acc loop independent
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(q, p, base_in, b11, b21, b12, b22)
                do s = 1, nums
                    do r = 1, numr
                        !$acc loop independent
                        do q = 1, numq
                            do p = 1, nump
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            end if
        end do
    end do
end subroutine exp_pqrs_sg

! -------------------------------------------------------------------------------
! 双精度 GPU 核心 (_dg)
! -------------------------------------------------------------------------------
subroutine exp_pqrs_dg(IOPPQ, IOPRS, sympq, symrs, popp, popr, numpq_f, numrs_f, numpq_c, numrs_c, iin, iout)
    use mod_symm
    implicit none
    integer, intent(in) :: IOPPQ, IOPRS, sympq, symrs, numpq_f, numrs_f, numpq_c, numrs_c
    integer, intent(in) :: popp(8), popr(8)
    real*8, intent(in)  :: iin(*)
    real*8, intent(out) :: iout(*)
    integer :: irreps, irrepr, irrepq, irrepp, nums, numr, numq, nump
    integer :: rs_off_c, rs_n1, rs_n2, pq_off_c, pq_n1, pq_n2
    integer :: s, r, q, p, i, base_in, b11, b21, b12, b22
    real*8  :: spq, srs
    integer, external :: isymoffso
    
    spq = 1.0d0; if (IOPPQ == 0) spq = -1.0d0
    srs = 1.0d0; if (IOPRS == 0) srs = -1.0d0

    !$acc parallel loop deviceptr(iout)
    do i = 1, numpq_f * numrs_f; iout(i) = 0.0d0; end do

    do irreps = 1, nirrep
        irrepr = dirprd(irreps, symrs)
        if (irrepr > irreps) cycle
        nums = popr(irreps); numr = popr(irrepr)
        if (nums == 0 .or. numr == 0) cycle
        rs_off_c = isymoffso(irreps, symrs, popr, popr, IOPRS)
        rs_n1 = isymoffso(irreps, symrs, popr, popr, 1); rs_n2 = isymoffso(irrepr, symrs, popr, popr, 1)

        do irrepq = 1, nirrep
            irrepp = dirprd(irrepq, sympq)
            if (irrepp > irrepq) cycle
            numq = popp(irrepq); nump = popp(irrepp)
            if (numq == 0 .or. nump == 0) cycle
            pq_off_c = isymoffso(irrepq, sympq, popp, popp, IOPPQ)
            pq_n1 = isymoffso(irrepq, sympq, popp, popp, 1); pq_n2 = isymoffso(irrepp, sympq, popp, popp, 1)

            if (irrepr == irreps .and. irrepp == irrepq) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 1, nums
                    !$acc loop private(q, p, base_in, b11, b21, b12, b22)
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr == irreps .and. irrepp < irrepq) then
                !$acc parallel loop deviceptr(iin, iout)
                do s = 1, nums
                    !$acc loop private(q, p, base_in, b11, b21, b12, b22)
                    do r = 1, s
                        if (IOPRS == 0 .and. r == s) cycle
                        do q = 1, numq
                            do p = 1, nump
                                if (IOPRS == 0) then; base_in = (rs_off_c - 1 + (s-1)*(s-2)/2 + r-1) * numpq_c
                                else;                 base_in = (rs_off_c - 1 + (s-1)*s/2 + r-1) * numpq_c; endif
                                base_in = base_in + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else if (irrepr < irreps .and. irrepp == irrepq) then
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(q, p, base_in, b11, b21, b12, b22)
                do s = 1, nums
                    do r = 1, numr
                        !$acc loop independent
                        do q = 1, numq
                            do p = 1, q
                                if (IOPPQ == 0 .and. p == q) cycle
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c
                                if (IOPPQ == 0) then; base_in = base_in + pq_off_c - 1 + (q-1)*(q-2)/2 + p
                                else;                 base_in = base_in + pq_off_c - 1 + (q-1)*q/2 + p; endif
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(2) deviceptr(iin, iout) private(q, p, base_in, b11, b21, b12, b22)
                do s = 1, nums
                    do r = 1, numr
                        !$acc loop independent
                        do q = 1, numq
                            do p = 1, nump
                                base_in = (rs_off_c - 1 + (s-1)*numr + r-1) * numpq_c + pq_off_c - 1 + (q-1)*nump + p
                                b11 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b21 = (rs_n1 - 1 + (s-1)*numr + r-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                b12 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n1 - 1 + (q-1)*nump + p
                                b22 = (rs_n2 - 1 + (r-1)*nums + s-1) * numpq_f + pq_n2 - 1 + (p-1)*numq + q
                                iout(b11) = iin(base_in); iout(b21) = iin(base_in) * spq
                                iout(b12) = iin(base_in) * srs; iout(b22) = iin(base_in) * spq * srs
                            end do
                        end do
                    end do
                end do
            end if
        end do
    end do
end subroutine exp_pqrs_dg
