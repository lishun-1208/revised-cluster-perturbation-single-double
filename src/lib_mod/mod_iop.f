! variable must be assignment if use libsrint:
!    prec, isd, igpu
      module mod_iop 
! 
      ! totgpu only use in readinp
      integer prec, isd, igpu, totgpu !for the precision and gpu
!            prec=8, isd=2, igpu=0, totgpu=0 for double precision calculations
      ! the directory for tmp file, default be './tmp/', must be complete routine 
      character(len=200) :: tmpdir = './tmp/'

      ! follow variable only use in readinp
      logical ldebug ! for debug
      integer cc_int    ! 1 for cfour; 2 for pyscf. default 2 
      integer maxiter   ! for control the max iteration step 
      integer cc_level  ! CC计算级别 0 CCSD 1 CCSD(T) 2 CC2
      integer save_tmp  ! 中间项控制参数
      real*8  ccsd_conv ! DIIS收敛标准
      real*8  energy_conv ! 能量收敛标准
      integer iop_nbo    !0 非nbo， 1 nbo基态， 2 nbo激发态 
      integer iop_soc    !0 旋轨耦合关闭 1 旋轨耦合开启 
      real*8  nbo_value  !本征值小于特定值的自然轨道被冻结
      integer eri_level  !0 all eri generated, 1 <abcd> not generated 2 <abcd> and <abci> not generated 
      integer eri_cfour  !0                    1                      2 
      logical labic
      integer iop_ref
      real*8  ee_conv    !wzf 
      integer estate_sym(8) !lyz
      
      ! follow variable not use in libsrint
      integer ioiabcd, ioiabci, iowmbej ! for decide whether to save integrals in memory or disk
      integer igpuabci ! for decide whether to calculate all the integrals in GPU except abcd
      integer igpuabcd ! for decide whether to save iabcd in gpu
      integer igpuall
      integer lcc_level  !gmg
      integer ioiabcdns, ioiabcins, iowmbejns, iowefamns !wzf
      integer tot_incore, totgpu_incore ! wzf  default 0
      character*200 tmpfile_path
      end module
        
