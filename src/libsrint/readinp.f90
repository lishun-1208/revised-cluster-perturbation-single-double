      subroutine readinp 
      use mod_ri
      use mod_iop
      use mod_orbit
      use mod_symm
      use mod_size
      implicit none 
      integer i, j, k, l1, l2
      integer mgb, mggb, ntmp
      logical yesno
      character*100 iname(60) 
!    ri=on or ri=off
!    memory=60  
!    gpumemory=10
!    prec=4 or prec=8 for single and double precision
!    gpu=off
!    debug=true or debug=false
!    cc_int=1    cc_int=1 for get integrals from cfour, cc_int=2 for get integrals from pyscf
!    maxiter = 20
!    ccsd_conv = 1.d-8
!    energy_conv = 1.d-7
!    cc_level=0, 1  ! 0 for CCSD and 1 for CCSD(T)
!    save_tmp=0 not save, save_tmp=1 save
!    ee_conv        !
!    estate_sym(8) !
      cc_int=2
      totgpu=0                      ! 0 不使用GPU， 1使用GPU
      mgb=40                        ! 内存， GB
      mggb=10                       ! 显存， GB
      icrsiz = 1024*1024*1024/4*40  !80GB   memory size 
      igpusize= 1024*1024*1024/4*10 !10GB   gpu size
      l_ri=.false.                  ! 是否使用DF
      ldebug=.false.                ! 是否检验模式
      maxiter=20                    ! 最大CC迭代
      ccsd_conv=1.d-8               ! DIIS收敛标准   
      energy_conv=1.d-7             ! 能量收敛标准
      cc_level=1                    ! CC计算级别
      save_tmp=1                    ! 中间项是否存储
      iop_nbo = 0                   ! 是否使用NBO，目前仅能使用2激发态。       
      iop_soc = 0
      nbo_value = 1.d-4
      eri_level = 0
      eri_cfour = 0
      labic = .true. 
      iop_ref = 0                   ! ref_state=rhf,uhf         
      ee_conv = 1.d-5  
        
!      write(6,*)'icrsiz', icrsiz
      iname=' '
      inquire(file='CCSOINP', exist=yesno)
      if(.not.yesno) then
        write(6,*)'CCSOINP file not exist in this file folder' 
        write(6,*)'use the default set with'
        write(6,*)'memory=40 GB, gpu memory=10 GB, &
                   ri=off, gpu=false, debug=false, maxiter=20'
        write(6,*)'ccsd(t) would be calculated and  &
                   tmp_file would be saved'
         return
      endif

      open(unit=30, file='CCSOINP', form='formatted')
      rewind 30
!      read(30,*, end=5000) (iname(i),i=1,100)
      i = 1
      do 
        read(30,*, end=5000) iname(i)
        i = i + 1
      enddo
5000  close(30)
      do i=1,100
         if(iname(i)/=' ') then
!           write(6,*)'i',i,' ', iname(i)
         else
            exit
         endif
      enddo

      l_ri=.false.
      do j = 1, i-1
         do k = 1, 80
            if(iname(j)(k:k)=='=') then
               l1 = k
            endif
            if(iname(j)(k:k)==' ') then
               l2 = k
               exit
            endif
         enddo
!         write(6,*)'iname', iname(j)
!         write(6,*)'iname',' ', iname(j)(1:l1-1)
         if(iname(j)(1:l1-1)=='cc_level') then
            read(iname(j)(l1+1:l2-1), "(I4)") cc_level
         endif
         if(iname(j)(1:l1-1)=='save_tmp') then
            read(iname(j)(l1+1:l2-1), "(F)") save_tmp
         endif
         if(iname(j)(1:l1-1)=='energy_conv') then
            read(iname(j)(l1+1:l2-1), "(I4)") ntmp
            energy_conv = 1.d-1**ntmp
         endif
         if(iname(j)(1:l1-1)=='ccsd_conv') then
            read(iname(j)(l1+1:l2-1), "(I4)") ntmp
            ccsd_conv = 1.d-1**ntmp 
         endif
         if(iname(j)(1:l1-1)=='ee_conv') then
            read(iname(j)(l1+1:l2-1), "(I4)") ntmp 
            ee_conv = 1.d-1**ntmp 
         endif
         if(iname(j)(1:l1-1)=='maxiter') then
            read(iname(j)(l1+1:l2-1), "(I4)") maxiter
!            write(6,*)'maxiter have been readed', maxiter 
         endif
         if(iname(j)(1:l1-1)=='gpumemory') then
            read(iname(j)(l1+1:l2-1), "(I4)") mggb
         endif
         if(iname(j)(1:l1-1)=='memory') then
            read(iname(j)(l1+1:l2-1), "(I4)") mgb
!            write(6,*)'Memory in GB-------------', mgb
         endif
         if(iname(j)(1:l1-1)=='prec') then
            read(iname(j)(l1+1:l2-1), "(I4)") prec
!            write(6,*)'Prec of data-------------', prec 
            if(prec==4) then 
!               write(6,*)'  Single precision data -'
                isd = 1 
            elseif(prec==8) then 
!               write(6,*)'  Double precision data -' 
                isd = 2
            endif
         endif
         if(iname(j)(1:l1-1)=='pysoc') then
            if(iname(j)(l1+1:l2-1)=='on') then
               iop_soc=1
            endif
         endif
         if(iname(j)(1:l1-1)=='gpu') then
            if(iname(j)(l1+1:l2-1)=='on') then
               totgpu=1
!            write(6,*)'Employ gpu acceleration--', totgpu
            endif 
         endif
         if(iname(j)(1:l1-1)=='nbo_state') then
            if(iname(j)(l1+1:l2-1)=='none') then
               iop_nbo=0
            elseif(iname(j)(l1+1:l2-1)=='ground') then
               iop_nbo=1
            elseif(iname(j)(l1+1:l2-1)=='excited') then
               iop_nbo=2
!               eri_level=1
            endif
         endif
         if(iname(j)(1:l1-1)=='nbo_value') then
            read(iname(j)(l1+1:l2-1), "(F)") nbo_value
         endif

         if(iname(j)(1:l1-1)=='ref_state') then
            if(iname(j)(l1+1:l2-1)=='rhf') then
               iop_ref=0
            elseif(iname(j)(l1+1:l2-1)=='uhf') then
               iop_ref=1
            endif
         endif

         if(iname(j)(1:l1-1)=='ri') then
            if(iname(j)(l1+1:l2-1)=='on') then
               l_ri=.true.
!            write(6,*)'Employ ri integrals------'
            endif 
         endif
         if(iname(j)(1:l1-1)=='debug') then
            if(iname(j)(l1+1:l2-1)=='true') then
               ldebug=.true.
!            write(6,*)'debug -------------------'
            endif 
         endif
         if(iname(j)(1:l1-1)=='eri_level') then
            read(iname(j)(l1+1:l2-1), "(I4)") eri_level
         endif
         if(iname(j)(1:l1-1)=='eri_cfour') then
            read(iname(j)(l1+1:l2-1), "(I4)") eri_cfour
         endif
         if(iname(j)(1:l1-1)=='cc_int') then
            read(iname(j)(l1+1:l2-1), "(I4)") cc_int 
         endif
         if(iname(j)(1:l1-1)=='cc_int') then
            read(iname(j)(l1+1:l2-1), "(I4)") cc_int 
         endif
         if(iname(j)(1:l1-1)=='save_abic') then
            if(iname(j)(l1+1:l2-1)=='false') then
               labic = .false.
            endif
         endif

      enddo
        if(ldebug) then 
            if(cc_level==0) then 
               write(6,*)'CCSD energy -------------' 
            elseif(cc_level==1) then 
               write(6,*)'CCSD(T) energy ----------'
            else
               write(6,*)'CC level not supportted -'
              stop
            endif
            if(save_tmp==1) then 
               write(6,*)'tmp file be saved--------' 
            elseif(cc_level==0) then 
               write(6,*)'tmp file be deleted -----'
            endif
             write(6,*)'energy_conv--------------', energy_conv
             write(6,*)'ccsd_conv----------------', ccsd_conv
             write(6,*)'Maxiter in ccsd----------', maxiter
             write(6,*)'Gpumemory in GB----------', mggb
            if(cc_int==1) then 
             write(6,*)'Get data in cfour--------', cc_int
            elseif(cc_int==2) then 
             write(6,*)'Get data in pyscf--------', cc_int
            endif
         endif
      icrsiz = 1024*1024*1024/4*mgb  !80GB   memory size 
      igpusize = 1024*1024*1024/4*mggb
      return
      end
