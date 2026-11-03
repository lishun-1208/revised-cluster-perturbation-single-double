      subroutine swapso(nsize, iin, ioffin, iout, ioffout) 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in) :: iin(*), iout(*) 
C
      if(prec==8.and.igpu==0) then 
         call dswap(nsize,  iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==0) then 
         call sswap(nsize,  iin, ioffin, iout, ioffout)
      elseif(prec==8.and.igpu==1) then 
         call dswapg(nsize,  iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==1) then 
         call sswapg(nsize,  iin, ioffin, iout, ioffout)
      endif
      return
      end
       
      subroutine dswapg(nsize, iin, ioffin, iout, ioffout)
      use cublas 
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*8, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublasdswap(nsize4, iin, ioffin4, iout, ioffout4)  
      return
      end
       
      subroutine sswapg(nsize, iin, ioffin, iout, ioffout)
      use cublas 
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublassswap(nsize4, iin, ioffin4, iout, ioffout4)  
      return
      end
       
