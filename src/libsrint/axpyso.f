      subroutine axpyso(nsize, fact, iin, ioffin, iout, ioffout)
      use cublas, only:cublassaxpy, cublasdaxpy 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in),managed :: iin(*), iout(*) 
      real*8, intent(in) :: fact
      real*4 fact4
      integer*4 ioffin4, ioffout4, nsize4 
C
      if(prec==8.and.igpu==0) then 
         call daxpy(nsize, fact, iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==0) then 
         fact4 = fact
         call saxpy(nsize, fact4, iin, ioffin, iout, ioffout)
      elseif(prec==8.and.igpu==1) then 
         call daxpyg(nsize, fact, iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==1) then 
         fact4 = fact
         ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
         call cublassaxpy(nsize4, fact4, iin, ioffin4, iout, ioffout4)
      endif
      return
      end
       
      subroutine daxpyg(nsize, fact, iin, ioffin, iout, ioffout)
      use cublas, only:cublassaxpy, cublasdaxpy 
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*8, intent(in), device :: iin(*), iout(*) 
      real*8, intent(in) :: fact
      real*4 fact4
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublasdaxpy(nsize4, fact, iin, ioffin4, iout, ioffout4)  
      return
      end
       
      subroutine saxpyg(nsize, fact, iin, ioffin, iout, ioffout)
      use cublas, only:cublassaxpy, cublasdaxpy 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in),managed :: iin(*), iout(*) 
      real*8, intent(in) :: fact
      real*4 fact4
      integer*4 ioffin4, ioffout4, nsize4 
C
      fact4=fact;ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
      call cublassaxpy(nsize4, fact4, iin, ioffin4, iout, ioffout4)  
      return
      end
       
