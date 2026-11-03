      subroutine scalso(nsize, fact, iin, ioffin) 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin
      real*4, intent(in) :: iin(*)
      real*8, intent(in) :: fact
      real*4 fact4 
C
      if(prec==8.and.igpu==0) then 
         call dscal(nsize, fact, iin, ioffin)  
      elseif(prec==4.and.igpu==0) then 
         fact4 = fact
         call sscal(nsize, fact4, iin, ioffin)
      elseif(prec==8.and.igpu==1) then 
         call dscalg(nsize, fact, iin, ioffin)  
      elseif(prec==4.and.igpu==1) then 
         call sscalg(nsize, fact, iin, ioffin)
      endif
      return
      end
       
      subroutine dscalg(nsize, fact, iin, ioffin)
      use cublas
      implicit none 
      integer, intent(in) :: nsize, ioffin
      real*8, intent(in),managed :: iin(*)
      real*8, intent(in) :: fact
      integer*4 ioffin4, ioffout4, nsize4 
C
      ioffin4=ioffin; nsize4=nsize
      call cublasdscal(nsize4, fact, iin, ioffin4)  
      return
      end
       
      subroutine sscalg(nsize, fact, iin, ioffin)
      use cublas
      implicit none 
      integer, intent(in) :: nsize, ioffin
      real*4, intent(in),managed :: iin(*)
      real*8, intent(in) :: fact
      real*4 fact4
      integer*4 ioffin4, ioffout4, nsize4 
C
      fact4=fact;ioffin4=ioffin; nsize4=nsize
      call cublassscal(nsize4, fact4, iin, ioffin4)  
      return
      end
       
