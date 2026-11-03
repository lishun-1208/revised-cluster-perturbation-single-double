      real*8 function dotso(nsize, iin, ioffin, iout, ioffout)
      use cublas, only:cublassdot, cublasddot 
      use mod_iop
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
      real*8, external :: ddot, ddotsog
      real*4, external :: sdot
C
      if(prec==8.and.igpu==0) then 
         dotso = ddot(nsize, iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==0) then 
         dotso = sdot(nsize, iin, ioffin, iout, ioffout)
      elseif(prec==8.and.igpu==1) then 
         ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
         dotso = ddotsog(nsize4, iin, ioffin4, iout, ioffout4)  
      elseif(prec==4.and.igpu==1) then 
         ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
         dotso = cublassdot(nsize4, iin, ioffin4, iout, ioffout4)
      endif
      return
      end
       
      real*8 function ddotsog(nsize, iin, ioffin, iout, ioffout)
      use cublas, only:cublassdot, cublasddot 
      implicit none 
      integer*4, intent(in) :: nsize, ioffin, ioffout 
      real*8, intent(in),managed :: iin(*), iout(*) 
      real*8, external :: ddot
      real*4, external :: sdot
C
      ddotsog = cublasddot(nsize, iin, ioffin, iout, ioffout) 
      return
      end 
