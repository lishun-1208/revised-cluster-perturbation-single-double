      double precision function nrm2so(nsize, iin, ioff) 
C
      use mod_iop
      implicit none 
C
      real*4, intent(in),managed :: iin(*) 
      integer, intent(in) :: nsize, ioff
      real*8, external :: dnrm2, snrm2g, dnrm2g
      real*4, external :: snrm2
C
      if(prec==4.and.igpu==0) then 
         nrm2so = snrm2(nsize, iin, ioff) 
      elseif(prec==8.and.igpu==0) then 
         nrm2so = dnrm2(nsize, iin, ioff) 
      elseif(prec==4.and.igpu==1) then 
         nrm2so = snrm2g(nsize, iin, ioff) 
      elseif(prec==8.and.igpu==1) then 
         nrm2so = dnrm2g(nsize, iin, ioff) 
      endif
      return 
      end 
      
      double precision function dnrm2g(nsize, iin, ioff) 
C
      use cublas
      implicit none 
C
      real*8, intent(in),managed :: iin(*) 
      integer, intent(in) :: nsize, ioff
C
      integer*4 nsize4, ioff4 
      nsize4 = nsize ; ioff4 = ioff
      dnrm2g = cublasdnrm2(nsize4, iin, ioff4) 
      return 
      end 
      
      double precision function snrm2g(nsize, iin, ioff) 
C
      use cublas
      implicit none 
C
      real*4, intent(in),managed :: iin(*) 
      integer, intent(in) :: nsize, ioff
C
      integer*4 nsize4, ioff4 
      nsize4 = nsize ; ioff4 = ioff
      snrm2g = cublassnrm2(nsize4, iin, ioff4) 
      return 
      end 
      
