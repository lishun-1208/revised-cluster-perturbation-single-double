      real*8 function dotso2(nsize, iin, ioffin, iout, ioffout)
      use cublas, only:cublassdot, cublasddot
      use mod_iop 
      implicit none 
      integer, intent(in) :: nsize, ioffin, ioffout 
      real*4, intent(in),managed :: iin(*), iout(*) 
      integer*4 ioffin4, ioffout4, nsize4 
      real*8, external :: ddot, ddotsog, sdotsog2, sdotso2
      real*4, external :: sdot
C
      if(prec==8.and.igpu==0) then 
         dotso2 = ddot(nsize, iin, ioffin, iout, ioffout)  
      elseif(prec==4.and.igpu==0) then 
         dotso2 = sdotso2(nsize, iin, ioffin, iout, ioffout)
      elseif(prec==8.and.igpu==1) then 
         ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
         dotso2 = ddotsog(nsize4, iin, ioffin4, iout, ioffout4)  
      elseif(prec==4.and.igpu==1) then 
         ioffin4=ioffin; ioffout4=ioffout; nsize4=nsize
         dotso2 = sdotsog2(nsize4, iin, ioffin4, iout, ioffout4)
      endif
      return
      end
C       
      real*8 function sdotsog2(nsize, iin, ioffin, iout, ioffout)
      implicit none
      integer*4, intent(in) :: nsize, ioffin, ioffout
      real*4, intent(in),device :: iin(*), iout(*)
      integer i
C
      sdotsog2=0.d0
!$acc parallel loop reduction(+:sdotsog2)
      do i = 1, nsize
         sdotsog2 = sdotsog2 + iin(i)*iout(i)
      enddo
!$acc end parallel loop
      return
      end
C
      real*8 function sdotso2(nsize, iin, ioffin, iout, ioffout)
      implicit none
      integer*4, intent(in) :: nsize, ioffin, ioffout
      real*4, intent(in) :: iin(*), iout(*)
      integer i
C
      sdotso2=0.d0
!$omp parallel do reduction(+:sdotso2)
      do i = 1, nsize
         sdotso2 = sdotso2 + iin(i)*iout(i)
      enddo
!$omp end parallel do 
      return
      end
 
