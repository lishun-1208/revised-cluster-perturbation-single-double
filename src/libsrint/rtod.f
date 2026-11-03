      subroutine rtod(iin, iout, length) 
C--- This program is to translate double precision data to single
C--- precision data.
      implicit none
      integer, intent(in) :: length
      double precision, intent(inout) :: iout(*)
      real*4, intent(inout) :: iin(*)
C
      integer i 
C
!$omp parallel do 
      do i = 1, length 
         iout(i) = dble(iin(i))    
      enddo
!$omp end parallel do
C  
      return
      end 
C
