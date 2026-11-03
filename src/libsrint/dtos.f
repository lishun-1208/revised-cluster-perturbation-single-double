      subroutine dtos(iin, iout, length) 
C--- This program is to translate double precision data to single
C--- precision data.
      implicit none
      integer, intent(in) :: length
      double precision, intent(in) :: iin(*)
      real*4, intent(out) :: iout(*)
C
      integer i 
C
!$omp parallel do 
      do i = 1, length 
         iout(i) = real(iin(i))    
      enddo
!$omp end parallel do 
C  
      return
      end 
C
