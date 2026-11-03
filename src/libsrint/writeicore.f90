  subroutine writeicore_nzero(iin, lth, outname) 
  implicit none 
  real*8, intent(in) :: iin(*) 
  integer, intent(in) :: lth
  character(*) :: outname
  integer i, iter 
   write(6,*) outname, lth
  iter = 0 
  do i = 1, lth
     if(dabs(iin(i))>1.d-10) then  
        write(6,*)'output', outname, i, iin(i)
        iter = iter + 1
        if(iter>20) exit 
     endif
  enddo

  return
  end
