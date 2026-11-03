C
C JJ07
C
      subroutine walltimeso(rt)
      real*8 rt 
      integer*4 istat, initime(8), timer(8)
      character*12 tstr(3)
      integer*4 rank
      save initime

C get current realtime and calculate current walltime
           call date_and_time(tstr(1), tstr(2), tstr(3), timer)
           call dat2sec(timer, rt)

      end subroutine
