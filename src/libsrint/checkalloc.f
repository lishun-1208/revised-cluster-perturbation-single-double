      subroutine checkalloc(istat, string)
      integer, intent(in) :: istat
      character*(*) :: string
      if (istat/=0) then
          write(6,*) 'allocating ', string, ' error!'
          stop
      endif
      end
