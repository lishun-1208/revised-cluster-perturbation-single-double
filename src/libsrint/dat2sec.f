      subroutine dat2sec(tim2, sec)
C     lets start with the most important Fortran command
      implicit none
C
C     MEH October 2010, Austin: should be safe for the next 1000 years or more  
C
          integer*4, dimension(8) :: tim2
          integer, dimension(8) :: tim
          integer :: i
          real*8    :: sec,numdays
          logical   :: leapyear

          tim(1)=tim2(1)
          tim(2)=tim2(2)
          tim(3)=tim2(3)
          tim(4)=tim2(4)
          tim(5)=tim2(5)
          tim(6)=tim2(6)
          tim(7)=tim2(7)
          tim(8)=tim2(8)

C  tim(1)      year
C  tim(2)      month
C  tim(3)      day of the month
C  tim(4)      time difference with utc in minutes
C  tim(5)      hour of the day
C  tim(6)      minutes of the hour
C  tim(7)      seconds of the minute
C  tim(8)      milliseconds of the second

C  calculate seconds elapsed since 01.01.2001 00:00

          sec=0.d0

C take care of years, get rid of the first 2000 years

          tim(1)=tim(1)-2000

C loop over all completed years since 01.01.2001 00:00

          do i=1,tim(1)-1

          leapyear=.false.
          if (mod(i,4).eq.0) leapyear=.true.
          if (mod(i,100).eq.0) leapyear=.false.
          if (mod(i,400).eq.0) leapyear=.true.
           
          if (leapyear) then
          sec=sec+i*60.d0*60.d0*24.d0*366.d0
          else
          sec=sec+i*60.d0*60.d0*24.d0*365.d0
          end if

          end do
 
C take care of months

C is the current year a leapyear ?

          leapyear=.false.
          if (mod(tim(1),4).eq.0) then 
          if ((mod(tim(1),100).ne.0).or.(mod(tim(1),400).eq.0)) then
          leapyear=.true.
          end if
          end if

C loop over all passed months in the current year

          do i=1,tim(2)-1

C get the number of days for the different months

          select case(i)
          case (2)
          numdays=28.d0
          if (leapyear) numdays=29.d0
          case (4)
          numdays=30.d0
          case (6)
          numdays=30.d0
          case (9)
          numdays=30.d0
          case (11)
          numdays=30.d0
          case default
          numdays=31.d0
          end select         

C add the days from the passed months

          sec=sec+numdays*60.d0*60.d0*24.d0
 
          end do

C now deal with the current month (its days, hours, and seconds)

          sec=sec
C add full hours
     &        +(tim(3)-1)*60.d0*60.d0*24.d0
C add the hours
     &        +tim(5)*60.d0*60.d0
C add the minutes
     &        +tim(6)*60.d0
C add the seconds
     &        +tim(7)
C add the milliseconds
     &        +tim(8)*0.001d0

      end subroutine
