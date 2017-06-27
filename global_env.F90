module global_env

   use iso_fortran_env
   
   integer, parameter :: DP=REAL64

   integer, parameter :: X=1, Y=2, Z=3
   integer, parameter :: fmain=11, fdbg=15, fdat1=21, fdat2=22, fdat3=23

   logical debug

   integer :: n_collisions

!  system parameters
   integer, dimension(3) :: box_size
   real(DP) :: k_on, k_off, k_3d, k_2d

contains

   function time_stamp()
   
      implicit none
      
      character(18) :: time_stamp
      character(8)  :: date
      character(10) :: time

      call date_and_time(date,time)
      
      time_stamp = date(5:6) // '/' // date(7:8) // '/' //  date(1:4) &
                             // '  ' // time(1:2) // ':' // time(3:4)

   end function time_stamp

end module global_env
