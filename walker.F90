module walker

   use global_env
   use calendar
   use cell

   integer, dimension(:,:), allocatable :: walker_pos
   integer, dimension(:),   allocatable :: walker_col
   
   integer :: n_walkers, max_walkers
   
 contains

   subroutine walker_init(max_walkers)

      implicit none
      
      integer, intent(in) :: max_walkers

      allocate(walker_pos(max_walkers,3))
      allocate(walker_col(max_walkers))
      
   end subroutine walker_init

   subroutine walker_predict_event(walker_id)
   !******************************************************
   ! Definition of action 
   !    1 = migration (-x)
   !    2 = migration (+x)
   !    3 = migration (-y)
   !    4 = migration (+y)
   !    5 = migration (-z)
   !    6 = migration (+z)
   !    7 = membrane association
   !    8 = membrane desociation
   !
   !   walker_action takes care of boundary conditions
   !******************************************************
   
      implicit none
      
      integer, intent(in) :: walker_id

      integer :: action
      real(DP) :: tau
      real(DP) :: g

      if (walker_pos(walker_id,Z)==box_size(Z))  then
         ! in front of the membrane
         call random_number(g)
         tau = -log(g)/(5.0_DP*k_3d+k_on)
         call random_number(g)
         if(k_on/(5.0_DP*k_3d+k_on)>g) then
            action = 7
         else
            call random_number(g)
            action = ceiling(g*5.0_DP)
         end if
      else if (walker_pos(walker_id,Z)>box_size(Z))  then
         ! on the membrane
         call random_number(g)
         tau = -log(g)/(4.0_DP*k_2d+k_off)
         call random_number(g)
         if(k_off/(4.0_DP*k_2d+k_off)>g) then
            action = 8
         else
            call random_number(g)
            action = ceiling(g*4.0_DP)
         end if  
      else
         ! in the cytosol
         call random_number(g)
         tau = -log(g)/(6.0_DP*k_3d)
         call random_number(g)
         action = ceiling(g*6.0_DP)
      end if

      call calendar_schedule_event(action,walker_id,0,time_current+tau)

   end subroutine walker_predict_event
   
   subroutine walker_action(walker_id,action)
   
      implicit none
      
      integer, intent(in) :: walker_id
      integer, intent(in) :: action
      integer, dimension(3) :: p ! position
      integer :: n
      logical :: occupied
      
      p = walker_pos(walker_id,1:3)
      
      select case (action)

      case(1)
         p(X) = p(X)-1
         if(p(X)<1) p(X)=box_size(X)  ! periodic

      case(2)
         p(X) = p(X)+1
         if(p(X)>box_size(X)) p(X)=1  ! periodic

      case(3)
        p(Y) = p(Y)-1
        if(p(Y)<1) p(Y)=box_size(Y)   ! periodic

      case(4)
         p(Y) = p(Y)+1
         if(p(Y)>box_size(Y)) p(Y)=1  ! periodic
         
      case(5)
         p(Z) = p(Z)-1
         if(p(Z)<1) p(Z)=1         ! reflection

      case(6)
         p(Z) = p(Z)+1
         if(debug .and. p(Z)>box_size(Z)) then
            write(fdbg,*) '**walker_action**: Unexpected action.'
            stop
         end if

      case(7)
         p(Z) = p(Z)+1
         if(debug .and. p(Z)/=box_size(Z)+1) then
            write(fdbg,*) '**walker_action**: Unexpected position.'
            stop
         end if
   
      case(8)
         p(Z) = p(Z)-1
         if(debug .and. p(Z)/=box_size(Z)) then
            write(fdbg,*) '**walker_action**: Dissociation failed.'
            stop
         end if
         
      end select
      
      call cell_find_walkers(cell_coordinates(p),cell_local_walkers)

      occupied=.false.
      n=1
      do while(cell_local_walkers(n)>0)
         if(cell_local_walkers(n) /= walker_id .and. all(p==walker_pos(cell_local_walkers(n),:))) then
            occupied=.true.
            exit
         end if
         n=n+1
      end do

      if(.not.occupied) then
         walker_pos(walker_id,:) = p
         call cell_crossing(walker_id,p)
      else
         n_collisions=n_collisions+1
      end if
      
   end subroutine walker_action
        
end module walker
      
