program pooltest

    use pool

    implicit none

    character(len = *), parameter :: fmt1 = "(A15,10(I3))"
    character(len = 6) :: option
    integer :: id
    n_walkers = 5
    max_walkers = 10

    call walker_id_pool_init()
    write(*,*) "walker_id_root:", walker_id_root
    write(*,*) "n_walkers:", n_walkers
    write(*,fmt1) "walker_id_pool:", walker_id_pool

    ! Get desired
10  write(*,*) "Would you like to add or remove a particle? Please &
                specify either 'add' or 'remove'."
    read(*,*) option

    select case (option)

        case('add')
            ! The assigned id should be 6

            ! Pull the next available id
            call walker_id_pool_pull(id)

            if (n_walkers .gt. max_walkers) then
                write(*,*) "Error: There is no more space for any more&
                            walkers."
                stop
            end if

            write(*,*) "Assigned ID:", id
            write(*,*) "walker_id_root:", walker_id_root
            write(*,*) "n_walkers:", n_walkers
            write(*,fmt1) "walker_id_pool:", walker_id_pool

            ! Check to see if the use wants to keep going
20          write(*,*) "Would you like to keep going? Please &
                        state 'yes' or 'no'."
            read(*,*) option

            select case (option)
                case('yes')
                    go to 10
                case('no')
                    stop
                case default
                    write(*,*) "Please choose either 'yes' &
                                or 'no'."
                    go to 20
            end select

        case('remove')

            if (n_walkers .le. 0) then
                write(*,*) "Error: There are no more particles &
                            in they system."
                stop
            end if

30          write(*,*) "Which particle would you like to remove?"
            read(*,*) id

            if (id .le. 0 .or. id .gt. max_walkers) then
                write(*,*) "Please select an appropriate ID."
                go to 30
            end if

            if (walker_id_pool(id) .gt. 0) then
                write(*,*) "That ID is unoccupied, so there is &
                            nothing to do."
                go to 40
            end if

            call walker_id_pool_push(id)
            write(*,*) "walker_id_root:", walker_id_root
            write(*,*) "n_walkers:", n_walkers
            write(*,fmt1) "walker_id_pool:", walker_id_pool

            ! Check to see if the use wants to keep going
40          write(*,*) "Would you like to keep going? Please &
                        state 'yes' or 'no'."
            read(*,*) option

            select case (option)
                case('yes')
                    go to 10
                case('no')
                    stop
                case default
                    write(*,*) "Please choose either 'yes' &
                                or 'no'."
                    go to 40
            end select

        case default
            write(*,*) "Please specify either 'add' or 'remove'."
            go to 10
    end select

end program pooltest
