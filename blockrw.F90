program blockrw

    use global_env
    use calendar
    use walker
    use cell

    implicit none

    integer :: max_events
    integer :: action, opt1, opt2

    real(DP) :: time, time_limit, time_stat
    real(DP) :: sigma, rho, xx

    integer, dimension(:,:), allocatable :: p0
    real(DP), dimension(:), allocatable :: rn
    integer, dimension(3) :: vertex, sides
    integer :: block_walkers

    integer :: i, j, k, l, m, n, part
    character(len = 10) :: arg
    debug = .false.

    if(iargc() > 0) then
        do i = 1, iargc()
            call getarg(i, arg)
            select case (trim(arg))
            case ('-debug')
                debug = .true.
            case ('-newseed')
                call newseed
            case default
                open(10, file = 'blockrw_err.log')
                write(10, *) "Unknown argument to program. Check your terminal command."
                close(10)
                stop
            end select
        end do
    end if

    open(10, file = 'block.dat')
    read(10, *) vertex
    read(10, *) sides
    close(10)

    open(10, file = 'config.dat')
    read(10, *) n_walkers
    read(10, *) box_size
    read(10, *) cell_size
    read(10, *) k_3d, k_2d
    read(10, *) k_on, k_off
    read(10, *) time_limit
    read(10, *) time_stat
    close(10)

    open(fmain, file = 'rwblock.log')
    write(fmain, *) '*** KMC FRAP Congested ***'
    write(fmain, *) 'Started at ', time_stamp()
    write(fmain, '(/a,i5)') 'Number of Random Walkers = ', n_walkers
    write(fmain, '(3(a,i4))') 'Box size = ', box_size(X), 'x', box_size(Y), 'x', box_size(Z)
    write(fmain, '(3(a,i4))') 'Cell size = ', cell_size(X), 'x', cell_size(Y), 'x', cell_size(Z)
    write(fmain, '(3(a,i4))') 'Block vertex = ', vertex(X), ', ', vertex(Y), ', ', vertex(Z)
    write(fmain, '(3(a,i4))') 'Block size = ', sides(X), ', ', sides(Y), ', ', sides(Z)
    write(fmain, '(2(a,f8.4),a)') 'Transition rate = ', k_3d, ' (3d), ', k_2d, ' (2d)'
    write(fmain, '(2(a,f8.4),a)') 'Transition rate = ', k_on, ' (on), ', k_off, ' (off)'
    write(fmain, '(a,e15.3)') 'Total simulation time = ', time_limit
    write(fmain, *)

    write(fmain, *) 'Output: onoff.dat'
    open(fdat1, file = 'onoff.dat')
    write(fmain, *) 'Output: diffusion.dat'
    open(fdat2, file = 'diffusion.dat')

    if (debug) then
        write(fmain, *) 'output: debug.dat'
        open(fdbg, file = 'debug.dat')
    end if

    ! Compute the number of walkers needed to define the blockage

    block_walkers = product(sides) - product(sides - 2)

    max_walkers = n_walkers + block_walkers
    max_events = 2 * max_walkers

    call calendar_init(max_walkers, max_events)

    call walker_init(max_walkers)

    ! Allocate the positions for the moving walkers
    allocate(rn(n_walkers))
    call random_number(rn)
    walker_pos(1 : n_walkers, X) = ceiling(rn * box_size(X))
    call random_number(rn)
    walker_pos(1 : n_walkers, Y) = ceiling(rn * box_size(Y))
    call random_number(rn)
    walker_pos(1 : n_walkers, Z) = ceiling(rn * box_size(Z))
    deallocate(rn)
    ! Allocate the positions for the blockage

    part = n_walkers + 1
    do i = 1, sides(X)
        do j = 1, sides(Z)
            walker_pos(part, X) = vertex(X) + i - 1
            walker_pos(part, Y) = vertex(Y)
            walker_pos(part, Z) = vertex(Z) + j - 1
            walker_pos(part + 1, X) = vertex(X) + i - 1
            walker_pos(part + 1, Y) = vertex(Y) + sides(Y) - 1
            walker_pos(part + 1, Z) = vertex(Z) + j - 1
            part = part + 2
        end do
    end do
    do k = 1, sides(Y) - 2
        do i = 1, sides(X)
            walker_pos(part, X) = vertex(X) + i - 1
            walker_pos(part, Y) = vertex(Y) + k
            walker_pos(part, Z) = vertex(Z)
            walker_pos(part + 1, X) = vertex(X) + i - 1
            walker_pos(part + 1, Y) = vertex(Y) + k
            walker_pos(part + 1, Z) = vertex(Z) + sides(Z) - 1
            part = part + 1
        end do
        do j = 2, sides(Z) - 1
            walker_pos(part, X) = vertex(X)
            walker_pos(part, Y) = vertex(Y) + k
            walker_pos(part, Z) = vertex(Z) + j - 1
            walker_pos(part + 1, X) = vertex(X) + sides(X) - 1
            walker_pos(part + 1, Y) = vertex(Y) + k
            walker_pos(part + 1, Z) = vertex(Z) + j - 1
            part = part + 2
        end do
    end do

    call cell_init(max_walkers, max_walkers, walker_pos)

    if (debug) call cell_test(max_walkers)

    do i = 1, n_walkers
        call walker_predict_event(i)
    end do

    do i = n_walkers + 1, max_walkers
        ! Using time_limit + i means that these events will never
        ! actually occur. This is the most important part of my code.
        call calendar_schedule_event(1, i, 0, time_limit + i)
    end do

    allocate(p0(n_walkers, 3))
    p0 = walker_pos(n_walkers, 3)

    call calendar_schedule_event(11, 0, 0, time_stat)
    call calendar_schedule_event(12, 0, 0, time_stat)

    n_collisions = 0
    time_current = 0.0_DP

    do while (time_current .le. time_limit)
        call calendar_find_event(action, opt1, opt2)
        if (action .le. 10) then
            call walker_action(opt1, action)
            call walker_predict_event(opt1)
        else if (action .le. 20) then
            select case (action)
                case(11)
                    sigma = count(walker_pos(n_walkers,Z) .gt. box_size(Z))
                    rho = count(walker_pos(n_walkers,Z) .eq. box_size(Z))
                    write(fdat1,'(f10.2,2e13.5)') time_current, rho, sigma
                    call calendar_schedule_event(11,0,0,time_current+time_stat)
                case(12)
                    xx = real(sum((walker_pos(n_walkers,3)-p0)**2),kind=DP)/real(n_walkers,kind=DP)
                    write(fdat2,'(f10.2,e13.5)') time_current, xx
                    call calendar_schedule_event(12,0,0,time_current+time_stat)
                case default
                    write(fmain,'(f10.2,a)') time_current, 'unknown event'
                    stop
            end select
        else 
            write(fmain,*) 'Unknown event'
            stop
        end if
    end do

    write(fmain,'(a,i10)') 'Number of collisions = ', n_collisions

    if (debug) call cell_test(max_walkers)

    close(fdat1)
    close(fdat2)

    open(fdat1, file = 'walkers.dat')
    write(fdat1,'(3i5)') (walker_pos(i, :), i = 1, n_walkers)
    close(fdat1)

    if (debug) close(fdbg)

    write(fmain, *) '*** Done at ', time_stamp(), ' ***'
    close(fmain)

end program blockrw
