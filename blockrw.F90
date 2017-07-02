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
    open(fdat2, file = 'diffusion.dat'

    if (debug) then
        write(fmain, *) 'output: debug.dat'
        open(fdbg, file = 'debug.dat')
    end if

    ! Compute the number of walkers needed to define the blockage

    block_walkers = product(sides) - product(sides - 2)

    max_walkers = n_walkers + block_walkers
    max_events = 2 * max_walkers

    call calender_init(max_walkers, max_events)

    call walker_init(max_walkers)

    ! Allocate the positions for the moving walkers
    allocate(rn(n_walkers))
    call random_number(rn)
    walker_pos(1 : n_walkers, X) = ceiling(rn * box_size(X))
    call random_number(rn)
    walker_pos(1 : n_walkers, Y) = ceiling(rn * box_size(Y))
    call random_number(rn)
    walker_pos(1 : n_walkers, Z) = ceiling(rn * box_size(Z))

    ! Allocate the positions for the blockage










