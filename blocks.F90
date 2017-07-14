module blocks

    use global_env
    use walker

    implicit none

    integer, dimension(3) :: vertex, sides
    integer :: block_walkers

contains

    subroutine blocks_alloc(n_walkers)

        integer, intent(in) :: n_walkers

        open(10, file = 'block.dat')
        read(10, *) vertex
        read(10, *) sides
        close(10)

        block_walkers = product(sides) - product(sides - 2)

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
    end subroutine

    subroutine block_schedule(n_walkers, time_limit)

        integer, intent(in) :: n_walkers, time_limit
        integer :: i

        do i = n_walkers +1, n_walkers + block_walkers
            call calendar_schedule_event(1, i, 0, time_limit + 1)
        end do

    end subroutine
