module pool

    implicit none

    integer, dimension(:), allocatable :: walker_id_pool
    integer :: n_walkers, max_walkers
    integer :: walker_id_root, i

contains

    subroutine walker_id_pool_init()

        allocate(walker_id_pool(max_walkers))

        walker_id_pool(1 : n_walkers) = 0
        do i = n_walkers + 1, max_walkers - 1
            walker_id_pool(i) = i + 1
        end do
        walker_id_pool(max_walkers) = -1

        walker_id_root = n_walkers + 1

    end subroutine

    subroutine walker_id_pool_push(id)

        ! Put an id back in the pool. This is the equivalent of
        ! removing a particle.

        integer, intent(in) :: id

        n_walkers = n_walkers - 1
        walker_id_pool(id) = walker_id_root
        walker_id_root = id

    end subroutine

    subroutine walker_id_pool_pull(id)

        ! Pull an id from the pool. This is the equivalent of adding
        ! a particle.

        integer, intent(out) :: id

        n_walkers = n_walkers + 1
        id = walker_id_root
        walker_id_root = walker_id_pool(id)
        walker_id_pool(id) = 0

    end subroutine

end module pool
