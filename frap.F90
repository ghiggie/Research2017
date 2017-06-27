module frap

    use global_env
    use walker
    use calendar
    use cell

    implicit none

    integer, dimension(:), allocatable :: frap_active

contains

    subroutine frap_bleach()

        integer :: i

        allocate(frap_active(n_walkers))

        do i = 1, n_walkers
            frap_active(i) = 1
        end do

        do i = 1, n_walkers
            if (walker_pos(i,Z) .eq. box_size(Z)) then
                frap_active(i) = 0
            end if
        end do

    end subroutine

    subroutine frap_measure(ratio)

        real(DP), intent(out) :: ratio

        integer :: num = 0, denom = 0
        integer :: i

        do i = 1, n_walkers
            if (walker_pos(i,Z) .eq. box_size(Z)) then
                denom = denom + 1
                if (frap_active(i) .eq. 1) then
                    num = num + 1
                end if
            end if
        end do

        ratio = num/denom

    end subroutine

end module
