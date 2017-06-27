module crowding

    use global_env

    implicit none

contains

    function distance(v1, v2)

        integer, dimension(3) :: v1, v2
        integer :: i
        real(DP) :: distance = 0, summ = 0

        do i = 1, 3
            summ = summ + (v1(i) - v2(i))**2
        end do

        distance = sqrt(summ)

    end function







end module
