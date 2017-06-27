subroutine newseed()

  implicit none

  integer :: i, j, n, nmax
  integer, dimension(8) :: val
  integer, dimension(64) :: f
  integer, dimension(:), allocatable :: seed

  call date_and_time(VALUES=val)

  call random_seed(SIZE=nmax)
  if(nmax>64) then
     write(11,*) 'Seed size is too large:', nmax
     stop
  end if
  allocate(seed(nmax))

  ! generate 64 integers based on current date and time
  call date_and_time(VALUES=val)
  f = (val(8)+val(7)*1000+val(6)*60000)

  call random_seed(GET=seed(1:nmax))

  seed(1:nmax) = seed(1:nmax)*f(1:nmax)

  call random_seed(PUT=seed(1:nmax))

end subroutine newseed

