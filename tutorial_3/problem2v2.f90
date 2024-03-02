subroutine swap(a,b)
    implicit none
    integer, intent(inout) :: a, b
    
    a = a + b
    b = a - b
    a = a - b
end subroutine swap

program main
    use omp_lib

    implicit none
    double precision :: array(100)
    integer :: a(100), thread_count = 1, pass, i, n=100
    character(100) :: numchar, name

    if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
       write(*, *) 'Error, one command line arguments is required, stopping the program'
       STOP
    end if
  
    call GET_COMMAND_ARGUMENT(0, name)
    write(*, *) name
    
    call GET_COMMAND_ARGUMENT(1, numchar)
    read(numchar, *) thread_count

    call random_number(array)
    a = nint(array*100)
    ! a = [1,5,3,7,0,2,4,9,6,8]

    write(*,'(100i6)') a

    !$omp parallel num_threads(thread_count) &
    !$omp default(none) private(i) shared(a,n,pass)
    do pass = 1, n           ! outer loop
        if (mod(pass,2) .eq. 0) then
            !$omp do
            do i = 2, n, 2            ! inner loop
                if (a(i-1) .gt. a(i)) call swap(a(i-1), a(i))
            end do
            !$omp end do
        else
            !$omp do
            do i = 2, n-1, 2          ! inner loop
                if (a(i) .gt. a(i+1)) call swap(a(i), a(i+1))
            end do
            !$omp end do
        end if
    end do
    !$omp end parallel
    write(*,*) 'the sorted array is:'
    write(*,'(100i6)') a
end program main