program main
    use omp_lib

    implicit none
    double precision :: array(16), sum = 0.0, PI = 4.0*atan(1.0)
    integer :: i, thread_count = 1
    character(100) :: numchar, name

    if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
       write(*, *) 'Error, one command line arguments is required, stopping the program'
       STOP
    end if
  
    call GET_COMMAND_ARGUMENT(0, name)
    write(*, *) name
    
    call GET_COMMAND_ARGUMENT(1, numchar)
    read(numchar, *) thread_count


    do i = 1, 16
        array(i) = sin(PI/i)
    end do
    
    !$omp parallel do num_threads(thread_count) reduction(-: sum)
    do i = 1, 16
        sum = sum - array(i)
    end do
    !$omp end parallel do

    write(*,*) 'the value of sum = ', sum
end program main