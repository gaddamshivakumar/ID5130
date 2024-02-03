module num_integration
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    double precision, parameter :: pi = 3.14159265358
    
contains

double precision function func(x)
    implicit none
    double precision, intent(in) :: x

    func = sin(x)/(2*x**3)
end function func

subroutine simpsons(n, a, b, thread_count, final_result)
    implicit none
    integer, intent(in) :: n, thread_count
    double precision, intent(in) :: a, b
    double precision, intent(out) :: final_result
    integer :: i
    double precision :: x, y, h, evensum = 0.0, oddsum = 0.0

    h = (b-a)/n
    final_result = func(a) + func(b)
    !$OMP PARALLEL DO NUM_THREADS(thread_count) &
    !$OMP PRIVATE(i,x,y) REDUCTION(+:evensum,oddsum)
    do i = 1, n-1
        ! my_rank = OMP_GET_THREAD_NUM()
        ! write(*, *) 'thread number is: ',my_rank
        x = a + h*i
        y = func(x)
        if (mod(i,2) .eq. 0) then
            evensum = evensum + y
        else
            oddsum = oddsum + y
        endif
    enddo
    !$OMP END PARALLEL DO
    final_result = (h/3)*(final_result + 4*oddsum + 2*evensum)
end subroutine simpsons
    
end module num_integration


!--------------------------- MAIN PROGRAM ------------------------------
program main
    use num_integration
    implicit none
    double precision :: a, b, final_result, exact_result
    integer :: n, thread_count = 1
    character(100) :: numchar
    character(100) :: name

    if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
       write(*, *) 'Error, one command line arguments is required, stopping the program'
       stop
    end if
  
    call GET_COMMAND_ARGUMENT(0, name)
    write(*, *) name
    
    call GET_COMMAND_ARGUMENT(1, numchar)
    read(numchar, *) thread_count
  
    n = 32
    a = 1.0
    b = PI
    final_result = 0.0
    exact_result = 0.198557

    call simpsons(n, a, b, thread_count, final_result)

    write(*, *) 'The final result is ', final_result
    write(*, *) 'The error is ', abs(exact_result-final_result)
end program main