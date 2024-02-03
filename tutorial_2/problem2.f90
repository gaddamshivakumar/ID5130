! numerical integration using trapezoidal rule -- parallel code
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

    subroutine trapezoidal_rule_critical(n, a, b, final_result)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: a, b
        double precision, intent(out) :: final_result
        double precision :: h, x
        integer :: i
#ifdef _OPENMP
        double precision :: total
        integer :: local_n
        double precision :: local_a, local_b
        integer :: my_rank, thread_count
#endif
    
        h = (b-a)/n
#ifdef _OPENMP
        my_rank = OMP_GET_THREAD_NUM()
        thread_count = OMP_GET_NUM_THREADS()
        local_n = n/thread_count;
        local_a = a + my_rank*local_n*h
        local_b = local_a + local_n*h
    
        total = (func(local_a) + func(local_b))/2.0
        do i = 1, local_n-1
           x = local_a + i*h
           total = total + func(x)
        end do
        total = total*h
        ! write(*, *) 'my rank is: ', my_rank
        
        !$OMP CRITICAL
        final_result = final_result + total
        !$OMP END CRITICAL
#else
        final_result = (func(a) + func(b))/2.0
        do i = 1, n-1
            x = a + i*h
            final_result = final_result + func(x)
        end do
        final_result = final_result*h
#endif
    end subroutine trapezoidal_rule_critical

    subroutine trapezoidal_rule_parfor(n, a, b, thread_count, final_result)
        implicit none
        integer, intent(in) :: n, thread_count
        double precision, intent(in) :: a, b
        double precision, intent(out) :: final_result
        double precision :: h, x, y
        integer :: i
    
        h = (b-a)/n
        final_result = (func(a) + func(b))/2.0
        !$OMP PARALLEL DO NUM_THREADS(thread_count) REDUCTION(+:final_result) PRIVATE(x,y,i)
        do i = 1, n-1
            x = a + i*h
            y = func(x)
            final_result = final_result + y
        end do
        !$OMP END PARALLEL DO
        final_result = final_result*h
    end subroutine trapezoidal_rule_parfor
end module num_integration


program main
    use num_integration
    implicit none
    double precision :: a, b, final_result, exact_result
    integer :: n
#ifdef _OPENMP
    character(100) :: numchar
    character(100) :: name
    integer :: thread_count = 1
#endif

#ifdef _OPENMP
    if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
       write(*, *) 'Error, one command line arguments is required, stopping the program'
       STOP
    end if
  
    call GET_COMMAND_ARGUMENT(0, name)
    write(*, *) name
    
    call GET_COMMAND_ARGUMENT(1, numchar)
    read(numchar, *) thread_count
#endif
  
    n = 32
    a = 1.0
    b = PI
    final_result = 0.0
    exact_result = 0.198557

#ifdef _OPENMP
    if ( .true. ) then
        !$OMP PARALLEL NUM_THREADS(thread_count)
        call trapezoidal_rule_critical(n, a, b, final_result)
        !$OMP END PARALLEL
    else
        call trapezoidal_rule_parfor(n, a, b, thread_count, final_result)
    endif
#else
    call trapezoidal_rule_critical(n, a, b, final_result)
#endif

    write(*, *) 'The final result is ', final_result
    write(*, *) 'The error is ', abs(exact_result-final_result)
end program main