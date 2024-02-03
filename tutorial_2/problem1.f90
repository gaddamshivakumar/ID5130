! numerical integration using trapezoidal rule -- serial code
module num_integration
    implicit none
    
    double precision, parameter :: pi = 3.14159265358
contains
    
    double precision function func(x)
        implicit none
        double precision, intent(in) :: x

        func = sin(x)/(2*x**3)
    end function func

    subroutine trapezoidal_rule(n, a, b, final_result)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: a, b
        double precision, intent(out) :: final_result
        double precision :: h, x
        integer :: i
    
        h = (b-a)/n
    
        final_result = (func(a) + func(b))/2.0
        do i = 1, n-1
           x = a + i*h
           final_result = final_result + func(x)
        end do
        final_result = final_result*h
    end subroutine trapezoidal_rule

end module num_integration


program main
    use num_integration
    implicit none
    character(100) :: name
    double precision :: a, b, final_result, exact_result
    integer :: n

    n = 32
    a = 1.0
    b = pi
    final_result = 0.0
    exact_result = 0.198557

    call GET_COMMAND_ARGUMENT(0,name)
    write(*, *) 'executable: ', name

    call trapezoidal_rule(n, a, b, final_result)
    write(*, *) 'The final result is ', final_result
    write(*, *) 'The error is ', abs(exact_result-final_result)
end program main