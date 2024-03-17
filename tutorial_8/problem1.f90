module procedures
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none    
contains
    
real(kind=dp) function func(x)
    implicit none
    real(kind=dp) :: x

    func = sin(x)/(2.0_dp*x**3)
end function func

real(kind=dp) function trapz_proc(la, lb, ln, h)
    implicit none
    integer, intent(in) :: ln
    real(kind=dp), intent(in) :: la, lb, h
    integer :: i
    real(kind=dp) :: total, x

    total = (func(la) + func(lb))/2.0_dp
    do i = 2, ln-1
        x = la + i*h
        total = total + func(x)
    end do
    trapz_proc = total
end function trapz_proc
end module procedures

program problem1
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use :: procedures
    implicit none
    include "mpif.h"
    real(kind_dp), parameter :: pi = 3.14159265358_dp
    integer :: i, myid, nprocs, mpierror, status(MPI_STATUS_SIZE)
    integer :: n, ln, iproc
    real(kind=dp) :: a, b, final_result, h, la, lb, lsum


    call MPI_INIT(mpierror)
 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

    n = 1024
    a = 0.0_dp
    b = pi
    final_result = 0.0_dp

    h = (b-a)/n
    ln = n/nprocs

    la = a + myid * ln*h
    lb = la + ln*h
    lsum = trapz_proc(la, lb, ln, h)

    if (myid.ne.0) then
        call MPI_SEND(lsum, 1, mpi_double, 0, 0, MPI_COMM_WORLD, mpierror)
    else
        final_result = lsum
        do iproc = 1, nprocs
            call MPI_RECV(lsum, 1, mpi_double, iproc, 0, MPI_COMM_WORLD, status, mpierror)
            final_result = final_result + lsum
        end do
    end if
    if (myid.eq.0) then
        print *, "The area is equal to", final_result
    end if
    call MPI_FINALIZE(mpierror)
end program problem1
