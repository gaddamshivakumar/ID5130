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
    do i = 1, ln-1
        x = la + i*h
        total = total + func(x)
    end do
    trapz_proc = total*h
end function trapz_proc

real(kind=dp) function simpsons(la, lb, ln, h)
    implicit none
    integer, intent(in) :: ln
    real(kind=dp), intent(in) :: la, lb, h
    integer :: i
    real(kind=dp) :: total, x

    total = 0.0_dp
    do i = 1, ln
        x = la + i*h
        if (mod(i,2).eq.0) then
            total = total + 2*func(x)
        else
            total = total + 4*func(x)
        end if
    end do
    simpsons = total
end function simpsons
end module procedures

program problem1
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use :: mpi
    use :: procedures

    implicit none
    real(kind=dp), parameter :: pi = 3.14159265358_dp, ana_val = 0.198557298811136
    integer :: i, myid, nprocs, mpierror, status(MPI_STATUS_SIZE)
    integer :: n, ln, iproc
    real(kind=dp) :: a, b, final_result, h, la, lb, lsum


    call MPI_INIT(mpierror)
 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

    if (myid.eq.0) then
        write(*,*) "enter the value of n(integer): "
        read(*,*) n
        write(*,*) "enter the values of a, b (real): "
        read(*,*) a, b
    end if

    if (nprocs.gt.1) then
        if (myid.eq.0) then
            do iproc = 1, nprocs-1
                call mpi_send(n, 1, mpi_integer, iproc, 1, mpi_comm_world, mpierror)
                call mpi_send(a, 1, mpi_double, iproc, 2, mpi_comm_world, mpierror)
                call mpi_send(b, 1, mpi_double, iproc, 3, mpi_comm_world, mpierror)
            end do
        else
            call mpi_recv(n, 1, mpi_integer, 0, 1, mpi_comm_world, status, mpierror)
            call mpi_recv(a, 1, mpi_double, 0, 2, mpi_comm_world, status, mpierror)
            call mpi_recv(b, 1, mpi_double, 0, 3, mpi_comm_world, status, mpierror)
        end if
    end if

    final_result = 0.0_dp

    h = (b-a)/n
    ln = n/nprocs

    la = a + myid * ln*h
    lb = la + ln*h
    lsum = simpsons(la, lb, ln, h)
    
    call mpi_reduce(lsum, final_result, 1, mpi_double, mpi_sum, 0, mpi_comm_world, mpierror)
    
    if (myid.eq.0) then
        final_result = h*(func(a)-func(b)+final_result)/3
        print *, "The area is equal to", final_result
        print *, "The error: ", abs(ana_val-final_result)*100/ana_val, "%"
    end if
    call MPI_FINALIZE(mpierror)
end program problem1
