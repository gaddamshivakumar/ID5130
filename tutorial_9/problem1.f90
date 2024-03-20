program problem1
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use :: mpi

    implicit none
    integer :: i, myid, nprocs, mpierror, status(MPI_STATUS_SIZE)
    integer :: iproc
    real(kind=dp) :: a, b, a_reduce, b_allreduce
    integer, allocatable :: array1(:), array2(:), array3(:), array_scattered(:)


    call MPI_INIT(mpierror)
 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

    ! testing mpi_reduce
    a = myid*1+1    
    call mpi_reduce(a, a_reduce, 1, mpi_double, mpi_sum, 0, mpi_comm_world, mpierror)
    
    if (myid.eq.0) then
        print *, a_reduce, nprocs*(nprocs+1)/2
    end if

    ! testing mpi_allreduce
    b = (myid+1.0_dp)*2.0_dp
    call mpi_allreduce(b, b_allreduce, 1, mpi_double, mpi_sum, mpi_comm_world, mpierror)
    print *, 'myid: ', myid, 'b = ', b, '   allreduce sum = ', b_allreduce

    ! testing mpi_scatter
    allocate(array1(nprocs*3), array2(nprocs*3),  array3(nprocs*3), array_scattered(3))
    if (myid==0) then
        array1 = [(i, i=1,nprocs*3)]
        write(*,'("myid: ",i2,"    array1 = ", *(i4))') myid, array1
    end if
    call mpi_scatter(array1,3,mpi_integer,array_scattered,3,mpi_integer,0,mpi_comm_world,mpierror)
    print *, 'myid: ', myid, 'scattered array1 = ', array_scattered

    ! tesing mpi_gather
    array_scattered = array_scattered + 1 + myid
    call mpi_gather(array_scattered,3,mpi_integer,array2,3,mpi_integer,0,mpi_comm_world,mpierror)
    if (myid==0) write(*,'("myid: ",i2,"    array2 = ", *(i4))') myid, array2

    ! testing mpi_allgather
    array_scattered = array_scattered - 1 - myid
    call mpi_allgather(array_scattered,3,mpi_integer,array3,3,mpi_integer,mpi_comm_world,mpierror)
    write(*,'("myid: ",i2,"    array3 = ", *(i4))') myid, array3

    call MPI_FINALIZE(mpierror)
end program problem1
