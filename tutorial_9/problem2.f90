 program problem2
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use :: mpi

    implicit none
    integer :: myid, nprocs, mpierror, status(MPI_STATUS_SIZE)
    integer :: i, iproc, nsize, prows
    real(kind=dp), allocatable :: amat(:,:), amat_rowblock(:,:), xvec(:), bvec_serial(:), bvec_mpi(:)&
    , bvec_process (:)

    call MPI_INIT(mpierror)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

    if (myid.eq.0) then
        write(*,*) "enter the size of the matrix (nsize):"
        read(*,*) nsize
    end if

    call mpi_bcast(nsize,1,mpi_integer,0,mpi_comm_world,mpierror)
    allocate(xvec(nsize))

    if (myid.eq.0) then    
        allocate(amat(nsize,nsize), bvec_serial(nsize), bvec_mpi(nsize))
        call random_number(amat)
        call random_number(xvec)
        amat = amat*10; xvec = xvec*10
        bvec_serial = matmul(amat,xvec)
        ! do i = 1, nsize
        !     print *, amat(i,:)
        ! end do
        ! print *, xvec
        ! print *, bvec_serial
        amat = transpose(amat)
    end if

    call mpi_barrier(mpi_comm_world,mpierror)
    call mpi_bcast(xvec,nsize,mpi_double,0,mpi_comm_world,mpierror)
    ! print *, 'myid: ', myid, 'x vector= ', xvec

    prows = nsize/nprocs

    call mpi_barrier(mpi_comm_world,mpierror)
    allocate(amat_rowblock(nsize,prows),bvec_process(prows))
    call mpi_scatter(amat,nsize*prows,mpi_double,amat_rowblock,nsize*prows,mpi_double,0,mpi_comm_world,mpierror)
    ! print *, 'myid: ', myid, 'amat rowblocks= ', amat_rowblock

    call mpi_barrier(mpi_comm_world,mpierror)
    bvec_process = matmul(transpose(amat_rowblock),xvec)
    call mpi_gather(bvec_process, prows, mpi_double, bvec_mpi, prows, mpi_double, 0, mpi_comm_world, mpierror)
    ! if (myid.eq.0) print *, 'printing mpi solution:', bvec_mpi-bvec_serial
    if (myid.eq.0) print *, 'max error: ', maxval(bvec_mpi-bvec_serial)

    call MPI_FINALIZE(mpierror)
end program problem2