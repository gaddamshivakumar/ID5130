program jacobi_serial
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer :: t1, t2, i, j, iter, n, iter_max = 500
    real(kind=dp), allocatable :: matA(:,:), vecb(:,:), vecx_old(:,:), vecx_new(:,:)

    ! read the dimension of the matrix from the user
    write(*, *) 'Please enter the size N:'
    read(*, *) n

    ! allocate the memory
    allocate(matA(1:n,1:n), vecb(1:n,1), vecx_old(1:n,1), vecx_new(1:n,1))

    ! initialize b-vector
    vecb = 0.0_dp; vecb(1,1) = 1

    ! initialize the initial guess for x-vector
    vecx_old = 0.0_dp

    ! initialize A-matrix
    do j = 1, n
        do i = 1, n
            if (i.eq.j) then
                matA(i,j) = i + j
            elseif (i.eq.1.and.j.eq.1) then
                matA(i,j) = 1
            elseif (i.eq.n.and.j.eq.2) then
                matA(i,j) = 2*n - 1
            else
                matA(i,j) = 1/n
            end if
        end do
    end do

    do iter = 1, iter_max

    end do
end program jacobi_serial