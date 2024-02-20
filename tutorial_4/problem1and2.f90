program mataddition_multiplication
    use, intrinsic ::  iso_fortran_env, only : DP => REAL64
    implicit none
    integer :: N, thread_count, i, j, k
    real(kind=dp) :: rnum, t1, t2
    real(kind=dp), allocatable :: matA(:,:), matB(:,:), matC(:,:), matD(:,:)

    ! read the dimension of the matrix from the user
    write(*, *) 'Please enter the number of rows and number of columns:'
    read(*, *) N
    ! read the no of threads for parallelization
    write(*, *) 'Please enter the number of threads:'
    read(*, *) thread_count

    !allocate the memory
    allocate(matA(1:N,1:N), matB(1:N,1:N), matC(1:N,1:N), matD(1:N,1:N))
    matA = 0.0_dp; matB = 0.0_dp; matC = 0.0_dp

    !populate the matrices
    do j = 1, N
        do i = 1, N
            call random_number(rnum)
            matA(i,j) = rnum * 10
            call random_number(rnum)
            matB(i,j) = rnum * 10
        end do
    end do

    ! write(*,*) ' '
    ! do i = 1, N
    ! write(*,'(*(f10.6))') matA(i,:)
    ! end do
    ! write(*,*) ' '
    ! do i = 1, N
    ! write(*,'(*(f10.6))') matB(i,:)
    ! end do

    ! matrix addition
    call cpu_time(t1)
    !$omp parallel do num_threads(thread_count)default(none) shared(matA, matB, matC, N) private(i,j)
    do j = 1, N
        do i = 1, N
            matC(i,j) = matA(i,j) + matB(i,j)
        end do
    end do
    !$omp end parallel do
    call cpu_time(t2)
    write(*,*) 'time taken: ',t2-t1

    ! write(*,*) ' '
    ! do i = 1, N
    ! write(*,'(*(f10.6))') matC(i,:)
    ! end do

    ! matrix multiplication
    call cpu_time(t1)
    !$omp parallel do num_threads(thread_count) default(none) shared(matA, matB, matD, N) private(i,j,k)
    do j = 1, N
        do i = 1, N
            do k = 1, N
                matD(i,j) = matD(i,j) + matA(i,k)*matB(k,j)
            end do
        end do
    end do
    !$omp end parallel do
    call cpu_time(t2)
    write(*,*) 'time taken: ',t2-t1
end program mataddition_multiplication