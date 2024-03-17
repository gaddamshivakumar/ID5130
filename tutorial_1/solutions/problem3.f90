subroutine multiplyMatrixVector(n, mat1, vec1, matvecProduct, compute_time)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: mat1(1:n, 1:n), vec1(1:n)
    double precision, dimension(1:n), intent(out) :: matvecProduct
    double precision, intent(out) :: compute_time
    integer :: i, j
    double precision :: t1, t2

    ! time at the beginning
    call cpu_time(t1)
  
    ! initialize elements of matProduct
    do i = 1, n
        matvecProduct(i) = 0.0
    end do
  
    ! multiply the matrix and the vector
    do j = 1, n
       do i = 1, n
             matvecProduct(i) = matvecProduct(i) + mat1(i, j) * vec1(j)
       end do
    end do
    
    ! time at the end 
    call cpu_time(t2)
    compute_time =  t2-t1
end subroutine multiplyMatrixVector

subroutine printMatrix(m, n, mat)
    implicit none
    integer, intent(in) :: m, n 
    double precision, dimension(1:m, 1:n), intent(in) :: mat
    integer :: i, j 
    
    write(*, *) 'Printing matrix:'
    
    do i = 1, m
       write(*, '(100f10.4)') (mat(i, j), j = 1, n)
    end do
end subroutine printMatrix

subroutine printVector(n, vec)
    implicit none
    integer, intent(in) :: n 
    double precision, dimension(1:n), intent(in) :: vec
    integer :: i
    
    write(*, *) 'Printing vector:'
    do i = 1, n
        write(*, '(100f10.4)') vec(i)
    end do
end subroutine printVector

program main
    implicit none
    double precision :: PI = 4.0*atan(1.0), rnum, compute_time
    double precision, allocatable :: matA(:,:), vecX(:), matvecAX(:)
    integer :: p, i, j
    logical :: printcond

    ! read the dimension of the matrix from the user
    write(*, *) 'Please enter the size N:'
    read(*, *) p

    ! check 
    if (p <= 0) then 
        write(*, *) 'Matrix dimension can only be positive'
        stop
    end if

    ! allocate memory
    allocate (matA(1:p,1:p), vecX(1:p), matvecAX(1:p))

    ! populate the matrix
    do j = 1, p
        do i = 1, p
            ! ii = i-1
            ! jj = j-1
            ! matA(i, j) = 0.5**(0.5*ii) * sin(ii*jj*PI/(p+1))
            call random_number(rnum)
            matA(i, j) = rnum
        end do
    end do
    ! populate the vector
    do i = 1, p
        ! ii = i-1
        ! jj = j-1
        ! matA(i, j) = 0.5**(0.5*ii) * sin(ii*jj*PI/(p+1))
        call random_number(rnum)
        vecX(i) = rnum
    end do

    if (p<=10) printcond = .true.
    if (printcond) then
        call printmatrix(p, p, matA)
        call printvector(p, vecX)
    endif

    call multiplyMatrixVector(p, matA, vecX, matvecAX, compute_time)
    if (printcond) call printvector(p, matvecAX)
    write(*,*) 'time taken for compute: ',compute_time,' sec'
    
end program main
  