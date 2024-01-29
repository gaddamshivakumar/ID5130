subroutine multiplyTwoMatrices(m, n, o, mat1, mat2, matProduct)
    implicit none
    integer, intent(in) :: m, n, o
    double precision, dimension(1:m, 1:n), intent(in) :: mat1
    double precision, dimension(1:n, 1:o), intent(in) :: mat2
    double precision, dimension(1:m, 1:o), intent(out) :: matProduct
    integer :: i, j, k
    double precision :: t1, t2 

    ! check the dimensions
    ! if (m /= n) then 
    ! write(*, *) 'Column size of mat1 is not equal to Row size of mat2'
    ! write(*, *) 'Cannot multiply the matrices, please choose m = n'
    ! stop
    ! end if

    ! time at the beginning
    call cpu_time(t1)

    ! initialize elements of matProduct
    do j = 1, o
    do i = 1, m 
        matProduct(i, j) = 0.0
    end do
    end do

    ! multiply the two matrices
    do j = 1, o
    do i = 1, m
        do k = 1, n
            matProduct(i, j) = matProduct(i, j) + mat1(i, k) * mat2(k, j)
        end do
    end do
    end do
    
    ! time at the end 
    call cpu_time(t2)

    write(*, *) 'Multiplication of matrices took ', t2-t1, ' seconds'
end subroutine multiplyTwoMatrices

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

subroutine transpose(m, n, mat, matT)
    implicit none
    double precision, dimension(m, n), intent(in) :: mat
    double precision, dimension(n, m), intent(out) :: matT
    integer :: m, n, i, j

    ! m = shape(mat(:,1))
    ! n = shape(mat(:,2))

    do j = 1, n
        do i = 1, m
            matT(j, i) = mat(i, j)
        end do
    end do
end subroutine

subroutine isequal(m, n, mat1, mat2, equality)
    implicit none
    double precision, dimension(1:m, 1:n), intent(in) :: mat1, mat2
    integer :: i, j, m, n
    logical, intent(out) :: equality

    equality = .true.
    do j = 1, n
        do i = 1, m
            if (mat1(i,j) .ne. mat2(i,j)) then
                equality = .false.
            end if
        end do
    end do

end subroutine

logical function fisequal(m, n, mat1, mat2)
    implicit none
    double precision, dimension(1:m, 1:n), intent(in) :: mat1, mat2
    integer :: i, j, m, n

    fisequal = .true.
    do j = 1, n
        do i = 1, m
            if (mat1(i,j) .ne. mat2(i,j)) then
                fisequal = .false.
            end if
        end do
    end do

end function fisequal


program main
    implicit none
    double precision :: PI = 4.0*atan(1.0), rnum
    double precision, dimension(:, :), allocatable :: matA, matB, matAB, matAT, matBT, matABT, matBTAT
    integer :: p, q, r, i, j, ii, jj
    logical :: printcond, equality

    ! read the dimension of the matrix from the user
    write(*, *) 'Please enter the number of rows and number of columns of A:'
    read(*, *) p, q
    write(*, *) 'Please eneter the number of columns of B:'
    read(*, *) r

    ! check 
    if (P <= 0 .or. q <= 0 .or. r <= 0) then 
        write(*, *) 'Matrix dimension can only be positive'
        stop
    end if

    ! allocate memory
    allocate (matA(1:p,1:q), matB(1:q,1:r), matAT(1:q,1:p), matBT(1:r,1:q))
    allocate (matAB(1:p,1:r), matABT(1:r,1:p), matBTAT(1:r,1:p))

    ! populate the matrices
    do j = 1, q
        do i = 1, p
            ii = i-1
            jj = j-1
            ! matA(i, j) = 0.5**(0.5*ii) * sin(ii*jj*PI/(p+1))
            call random_number(rnum)
            matA(i, j) = rnum
        end do
    end do
    do j = 1, r
        do i = 1, q
            ii = i-1
            jj = j-1
            ! matB(i, j) = 0.5**(0.5*ii) * cos(ii*jj*PI/(q+1))
            call random_number(rnum)
            matB(i, j) = rnum
        end do
    end do

    if (p<=6 .and. q<=6 .and. r<=6) printcond = .true.

    if (printcond) then
        call printmatrix(p, q, matA)
        call printmatrix(q, r, matB)
    end if

    call transpose(p, q, matA, matAT)
    if (printcond) call printmatrix(q, p, matAT)

    call transpose(q, r, matB, matBT)
    if (printcond) call printmatrix(r, q, matBT)

    call multiplyTwoMatrices(p, q, r, matA, matB, matAB)
    if (printcond) call printmatrix(p, r, matAB)

    call transpose(p, r, matAB, matABT)
    if (printcond) call printmatrix(r, p, matABT)

    call multiplyTwoMatrices(r, q, p, matBT, matAT, matBTAT)
    if (printcond) call printmatrix(r, p, matBTAT)

    call isequal(r, p, matABT, matBTAT, equality)
    ! equality = fisequal(r, p, matABT, matBTAT) !not working
    write(*, *) ' ---- ---- ---- '
    if (equality) then
        write(*, *) 'transpose(AB) = transpose(B)*transpose(A) -- verified!'
    else
        write(*, *) 'something went wrong!'
    endif
    write(*, *) ' ---- ---- ---- '
    
end program main