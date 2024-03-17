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

subroutine addTwoMatrices(m, n, mat1, mat2, matSum)
    implicit none
    integer, intent(in) :: m, n
    double precision, dimension(1:m, 1:n), intent(in) :: mat1, mat2
    double precision, dimension(1:m, 1:n), intent(out) :: matSum
    integer :: i, j
    double precision :: t1, t2 
  
    ! time at the beginning
    call cpu_time(t1)
  
    ! add the two matrices
    do j = 1, n
       do i = 1, m
          matSum(i, j) = mat1(i, j) + mat2(i, j)
       end do
    end do
  
    ! time at the end 
    call cpu_time(t2)
    write(*, *) 'Addition of matrices took ', t2-t1, ' seconds'
end subroutine addTwoMatrices

subroutine transpose(m, n, mat, matT)
    implicit none
    double precision, dimension(m, n), intent(in) :: mat
    double precision, dimension(n, m), intent(out) :: matT
    integer :: m, n, i, j

    do j = 1, n
        do i = 1, m
            matT(j, i) = mat(i, j)
        end do
    end do
end subroutine

subroutine issymmetric(issym, n, mat)
    implicit none
    double precision, dimension(1:n, 1:n), intent(in) :: mat
    logical :: issym
    integer :: i, j, n

    issym = .true.
    do j = 1, n
        do i = j+1, n
            if ( mat(i,j) .ne. mat(j,i) ) then
                issym = .false.
            end if
        enddo
    enddo
end subroutine
            

program main
    implicit none
    double precision :: PI = 4.0*atan(1.0), rnum
    double precision, dimension(:, :), allocatable :: matA, matAT, matApAT
    integer :: p, i, j
    logical :: printcond, issym

    ! read the dimension of the matrix from the user
    write(*, *) 'Please enter the size of the square matrix "A":'
    read(*, *) p

    ! check 
    if (p <= 0) then 
        write(*, *) 'Matrix dimension can only be positive'
        stop
    end if

    ! allocate memory
    allocate (matA(1:p,1:p), matAT(1:p,1:p), matApAT(1:p,1:p))

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

    if (p<=10) printcond = .true.
    if (printcond) call printmatrix(p, p, matA)

    call transpose(p, p, matA, matAT)
    if (printcond) call printmatrix(p, p, matAT)

    call addTwoMatrices(p, p, matA, matAT, matApAT)
    if (printcond) call printmatrix(p, p, matApAT)

    call issymmetric(issym, p, matApAT)
    write(*, *) ' ---- ---- ---- '
    if (issym) then
        write(*, *) 'A + transpose(A) is symmetric -- verified!'
    else
        write(*, *) 'something went wrong!'
    endif
    write(*, *) ' ---- ---- ---- '

end program main