program test
    implicit none
    real(kind=8) :: a(3,2), i

    a(:,1) = [1.0, 2.0, 3.0]
    a(:,2) = [4.0, 5.0, 6.0]

    ! do i = 1,3
    !     print *, a(i,:)
    ! end do

    print *, a(2,-1)

    

end program