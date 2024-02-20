module procedures
    use, intrinsic :: iso_fortran_env, only : DP => REAL64
    implicit none
        
contains
    
real(kind=8) function fun_u(x)
    implicit none
    real(kind=dp) :: x

    fun_u = 7 - x*tan(x)
end function fun_u

real(kind=8) function fun_du_ana(x)
    implicit none
    real(kind=dp) :: x

    fun_du_ana = -tan(x) - x/cos(x)**2
end function fun_du_ana

subroutine fun_du_1(du_1, xpoints, dx, thread_count, npoints)
    implicit none
    integer, intent(in) :: thread_count, npoints
    real(kind=dp), intent(in) :: xpoints(*), dx
    real(kind=dp), intent(out) :: du_1(*)
    integer :: ipoint

    du_1(1) = (fun_u(xpoints(2))-fun_u(xpoints(1)))/dx
    !$omp parallel do num_threads(thread_count)
    do ipoint = 2, npoints
        du_1(ipoint) = (fun_u(xpoints(ipoint))-fun_u(xpoints(ipoint-1)))/dx
    end do
    !$omp end parallel do
end subroutine fun_du_1

subroutine fun_du_2(du_2, xpoints, dx, thread_count, npoints)
    implicit none
    integer, intent(in) :: thread_count, npoints
    real(kind=dp), intent(in) :: xpoints(*), dx
    real(kind=dp), intent(out) :: du_2(*)
    integer :: ipoint

    du_2(1) = (fun_u(xpoints(2))-fun_u(xpoints(1)))/dx
    du_2(npoints) = (fun_u(xpoints(npoints))-fun_u(xpoints(npoints-1)))/dx
    !$omp parallel do num_threads(thread_count)
    do ipoint = 2, npoints-1
        du_2(ipoint) = (fun_u(xpoints(ipoint+1))-fun_u(xpoints(ipoint-1)))/(2*dx)
    end do
    !$omp end parallel do
end subroutine fun_du_2

subroutine fun_du_4(du_4, xpoints, dx, thread_count, npoints)
    implicit none
    integer, intent(in) :: thread_count, npoints
    real(kind=dp), intent(in) :: xpoints(*), dx
    real(kind=dp), intent(out) :: du_4(*)
    integer :: ipoint

    du_4(1) = (fun_u(xpoints(2))-fun_u(xpoints(1)))/dx
    du_4(2) = (fun_u(xpoints(3))-fun_u(xpoints(1)))/(2*dx)
    du_4(npoints-1) = (fun_u(xpoints(npoints))-fun_u(xpoints(npoints-2)))/(2*dx)
    du_4(npoints) = (fun_u(xpoints(npoints))-fun_u(xpoints(npoints-1)))/dx
    !$omp parallel do num_threads(thread_count)
    do ipoint = 3, npoints-2
        du_4(ipoint) = (fun_u(xpoints(ipoint-2))-8*fun_u(xpoints(ipoint-1))+8*fun_u(xpoints(ipoint+1))&
        -fun_u(xpoints(ipoint+2)))/(12*dx)
    end do
    !$omp end parallel do
end subroutine fun_du_4
end module procedures



program fderivative
    use, intrinsic :: iso_fortran_env, only : DP => REAL64
    use :: procedures
    implicit none
    integer :: i, npoints, thread_count
    real(kind=dp), parameter :: xleft=-1, xright=1
    real(kind=dp) :: dx
    real(kind=dp), allocatable :: xpoints(:), du_ana(:), du_1(:), du_2(:), du_4(:)
    ! real(kind=dp) :: fun_u, fun_du_ana

    ! read the gridsize dx
    write(*,*) "enter the grid size 'dx':"
    read(*,*) dx
    ! read the number of therads
    write(*,*) "enter the number of threads:"
    read(*,*) thread_count

    npoints = int((xright-xleft)/dx) + 1

    ! allocate the memory
    allocate(xpoints(1:npoints), du_ana(1:npoints), du_1(1:npoints), du_2(1:npoints), du_4(1:npoints))
    xpoints = 0.0_dp; du_ana = 0.0_dp; du_1 = 0.0_dp; du_2 = 0.0_dp; du_4 = 0.0_dp

    xpoints = [((xleft+(i-1)*dx), i = 1,npoints)]
    du_ana = [(fun_du_ana(xpoints(i)), i = 1,npoints)]
    call fun_du_1(du_1, xpoints, dx, thread_count, npoints)
    call fun_du_2(du_2, xpoints, dx, thread_count, npoints)
    call fun_du_4(du_4, xpoints, dx, thread_count, npoints)

    print '(*(f6.2))', xpoints(1:25)
    print '(*(f6.2))', du_ana(1:25)
    print '(*(f6.2))', abs((du_ana(1:25) - du_1(1:25))/du_ana(1:25))*100
    print '(*(f6.2))', abs((du_ana(1:25) - du_2(1:25))/du_ana(1:25))*100
    print '(*(f6.2))', abs((du_ana(1:25) - du_4(1:25))/du_ana(1:25))*100

end program fderivative