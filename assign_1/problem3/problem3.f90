module procedures
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use :: omp_lib
    implicit none
contains

! function (analytical) values of phi -----------------------------------------------
real(kind=dp) function fun_phi(x, y)
    implicit none
    real(kind=dp), intent(in) :: x, y

    fun_phi = (x**2 - 1.0_dp)*(y**2 - 1.0_dp)
end function fun_phi

! function values of q -------------------------------------------------
elemental real(kind=dp) function fun_q(x, y)
    implicit none
    real(kind=dp), intent(in) :: x, y

    fun_q = 2.0_dp * (2.0_dp - x**2 - y**2)
end function fun_q

! Guass-Seidal serial --------------------------------------------------
subroutine gs_serial(npoints, delta, tolerance, xpoints, ypoints, phi, niters)
    implicit none
    integer, intent(in) :: npoints
    integer, intent(out) :: niters
    real(kind=dp), intent(in) :: delta, tolerance, xpoints(:), ypoints(:)
    real(kind=dp), intent(inout) :: phi(:,:)
    integer :: i, j, ix, jy
    real(kind=dp) :: error, phi_ana
    logical :: next_iter = .true.

    niters = 0
    do while(next_iter)     ! outer iteration loop
        next_iter = .false.
        ! loop over the grid points
        do jy = 2, npoints-1
            do ix = 2, npoints-1
                i = (npoints+1) - jy
                j = ix
                phi_ana = fun_phi(xpoints(ix),ypoints(jy))
                phi(i,j) = 0.25_dp*(phi(i,j+1)+phi(i,j-1)+phi(i-1,j)+phi(i+1,j)) + &
                0.25_dp*(delta**2)*fun_q(xpoints(ix),ypoints(jy))
                error = abs((phi(i,j)-phi_ana)/phi_ana)
                if (error > tolerance) then
                    next_iter = .true.
                end if
            end do
        end do
        niters = niters + 1
    end do

end subroutine gs_serial

! Gauss-Seidal parallel : diagonal approach ----------------------------
subroutine gs_diagonal(thread_count, npoints, delta, tolerance, xpoints, ypoints, phi, niters)
    implicit none
    integer, intent(in) :: npoints, thread_count
    integer, intent(out) :: niters
    real(kind=dp), intent(in) :: delta, tolerance, xpoints(:), ypoints(:)
    real(kind=dp), intent(inout) :: phi(:,:)
    integer :: i, j, ix, jy, ix_start, ix_end, ndiag, diag
    real(kind=dp) :: error, phi_ana
    logical :: next_iter = .true.

    ndiag = 2*npoints-1 ! no. of diagonals

    niters = 0
    do while(next_iter)     ! outer iteration loop
        next_iter = .false.
        ! loop over the diagonals
        !$omp parallel num_threads(thread_count) default(shared) private(ix,jy,i,j,phi_ana,error)
        do diag = 2, ndiag-1
            if (diag<=npoints) then
                ix_start = 1; ix_end = diag
            else
                ix_start = diag - npoints + 1; ix_end = npoints
            end if
            ! the parallel loops
            !$omp do
            do ix = ix_start+1, ix_end-1
                jy = diag-ix+1
                i = (npoints+1) - jy
                j = ix
                phi_ana = fun_phi(xpoints(ix),ypoints(jy))
                phi(i,j) = 0.25_dp*(phi(i,j+1)+phi(i,j-1)+phi(i-1,j)+phi(i+1,j)) + &
                0.25_dp*(delta**2)*fun_q(xpoints(ix),ypoints(jy))
                error = abs((phi(i,j)-phi_ana)/phi_ana)
                if (error > tolerance) then
                    next_iter = .true.
                end if
            end do
            !$omp end do
        end do
        !$omp end parallel
        niters = niters + 1
    end do
end subroutine gs_diagonal

! Gauss-Seidal parallel : red-black approach ---------------------------
subroutine gs_redblack(thread_count, npoints, delta, tolerance, xpoints, ypoints, phi, niters)
    implicit none
    integer, intent(in) :: npoints, thread_count
    integer, intent(out) :: niters
    real(kind=dp), intent(in) :: delta, tolerance, xpoints(:), ypoints(:)
    real(kind=dp), intent(inout) :: phi(:,:)
    integer :: i, j, ix, jy
    real(kind=dp) :: error, phi_ana
    logical :: next_iter = .true.

    niters = 0
    do while(next_iter)     ! outer iteration loop
        next_iter = .false.
        !$omp parallel num_threads(thread_count) default(shared) private(ix,jy,i,j,phi_ana,error)

        ! calculation of red grid points -- parallel loops
        !$omp do
        do jy = 2, npoints-1
            do ix = 2, npoints-1
                i = (npoints+1) - jy
                j = ix
                if (mod(ix+jy,2)==1) then
                    phi_ana = fun_phi(xpoints(ix),ypoints(jy))
                    phi(i,j) = 0.25_dp*(phi(i,j+1)+phi(i,j-1)+phi(i-1,j)+phi(i+1,j)) + &
                    0.25_dp*(delta**2)*fun_q(xpoints(ix),ypoints(jy))
                    error = abs((phi(i,j)-phi_ana)/phi_ana)
                    if (error > tolerance) then
                        next_iter = .true.
                    end if
                end if
            end do
        end do
        !$omp end do

        ! calculation of black grid points -- parallel loops
        !$omp do
        do jy = 2, npoints-1
            do ix = 2, npoints-1
                i = (npoints+1) - jy
                j = ix
                if (mod(ix+jy,2)==0) then
                    phi_ana = fun_phi(xpoints(ix),ypoints(jy))
                    phi(i,j) = 0.25_dp*(phi(i,j+1)+phi(i,j-1)+phi(i-1,j)+phi(i+1,j)) + &
                    0.25_dp*(delta**2)*fun_q(xpoints(ix),ypoints(jy))
                    error = abs((phi(i,j)-phi_ana)/phi_ana)
                    if (error > tolerance) then
                        next_iter = .true.
                    end if
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel
        niters = niters + 1
    end do

end subroutine gs_redblack

! print any 2nd rank matrix --------------------------------------------
subroutine print_matrix(mat)
    implicit none
    real(kind=dp), intent(in) :: mat(:,:)
    integer :: i, size(2)

    size = shape(mat)
    do i = 1, size(1)
        write(*, '(*(f12.8))') mat(i,:)
    end do
end subroutine print_matrix
end module
!--------------------------------------------------------- end of module



! ----------------------------------------------------------------------
!                              MAIN PROGRAM 
!-----------------------------------------------------------------------
program problem3
    use, intrinsic :: iso_fortran_env, only : dp => real64, error_unit
    use :: procedures
    implicit none
    integer :: unit_nr, istat
    character(len=1024) :: msg
    integer, parameter :: thread_count=16
    real(kind=dp), parameter :: delta=0.005_dp, xleft=-1.0_dp, xright=1.0_dp, yleft=-1.0_dp, &
    yright=1.0_dp, tolerance = 0.01_dp
    integer :: i, j, ix, jy, npoints, niters_gss, niters_gsd, niters_gsrb
    real(kind=dp) :: tstart, tend
    real(kind=dp), allocatable :: xpoints(:), ypoints(:) 
    real(kind=dp), allocatable :: phi_ana(:,:), phi_gs_serial(:,:), phi_gs_diag(:,:), phi_gs_rb(:,:),&
    dummy(:,:)

    npoints = int((xright-xleft)/delta) + 1

    ! allocate the memory
    allocate(xpoints(1:npoints), ypoints(1:npoints))
    allocate(phi_ana(1:npoints,1:npoints), phi_gs_serial(1:npoints,1:npoints))
    allocate(phi_gs_diag(1:npoints,1:npoints), phi_gs_rb(1:npoints,1:npoints))
    allocate(dummy(1:npoints,1:npoints))

    ! initialize the grid points
    xpoints = [(xleft+(i-1)*delta, i = 1, npoints)]
    ypoints = [(yleft+(i-1)*delta, i = 1, npoints)]

    ! initialize/incorporate the boundary conditions
    phi_ana = 0.0_dp
    phi_gs_serial = 0.0_dp
    phi_gs_diag = 0.0_dp    
    phi_gs_rb = 0.0_dp

    ! print *, npoints
    ! write(*,'(*(f10.4))') xpoints
    ! write(*,'(*(f10.4))') ypoints

    ! analytical solution
    do jy = 2, npoints-1
        do ix = 2, npoints-1
            i = (npoints+1) - jy
            j = ix
            phi_ana(i,j) = fun_phi(xpoints(ix),ypoints(jy))
        end do
    end do
    ! call print_matrix(phi_ana)

    ! call cpu_time(tstart)
    ! do i = 1, 1000000  ! for calculation of time
    ! call gs_serial(npoints, delta, tolerance, xpoints, ypoints, phi_gs_serial, niters_gss)
    ! end do
    ! call cpu_time(tend)
    ! print *, 'GS-serial time: ',tend-tstart

    ! call cpu_time(tstart)
    ! do i = 1, 2  ! for calculation of time
    ! call gs_diagonal(thread_count, npoints, delta, tolerance, xpoints, ypoints, phi_gs_diag, niters_gsd)
    ! end do
    ! call cpu_time(tend)
    ! print *, 'GS-diag time: ',tend-tstart

    call cpu_time(tstart)
    do i = 1, 2  ! for calculation of time
    call gs_redblack(thread_count, npoints, delta, tolerance, xpoints, ypoints, phi_gs_rb, niters_gsrb)
    end do
    call cpu_time(tend)
    print *, 'GS-redblack time: ',tend-tstart
    

    stop
    ! ------------------------------------------------------------------
    ! Writing data to matlab stuct file
    ! ------------------------------------------------------------------
    open (newunit=unit_nr, file='output_n21.m', access='sequential', action='write', status='new', &
    form='formatted', iostat=istat, iomsg=msg)
    if (istat /= 0) then
        write (unit=error_unit, fmt='(2A)') 'error: ', trim(msg) !2A is two strings
        stop 1
    end if
    ! number of points
    write (unit=unit_nr, fmt=*) 'data.npoints = ',npoints,';'
    ! xpoints
    write (unit=unit_nr, fmt=*) 'data.xpoints = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) xpoints(i),';'
    end do
    write (unit=unit_nr, fmt=*) '];'
    ! ypoints
    write (unit=unit_nr, fmt=*) 'data.ypoints = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) ypoints(i),';'
    end do
    write (unit=unit_nr, fmt=*) '];'
    ! analytical solution
    write (unit=unit_nr, fmt=*) 'data.phi_ana = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) phi_ana(i, :),';'
    end do
    write (unit=unit_nr, fmt=*) '];'
    ! Gauss seidal serial solution, no of iterations
    write (unit=unit_nr, fmt=*) 'data.niters_gss = ',niters_gss,';'
    write (unit=unit_nr, fmt=*) 'data.phi_gs_serial = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) phi_gs_serial(i, :),';'
    end do
    write (unit=unit_nr, fmt=*) '];'
    ! Gauss seidal diagonal parallel solution, no of iterations
    write (unit=unit_nr, fmt=*) 'data.niters_gss = ',niters_gsd,';'
    write (unit=unit_nr, fmt=*) 'data.phi_gs_diag = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) phi_gs_diag(i, :),';'
    end do
    write (unit=unit_nr, fmt=*) '];'
    ! Gauss seidal redblack parallel solution, no of iterations
    write (unit=unit_nr, fmt=*) 'data.niters_gss = ',niters_gsrb,';'
    write (unit=unit_nr, fmt=*) 'data.phi_gs_rb = ['
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) phi_gs_rb(i, :),';'
    end do
    write (unit=unit_nr, fmt=*) '];'

    close (unit=unit_nr)
    ! ------------------------------------------------------------------
end program problem3