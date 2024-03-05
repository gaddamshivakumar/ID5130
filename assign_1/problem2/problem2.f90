module procedures
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use omp_lib
    implicit none
contains

! function values ------------------------------------------------------
elemental real(kind=dp) function fun(x)
    implicit none
    real(kind=dp), intent(in) :: x

    fun = sin(5*x)
end function fun

! analytical values of derivative of the function ----------------------
elemental real(kind=dp) function dfun_ana(x)
    implicit none
    real(kind=dp), intent(in) :: x

    dfun_ana = 5*cos(5*x)
end function dfun_ana

! construction of system of linear equations from pade's scheme --------
subroutine pade_scheme(npoints, h, xpoints, amatrix, bvector)
    implicit none
    integer, intent(in) :: npoints
    real(kind=dp), intent(in) :: h, xpoints(:)
    real(kind=dp), intent(inout) ::  amatrix(:,:), bvector(:)
    integer :: i

    ! Construct the A-matrix
    amatrix(1,1) = 1.0_dp; amatrix(1,2) = 2.0_dp
    amatrix(npoints, npoints-1) = 2.0_dp; amatrix(npoints, npoints) = 1.0_dp
    do i = 2, npoints-1
        amatrix(i, i-1) = 1.0_dp
        amatrix(i, i) = 4.0_dp
        amatrix(i, i+1) = 1.0_dp
    end do
    ! Construct the b-vector
    bvector(1) = (1/h)*(-5.0_dp/2.0_dp*fun(xpoints(1)) + 2.0_dp*fun(xpoints(2)) + &
    1.0_dp/2.0_dp*fun(xpoints(3)))
    bvector(npoints) = (1/h)*(-1.0_dp/2.0_dp*fun(xpoints(npoints-2)) - &
    2.0_dp*fun(xpoints(npoints-1)) + 5.0_dp/2.0_dp*fun(xpoints(npoints)))
    do i = 2, npoints-1
        bvector(i) = (3.0_dp/h)*(fun(xpoints(i+1)) - fun(xpoints(i-1)))
    end do
end subroutine pade_scheme


! serial - lu decomposition --------------------------------------------
subroutine serial_ludecomp(npoints, amatrix, bvector, df_sol_s)
    implicit none
    integer, intent(in) :: npoints
    real(kind=dp), intent(in) :: amatrix(:,:), bvector(:)
    real(kind=dp), intent(out) :: df_sol_s(:)
    integer :: i
    real(kind=dp) :: lvec(npoints), uvec(npoints), zvec(npoints)!, lmat(npoints,npoints), umat(npoints,npoints)

    ! calculate l-elements and u-elements
    uvec(1) = amatrix(1,1)
    lvec(1) = 0.0_dp
    do i = 2, npoints
        lvec(i) = amatrix(i,i-1)/uvec(i-1)
        uvec(i) = amatrix(i,i) - lvec(i)*amatrix(i-1,i)
    end do
    ! solving for z such that Lz=b
    zvec(1) = bvector(1)
    do i = 2, npoints
        zvec(i) = bvector(i) - lvec(i)*zvec(i-1)
    end do
    ! solving for x such that Ux=z
    df_sol_s(npoints) = zvec(npoints)/uvec(npoints)
    do i = npoints-1, 1, -1
        df_sol_s(i) = (zvec(i)-amatrix(i,i+1)*df_sol_s(i+1))/uvec(i)
    end do
end subroutine serial_ludecomp

! parallel -- recursive doubling ---------------------------------------
subroutine recursive_doubling(thread_count, npoints, amatrix, bvector, df_sol_p)
    implicit none
    integer, intent(in) :: thread_count, npoints
    real(kind=dp), intent(in) :: amatrix(:,:), bvector(:)
    real(kind=dp), intent(out) :: df_sol_p(:)
    integer :: ieq, istep, nsteps, k, k_1
    real(kind=dp) :: alpha, beta, diag_a(npoints,2), diag_b(npoints,2), diag_c(npoints,2), &
    vec_y(npoints,2)

    ! initialize a, b, c diagonals
    diag_a = 0.0_dp; diag_b = 0.0_dp; diag_c = 0.0_dp; vec_y = 0.0_dp
    diag_a(npoints,1) = amatrix(npoints,npoints-1)
    diag_c(1,1) = amatrix(1,2)
    diag_b(1,1) = amatrix(1,1); diag_b(npoints,1) = amatrix(npoints,npoints)
    vec_y(:,1) = bvector(:)
    do ieq = 2, npoints-1
        diag_a(ieq,1) = amatrix(ieq,ieq-1)
        diag_b(ieq,1) = amatrix(ieq,ieq)
        diag_c(ieq,1) = amatrix(ieq,ieq+1)
    end do

    ! elimination phase: compute a, b, c, y
    nsteps = int(log(npoints*1.0_dp)/log(2.0_dp)) + 1 ! log(npoints)(2)
    !if npoints is not the powers of 2, nsteps would require one more step.

    !$omp parallel num_threads(thread_count) default(shared) private(ieq,alpha,beta)
    ! print *, 'no of threads: ', omp_get_num_threads()
    do istep = 1, nsteps
        k = mod(istep,2)+1
        if (k.eq.2) then
            k_1 = 1
        elseif (k.eq.1) then
            k_1 = 2
        end if
        !$omp do
        do ieq = 1, npoints ! the parallel loop
            ! print *, 'rank of the thread: ', omp_get_thread_num(), 'eq: ', ieq-1
            if ((2**(istep-1)+1.le.ieq).and.(ieq.le.npoints)) then
                alpha = -diag_a(ieq,k_1)/diag_b(ieq-2**(istep-1),k_1)
            else
                alpha = 0.0_dp
            end if
            if ((1.le.ieq).and.(ieq.le.(npoints-2**(istep-1)))) then
                beta = -diag_c(ieq,k_1)/(diag_b(ieq+2**(istep-1),k_1))
            else
                beta = 0.0_dp
            end if
            diag_a(ieq,k) = alpha*diag_a(ieq-2**(istep-1),k_1)
            diag_c(ieq,k) = beta*diag_c(ieq+2**(istep-1),k_1)
            diag_b(ieq,k) = alpha*diag_c(ieq-2**(istep-1),k_1) + diag_b(ieq,k_1) + beta*diag_a(ieq+2**(istep-1),k_1)
            vec_y(ieq,k) = alpha*vec_y(ieq-2**(istep-1),k_1) + vec_y(ieq,k_1) + beta*vec_y(ieq+2**(istep-1),k_1)
        end do
        !$omp end do
    end do
    !$omp end parallel

    ! solution phase: compute x (df_sol_p)
    do ieq = 1, npoints
        df_sol_p(ieq) = vec_y(ieq,k)/diag_b(ieq,k)
    end do
end subroutine recursive_doubling


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
end module procedures
!--------------------------------------------------------- end of module



! ----------------------------------------------------------------------
!                              MAIN PROGRAM 
!-----------------------------------------------------------------------
program problem2
    use, intrinsic :: iso_fortran_env, only : dp => real64, error_unit
    use :: procedures
    implicit none
    integer :: unit_nr, istat, i
    character(len=1024) :: msg

    integer, parameter :: npoints=1000, thread_count=1
    real, parameter :: xleft=0.0_dp, xright=3.0_dp
    real(kind=dp) :: tstart, tend, h
    real(kind=dp), dimension(npoints) :: xpoints, bvector, df_ana, df_sol_s, df_sol_p
    real(kind=dp), dimension(npoints,npoints) :: amatrix

    ! initialize the grid points
    h = (xright-xleft)/(npoints-1)
    xpoints = [(xleft+(i-1)*h, i = 1,npoints)]
    df_ana = dfun_ana(xpoints)

    ! print *, xpoints
    ! print *, fun(xpoints)
    ! print *, df_ana

    ! obtain the A-matrix and b-vector
    amatrix=0.0_dp; bvector=0.0_dp
    call pade_scheme(npoints, h, xpoints, amatrix, bvector)
    ! call print_matrix(amatrix)
    ! print *, bvector

    ! serial program
    call cpu_time(tstart)
    call serial_ludecomp(npoints, amatrix, bvector, df_sol_s)
    ! print *, df_sol_s
    call cpu_time(tend)
    print *, 'serial time: ',tend-tstart

    ! parallel program - threads
    call cpu_time(tstart)
    ! do i = 1, 100000  ! for calculation of time
    call recursive_doubling(thread_count, npoints, amatrix, bvector, df_sol_p)
    ! print *, df_sol_p
    ! end do
    call cpu_time(tend)
    print *, 'parallel time: ',tend-tstart

    print *, 'max error (analytical vs serial): ', maxval(abs(df_ana(10:90)-df_sol_s(10:90)))
    print *, 'max error (analytical vs parallel): ', maxval(abs(df_ana(10:90)-df_sol_p(10:90)))
    print *, 'max error (serial vs parallel): ', maxval(abs(df_sol_s-df_sol_p))

    ! write the output data
    open (newunit=unit_nr, file='output_n25.csv', access='sequential', action='write', status='new', &
    form='formatted', iostat=istat, iomsg=msg)
    if (istat /= 0) then
        write (unit=error_unit, fmt='(2A)') 'error: ', trim(msg) !2A is two strings
        stop 1
    end if
    do i = 1, npoints
        write (unit=unit_nr, fmt=*) xpoints(i), df_ana(i), df_sol_s(i), df_sol_p(i)
    end do
    close (unit=unit_nr)
end program problem2