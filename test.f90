program test
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer :: nsize = 7, nprocs = 3
    real(kind=dp) :: nblocks
    nblocks = floor(real(nsize)/nprocs)
    print *, nblocks

end program