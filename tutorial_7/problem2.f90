program main
    implicit none
    include "mpif.h"
    integer :: i, myid, size, mpierror, tag, status(MPI_STATUS_SIZE)
    character(len=50) :: message_send, message_recv

    call MPI_INIT(mpierror)
 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)
    
end program main