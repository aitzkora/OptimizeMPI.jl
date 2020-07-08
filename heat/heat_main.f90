program heat_main
    use heat_solve
    use mpi
    use iso_c_binding, only: c_int32_t, c_double
    implicit none
    integer(c_int32_t) :: n,  proc,  ierr
    integer :: rank_w, size_w, narg, i
    logical :: verbose
    integer :: snapshot_size, snapshot_step
    real(kind=c_double) :: t_final, h, error
    real(kind=c_double), allocatable :: u_target(:,:), x(:), xx(:,:), yy(:,:), boundary(:)
    character(len=32)::param, variable_format


    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)

    verbose = (rank_w == 0)

    narg = command_argument_count()
    if (narg < 3) then
        call usage()
    end if

    call get_command_argument( 1, param )
    read (param,*) n 

    call get_command_argument( 2, param )
    read (param,*) t_final

    call get_command_argument( 3, param )
    read (param,*) proc

    if (proc * proc /= size_w)   then
        if (verbose) print *, 'the total number of processus differs from the product proc * proc ', & 
            & proc, ' x ' , proc , '/= ', size_w
        call MPI_FINALIZE( ierr )
        stop
    endif
   
    h = 1.d0 / (n+1) 
    x = [ (h*i, i=1, n) ]
    xx = spread(x, 2, n)
    yy = reshape(spread( x, 1, n), [n,n])
    u_target = exp( xx ) * sin( yy )
    allocate ( boundary(4*n+4) )
    call random_number( boundary )
    call compute_error(n, proc, t_final, boundary, u_target, error) 
    if (rank_w == 0) then
        print *, "error = ", error
    end if
    call MPI_FINALIZE( ierr )

contains

    subroutine usage()
        implicit none
        character(len=255)::name

        call get_command_argument(0,name)
        print *, 'Usage: mpirun -np proc * proc ', TRIM(name), ' n t_final proc '
        print *, '    n   number of discretisation points in each direction'
        print *, '    t_final  total duration of the simulation'
        print *, '    proc  number of MPI processes in each direction'
        stop
    end subroutine usage

end program heat_main
