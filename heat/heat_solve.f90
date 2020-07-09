module heat_solve
    use iso_c_binding, only: c_int32_t, c_double
    use iso_fortran_env
    use mat_utils, only : print_mat
    use mpi
    use heat
    use communications
    implicit none
    public :: compute_error, set_boundary
 
contains

  subroutine compute_error(n, proc, t_final, boundary, u_target, error) bind( C, name="compute_error" )

    integer(c_int32_t), intent(in) :: n ! size of the global square matrix
    integer(c_int32_t), intent(in) :: proc ! nb of processes (in each dimensions)
    real(c_double), intent(in) :: t_final
    real(c_double), intent(in) :: u_target(1:n, 1:n)
    real(c_double), intent(in) :: boundary(4*n+4)
    real(c_double), intent(out) :: error

    integer(c_int32_t) :: i, n_loc, ierr
    integer(c_int32_t) :: rank_w, size_w, rank_2D, comm2D, type_row
    integer(c_int32_t), parameter :: ndims = 2, North = 1, South = 2, East = 3, West = 4
    logical :: is_master, reorder = .true.
    integer(c_int32_t), dimension(4) :: neighbour
    real(c_double) :: h, dt, error_loc, t
    real(c_double), allocatable :: u_in(:,:), u_out(:,:), sol_space(:,:)

    integer(c_int32_t), dimension(ndims) :: dims , coords
    logical, dimension(ndims) :: periods = .false.


    call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)

    is_master = (rank_w == 0)

    n_loc = 2 + n / proc

    h = 1.d0 / n
    dt = h ** 2 / 4.d0

    ! construction of the cartesion topology
    dims(1) = proc
    dims(2) = proc

    call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, periods, reorder, comm2D, ierr )
    call MPI_COMM_RANK( comm2D, rank_2D, ierr )

    !! Fetch the processus coordinates in the 2-D grid
    call MPI_CART_COORDS( comm2D, rank_2D, ndims, coords, ierr )

    !! Creation of a non-contiguous in memory column type
    !! to address Fortran storage: no stride of n_loc
    call MPI_TYPE_VECTOR( n_loc - 2, 1, n_loc, MPI_DOUBLE_PRECISION, type_row, ierr )
    call MPI_TYPE_COMMIT( type_row, ierr )

    !! fetching the neighbor ranks
    call MPI_CART_SHIFT(COMM2D, 0, 1, neighbour( North ), neighbour( South ), ierr )
    call MPI_CART_SHIFT(COMM2D, 1, 1, neighbour( West ), neighbour( East ), ierr )

    allocate( u_in(n_loc, n_loc) ) 
    allocate( u_out(n_loc, n_loc) )

    call set_boundary( coords, proc, u_in, boundary )
    call set_boundary( coords, proc, u_out, boundary )
    t = 0.d0
    do 
      if (t > t_final ) exit
      call heat_kernel( u_in, u_out )
      call ghosts_swap( comm2D, type_row, neighbour, u_in )
      u_in = u_out
      t = t + dt
    end do

    ! We gather the solution on process 0
    allocate ( sol_space(n, n) )
    call gather_solution( sol_space, n, u_in, ndims, comm2D, is_master )
    if (is_master) then
      error = sum( (sol_space-u_target)**2)
    end if
    call MPI_Bcast(error, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    deallocate( sol_space )

    deallocate( u_in )
    deallocate( u_out )

    call MPI_TYPE_FREE( type_row, ierr )

  end subroutine compute_error

  subroutine set_boundary(coo, proc, u, boundary)

    integer, intent(in) :: proc
    integer, intent(in) :: coo(2)
    real(c_double), intent(out) :: u(:,:)
    real(c_double), intent(in) :: boundary(:)
    integer :: n_loc, offset, offset_loc, n, i

    n_loc = size(u, 1)
    n  = (size(boundary, 1) - 4 )/ 4
    u = 0.d0
    print *, "coo = ", coo

    ! upper line :  range [1:n+2]
    if (coo(1) == 0)  then
       offset = coo(2)*(n_loc-2)
       u(1, 1:n_loc) = boundary(offset + 1: offset + n_loc )
    end if

    ! lower line : range -> [3n+3:4n+4]
    if (coo(1) == (proc - 1)) then 
        offset = 3*n+2 + coo(2) * (n_loc-2) 
        u(size( u, 1 ),  1 : n_loc  ) = boundary( offset + 1: offset + n_loc )
    endif

    !! left column : range -> [n+3:3n+2:2]  
    if (coo(2) == 0) then 
      offset = n+2 + coo(1) * (2*n_loc -4) 
      if (coo(1) /= 0) then
        u(1:n_loc, 1) = boundary( offset -1: offset + 2*n_loc-2 : 2)
      elseif (coo(1) == (proc -1) ) then
        u(1:n_loc-1, 1) = boundary( offset + 1 : offset + 2*n_loc-2 : 2)
      else
        u(2:n_loc, 1) = boundary( offset + 1 : offset + 2*n_loc-2 : 2)
      end if
    end if
    !! right column : range -> [n+4:3n+2:2]
    if (coo(2) == (proc - 1)) then 
        offset = n+3 + coo(1) * (2*n_loc -4) 
        if (coo(1) == (proc - 1 ) ) then
          u(1:n_loc-1, size(u, 2) ) = boundary( offset - 1 : offset + 2*n_loc-4 : 2)
        else if (coo(1) /= 0) then
          u(1:n_loc, size( u, 2 ) ) = boundary( offset - 1 : offset + 2*n_loc-2 : 2)
        else 
          u(2:n_loc, size( u, 2 ) ) = boundary( offset + 1 : offset + 2*n_loc-2 : 2)
        end if
    end if

  end subroutine set_boundary

end module heat_solve
