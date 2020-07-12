module heat_solve
    use iso_c_binding, only: c_int32_t, c_double
    use iso_fortran_env
    use mat_utils, only : print_mat
    use mpi
    use heat
    use communications
    implicit none
    private 
    public :: compute_error, set_boundary, compute_f_p_t
 
contains

  subroutine compute_error(n, proc, t_final, boundary, u_target, error, gradient) bind( C, name="compute_error" )

    integer(c_int32_t), intent(in) :: n ! size of the global square matrix
    integer(c_int32_t), intent(in) :: proc ! nb of processes (in each dimensions)
    real(c_double), intent(in) :: t_final
    real(c_double), intent(in) :: u_target(1:n, 1:n)
    real(c_double), intent(in) :: boundary(4*n+4)
    real(c_double), intent(out) :: error
    real(c_double), optional, intent(out) :: gradient(4*n+4)

    integer(c_int32_t) :: n_loc, ierr, p
    integer(c_int32_t) :: rank_w, size_w, rank_2D, comm2D, type_row
    integer(c_int32_t), parameter :: ndims = 2, North = 1, South = 2, East = 3, West = 4
    logical :: is_master, reorder = .true.
    integer(c_int32_t), dimension(4) :: neighbour
    real(c_double) :: h, dt, error_loc, t
    real(c_double), allocatable :: u_in(:,:), u_out(:,:), target_loc(:,:), target_ext(:,:)
    real(c_double), allocatable :: lambda(:,:), lambda_tmp(:,:), f_p(:), zeros(:)

    integer(c_int32_t):: dims(ndims) , coords(ndims)
    logical :: periods(ndims) = [.false., .false.], with_gradient

    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)

    is_master = (rank_w == 0)

    n_loc = 2 + n / proc
    p = 4 * n + 4

    h = 1.d0 / (n+1)
    dt = h ** 2 / 4.d0 
   
    allocate(target_loc(n_loc, n_loc))

    allocate(target_ext(n+2,n+2))
    target_ext = 0.d0 
    target_ext(2:n+1,2:n+1) = u_target
 
    if (present(gradient)) then
        with_gradient = .true.
        allocate(lambda(n_loc, n_loc))
        allocate(f_p(p))
        gradient = 0.d0
        zeros = spread(0.d0, 1, p)
    else
        with_gradient = .false.
    end if

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
    u_in = 0.d0
    u_out = 0.d0
    call set_boundary( coords, proc, u_in, boundary )
    u_out = u_in
    t = 0.d0
    ! forward phase
    do 
      if (t > t_final ) exit
      u_out = heat_kernel( u_in )
      u_in = u_in - dt / h**2 * u_out
      call ghosts_swap( comm2D, type_row, neighbour, u_in )
      t = t + dt
    end do

    associate (of_x => coords(1) * (n_loc-2), of_y => coords(2) * (n_loc-2))
      target_loc = target_ext(of_x + 1 : of_x + n_loc, of_y + 1 : of_y + n_loc) 
    end associate
    error = 0.d0 
    error_loc = sum ((u_in(2:n_loc-1,2:n_loc-1) - target_loc(2:n_loc-1,2:n_loc-1))**2)
    call MPI_ALLREDUCE(error_loc, error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if ( with_gradient ) then
      ! backward  phase
      lambda= 2.d0 * (u_in - target_loc)
      call set_boundary( coords, proc, lambda, zeros)
      t = t_final 
      do
        if (t <= 0) exit
        ! ∇f += ∂ₚFᵀ(λ)
        call compute_f_p_t(coords, proc, lambda(2:n_loc-1, 2:n_loc-1), f_p)
        gradient = gradient + f_p  
        ! λₙ₊₁ = ∂ᵤF(λₙ)
        lambda_tmp = heat_kernel( lambda )
        lambda = lambda - dt * (n+1)**2 * lambda_tmp
        call ghosts_swap( comm2D, type_row, neighbour, lambda )
        call set_boundary(coords, proc, lambda, zeros)
        t = t - dt
      end do
      deallocate (f_p)
      deallocate (lambda)
      deallocate (lambda_tmp)
    end if

    deallocate( u_in )
    deallocate( u_out )
     
    call MPI_TYPE_FREE( type_row, ierr )

  end subroutine compute_error

  ! beware, this function claims for putting zeroes on the boundary not inside!
  subroutine set_boundary(coo, proc, u, boundary)

    integer, intent(in) :: proc
    integer, intent(in) :: coo(2)
    real(c_double), intent(out) :: u(:,:)
    real(c_double), intent(in) :: boundary(:)
    integer :: n_loc, offset, n

    n_loc = size(u, 1)
    n  = (size(boundary, 1) - 4 )/ 4

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

  subroutine compute_f_p_t(coo, proc, lambda_loc, f_p)

    integer, intent(in) :: coo(2)
    integer, intent(in) :: proc
    real(c_double), intent(in) :: lambda_loc(:,:)
    real(c_double), intent(out) :: f_p(:)
    real(c_double), allocatable :: f_p_loc(:)
    real(c_double) :: dt, h
    integer :: n_loc, offset, n, ierr, p

    n_loc = size(lambda_loc, 1) + 2
    n = (n_loc - 2) * proc
    p = 4 * n + 4
    allocate (f_p_loc(p))
    f_p = 0.d0  
    f_p_loc = 0.d0

    ! upper line 
    if (coo(1) == 0)  then
      offset = coo(2)*(n_loc-2)
      associate ( f => f_p_loc(offset + 2 : offset + n_loc -1) )
        f = f - lambda_loc(1, :)
      end associate
    end if

    ! lower line 
    if (coo(1) == (proc - 1)) then 
      offset = 3*n+2 + coo(2) * (n_loc-2) 
      associate ( f => f_p_loc(offset + 2 : offset + n_loc -1))
        f = f - lambda_loc(n_loc-2, :)
    end associate
    endif

    ! left column 
    if (coo(2) == 0) then 
      offset = coo(1) * (2 *n_loc - 4) + n + 1 
      associate ( f => f_p_loc(offset + 2 : offset + 2*n_loc-4 : 2) )
        f  = f - lambda_loc(:, 1)
      end associate
    end if
    ! right column 
    if (coo(2) == (proc - 1)) then 
      offset = coo(1) * (2*n_loc - 4) + n + 2
      associate ( f => f_p_loc(offset + 2 : offset + 2*n_loc -4 :2) )
        f = f - lambda_loc( : , n_loc - 2)
      end associate
    end if


    call MPI_ALLREDUCE(f_p_loc, f_p, 4*n+4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    h = 1.d0 / (n + 1)
    dt = h ** 2 / 4.d0
    f_p = -dt * f_p  * (n+1)**2

    deallocate(f_p_loc)

  end subroutine compute_f_p_t


end module heat_solve
