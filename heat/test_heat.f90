module test_heat
    use heat_solve
    use asserts
    implicit none

contains

    subroutine test_set_boundary
      integer(c_int32_t) :: i, n_loc, ierr, proc, n
      integer(c_int32_t) :: rank_w, size_w, rank_2D, comm2D
      integer(c_int32_t), parameter :: ndims = 2, North = 1, South = 2, East = 3, West = 4
      logical :: is_master, reorder = .true.
      integer(c_int32_t) :: neighbour(4), coords(2)
      real(c_double) , allocatable :: boundary(:), u(:,:)

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)
      proc = 2 
      n = 10
      n_loc = 2 + n / proc
      allocate( u(n_loc, n_loc) )

      boundary = [(1.d0 * i, i =1,44) ]

      call assert_equals(size_w, 4)
      call MPI_CART_CREATE( MPI_COMM_WORLD, 2, [proc,proc], [.false.,.false.], .true., comm2D, ierr )
      call MPI_COMM_RANK( comm2D, rank_2D, ierr )

      !! Fetch the processus coordinates in the 2-D grid
      call MPI_CART_COORDS( comm2D, rank_2D, ndims, coords, ierr )

      call set_boundary(coords, proc, u, boundary)

      deallocate (u)
      call MPI_FINALIZE( ierr )
      

    end subroutine test_set_boundary




end module test_heat
