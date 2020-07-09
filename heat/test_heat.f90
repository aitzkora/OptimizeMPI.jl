module test_heat
    use heat_solve
    use asserts
    use mat_utils
    implicit none

contains

    subroutine test_set_boundary
      integer(c_int32_t) :: i, n_loc, ierr, proc, n
      integer(c_int32_t) :: rank_w, size_w, rank_2D, comm2D
      integer(c_int32_t), parameter :: ndims = 2, North = 1, South = 2, East = 3, West = 4
      logical :: is_master, reorder = .true.
      integer(c_int32_t) :: neighbour(4), coords(2)
      real(c_double) , allocatable :: boundary(:), u(:,:), u_check(:,:)

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)
      proc = 3
      n = 6
      n_loc = 2 + n / proc
      allocate( u(n_loc, n_loc) )

      boundary = [(1.d0 * i, i =1,4*n+4) ]

      call assert_equals(size_w, 9)
      call MPI_CART_CREATE( MPI_COMM_WORLD, 2, [proc,proc], [.false.,.false.], .true., comm2D, ierr )
      call MPI_COMM_RANK( comm2D, rank_2D, ierr )

      !! Fetch the processus coordinates in the 2-D grid
      call MPI_CART_COORDS( comm2D, rank_2D, ndims, coords, ierr )

      call set_boundary(coords, proc, u, boundary)

      select case(rank_w) 
        case(0) 
            u_check = 1.d0 * reshape([[1 ,2,3,4],&
                                      [9 ,0,0,0],&
                                      [11,0,0,0],&
                                      [13,0,0,0]],[4,4])
        case(1) 
            u_check = 1.d0 * reshape([[3 ,4,5,6],&
                                      [0 ,0,0,0],&
                                      [00,0,0,0],&
                                      [00,0,0,0]],[4,4])
        case(2) 
            u_check = 1.d0 * reshape([[5,6,7, 8],&
                                      [0,0,0,10],&
                                      [0,0,0,12],&
                                      [0,0,0,14]],[4,4])
        case(3) 
            u_check = 1.d0 * reshape([[11,0,0,0],&
                                      [13,0,0,0],&
                                      [15,0,0,0],&
                                      [17,0,0,0]],[4,4])
        case(4) 
            u_check = reshape(spread(0.d0, 1, 16), [4,4])
        case(5) 
            u_check = 1.d0 * reshape([[0,0,0,12],&
                                      [0,0,0,14],&
                                      [0,0,0,16],&
                                      [0,0,0,18]],[4,4])
        case(6) 
            u_check = 1.d0 * reshape([[15, 0, 0, 0],&
                                      [17, 0, 0, 0],&
                                      [19, 0, 0, 0],&
                                      [21,22,23,24]],[4,4])
        case(7) 
            u_check = 1.d0 * reshape([[0, 0, 0, 0],&
                                      [0, 0, 0, 0],&
                                      [0, 0, 0, 0],&
                                      [23,24,25,26]],[4,4])
        case(8) 
            u_check = 1.d0 * reshape([[0, 0, 0, 16],&
                                      [0, 0, 0, 18],&
                                      [0, 0, 0, 20],&
                                      [25,26,27,28]],[4,4])
 
        case default
          stop "bad number proc"
      end select
      u_check = transpose(u_check)
      call print_mat( u )
      call assert_equals( u, u_check)
      deallocate (u)
      call MPI_FINALIZE( ierr )

    end subroutine test_set_boundary

    subroutine test_compute_error()
      real(c_double) :: h, error, t_final
      real(c_double) , allocatable :: boundary(:), u_target(:,:), x(:), xx(:,:), yy(:,:)
      integer(c_int32_t) :: rank_w, size_w, i, proc, n, ierr

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)
      
      t_final = 1.d0
      n = 6
      proc = 3
      h = 1.d0 / (n+1) 
      x = [ (h*i, i=1, n) ]
      xx = spread(x, 2, n)
      yy = reshape(spread( x, 1, n), [n,n])
      u_target = exp( xx ) * sin( yy )
      boundary = [ (1.d0*i/(4*n+4),i=1,4*n+4)]
      call compute_error(n, proc, t_final, boundary, u_target, error) 
      if (rank_w == 0) then
          call assert_equals_d(7.0749927571023630d0, error)
      end if
      call MPI_FINALIZE( ierr )
    end subroutine test_compute_error

end module test_heat
