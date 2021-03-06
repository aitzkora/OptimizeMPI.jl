module test_heat
    use heat_solve
    use asserts
    use mat_utils
    use mpi
    implicit none

contains

    subroutine test_set_boundary
      integer(c_int32_t) :: i, n_loc, ierr, proc, n
      integer(c_int32_t) :: rank_w, size_w, rank_2D, comm2D
      integer(c_int32_t), parameter :: ndims = 2
      integer(c_int32_t) :: coords(ndims)
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
      call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, [proc,proc], [.false.,.false.], .true., comm2D, ierr )
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
          call assert_equals_d(7.0749927571023630d0, error , 100.d0 * epsilon(1.d0))
      end if
      call MPI_FINALIZE( ierr )
    end subroutine test_compute_error

    subroutine test_compute_gradient()
      real(c_double) :: h, error, t_final
      real(c_double) , allocatable :: x(:), xx(:,:), yy(:,:)
      real(c_double) , allocatable :: boundary(:), u_target(:,:), grad(:), grad_check(:)
      integer(c_int32_t) :: rank_w, size_w, i, proc, n, ierr, p

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)
      
      t_final = 1.d0
      n = 6
      proc = 3
      h = 1.d0 / (n+1) 
      p = 4*n+ 4
      x = [ (h*i, i=1, n) ]
      xx = spread(x, 2, n)
      yy = reshape(spread( x, 1, n), [n,n])
      u_target = exp( xx ) * sin( yy )
      boundary = [ (1.d0*i/p,i=1,p)]
      grad = spread(0.d0, 1, p)
      call compute_error(n, proc, t_final, boundary, u_target, error, grad) 
      grad_check = [ 0.0d0, -0.09842815768006086d0, -0.4374761608700516d0, -0.8300004174911652d0, &
                     -1.1406120122679873d0, -1.2351524529267455d0, -0.9449325930916592d0, 0.0d0, -0.09842815768006086d0, &
                     -0.9449325930916592d0, -0.12776638217426417d0, -1.462891498489563d0, -0.10490816071833992d0, & 
                     -1.7652324922933025d0, -0.038871413460823526d0, -1.9083984256932705d0, 0.05531777209912303d0, &
                     -1.841839556884301d0, 0.1271842280918711d0, -1.3711354111606837d0, 0.0d0, 0.1271842280918711d0, &
                     -0.28458767577979044d0, -0.8661527030471672d0, -1.3900462605967137d0, -1.6513420741741534d0, &
                     -1.3711354111606837d0, 0.0d0 ]

      !if (rank_w ==  0) then 
      !    print '(8(1xf9.5)/6((1xf9.5)61xf9.5/)8(1xf9.5))', grad_check-grad
      !end if
      if (rank_w == 0) then
          call assert_equals_d_vec(grad_check, grad , 100.d0 * epsilon(1.d0))
      end if
      call MPI_FINALIZE( ierr )
    end subroutine test_compute_gradient


    subroutine test_f_p()
       
       integer(c_int32_t), parameter :: ndims = 2
       integer(c_int32_t) :: n, i, j, n_loc , proc, ierr, rank_2D, comm2D, coords(2), offset_1, offset_2, rank_w , size_w
       real(c_double), allocatable :: f_p_check(:) , f_p(:), lambda(:,:), lambda_loc(:,:)

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, rank_w, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, size_w, ierr)
       proc = 3
       n = 6
       n_loc = 2 + n / proc
       call MPI_CART_CREATE( MPI_COMM_WORLD, 2, [proc,proc], [.false.,.false.], .true., comm2D, ierr )
       call MPI_COMM_RANK( comm2D, rank_2D, ierr )

       !! Fetch the processus coordinates in the 2-D grid
       call MPI_CART_COORDS( comm2D, rank_2D, ndims, coords, ierr )

       ! computed by the julia code in (Matrix(∂ₚF(n))'*[ 1. * i / 36. for i=1:36])'
       f_p_check = [ 0.d0, 0.006944444444444443d0, 0.048611111111111105d0, 0.09027777777777776d0, 0.13194444444444442d0, &
                    0.17361111111111108d0, 0.21527777777777776d0, 0.d0, 0.006944444444444443d0, 0.21527777777777776d0, &
                    0.013888888888888886d0, 0.22222222222222218d0, 0.02083333333333333d0, 0.22916666666666663d0, & 
                    0.027777777777777773d0, 0.23611111111111108d0,0.03472222222222222d0, 0.24305555555555552d0, &
                    0.04166666666666666d0,0.24999999999999997d0,0.d0, 0.04166666666666666d0,0.08333333333333331d0,& 
                    0.12499999999999999d0, 0.16666666666666663d0,0.20833333333333331d0, 0.24999999999999997d0, 0.d0 ]
     
       allocate(f_p, mold=f_p_check)
       lambda = reshape([ (1.d0 * i / (n * n), i=1, n * n) ], [n, n]) 
       offset_1 = coords(1) * (n_loc - 2)
       offset_2 = coords(2) * (n_loc - 2)
       lambda_loc = lambda(offset_1 + 1 : offset_1 + n_loc-2 , offset_2 + 1 : offset_2 + n_loc-2 )
       call compute_f_p_t(coords, proc, lambda_loc, f_p)

       call assert_equals( f_p, f_p_check, epsilon(1.d0))

       call MPI_FINALIZE( ierr )

    end subroutine test_f_p

end module test_heat
