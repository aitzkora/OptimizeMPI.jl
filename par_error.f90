module par_error
    use mpi
    use iso_c_binding, only: c_int32_t, c_double
    public:: compute_error
contains

    ! the calling subroutine
    ! compute in a // way : 
    !    - the objective function 0.5 * | x -c |**2
    !    - the gradient x - c  

    subroutine compute_error(n, x, c, f, df) bind(C, name="compute_error")
        implicit none
        integer(c_int32_t), intent(in) :: n  ! size of the array
        real(c_double), intent(in) :: x(n), c(n) !
        real(c_double), intent(out) :: f
        real(c_double), intent(out) :: df(n)

        real(c_double) ::  f_loc
        integer(c_int32_t) :: code
        
        f = 0.d0
        ! computing objective
        f_loc = 0.5d0 * sum( (x - c)**2 )
        call MPI_ALLREDUCE( f_loc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, code )
        
        ! computing gradient and ...
        df = x - c
    end subroutine compute_error

end module par_error
