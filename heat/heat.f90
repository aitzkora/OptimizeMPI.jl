module heat

  use iso_c_binding, only: c_int32_t, c_double
  public :: heat_kernel

contains

  subroutine heat_kernel(u_in, u_out)
    implicit none
    real(c_double), dimension(:, :), intent(in) :: u_in
    real(c_double), dimension(:, :), intent(out) :: u_out
    integer(c_int32_t) :: i, j
  
    do j = 2, size( u_in, 2) - 1
      do i = 2, size( u_in, 1) - 1
        u_out(i,j) = 4.d0 * u_in(i,j) - u_in(i - 1, j) - u_in(i + 1, j) & 
                          - u_in(i, j - 1) - u_in(i, j + 1)
      end do
    end do
  
  end subroutine heat_kernel

end module heat
