module mat_utils
   use iso_c_binding
   use iso_fortran_env
   public :: print_mat
contains
 subroutine print_mat( u )
   implicit none
   character(len=48) :: variable_format
   real(c_double), intent(in) :: u(:,:)
   write (variable_format, '(ai0ai0a)') "(", size( u, 1 ), "(", size( u, 2 ), "(1xf8.5)/))"
   write (output_unit, variable_format) transpose(u)
 end subroutine print_mat

end module
