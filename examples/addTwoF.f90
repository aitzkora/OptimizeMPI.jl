module example
    use iso_c_binding
contains
  subroutine addTwoF(x) bind(C, name ="addTwoF")
    integer(c_int), intent (inout) :: x
    x = x + 2
  end subroutine 
end module
