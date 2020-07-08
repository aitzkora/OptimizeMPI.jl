module asserts

  implicit none
  public :: assert_equals

  interface assert_equals

    module procedure :: assert_equals_d
    module procedure :: assert_equals_vec
    module procedure :: assert_equals_integer_vec
    module procedure :: assert_equals_integer
    module procedure :: assert_equals_string
    module procedure :: assert_equals_logical

  end interface


contains

! assert subroutine to test equality between integer vectors
  subroutine assert_equals_integer_vec(x, y)
    integer, intent (in) :: x(:), y(:)

    if (any((x/=y))) then
      print *, abs(x-y)
      stop -1
    end if

  end subroutine

! assert subroutine to test equality between integer
  subroutine assert_equals_integer(x, y)
    integer, intent (in) :: x, y

    if (x/=y) then
      print *, abs(x-y)
      stop -1
    end if

  end subroutine

! assert subroutine to test equality between strings
  subroutine assert_equals_string(x, y)
    character (*), intent (in) :: x, y

    if (x/=y) then
      print *, x // '/=' // y
      stop -1
    end if

  end subroutine

! assert subroutine to test equality between logical
  subroutine assert_equals_logical(x, y)
    logical, intent (in) :: x, y

    if (x .neqv. y) then
      print *, x, '.neqv.', y
      stop -1
    end if

  end subroutine

! assert subroutine to test equality between reals
  subroutine assert_equals_d(x, y, prec)
    real , intent (in) :: x, y
    real , optional :: prec
    real :: eps = epsilon(x)

    if (present(prec)) then
      eps = prec
    end if
    if (abs(x-y)>eps) then
      print *, abs(x-y)
      stop -1
    end if

  end subroutine


  subroutine assert_equals_vec(x, y, prec)
    real , intent (in) :: x(:), y(:), prec

    if (sum(abs(x-y))>prec) then
      print *, abs(x-y)
      stop -1
    end if

  end subroutine

end module
