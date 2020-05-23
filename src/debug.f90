module mod_soild_debug
  use mod_soild_util

  private
  public :: soild_test_set_myrank
  public :: soild_write_header
  public :: soild_debug_header
  public :: soild_debug_int
  public :: soild_debug_real
  public :: soild_debug_char
  public :: soild_debug_logical
  public :: soild_debug_time
  public :: soild_plot_time
  public :: soild_plot_solver

  integer(kint), save :: myrank = 0
  integer(kint), parameter :: flog = 30

contains

  subroutine soild_test_set_myrank(set)
    implicit none
    integer(kint) :: set
    character :: output_dir*100
    myrank = set

    !output_dir = "log/"
    !call system('if [ ! -d log ]; then (echo "** create log"; mkdir -p log); fi')

    !if(myrank == 0)then
    !  open(flog, file = "./log/analysis.log", status = "replace")
    !  write(flog, "(a)")"* s-version FEM analysis"
    !endif
  end subroutine soild_test_set_myrank

  subroutine soild_write_header(header)
    implicit none
    character(*) :: header

    if(myrank == 0) write(*,"(a)")"* "//trim(header)
  end subroutine soild_write_header

  subroutine soild_debug_header(header)
    implicit none
    character(*) :: header

    if(myrank == 0)then
      write(*,"(a)")"** soild debug: "//trim(header)
    endif
  end subroutine soild_debug_header

  subroutine soild_debug_int(header, n)
    implicit none
    integer(kint) :: n
    character(*) :: header

    if(myrank == 0) write(*,"(a,i12)")"** soild debug: "//trim(header)//": ", n
  end subroutine soild_debug_int

  subroutine soild_debug_real(header, r)
    implicit none
    real(kdouble) :: r
    character(*) :: header

    if(myrank == 0) write(*,"(a,1pe12.5)")"** soild debug: "//trim(header)//": ", r
  end subroutine soild_debug_real

  subroutine soild_debug_char(header, char)
    implicit none
    character(*) :: header, char

    if(myrank == 0) write(*,"(a,a)")"** soild debug: "//trim(header)//": ", trim(char)
  end subroutine soild_debug_char

  subroutine soild_debug_logical(header, l)
    implicit none
    logical :: l
    character(*) :: header

    if(myrank == 0) write(*,"(a,l)")"** soild debug: "//trim(header)//": ", l
  end subroutine soild_debug_logical

  subroutine soild_debug_time(step, time)
    implicit none
    integer(kint) :: step
    real(kdouble) :: time
    if(myrank == 0)then
      write(*,"(a,i8,1pe12.5)")"* current time step: ", step, time
      !write(flog,"(a,i8,1pe12.5)")"* current time step: ", step, time
    endif
  end subroutine soild_debug_time

  subroutine soild_plot_time(header, time)
    implicit none
    real(kdouble) :: time
    character(*) :: header

    if(myrank == 0)then
      write(*,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
      !write(flog,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
    endif
  end subroutine soild_plot_time

  subroutine soild_plot_solver(iter, resid)
    implicit none
    integer(kint) :: iter
    real(kdouble) :: resid

    !if(myrank == 0)then
    !  write(flog,"(a,i8)")     "  - monolis converge iter    :", iter
    !  write(flog,"(a,1pe10.3)")"  - monolis converge residual:", resid
    !endif
  end subroutine soild_plot_solver

end module mod_soild_debug