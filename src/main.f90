
program main
  use mod_soild_util
  implicit none
  type(meshdef) :: mesh
  real(kdouble) :: t1, t2, t3

  call monolis_initialize(monolis)
  !call soild_test_set_myrank(monolis%COM%myrank)
  t1 = monolis_get_time()

  !> FEM part
  !call soild_write_header("Solid FEM")
  !call soild_input_param(mesh)
  !call soild_input_mesh(mesh)

  t2 = monolis_get_time()
  !call soild_plot_time("input", t2 - t1)

  !call soild_linear_static(mesh)

  t3 = monolis_get_time()
  !call soild_plot_time("total ", t3 - t1)

  call monolis_finalize(monolis)
end program main
