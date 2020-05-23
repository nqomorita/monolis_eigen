
program main
  use mod_soild_util
  use mod_soild_io
  use mod_soild_analysis
  use mod_soild_debug
  implicit none
  type(meshdef) :: mesh
  type(paramdef) :: param
  type(vardef) :: var
  real(kdouble) :: t1, t2, t3

  call monolis_initialize(monolis)
  call soild_test_set_myrank(monolis%COM%myrank)
  t1 = monolis_get_time()

  !> FEM part
  call soild_write_header("Solid FEM")
  call soild_input_param(param)
  call soild_input_mesh(mesh, param)

  t2 = monolis_get_time()
  call soild_plot_time("input", t2 - t1)

  call solid_linear_static(mesh, param, var)

  t3 = monolis_get_time()
  call soild_plot_time("total ", t3 - t1)

  call monolis_finalize(monolis)
end program main
