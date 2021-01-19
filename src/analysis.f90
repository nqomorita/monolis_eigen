module mod_soild_analysis
  use mod_soild_util
  use mod_soild_io
  use mod_soild_matrix
  use mod_soild_solver
  use mod_soild_debug

contains

  subroutine solid_eigen(mesh, param, var)
    use mod_monolis_util
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    real(kdouble) :: t1, t2, t3, t4, t5

    call soild_write_header("solid_eigen")

    t1 = monolis_get_time_sync()

    call init_mesh(mesh, param, var)
    call init_matrix(mesh)

    t2 = monolis_get_time_sync()
    call soild_plot_time("nonzero detection", t2 - t1)

    call get_stiff_matrix(mesh, var, param)
    call get_mass_matrix(mesh, var, param)
    call monolis_mass_scaling_fw(monolis%PRM, monolis%COM, monolis%MAT, var%mass)

    call bound_condition(mesh, param, var)

    t3 = monolis_get_time_sync()
    call soild_plot_time("matrix generation", t3 - t2)

    call solver(mesh, param, var)

    t4 = monolis_get_time_sync()
    call soild_plot_time("solver", t4 - t3)

    call monolis_mass_scaling_bk(mesh, param, var, var%mass)
    call outout_res(mesh, param, var)
    call finalize_mesh(mesh, var)

    t5 = monolis_get_time_sync()
    call soild_plot_time("output", t5 - t4)
  end subroutine solid_eigen

end module mod_soild_analysis
