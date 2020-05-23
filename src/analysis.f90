module mod_soild_analysis
  use mod_soild_util
  use mod_soild_io
  use mod_soild_debug

contains

  subroutine solid_linear_static(mesh, param, var)
    use mod_monolis_util
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    type(matdef) :: mat
    real(kdouble) :: t1, t2, t3, t4, t5, t6

    call soild_write_header("solid_linear_static")

    param%cur_time_step  = 0
    param%cur_nrstep = 1
    t1 = monolis_get_time_sync()

    call soild_debug_time(param%cur_time_step, 0.0d0)
    call init_mesh(mesh)
    call init_matrix(mesh, mat)

    t2 = monolis_get_time_sync()
    call soild_plot_time("nonzero detection", t2 - t1)

    !call get_stiff_matrix(mesh, mat)
    !call load_condition(mesh)
    !call get_RHS(mesh, mat)
    !call bound_condition(mesh, meshL, mat)

    t3 = monolis_get_time_sync()
    call soild_plot_time("matrix generation", t3 - t2)

    !call solver(mesh, mat)

    t4 = monolis_get_time_sync()
    call soild_plot_time("solver", t4 - t3)

    !call get_reaction_force(mesh, meshL)
    !call stress_update(mesh, meshL)
    !call delta_u_update(mesh)
    !call u_update(mesh)

    t5 = monolis_get_time_sync()
    call soild_plot_time("stress calculation", t5 - t4)

    !call outout_res(mesh)
    !call finalize_mesh(mesh)

    t6 = monolis_get_time_sync()
    call soild_plot_time("output", t6 - t5)
  end subroutine solid_linear_static

end module mod_soild_analysis
