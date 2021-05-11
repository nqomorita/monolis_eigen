module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug

contains

  subroutine solver(mesh, param, var)
    use mod_monolis_util
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var

    call soild_debug_header("solver")

    call monolis_param_set_method(monolis, monolis_iter_CG)
    call monolis_param_set_precond(monolis, monolis_prec_SOR)
    call monolis_param_set_maxiter(monolis, 100000)
    call monolis_param_set_tol(monolis, 1.0d-8)
    call monolis_param_set_is_scaling(monolis, .false.)
    call monolis_param_set_is_reordering(monolis, .false.)
    call monolis_param_set_is_debug(monolis, .false.)
    call monolis_param_set_show_time(monolis, .false.)
    call monolis_param_set_show_iterlog(monolis, .false.)
    call monolis_param_set_show_summary(monolis, .false.)

    !call monolis_eigen_inverted_lobpcg(monolis, &
    !  & param%n_get_eigen, param%thresh, 100, var%u)

    call monolis_eigen_inverted_standard_lanczos(monolis, &
      & param%n_get_eigen, param%thresh, 200, var%val, var%vec, var%is_bc)
  end subroutine solver
end module mod_soild_solver
