module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug

contains

  subroutine solver(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var

    call soild_debug_header("solver")

    call monolis_set_method (monolis, monolis_iter_CG)
    !call monolis_set_precond(monolis, monolis_prec_DIAG)
    call monolis_set_precond(monolis, monolis_prec_MUMPS)
    call monolis_set_maxiter(monolis, 100000)
    call monolis_set_tolerance(monolis, 1.0d-8)
    call monolis_show_timelog (monolis, .false.)
    call monolis_show_iterlog (monolis, .false.)
    call monolis_show_summary (monolis, .true.)

    !call monolis_eigen_standard_lobpcg(monolis, &
    !  & param%n_get_eigen, param%thresh, 10, var%val, var%vec, var%is_bc)

    call monolis_eigen_inverted_standard_lanczos_R(monolis, monoCOM, &
      & param%n_get_eigen, param%thresh, 200, var%val, var%vec, var%is_bc)

  end subroutine solver
end module mod_soild_solver
