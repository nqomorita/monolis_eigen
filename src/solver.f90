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

    monolis%PRM%maxiter = 10000
    monolis%PRM%precond = monolis_prec_MUMPS
    monolis%PRM%tol = 1.0d-8
    monolis%PRM%is_scaling = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x = .true.
    monolis%PRM%is_debug = .false.
    monolis%PRM%show_summary = .true.
    monolis%PRM%show_time = .false.
    monolis%PRM%show_iterlog = .false.

    !call monolis_eigen_inverted_lobpcg(monolis, &
    !  & param%n_get_eigen, param%thresh, 100, var%u)

    call monolis_eigen_inverted_standard_lanczos(monolis, &
      & param%n_get_eigen, param%thresh, 200, var%val, var%vec, var%is_bc)
  end subroutine solver
end module mod_soild_solver
