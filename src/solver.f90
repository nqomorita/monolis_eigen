module mod_soild_solver
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine solver(mesh, var)
    use mod_monolis_util
    use mod_monolis_solve
    use mod_soild_debug
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    call soild_debug_header("solver")

    monolis%PRM%maxiter = 100000
    monolis%PRM%precond = 1
    monolis%PRM%tol = 1.0d-8
    monolis%PRM%is_scaling = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x = .true.
    monolis%PRM%is_debug = .false.
    monolis%PRM%show_summary = .true.
    monolis%PRM%show_time = .true.
    monolis%PRM%show_iterlog = .false.

    call monolis_solve(monolis, var%B, var%X)
    call soild_plot_solver(monolis%PRM%curiter, monolis%PRM%curresid)

    if(monolis%PRM%curresid > monolis%PRM%tol)then
      if(monolis%COM%myrank == 0) write(*,"(a)") "*** ERROR: monolis solver is not converge"
      stop
    endif
  end subroutine solver

  function is_convergence(mesh, var, step)
    use mod_monolis_util
    use mod_monolis_linalg
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: step
    real(kdouble), save :: b0nrm
    real(kdouble) :: bnrm, rnrm, rnrmmax
    logical :: is_convergence

    is_convergence = .false.

    bnrm = 0.0d0
    rnrm = 0.0d0
    call monolis_inner_product_R(monolis%COM, mesh%nnode, ndof, var%B, var%B, bnrm)
    bnrm = bnrm

    if(step == 1)then
      b0nrm = bnrm
      write(*,"(a,1pe12.5)")"  ** NR        b0nrm: ", dsqrt(b0nrm)
    else
      rnrm    = dsqrt(bnrm/b0nrm)
      rnrmmax = dabs(maxval(var%B))
      write(*,"(a,1pe12.5,a,1pe12.5)")"  ** NR     residual: ", rnrm, ", ", rnrmmax
      if(rnrm < 1.0d-6 .or. rnrmmax < 1.0d-8) is_convergence = .true.
    endif
  end function is_convergence

end module mod_soild_solver