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

  subroutine is_convergence(mesh, var)
    use mod_monolis_util
    use mod_monolis_linalg
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    real(kdouble), save :: b0nrm
    real(kdouble) :: bnrm, rnrm, rnrmmax

    bnrm = 0.0d0
    call monolis_inner_product_R(monolis%COM, mesh%nnode, ndof, var%B, var%B, bnrm)
    bnrm = dsqrt(bnrm)
    !if(mesh%cur_nrstep == 1)then
    !  b0nrm = bnrm
    !  write(*,"(a,1pe12.5)")"  ** NR        b0nrm: ", b0nrm
    !else
    !  rnrm    = bnrm/b0nrm
    !  rnrmmax = dabs(maxval(mat%B))
    !  write(*,"(a,1pe12.5,a,1pe12.5)")"  ** NR     residual: ", rnrm, ", ", rnrmmax
    !endif
  end subroutine is_convergence

end module mod_soild_solver