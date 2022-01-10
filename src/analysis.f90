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

    call output_mm(var%mass)

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

  subroutine output_mm(mass)
    integer(4) :: n, nz, i, j, jS, jE, in, k1, k2
    real(8) :: val, mass(:)

    n  = monolis%MAT%N*monolis%MAT%NDOF
    nz = monolis%MAT%NDOF*monolis%MAT%NDOF*monolis%MAT%NZ
    open(20, file = "A.mtx", status = "replace")
      write(20,"(a)")"%%MatrixMarket matrix coordinate real symmetric"
      write(20,*) n, n, nz

      do i = 1, monolis%MAT%N
        jS = monolis%MAT%index(i-1) + 1
        jE = monolis%MAT%index(i)
        do j = jS, jE
          in = monolis%MAT%item(j)
          do k1 = 1, 3
          do k2 = 1, 3
            val = monolis%MAT%A(9*(j-1) + 3*(k1-1) + k2)
            write(20,"(i0,x,i0,x,1pe22.14)") 3*(i-1)+k1, 3*(in-1)+k2, val
          enddo
          enddo
        enddo
      enddo
    close(20)

    open(20, file = "B.dat", status = "replace")
      write(20,*) n
      do i = 1, n
        write(20,"(1pe22.14)") mass(i)
      enddo
    close(20)
  end subroutine output_mm

end module mod_soild_analysis
