module mod_soild_update
  use mod_soild_util
  use mod_soild_debug
  use mod_soild_c3d8

contains

  subroutine delta_u_update(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i

    call soild_debug_header("delta_u_update")

    do i = 1, mesh%nnode
      var%du(3*i-2) = var%du(3*i-2) + var%X(3*i-2)
      var%du(3*i-1) = var%du(3*i-1) + var%X(3*i-1)
      var%du(3*i  ) = var%du(3*i  ) + var%X(3*i  )
    enddo
  end subroutine delta_u_update

  subroutine u_update(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i

    call soild_debug_header("u_update")

    do i = 1, mesh%nnode
      var%u(3*i-2) = var%u(3*i-2) + var%du(3*i-2)
      var%u(3*i-1) = var%u(3*i-1) + var%du(3*i-1)
      var%u(3*i  ) = var%u(3*i  ) + var%du(3*i  )
    enddo
  end subroutine u_update

  subroutine stress_update(mesh, var, param)
    use mod_soild_matrix
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, j, in, icel
    real(kdouble) :: func(8,8), inv(8,8), tmp
    real(kdouble) :: nstrain(8,6), nstress(8,6)
    real(kdouble) :: estrain(6),   estress(6)
    real(kdouble) :: q(24), r(3)
    integer(kint), allocatable :: inode(:)

    call soild_debug_header("stress_update")

    do i = 1, mesh%nnode
      do j = 1, 6
        var%nstrain(j,i) = 0.0d0
        var%nstress(j,i) = 0.0d0
      enddo
    enddo

    allocate(inode(mesh%nnode), source = 0)

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func(i,:))
    enddo
    call monolis_get_inverse_matrix(8, func, inv)

    var%q = 0.0d0

    do icel = 1, mesh%nelem
      call C3D8_update(mesh, var, param, icel, q)
      call C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)

      do i = 1, 8
        in = mesh%elem(i,icel)
        inode(in) = inode(in) + 1
        do j = 1, 6
          var%nstrain(j,in) = var%nstrain(j,in) + nstrain(i,j)
          var%nstress(j,in) = var%nstress(j,in) + nstress(i,j)
        enddo
        var%q(3*in-2) = var%q(3*in-2) + q(3*i-2)
        var%q(3*in-1) = var%q(3*in-1) + q(3*i-1)
        var%q(3*in  ) = var%q(3*in  ) + q(3*i  )
      enddo

      do j = 1, 6
        var%estrain(j,icel) = estrain(j)
        var%estress(j,icel) = estress(j)
      enddo
    enddo

    do i = 1, mesh%nnode
      tmp = 1.0d0/dble(inode(i))
      do j = 1, 6
        var%nstrain(j,i) = var%nstrain(j,i) * tmp
        var%nstress(j,i) = var%nstress(j,i) * tmp
      enddo
    enddo

    do i = 1, mesh%nnode
      call get_mises(var%nstress(1:6,i), var%nmises(i))
    enddo
    do i = 1, mesh%nelem
      call get_mises(var%estress(1:6,i), var%emises(i))
    enddo
  end subroutine stress_update

end module mod_soild_update