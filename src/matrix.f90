module mod_soild_matrix
  use mod_soild_util
  use mod_soild_debug
  use mod_soild_c3d8

contains

  subroutine get_stiff_matrix(mesh, var, param)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24), x(3,8)

    call soild_debug_header("get_stiff_matrix")

    do icel = 1, mesh%nelem
      call get_element_node_id(icel, mesh%elem, elem)
      call get_element_node(elem, mesh%node, x)
      call C3D8_stiff(mesh, var, param, icel, x, stiff)
      call monolis_add_matrix_to_sparse_matrix(monolis, 8, elem, stiff)
    enddo
  end subroutine get_stiff_matrix

  subroutine get_mass_matrix(mesh, var, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, icel, in
    integer(kint) :: elem(8)
    real(kdouble) :: mass(24,24), x(3,8)

    call soild_debug_header("get_mass_matrix")

    var%mass = 0.0d0

    do icel = 1,mesh%nelem
      call get_element_node_id(icel, mesh%elem, elem)
      call get_element_node(elem, mesh%node, x)
      call C3D8_mass(mesh, var, param, icel, x, mass)
      do i = 1, 8
        in = elem(i)
        var%mass(3*in-2) = var%mass(3*in-2) + mass(3*i-2,3*i-2)
        var%mass(3*in-1) = var%mass(3*in-1) + mass(3*i-1,3*i-1)
        var%mass(3*in  ) = var%mass(3*in  ) + mass(3*i  ,3*i  )
      enddo
    enddo
  end subroutine get_mass_matrix

  subroutine bound_condition(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, in, dof, nb
    real(kdouble) :: val, tmp
    integer(kint), allocatable :: indexR(:), itemR(:), permA(:)
    real(kdouble), allocatable :: b(:)  !> solution vector of Ax = b

    call soild_debug_header("bound_condition")
    allocate(b(3*mesh%nnode), source = 0.0d0)

    do nb = 1, param%nbound
      in  = param%ibound(1, nb)
      dof = param%ibound(2, nb)
      if(ndof < dof) stop "*** error: 3 < dof"
      call monolis_set_Dirichlet_bc(monolis, b, in, dof, 0.0d0)
      var%is_bc(3*(in-1)+dof) = .true.
    enddo
  end subroutine bound_condition
end module mod_soild_matrix