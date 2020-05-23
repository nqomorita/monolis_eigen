module mod_soild_matrix
  use mod_soild_util
  use mod_soild_debug
  !use mod_soild_c3d8

contains

  subroutine get_stiff_matrix(mesh, mat)
    implicit none
    type(meshdef) :: mesh
    type(matdef) :: mat
    integer(kint) :: i, icel
    integer(kint) :: elem(8)
    real(kdouble) :: stiff(24,24)

    call soild_debug_header("get_stiff_matrix")

    mat%A = 0.0d0
    do icel = 1,mesh%nelem
      do i = 1, 8
        elem(i) = mesh%elem(i, icel)
      enddo
      !call C3D8_stiff(mesh, icel, elem, stiff)
      call monolis_sparse_matrix_assemble(mat%index, mat%item, mat%A, 8, 3, elem, elem, stiff)
    enddo
  end subroutine get_stiff_matrix

  subroutine load_condition(var, param)
    implicit none
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, in, dof
    real(kdouble) :: val

    call soild_debug_header("load_condition")
    var%f = 0.0d0
    do i = 1, param%ncload
      in  = param%icload(1, i)
      if(in == -1) cycle

      dof = param%icload(2, i)
      val = param%cload(i)
      if(ndof < dof) stop "*** error: 3 < dof"
      var%f(ndof*(in-1) + dof) = val
    enddo
  end subroutine load_condition

  subroutine get_RHS(mesh, var, mat)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(matdef) :: mat
    integer(kint) :: i

    call soild_debug_header("get_RHS")

    do i = 1, mesh%nnode
      mat%B(3*i-2) = var%f(3*i-2) - var%q(3*i-2)
      mat%B(3*i-1) = var%f(3*i-1) - var%q(3*i-1)
      mat%B(3*i  ) = var%f(3*i  ) - var%q(3*i  )
    enddo
  end subroutine get_RHS

  subroutine bound_condition(mesh, param, mat)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(matdef) :: mat
    integer(kint) :: i, in, k, kS, kE, idof, nb
    integer(kint), allocatable :: indexR(:), itemR(:), permA(:)
    real(kdouble) :: val

    call soild_debug_header("bound_condition")

    call monolis_get_CRR_format(mesh%nnode, mat%index, mat%item, indexR, itemR, permA)

    do nb = 1, param%nbound
      in  = param%ibound(1, nb)
      if(in == -1) cycle

      idof = param%ibound(2, nb)
      val = param%bound(nb)

      if(idof < 0 .or. 3 < idof) stop "*** error: 3 < dof"
      if(idof == 0)then
        kS = 1; kE = 3
      else
        kS = idof; kE = idof
      endif

      do k = kS, kE
        call monolis_sparse_matrix_add_bc(mat%index, mat%item, mat%A, mat%B, &
          & indexR, itemR, permA, 3, in, k, val)
      enddo
    enddo
  end subroutine bound_condition
end module mod_soild_matrix