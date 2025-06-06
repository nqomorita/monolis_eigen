module mod_soild_util
  use mod_monolis

  integer(kint), parameter :: ndof = 3
  integer(kint), save :: comm_size = 1
  integer(kint), save :: myrank = 0
  logical, save :: isdebug = .true.

  type meshdef
    integer(kint) :: nnode
    integer(kint) :: nelem
    integer(kint) :: nbase_func
    integer(kint), allocatable :: nid(:)
    integer(kint), allocatable :: eid(:)
    real(kdouble), allocatable :: node(:,:)
    integer(kint), allocatable :: elem(:,:)
  end type meshdef

  type paramdef
    integer(kint) :: n_get_eigen
    real(kdouble) :: thresh

    integer(kint) :: elem_type

    !> for boundary condition
    integer(kint) :: nbound
    integer(kint), allocatable :: ibound(:,:)
    real(kdouble), allocatable :: bound(:)

    !> for material property
    real(kdouble) :: E, mu, rho
  end type paramdef

  type vardef
    !> for analysis
    real(kdouble), allocatable :: val(:)
    real(kdouble), allocatable :: vec(:,:)
    real(kdouble), allocatable :: mass(:)
    logical, allocatable :: is_bc(:)
  end type vardef

  type(monolis_structure) :: monolis
  type(monolis_com) :: monoCOM

contains

  subroutine init_mesh(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var

    allocate(var%val(param%n_get_eigen), source = 0.0d0)
    allocate(var%vec(3*mesh%nnode, param%n_get_eigen), source = 0.0d0)
    allocate(var%mass(3*mesh%nnode), source = 0.0d0)
    allocate(var%is_bc(3*mesh%nnode), source = .false.)
  end subroutine init_mesh

  subroutine init_matrix(mesh)
    implicit none
    type(meshdef) :: mesh

    call monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, mesh%nnode, 8, ndof, mesh%nelem, mesh%elem)
  end subroutine init_matrix

  subroutine monolis_mass_scaling_fw(monoPRM, monoCOM, monoMAT, mass)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NDOF2
    integer(kint) :: inod
    integer(kint) :: i, j, jS, jE, in, k, l
    real(kdouble) :: mass(:)
    real(kdouble), pointer :: A(:)
    real(kdouble), allocatable :: diag(:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: t1, t2

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    allocate(diag(NDOF*NP), source = 0.0d0)

    do i = 1, 3*NP
      diag(i) = 1.0d0 / dsqrt(dabs(mass(i)))
    enddo

!    call monolis_update_R(monoCOM, NDOF, diag, tcomm)

    do i = 1, NP
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        do k = 1, NDOF
          do l = 1, NDOF
            A(NDOF2*(j-1) + NDOF*(k-1) + l) = &
          & A(NDOF2*(j-1) + NDOF*(k-1) + l)*diag(NDOF*(i-1)+k)*diag(NDOF*(in-1)+l)
          enddo
        enddo
      enddo
    enddo
  end subroutine monolis_mass_scaling_fw

  subroutine monolis_mass_scaling_bk(mesh, param, var, mass)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, j
    real(kdouble) :: mass(:)
    real(kdouble), allocatable :: diag(:)

    allocate(diag(mesh%nnode*ndof), source = 0.0d0)

    do i = 1, mesh%nnode*ndof
      diag(i) = 1.0d0 / dsqrt(dabs(mass(i)))
    enddo

    do i = 1, param%n_get_eigen
      do j = 1, mesh%nnode*ndof
        var%vec(j,i) = var%vec(j,i) * diag(j)
      enddo
    enddo
  end subroutine monolis_mass_scaling_bk

  subroutine finalize_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

  end subroutine finalize_mesh

  subroutine get_element_node_id(eid, elem, elemid)
    implicit none
    integer(kint) :: i, eid, elem(:,:), elemid(:)
    do i = 1, 8
      elemid(i) = elem(i,eid)
    enddo
  end subroutine get_element_node_id

  subroutine get_element_node(elem, node, x)
    implicit none
    integer(kint) :: i, in, j, elem(:)
    real(kdouble) :: node(:,:), x(:,:)

    do i = 1, 8
      in = elem(i)
      do j = 1, ndof
        x(1,i) = node(1,in)
        x(2,i) = node(2,in)
        x(3,i) = node(3,in)
      enddo
    enddo
  end subroutine get_element_node
end module mod_soild_util