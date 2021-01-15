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
    !> for time step loop
    integer(kint) :: cur_time_step

    !> for boundary condition
    integer(kint) :: nbound
    integer(kint), allocatable :: ibound(:,:)
    real(kdouble), allocatable :: bound(:)

    !> for material property
    real(kdouble) :: E, mu, rho
  end type paramdef

  type vardef
    !> for analysis
    real(kdouble), allocatable :: x(:)  !> solution vector of Ax = b
    real(kdouble), allocatable :: b(:)  !> solution vector of Ax = b
    real(kdouble), allocatable :: u(:)
  end type vardef

  type(monolis_structure) :: monolis

contains

  subroutine init_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

    allocate(var%u (3*mesh%nnode), source = 0.0d0)
    allocate(var%x (3*mesh%nnode), source = 0.0d0)
    allocate(var%b (3*mesh%nnode), source = 0.0d0)
  end subroutine init_mesh

  subroutine init_matrix(mesh)
    implicit none
    type(meshdef) :: mesh

    call monolis_get_nonzero_pattern(monolis, mesh%nnode, 8, ndof, mesh%nelem, mesh%elem)
    !monolis%MAT%N = mesh%nnode
  end subroutine init_matrix

  subroutine finalize_mesh(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var

  end subroutine finalize_mesh

end module mod_soild_util