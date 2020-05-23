module mod_soild_util
  use mod_monolis

  integer(kint), parameter :: ndof = 3
  logical, save :: isNLGeom = .false.
  logical, save :: isdebug = .true.

  type gaussdef
    real(kdouble) :: strain(6) = 0.0d0
    real(kdouble) :: stress(6) = 0.0d0
  end type gaussdef

  type matdef
    integer(kint) :: n, ndof
    integer(kint), pointer :: item(:)
    integer(kint), pointer :: index(:)
    real(kdouble), pointer :: A(:)
    real(kdouble), pointer :: B(:)
    real(kdouble), pointer :: X(:)
  end type matdef

  type meshdef
    integer(kint) :: nnode
    integer(kint) :: nelem
    integer(kint), pointer :: global_nid(:)
    integer(kint), pointer :: global_eid(:)
    real(kdouble), pointer :: node(:,:)
    integer(kint), pointer :: elem(:,:)
  end type meshdef

  type paramdef
    !> for time step loop
    integer(kint) :: cur_time_step
    integer(kint) :: max_time_step
    real(kdouble) :: delta_t

    !> for NR loop
    integer(kint) :: cur_nrstep
    integer(kint) :: max_nrstep

    !> for boundary condition
    integer(kint) :: nbound
    integer(kint), pointer :: ibound(:,:)
    real(kdouble), pointer :: bound(:)
    integer(kint), pointer :: is_bound(:)

    integer(kint) :: ncload
    integer(kint), pointer :: icload(:,:)
    real(kdouble), pointer :: cload(:)

    !> for material property
    real(kdouble) :: E, mu, rho
  end type paramdef

  type varmdef
    !> for analysis
    real(kdouble), pointer :: x(:)  !> solution vector of Ax = b
    real(kdouble), pointer :: x_raw(:)  !> solution vector of Ax = b
    real(kdouble), pointer :: u(:)  !> displacement
    real(kdouble), pointer :: du(:) !> delta displacement
    real(kdouble), pointer :: q(:)  !> internal force
    real(kdouble), pointer :: f(:)  !> external force
    real(kdouble), pointer :: f_reaction(:) !> reaction force
    real(kdouble), pointer :: f_reaction_raw(:) !> reaction force
    real(kdouble), pointer :: traction(:) !> reaction force
    real(kdouble), pointer :: g(:)  !> vector for Newmark-Beta
    real(kdouble), pointer :: a(:)  !> acceleration
    real(kdouble), pointer :: a_prev(:)  !> previous acceleration
    real(kdouble), pointer :: v(:)  !> velosity
    real(kdouble), pointer :: v_prev(:)  !> previous velosity

    !> for results
    type(gaussdef), pointer :: gauss(:,:)
    real(kdouble), pointer :: nstrain(:,:)
    real(kdouble), pointer :: nstress(:,:)
    real(kdouble), pointer :: nmises(:)
    real(kdouble), pointer :: estrain(:,:)
    real(kdouble), pointer :: estress(:,:)
    real(kdouble), pointer :: emises(:)
  end type varmdef

  type(monolis_structure) :: monolis

contains

  subroutine init_mesh(mesh)
    implicit none
    type(meshdef) :: mesh

    !allocate(mesh%gauss(8, mesh%nelem))
    !allocate(mesh%nstrain(6, mesh%nnode), source = 0.0d0)
    !allocate(mesh%nstress(6, mesh%nnode), source = 0.0d0)
    !allocate(mesh%nmises (mesh%nnode), source = 0.0d0)
    !allocate(mesh%estrain(6, mesh%nelem), source = 0.0d0)
    !allocate(mesh%estress(6, mesh%nelem), source = 0.0d0)
    !allocate(mesh%emises (mesh%nelem), source = 0.0d0)
    !allocate(mesh%u (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%du(3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%q (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%f (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%f_reaction (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%f_reaction_raw (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%traction (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%g (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%a (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%a_prev (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%v (3*mesh%nnode), source = 0.0d0)
    !allocate(mesh%v_prev (3*mesh%nnode), source = 0.0d0)
  end subroutine init_mesh

  subroutine init_matrix_soild(mesh, mat)
    use iso_c_binding
    implicit none
    type(meshdef) :: mesh
    type(matdef) :: mat
    integer(kint) :: i, j, nnode, nz, jS, jE
    integer(c_int), pointer :: index(:), item(:)

    call monolis_get_mesh_to_nodal(mesh%nnode, mesh%nelem, 8, mesh%elem, index, item)

    nnode = mesh%nnode
    mat%n = nnode
    mat%ndof = ndof
    allocate(mat%x(ndof*nnode), source = 0.0d0)
    allocate(mat%b(ndof*nnode), source = 0.0d0)
    allocate(mat%index(0:nnode), source = 0)
    do i = 1, nnode
      mat%index(i) = index(i+1) + i
    enddo

    nz = mat%index(nnode)
    allocate(mat%A(ndof*ndof*nz), source = 0.0d0)
    allocate(mat%item(nz), source = 0)
    do i = 1, nnode
      jS = mat%index(i-1) + 1
      jE = mat%index(i)
      mat%item(jS) = i
      do j = jS+1, jE
        mat%item(j) = item(j-i) + 1
      enddo
      call monolis_qsort_int(mat%item(jS:jE), 1, jE - jS + 1)
    enddo
  end subroutine init_matrix_soild

  subroutine finalize_mesh(mesh)
    implicit none
    type(meshdef) :: mesh

    !if(associated(mesh%node)) deallocate(mesh%node)
    !if(associated(mesh%elem)) deallocate(mesh%elem)
    !if(associated(mesh%ibound)) deallocate(mesh%ibound)
    !if(associated(mesh%bound)) deallocate(mesh%bound)
    !if(associated(mesh%is_bound)) deallocate(mesh%is_bound)
    !if(associated(mesh%icload)) deallocate(mesh%icload)
    !if(associated(mesh%cload)) deallocate(mesh%cload)
    !if(associated(mesh%gauss)) deallocate(mesh%gauss)
    !if(associated(mesh%nstrain)) deallocate(mesh%nstrain)
    !if(associated(mesh%nstress)) deallocate(mesh%nstress)
    !if(associated(mesh%nmises)) deallocate(mesh%nmises)
    !if(associated(mesh%estrain)) deallocate(mesh%estrain)
    !if(associated(mesh%estress)) deallocate(mesh%estress)
    !if(associated(mesh%emises)) deallocate(mesh%emises)
    !if(associated(mesh%u)) deallocate(mesh%u)
    !if(associated(mesh%du)) deallocate(mesh%du)
    !if(associated(mesh%q)) deallocate(mesh%q)
    !if(associated(mesh%f)) deallocate(mesh%f)
    !if(associated(mesh%g)) deallocate(mesh%g)
    !if(associated(mesh%a)) deallocate(mesh%a)
    !if(associated(mesh%a_prev)) deallocate(mesh%a_prev)
    !if(associated(mesh%v)) deallocate(mesh%v)
    !if(associated(mesh%v_prev)) deallocate(mesh%v_prev)
  end subroutine finalize_mesh

end module mod_soild_util