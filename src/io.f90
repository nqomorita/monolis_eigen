module mod_soild_io
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine soild_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i

    open(10, file="input.dat", status='old')
      read(10,*) param%n_get_eigen
      read(10,*) param%thresh
      read(10,*) param%E
      read(10,*) param%mu
      read(10,*) param%rho
    close(10)
  end subroutine soild_input_param

  subroutine soild_input_mesh(mesh, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(kint) :: i, in, ndof
    integer(kint), allocatable :: nid(:), perm(:)
    character :: cnum*5, fname*100

    call soild_debug_header("soild_input_mesh")

    if(comm_size > 1)then
      call modify_finename("node", fname)
      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node)

      call modify_finename("elem", fname)
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)

      call modify_finename("node.id", fname)
      call monolis_input_id(fname, mesh%nid)

      call modify_finename("elem.id", fname)
      call monolis_input_id(fname, mesh%eid)
    else
      call modify_finename("node", fname)
      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node, mesh%nid)

      call modify_finename("elem", fname)
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem, mesh%eid)

      call global_to_local_elem(mesh%nnode, mesh%nid, mesh%nelem, mesh%elem, mesh%nbase_func)
    endif

    fname = "bc.dat"
    call monolis_input_condition(fname, param%nbound, ndof, param%ibound, param%bound)

    call global_to_local_conditoin(mesh%nnode, mesh%nid, param%nbound, param%ibound)

    call soild_debug_int("nnode", mesh%nnode)
    call soild_debug_int("nnode", mesh%nelem)
    call soild_debug_int("nbound", param%nbound)
  end subroutine soild_input_mesh

  subroutine modify_finename(fname_in, fname)
    implicit none
    character(*) :: fname_in
    character :: fname*100, cnum*5, output_dir*8

    if(comm_size > 1)then
       output_dir = "parted/"
       write(cnum,"(i0)") myrank
      fname = trim(output_dir)//trim(fname_in)//"."//trim(cnum)
    else
      fname = trim(fname_in)//".dat"
    endif
  end subroutine modify_finename

  subroutine global_to_local_elem(nnode, nid, nelem, e, nenode)
    implicit none
    integer(kint) :: i, in, j, id, nenode
    integer(kint) :: nnode, nid(:)
    integer(kint) :: nelem, e(:,:)
    integer(kint), allocatable :: perm(:), temp(:)

    allocate(temp(nnode), source = 0)
    allocate(perm(nnode), source = 0)
    do i = 1, nnode
      temp(i) = nid(i)
      perm(i) = i
    enddo
    call monolis_qsort_int_with_perm(temp, 1, nnode, perm)

    do i = 1, nelem
      do j = 1, nenode
        in = e(j,i)
        call monolis_bsearch_int(temp, 1, nnode, in, id)
        if(id == -1)then
          e(j,i) = -1
        else
          e(j,i) = perm(id)
        endif
      enddo
    enddo
  end subroutine global_to_local_elem

  subroutine global_to_local_conditoin(nnode, nid, nb, b)
    implicit none
    integer(kint) :: i, in, j, id
    integer(kint) :: imax, imin, nb
    integer(kint) :: nnode, nid(:)
    integer(kint) :: b(:,:)
    integer(kint), allocatable :: perm(:)

    allocate(perm(nnode), source = 0)
    do i = 1, nnode
      perm(i) = i
    enddo
    call monolis_qsort_int_with_perm(nid, 1, nnode, perm)

    do i = 1, nb
      in = b(1,i)
      call monolis_bsearch_int(nid, 1, nnode, in, id)
      if(id == -1)then
        b(1,i) = -1
      else
        b(1,i) = perm(id)
      endif
    enddo
  end subroutine global_to_local_conditoin

  subroutine outout_res(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, id, nnode, nelem
    character :: cstep*5, cnum*5, output_dir*100

    call soild_debug_header("outout_res")

    nnode = mesh%nnode
    nelem = mesh%nelem

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    call convert_to_real(mesh, param, var)

    open(20, file=trim(output_dir)//'result.eigen_value.dat', status='replace')
    write(*,"(a)")" Mode No    Freq. [Hz]"
    do id = 1, param%n_get_eigen
      write(20,"(1pe12.5)") var%val(id)
      write(*,"(i8,1pe14.5)") id, var%val(id)
    enddo
    close(20)

    do id = 1, param%n_get_eigen
      write(cstep,"(i5.5)")id

      if(monolis%COM%myrank == 0)then
      open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.pvtu', status='replace')
        write(20,"(a)")'<?xml version="1.0"?>'
        write(20,"(a)")'<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">'
        write(20,"(a)")'<PUnstructuredGrid>'
        write(20,"(a)")'<PPoints>'
        write(20,"(a)")'<PDataArray type="Float32" NumberOfComponents="3"/>'
        write(20,"(a)")'</PPoints>'
        write(20,"(a)")'<PCells>'
        write(20,"(a)")'<PDataArray type="Int32" Name="connectivity" format="appended"/>'
        write(20,"(a)")'<PDataArray type="Int32" Name="offsets" format="appended"/>'
        write(20,"(a)")'<PDataArray type="Int32" Name="types" format="appended"/>'
        write(20,"(a)")'</PCells>'
        write(20,"(a)")'<PPointData>'
        write(20,"(a)")'<PDataArray type="Float32" Name="disp" NumberOfComponents="3" format="appended"/>'
        write(20,"(a)")'</PPointData>'
        write(20,"(a)")'<PCellData>'
        write(20,"(a)")'</PCellData>'
        do i = 0, monolis%COM%commsize - 1
          write(cnum,"(i0)") i
          write(20,"(a)")'<Piece Source="./result.'//trim(cstep)//'.'//trim(cnum)//'.vtu"/>'
        enddo
        write(20,"(a)")'</PUnstructuredGrid>'
        write(20,"(a)")'</VTKFile>'
      close(20)
      endif

      write(cnum,"(i0)")monolis%COM%myrank

      open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.'//trim(cnum)//'.vtu', status='replace')
        write(20,"(a)")'<?xml version="1.0"?>'
        write(20,"(a)")'<VTKFile type="UnstructuredGrid" version="1.0">'
        write(20,"(a)")'<UnstructuredGrid>'
        write(20,"(a,i0,a,i0,a)")'<Piece NumberOfPoints="', nnode, '" NumberOfCells="', nelem, '">'
        write(20,"(a)")'<Points>'
        write(20,"(a)")'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        do i = 1, nnode
          write(20,"(1p3e20.12)")mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
        enddo

        write(20,"(a)")'</DataArray>'
        write(20,"(a)")'</Points>'
        write(20,"(a)")'<Cells>'
        write(20,"(a)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
        do i = 1, nelem
          write(20,"(8i8)")mesh%elem(1,i)-1, mesh%elem(2,i)-1, mesh%elem(3,i)-1, mesh%elem(4,i)-1, &
                           mesh%elem(5,i)-1, mesh%elem(6,i)-1, mesh%elem(7,i)-1, mesh%elem(8,i)-1
        enddo

        write(20,"(a)")'</DataArray>'
        write(20,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, nelem
          write(20,"(x,i0,$)")8*i
        enddo
        write(20,*)""

        write(20,"(a)")'</DataArray>'
        write(20,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
        do i = 1, nelem
          write(20,"(i3,$)")12
        enddo
        write(20,*)""

        write(20,"(a)")'</DataArray>'
        write(20,"(a)")'</Cells>'

        write(20,"(a)")'<PointData>'
        write(20,"(a)")'<DataArray type="Float32" Name="disp" NumberOfComponents="3" format="ascii">'
        do i = 1, nnode
          write(20,"(1p3e12.4)")var%vec(3*i-2,id), var%vec(3*i-1,id), var%vec(3*i,id)
        enddo
        write(20,"(a)")'</DataArray>'
        write(20,"(a)")'</PointData>'

        write(20,"(a)")'<CellData>'
        write(20,"(a)")'</CellData>'

        write(20,"(a)")'</Piece>'
        write(20,"(a)")'</UnstructuredGrid>'
        write(20,"(a)")'</VTKFile>'
      close(20)
    enddo
  end subroutine outout_res

  subroutine convert_to_real(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, j, nnode, nelem
    real(kind=kdouble) :: thr

    nnode = mesh%nnode
    nelem = mesh%nelem
    thr = 1.0d-30

    do i = 1, param%n_get_eigen
      do j = 1, 3*nnode
        if(dabs(var%vec(j,i)) < thr) var%vec(j,i) = 0.0d0
      enddo
    enddo
  end subroutine convert_to_real
end module mod_soild_io