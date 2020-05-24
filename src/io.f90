module mod_soild_io
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine soild_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i

    open(10, file="input.dat", status='old')
      read(10,*) i
      if(i == 1) isNLGeom = .true.
      read(10,*) param%max_nrstep
      read(10,*) param%E
      read(10,*) param%mu
      read(10,*) param%rho
    close(10)
  end subroutine soild_input_param

  subroutine soild_input_mesh(mesh, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(kint) :: i, in
    integer(kint), allocatable :: nid(:), perm(:)
    character :: cnum*5, fname*100

    call soild_debug_header("soild_input_mesh")

    if(comm_size > 1)then
      call modify_finename("node", fname)
      call monolis_input_mesh_node(fname, mesh%nnode_in, mesh%nnode, mesh%node)

      call modify_finename("elem", fname)
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)

      call modify_finename("node.id", fname)
      call monolis_input_id(fname, mesh%nid)

      call modify_finename("elem.id", fname)
      call monolis_input_id(fname, mesh%eid)
    else
      call modify_finename("node", fname)
      call monolis_input_mesh_node(fname, mesh%nnode_in, mesh%nnode, mesh%node, mesh%nid)

      call modify_finename("elem", fname)
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem, mesh%eid)
    endif

    fname = "bc.dat"
    call monolis_input_condition(fname, param%nbound, param%ibound, param%bound)

    fname = "load.dat"
    call monolis_input_condition(fname, param%ncload, param%icload, param%cload)

    call global_to_local(mesh%nnode, mesh%nid, mesh%nelem, mesh%elem, mesh%nbase_func, &
      param%nbound, param%ibound, param%ncload, param%icload)

    call soild_debug_int("nnode", mesh%nnode)
    call soild_debug_int("nnode", mesh%nelem)
    call soild_debug_int("nbound", param%nbound)
    call soild_debug_int("ncload", param%ncload)
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

  subroutine global_to_local(nnode, nid, nelem, e, nenode, nb, b, nc, c)
    implicit none
    integer(kint) :: i, in, j, id, nenode
    integer(kint) :: imax, imin, nb, nc
    integer(kint) :: nnode, nid(:)
    integer(kint) :: nelem, e(:,:), b(:,:), c(:,:)
    integer(kint), allocatable :: perm(:)

    allocate(perm(nnode), source = 0)
    do i = 1, nnode
      perm(i) = i
    enddo
    call monolis_qsort_int_with_perm(nid, 1, nnode, perm)

    do i = 1, nelem
      do j = 1, nenode
        in = e(j,i)
        call monolis_bsearch_int(nid, 1, nnode, in, id)
        if(id == -1)then
          e(j,i) = -1
        else
          e(j,i) = perm(id)
        endif
      enddo
    enddo

    do i = 1, nb
      in = b(1,i)
      call monolis_bsearch_int(nid, 1, nnode, in, id)
      if(id == -1)then
        b(1,i) = -1
      else
        b(1,i) = perm(id)
      endif
    enddo

    do i = 1, nc
      in = c(1,i)
      call monolis_bsearch_int(nid, 1, nnode, in, id)
      if(id == -1)then
        c(1,i) = -1
      else
        c(1,i) = perm(id)
      endif
    enddo
  end subroutine global_to_local

  subroutine outout_res(mesh, param, var)
    implicit none
    type(paramdef) :: param
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, id, nnode, nelem
    character :: cstep*5, cnum*5, output_dir*100

    call soild_debug_header("outout_res")

    nnode = mesh%nnode
    nelem = mesh%nelem

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    call convert_to_real(mesh, var)

    write(cstep,"(i5.5)")param%cur_time_step

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
      write(20,"(a)")'<PDataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nstress" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nmises" NumberOfComponents="1" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="appended"/>'
      write(20,"(a)")'</PPointData>'
      write(20,"(a)")'<PCellData>'
      write(20,"(a)")'<PDataArray type="Float32" Name="estrain" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="estress" NumberOfComponents="6" format="appended"/>'
      write(20,"(a)")'<PDataArray type="Float32" Name="emises" NumberOfComponents="1" format="appended"/>'
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
        write(20,"(1p3e12.4)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nstrain(1,i), var%nstrain(2,i), var%nstrain(3,i), &
                            & var%nstrain(4,i), var%nstrain(5,i), var%nstrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nstress" NumberOfComponents="6" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nstress(1,i), var%nstress(2,i), var%nstress(3,i), &
                            & var%nstress(4,i), var%nstress(5,i), var%nstress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nmises" NumberOfComponents="1" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%nmises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="ascii">'
      do i = 1, nnode
        write(20,"(1p6e12.4)")var%f_reaction(3*i-2), var%f_reaction(3*i-1), var%f_reaction(3*i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</PointData>'

      write(20,"(a)")'<CellData>'
      write(20,"(a)")'<DataArray type="Float32" Name="estrain" NumberOfComponents="6" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%estrain(1,i), var%estrain(2,i), var%estrain(3,i), &
                            & var%estrain(4,i), var%estrain(5,i), var%estrain(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="estress" NumberOfComponents="6" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%estress(1,i), var%estress(2,i), var%estress(3,i), &
                            & var%estress(4,i), var%estress(5,i), var%estress(6,i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'<DataArray type="Float32" Name="emises" NumberOfComponents="1" format="ascii">'
      do i = 1, nelem
        write(20,"(1p6e12.4)")var%emises(i)
      enddo
      write(20,"(a)")'</DataArray>'
      write(20,"(a)")'</CellData>'

      write(20,"(a)")'</Piece>'
      write(20,"(a)")'</UnstructuredGrid>'
      write(20,"(a)")'</VTKFile>'
    close(20)
  end subroutine outout_res

  subroutine convert_to_real(mesh, var)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    integer(kint) :: i, j, nnode, nelem
    real(kind=kdouble) :: thr

    nnode = mesh%nnode
    nelem = mesh%nelem
    thr = 1.0d-30

    do i = 1, nnode
      if(var%nmises(i) < thr) var%nmises(i) = 0.0d0
    enddo

    do i = 1, nnode
      do j = 1, 3
        if(dabs(var%u    (3*i-3+j)) < thr) var%u    (3*i-3+j) = 0.0d0
        if(dabs(var%f_reaction(3*i-3+j)) < thr) var%f_reaction(3*i-3+j) = 0.0d0
      enddo
    enddo

    do i = 1, nnode
      do j = 1, 6
        if(dabs(var%nstrain(j,i)) < thr) var%nstrain(j,i) = 0.0d0
        if(dabs(var%nstress(j,i)) < thr) var%nstress(j,i) = 0.0d0
      enddo
    enddo

    do i = 1, nelem
      if(dabs(var%emises(i)) < thr) var%emises(i) = 0.0d0
    enddo

    do i = 1, nelem
      do j = 1, 6
        if(dabs(var%estrain(j,i)) < thr) var%estrain(j,i) = 0.0d0
        if(dabs(var%estress(j,i)) < thr) var%estress(j,i) = 0.0d0
      enddo
    enddo
  end subroutine convert_to_real
end module mod_soild_io