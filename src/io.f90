module mod_soild_io
  use mod_soild_util
  use mod_soild_debug
contains

  subroutine soild_input_param(param)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i

    open(10, file="input.txt", status='old')
      read(10,*) param%n_get_eigen
      read(10,*) param%thresh
      read(10,*) param%E
      read(10,*) param%mu
      read(10,*) param%rho
      read(10,*) param%elem_type
    close(10)

    if(param%elem_type < 1 .or. 2 < param%elem_type) param%elem_type = 0
  end subroutine soild_input_param

  subroutine soild_input_mesh(mesh, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(kint) :: i, in, ndof
    integer(kint), allocatable :: nid(:), perm(:)
    character :: cnum*5, fname*100

    call soild_debug_header("soild_input_mesh")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.txt")
    call monolis_input_node(fname, mesh%nnode, mesh%node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.txt")
    call monolis_input_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "bc.txt")
    call monolis_input_bc_R(fname, param%nbound, ndof, param%ibound, param%bound)

    call soild_debug_int("nnode", mesh%nnode)
    call soild_debug_int("nnode", mesh%nelem)
    call soild_debug_int("nbound", param%nbound)
  end subroutine soild_input_mesh

  subroutine outout_res(mesh, param, var)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(kint) :: i, id, nnode, nelem
    real(kdouble) :: angle, freq
    real(kdouble), parameter :: PI = 3.14159265359d0
    character :: cstep*5, cnum*5, output_dir*100

    call soild_debug_header("outout_res")

    nnode = mesh%nnode
    nelem = mesh%nelem

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    call convert_to_real(mesh, param, var)

    open(20, file=trim(output_dir)//'result.eigen_value.txt', status='replace')
    write(20,"(a)")" Mode No    Freq. [Hz]"
    write(*,"(a)")" Mode No    Freq. [Hz]"
    do id = 1, param%n_get_eigen
      angle = dsqrt(var%val(id))
      freq  = angle*0.5d0/PI
      write(20,"(i8,1pe14.5)") id, freq
      write(*,"(i8,1pe14.5)") id, freq
    enddo
    close(20)

    do id = 1, param%n_get_eigen
      write(cstep,"(i5.5)")id

      !if(monolis%COM%myrank == 0)then
      !open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.pvtu', status='replace')
      !  write(20,"(a)")'<?xml version="1.0"?>'
      !  write(20,"(a)")'<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">'
      !  write(20,"(a)")'<PUnstructuredGrid>'
      !  write(20,"(a)")'<PPoints>'
      !  write(20,"(a)")'<PDataArray type="Float32" NumberOfComponents="3"/>'
      !  write(20,"(a)")'</PPoints>'
      !  write(20,"(a)")'<PCells>'
      !  write(20,"(a)")'<PDataArray type="Int32" Name="connectivity" format="appended"/>'
      !  write(20,"(a)")'<PDataArray type="Int32" Name="offsets" format="appended"/>'
      !  write(20,"(a)")'<PDataArray type="Int32" Name="types" format="appended"/>'
      !  write(20,"(a)")'</PCells>'
      !  write(20,"(a)")'<PPointData>'
      !  write(20,"(a)")'<PDataArray type="Float32" Name="disp" NumberOfComponents="3" format="appended"/>'
      !  write(20,"(a)")'</PPointData>'
      !  write(20,"(a)")'<PCellData>'
      !  write(20,"(a)")'</PCellData>'
      !  do i = 0, monolis%COM%commsize - 1
      !    write(cnum,"(i0)") i
      !    write(20,"(a)")'<Piece Source="./result.'//trim(cstep)//'.'//trim(cnum)//'.vtu"/>'
      !  enddo
      !  write(20,"(a)")'</PUnstructuredGrid>'
      !  write(20,"(a)")'</VTKFile>'
      !close(20)
      !endif

      write(cnum,"(i0)")monolis_mpi_get_global_my_rank()

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