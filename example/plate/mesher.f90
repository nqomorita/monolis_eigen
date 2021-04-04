program main
  implicit none
  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8
  integer(kint) :: i, j, k, in, ierr, iid
  integer(kint) :: mysize, myrank
  integer(kint) :: NNODE, NELEM
  integer(kint) :: NX, NY, NZ
  integer(kint), allocatable :: elemNode(:,:)
  integer(kint) :: EX, EY, EZ, iarg
  integer(kint) :: n_sample, sample_id
  integer(kint) :: i_mesh_thetax, i_mesh_thetay, i_mesh_thetaz, i_mesh_E, i_mesh_nu
  integer(kint) :: i_mesh_part
  integer(kint) :: ip1, ip2, ip3, ip4, ip5, ip6, ip7, ip8, ip9, ip10, ip11, ip12, ip13, ip14
  integer(kint) :: mesh_param(14)
  integer(kint) :: mesh_part(3,15)
  real(kdouble) :: eps
  real(kdouble) :: XX, YY, ZZ
  real(kdouble) :: TX, TY, TZ, LX, LY, LZ
  real(kdouble) :: x1, x2, y1, y2, z1, z2, ds, dc
  real(kdouble) :: mesh_thera(2,3)
  real(kdouble) :: mesh_E(2)
  real(kdouble) :: mesh_nu(2)
  real(kdouble) :: mesh_node(3,8)
  real(kdouble), allocatable :: coord(:,:), tmp_coord(:,:)
  character :: temp*100

  iarg = iargc()
  if(iarg == 6)then
    call getarg(1, temp)
    read(temp,*) LX
    call getarg(2, temp)
    read(temp,*) LY
    call getarg(3, temp)
    read(temp,*) LZ
    call getarg(4, temp)
    read(temp,*) EX
    call getarg(5, temp)
    read(temp,*) EY
    call getarg(6, temp)
    read(temp,*) EZ
  else
    EX = 1
    EY = 1
    EZ = 1
  endif

  NELEM =  EX   * EY   * EZ
  NNODE = (EX+1)*(EY+1)*(EZ+1)
  NX = EX + 1
  NY = EY + 1
  NZ = EZ + 1

  write(*,"(a,i8)")"NNODE   :", NNODE
  write(*,"(a,i8)")"NELEM   :", NELEM
  write(*,"(a,1pe12.5)")"Length X:", LX
  write(*,"(a,1pe12.5)")"Length Y:", LY
  write(*,"(a,1pe12.5)")"Length Z:", LZ

  TX = LX/dble(EX)
  TY = LY/dble(EY)
  TZ = LZ/dble(EZ)

  allocate(elemNode (8,NELEM))
  allocate(coord    (3,NNODE))
  allocate(tmp_coord(3,NNODE))
  elemNode = 0
  elemNode = 0.0d0
  tmp_coord = 0.0d0

  in = 1
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        coord(1,in) = dble(i-1)*TX
        coord(2,in) = dble(j-1)*TY
        coord(3,in) = dble(k-1)*TZ
        in = in + 1
      enddo
    enddo
  enddo

  in = 1
  do k = 1, EZ
    do j = 1, EY
      do i = 1, EX
        elemNode(1,in) = (k-1)*NX*NY + (j-1)*NX + i
        elemNode(2,in) = (k-1)*NX*NY + (j-1)*NX + i + 1
        elemNode(3,in) = (k-1)*NX*NY + (j  )*NX + i + 1
        elemNode(4,in) = (k-1)*NX*NY + (j  )*NX + i
        elemNode(5,in) = (k  )*NX*NY + (j-1)*NX + i
        elemNode(6,in) = (k  )*NX*NY + (j-1)*NX + i + 1
        elemNode(7,in) = (k  )*NX*NY + (j  )*NX + i + 1
        elemNode(8,in) = (k  )*NX*NY + (j  )*NX + i
        in = in + 1
      enddo
    enddo
  enddo

  open(30, file="node.dat",status='replace')
    !node
    write(30,"(i0)")NNODE
    do i = 1, NNODE
      write(30,"(1pe12.5,a,1pe12.5,a,1pe12.5)")coord(1,i),",",coord(2,i),",",coord(3,i)
    enddo
  close(30)

  open(30, file="elem.dat",status='replace')
    !element
    write(30,"(i0, i2)")NELEM, 8
    do i = 1, NELEM
      write(30,"(i8,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8)")elemNode(1,i),",",elemNode(2,i),",", &
    & elemNode(3,i),",",elemNode(4,i),",",elemNode(5,i),",",elemNode(6,i),",",elemNode(7,i),",",elemNode(8,i)
    enddo
  close(30)

  open(30, file="bc.dat",status='replace')
    !ngrp
    ip1 = 0
    do k = 1, NZ
      do j = 1, NY
        do i = 1, NX
          if(i == 1)then
            ip1 = ip1 + 1
          endif
        enddo
      enddo
    enddo

    write(30,"(i0,x,i0)")3*ip1, 3
    in = 1
    do k = 1, NZ
      do j = 1, NY
        do i = 1, NX
          if(i == 1)then
            write(30,"(i8,a)") in, ", 1, 0.0"
            write(30,"(i8,a)") in, ", 2, 0.0"
            write(30,"(i8,a)") in, ", 3, 0.0"
          endif

          in = in + 1
        enddo
      enddo
    enddo
  close(30)
end program main
