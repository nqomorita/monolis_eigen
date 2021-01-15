module mod_soild_c3d8
  use mod_soild_util
contains

  subroutine C3D8_stiff(mesh, var, param, icel, elem, stiff)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x0(3,8), u(3,8), stiff(24,24)
    real(kdouble) :: r(3), wg, det
    real(kdouble) :: B(6,24), D(6,6), dndx(8,3)

    wg     = 1.0d0
    stiff  = 0.0d0

    do i = 1, 8
      in = elem(i)
      x0(1,i) = mesh%node(1,in)
      x0(2,i) = mesh%node(2,in)
      x0(3,i) = mesh%node(3,in)
    enddo

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x0, r, dndx, det)
      call C3D8_Bmat(dndx, B)
      call C3D8_Dmat(param%E, param%mu, D)
      call C3D8_Kmat(D, B, wg, det, stiff)
    enddo
  end subroutine C3D8_stiff

  subroutine C3D8_Bmat(dndx, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: B(6,24), dndx(8,3), dudx(3,3)

    B = 0.0d0
    do i = 1, 8
      i1 = 3*i-2
      i2 = 3*i-1
      i3 = 3*i
      B(1,i1) = dndx(i,1)
      B(2,i2) = dndx(i,2)
      B(3,i3) = dndx(i,3)
      B(4,i1) = dndx(i,2)
      B(4,i2) = dndx(i,1)
      B(5,i2) = dndx(i,3)
      B(5,i3) = dndx(i,2)
      B(6,i1) = dndx(i,3)
      B(6,i3) = dndx(i,1)
    enddo
  end subroutine C3D8_Bmat

  subroutine C3D8_Dmat(E, mu, D)
    implicit none
    real(kdouble) :: D(6,6), E, mu, g

    D = 0.0d0
    g = E / ((1.0d0+mu) * (1.0d0-2.0d0*mu))

    D(1,1) = g*(1.0d0-mu)
    D(1,2) = g*mu
    D(1,3) = g*mu
    D(2,1) = g*mu
    D(2,2) = g*(1.0d0-mu)
    D(2,3) = g*mu
    D(3,1) = g*mu
    D(3,2) = g*mu
    D(3,3) = g*(1.0d0-mu)
    D(4,4) = 0.5d0*g*(1.0d0-2.0d0*mu)
    D(5,5) = 0.5d0*g*(1.0d0-2.0d0*mu)
    D(6,6) = 0.5d0*g*(1.0d0-2.0d0*mu)
  end subroutine C3D8_Dmat

  subroutine C3D8_Kmat(D, B, wg, det, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det

    DB = matmul(D, B)

    do i = 1, 24
      do j = 1, 24
        do k = 1, 6
          stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C3D8_Kmat

end module mod_soild_c3d8
