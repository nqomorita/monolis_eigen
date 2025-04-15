module mod_soild_c3d8
  use mod_soild_util

contains

  subroutine C3D8_stiff(mesh, var, param, icel, x, stiff)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x(3,8), stiff(24,24)
    real(kdouble) :: r(3), wg, det
    real(kdouble) :: B(6,24), D(6,6), dndx(8,3)

    wg    = 1.0d0
    stiff = 0.0d0
    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x, r, dndx, det)
      call C3D8_Bmat(dndx, B)
      call C3D8_Dmat(param%E, param%mu, D)
      call C3D8_Kmat(D, B, wg, det, stiff)
    enddo
  end subroutine C3D8_stiff

  subroutine C3D8IC_stiff(mesh, var, param, icel, x, stiff)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x(3,8), stiff(24,24), deriv(8,3), stmp(33,33)
    real(kdouble) :: r(3), wg, det, xj(3,3), inv(3,3), S(24,9)
    real(kdouble) :: B(6,33), D(6,6), dndx(11,3), K22(9,9), K22_inv(9,9)

    wg    = 1.0d0
    stiff = 0.0d0
    stmp  = 0.0d0

    !> inverse of Jacobian at element center
    r = 0.0d0
    call monolis_C3D8_shapefunc_deriv(r, deriv)
    xj = matmul(x, deriv)
    call C3D8_get_inverse_matrix(xj, inv, det)
    inv = inv*det

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x, r, dndx(1:8,:), det)
      dndx( 9,:) = -2.0d0*r(1)*inv(1,:)/det
      dndx(10,:) = -2.0d0*r(2)*inv(2,:)/det
      dndx(11,:) = -2.0d0*r(3)*inv(3,:)/det

      call C3D8_Bmat_IC(dndx, B)
      call C3D8_Dmat(param%E, param%mu, D)
      call C3D8_Kmat_IC(D, B, wg, det, stmp)
    enddo

    K22 = stmp(25:33,25:33)
    call monolis_get_inverse_matrix_R(9, K22, K22_inv)
    S = matmul(stmp(1:24,25:33), K22_inv)
    stiff = stmp(1:24,1:24) - matmul(S, stmp(25:33,1:24))
  end subroutine C3D8IC_stiff

  subroutine C3D8_get_inverse_matrix(xj, inv, det)
    implicit none
    real(kdouble) :: xj(3,3), inv(3,3), det, detinv

    det = xj(1,1) * xj(2,2) * xj(3,3) &
        + xj(2,1) * xj(3,2) * xj(1,3) &
        + xj(3,1) * xj(1,2) * xj(2,3) &
        - xj(3,1) * xj(2,2) * xj(1,3) &
        - xj(2,1) * xj(1,2) * xj(3,3) &
        - xj(1,1) * xj(3,2) * xj(2,3)

    if(det < 0.0d0) stop "determinant < 0.0"

    detinv = 1.0d0/det
    inv(1,1) = detinv * ( xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))
    inv(1,2) = detinv * (-xj(1,2)*xj(3,3) + xj(3,2)*xj(1,3))
    inv(1,3) = detinv * ( xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
    inv(2,1) = detinv * (-xj(2,1)*xj(3,3) + xj(3,1)*xj(2,3))
    inv(2,2) = detinv * ( xj(1,1)*xj(3,3) - xj(3,1)*xj(1,3))
    inv(2,3) = detinv * (-xj(1,1)*xj(2,3) + xj(2,1)*xj(1,3))
    inv(3,1) = detinv * ( xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2))
    inv(3,2) = detinv * (-xj(1,1)*xj(3,2) + xj(3,1)*xj(1,2))
    inv(3,3) = detinv * ( xj(1,1)*xj(2,2) - xj(2,1)*xj(1,2))
  end subroutine C3D8_get_inverse_matrix

  subroutine C3D8_mass(mesh, var, param, icel, x, mass)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    integer(kint) :: elem(8)
    real(kdouble) :: x(3,8), mass(24,24)
    real(kdouble) :: r(3), func(8), wg, det
    real(kdouble) :: dndx(8,3)

    wg = 1.0d0
    mass = 0.0d0
    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func)
      call monolis_C3D8_get_global_deriv(x, r, dndx, det)
      call C3D8_Mmat(param, func, wg, det, mass)
    enddo

    call get_lumped_mass(8, 3, mass)
  end subroutine C3D8_mass

  subroutine get_lumped_mass(nn, ndof, mass)
    implicit none
    integer(kint) :: i, j, nn, ndof
    real(kdouble) :: lumped(nn*ndof), mass(:,:)
    real(kdouble) :: diag_mass, total_mass

    total_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      do j = 1, nn*ndof, ndof
        total_mass = total_mass + mass(j,i)
      enddo
    enddo

    diag_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      diag_mass = diag_mass + mass(i,i)
    enddo

    lumped = 0.0d0
    diag_mass = 1.0d0/diag_mass
    do i = 1, nn*ndof
      lumped(i) = lumped(i) + mass(i,i)*total_mass*diag_mass
    enddo

    mass = 0.0d0
    do i = 1, nn*ndof
      mass(i,i) = lumped(i)
    enddo
  end subroutine get_lumped_mass

  subroutine C3D8_Mmat(param, func, wg, det, mass)
    implicit none
    type(paramdef) :: param
    integer(kint) :: i, j, k
    real(kdouble) :: mass(24,24), D(3,3), DN(3,24), wg, det, rho
    real(kdouble) :: func(8), N(3,24)

    rho = param%rho

    N = 0.0d0
    do j = 1, 8
      N(1,3*j-2) = func(j)
      N(2,3*j-1) = func(j)
      N(3,3*j  ) = func(j)
    enddo

    D = 0.0d0
    D(1,1) = rho
    D(2,2) = rho
    D(3,3) = rho

    DN = matmul(D, N)

    do i = 1, 24
      do j = 1, 24
        do k = 1, 3
          mass(j,i) = mass(j,i) + N(k,j)*DN(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C3D8_Mmat

  subroutine C3D8_volume(node, volume)
    implicit none
    integer(kint) :: i, j
    real(kdouble) :: node(3,8), mass(8,8), volume
    real(kdouble) :: r(3), func(8), wg, det
    real(kdouble) :: dndx(8,3)

    wg = 1.0d0
    mass = 0.0d0
    volume = 0.0d0

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_shapefunc(r, func)
      call monolis_C3D8_get_global_deriv(node, r, dndx, det)
      call C3D8_Vmat(func, wg, det, mass)
    enddo

    do i = 1, 8
      do j = 1, 8
        volume = volume + mass(j,i)
      enddo
    enddo
  end subroutine C3D8_volume

  subroutine C3D8_Vmat(func, wg, det, mass)
    implicit none
    integer(kint) :: i, j
    real(kdouble) :: mass(8,8), wg, det
    real(kdouble) :: func(8)

    do i = 1, 8
      do j = 1, 8
        mass(j,i) = mass(j,i) + func(j)*func(i)*wg*det
      enddo
    enddo
  end subroutine C3D8_Vmat

  subroutine C3D8_Bmat(dndx, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: B(6,24), dndx(8,3)

    B = 0.0d0
    do i = 1,8
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

  subroutine C3D8_Bmat_IC(dndx, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: B(6,33), dndx(11,3)

    B = 0.0d0
    do i = 1, 11
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
  end subroutine C3D8_Bmat_IC

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

  subroutine C3D8_Kmat_IC(D, B, wg, det, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(33,33), D(6,6), B(6,33), DB(6,33), wg, det

    DB = matmul(D, B)

    do i = 1, 33
      do j = 1, 33
        do k = 1, 6
          stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo
  end subroutine C3D8_Kmat_IC
end module mod_soild_c3d8
