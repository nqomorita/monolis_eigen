program main
  implicit none
  real(8), parameter :: E = 100000.0d0
  real(8), parameter :: mu = 0.3d0
  real(8), parameter :: P = 1.0d0
  real(8), parameter :: L = 1.0d0
  integer(4) :: i, in, j, k
  integer(4) :: nnode, nelem
  integer(4), allocatable :: nid(:), elem(:,:), perm(:)
  real(8) :: l2, l2t, u2, u2t, nodet(3,8), ut(3,8), diff, gspt(3), wt, det
  real(8) :: func(8), sim(3), theo(3), x(3), r(3), s(3), t2
  real(8), allocatable :: node(:,:), u(:,:)
  !character :: cid1*4, cid2*4

  !> for gauss integration
  integer(4), parameter :: intp = 8
  real(8), allocatable :: gsp(:,:)
  real(8), allocatable :: w(:)

  open(10, file="node.dat", status='old')
    read(10,*) nnode
    allocate(node(3, nnode), Source = 0.0d0)
    allocate(nid(nnode), Source = 0)
    do i = 1, nnode
      read(10,*)nid(i), node(1,i), node(2,i), node(3,i)
    enddo
  close(10)

  open(10, file="elem.dat", status='old')
    read(10,*) nelem
    allocate(elem(8, nelem), Source = 0)
    do i = 1, nelem
      read(10,*)in, elem(1,i), elem(2,i), elem(3,i), elem(4,i), &
                  & elem(5,i), elem(6,i), elem(7,i), elem(8,i)
    enddo
  close(10)

  open(10, file="u.dat", status='old')
    read(10,*) in
    if(in /= nnode) stop "error"
    allocate(u(3, nnode), Source = 0.0d0)
    do i = 1, nnode
      read(10,*)in, u(1,i), u(2,i), u(3,i)
    enddo
  close(10)

  call global_to_local(nnode, nid, nelem, elem, 8, perm)
  deallocate(perm)

  call set_C3D8_integral_point_kgl(intp)

  l2 = 0.0d0
  u2 = 0.0d0
  do i = 1, nelem
    do j = 1, 8
      nodet(1,j) = node(1,elem(j,i))
      nodet(2,j) = node(2,elem(j,i))
      nodet(3,j) = node(3,elem(j,i))
      ut(1,j) = u(1,elem(j,i))
      ut(2,j) = u(2,elem(j,i))
      ut(3,j) = u(3,elem(j,i))
    enddo

    l2t = 0.0d0
    u2t = 0.0d0
    do j = 1, intp**3
      gspt = gsp(:,j)
      wt = w(j)
      call C3D8_get_global_deriv(nodet, gspt, det)

      !> simulated value
      call C3D8_shapefunc(gspt, func)
      sim = matmul(ut, func)

      !> theoreticla value
      call C3D8_get_global_pos(nodet, gspt, x)
      call get_spherical_coordinate(x, r)
      call get_theoretical_value(r, s)
      call get_cartesian_disp(r, s, x, theo)

      !> diff
      t2 = get_l2(theo - sim)
      l2t = l2t + t2*wt*det

      t2 = get_l2(theo)
      u2t = u2t + t2*wt*det
    enddo
    l2 = l2 + l2t
    u2 = u2 + u2t
  enddo
  write(*,"(a,1pe12.5)") "l2   ", dsqrt(l2)
  write(*,"(a,1pe12.5)") "u2   ", dsqrt(u2)
  write(*,"(a,1pe12.5)") "error", dsqrt(l2/u2)

contains

  subroutine get_cartesian_disp(r, s, x, u)
    implicit none
    real(8) :: r(3), s(3), u(3), x(3)

    u(1) = (s(1)*dsin(r(3)) + s(3)*dcos(r(3)))*dcos(r(2)) - P/E*mu*x(1)
    u(2) = (s(1)*dsin(r(3)) + s(3)*dcos(r(3)))*dsin(r(2)) - P/E*mu*x(2)
    u(3) =  s(1)*dcos(r(3)) - s(3)*dsin(r(3))             + P/E*x(3)
  end subroutine get_cartesian_disp

  function get_l2(a)
    implicit none
    real(8) :: get_l2, a(3)
    get_l2 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
   ! get_l2 = dsqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
  end function get_l2

  subroutine C3D8_get_global_pos(node, r, pos)
    implicit none
    real(8) :: node(3,8), r(3), pos(3)
    real(8) :: func(8)

    call C3D8_shapefunc(r, func)
    pos = matmul(node, func)
  end subroutine C3D8_get_global_pos

  subroutine get_theoretical_value(r_in, u)
    implicit none
    real(8) :: A, B, C, lambda, rho
    real(8) :: r_in(3), u(3), r, theta, c1, c2, c3

    lambda = mu*E/((1.0d0 + mu)*(1.0d0 - 2.0d0*mu))
    rho = E/(2.0d0 + 2.0d0*mu)

    c3 = P/(8.0d0*rho)
    A =-L*L*L*c3*(13.0d0 - 10.0d0*mu)/(7.0d0 - 5.0d0*mu)
    B = L*L*L*L*L*c3*1.0d0 /(7.0d0 - 5.0d0*mu)
    C = L*L*L*c3*( 5.0d0 - 10.0d0*mu)/(7.0d0 - 5.0d0*mu)

    u = 0.0d0
    r = r_in(1)
    theta = r_in(3)
    c1 = (5.0d0 - 4.0d0*mu)/(1.0d0 - 2.0d0*mu)*C/(r*r) - 9.0d0*B/(r*r*r*r)
    u(1) = - A/(r*r) - 3.0d0*B/(r*r*r*r) + c1*dcos(2.0d0*theta)
    c2 = - 2.0d0*C/(r*r) - 6.0d0*B/(r*r*r*r)
    u(2) = 0.0d0
    u(3) = c2*dsin(2.0d0*theta)
  end subroutine get_theoretical_value

  subroutine get_spherical_coordinate(x, r)
    implicit none
    real(8) :: x(3), r(3), l2

    r(1) = dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))

    if(x(1) == 0.0d0)then
      r(2) = 0.5d0*3.141592653589d0
    else
      r(2) = datan(x(2)/x(1))
    endif

    if(x(3) == 0.0d0)then
      r(3) = 0.0d0
    else
      l2 = dsqrt(x(1)*x(1) + x(2)*x(2))
      r(3) = datan(l2/x(3))
    endif
  end subroutine get_spherical_coordinate

  subroutine global_to_local(nnode, nid, nelem, e, nenode, perm)
    implicit none
    integer(4) :: i, in, j, nenode
    integer(4) :: imax, imin
    integer(4) :: nnode, nid(:)
    integer(4) :: nelem, e(:,:)
    integer(4), allocatable :: temp(:)
    integer(4), optional, allocatable :: perm(:)

    imax = maxval(nid)
    imin = minval(nid)
    allocate(temp(imin:imax))
    temp = 0

    in = 1
    do i = 1, nnode
      temp(nid(i)) = in
      in = in + 1
    enddo

    do i = 1, nelem
      do j = 1, nenode
        in = e(j,i)
        e(j,i) = temp(in)
      enddo
    enddo

    if(present(perm))then
      allocate(perm(imin:imax))
      do i = imin, imax
        perm(i) = temp(i)
      enddo
    endif
    deallocate(temp)
  end subroutine global_to_local

  subroutine set_C3D8_integral_point_kgl(intp)
    implicit none
    integer(4) :: intp, j, k, l
    real(8), allocatable :: r(:,:)

    allocate(gsp(3,intp**3))
    allocate(w(intp**3))
    gsp = 0.0d0
    w = 0.0d0

    r = get_integ_points_and_weights(intp)

    do j = 1, intp
      do k = 1, intp
        do l = 1, intp
          gsp(1,l + (k-1)*intp + (j-1)*intp**2) = r(1,l)
          gsp(2,l + (k-1)*intp + (j-1)*intp**2) = r(1,k)
          gsp(3,l + (k-1)*intp + (j-1)*intp**2) = r(1,j)
        enddo
      enddo
    enddo
    do j = 1, intp
      do k = 1, intp
        do l = 1, intp
          w(l + (k-1)*intp + (j-1)*intp**2) = r(2,j)*r(2,k)*r(2,l)
        enddo
      enddo
    enddo
  end subroutine set_C3D8_integral_point_kgl

  function get_integ_points_and_weights(n)
    implicit none
    integer(4) :: n, i, k, iter
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    real(8) :: get_integ_points_and_weights(2, n)
    real(8) :: x, f, df, dx
    real(8), allocatable :: p0(:), p1(:), tmp(:)

    p0 = [1.d0]
    p1 = [1.d0, 0.d0]

    do k = 2, n
       tmp = (2.0d0*dble(k) - 1)*[p1,0.d0] - (k - 1.0d0)*[0.d0, 0.d0, p0]
       tmp = tmp/dble(k)
       p0 = p1
       p1 = tmp
    enddo

    do i = 1, n
      x = dcos(pi*(i - 0.25d0)/(n + 0.5d0))
      do iter = 1, 10
        f = p1(1)
        df = 0.0d0
        do k = 2, size(p1)
          df = f + x*df
          f  = p1(k) + x * f
        enddo
        dx =  f / df
        x = x - dx
        if (abs(dx) < 1.d2*epsilon(dx)) exit
      enddo

      get_integ_points_and_weights(1,i) = x
      get_integ_points_and_weights(2,i) = 2.0d0/((1.0d0 - x*x)*df*df)
    enddo
  end function get_integ_points_and_weights

  subroutine C3D8_get_global_deriv(node, r, det)
    implicit none
    real(8) :: node(3,8), r(3), deriv(8,3) !dndx(8,3)
    real(8) :: xj(3,3), inv(3,3), det

    call C3D8_shapefunc_deriv(r, deriv)
    xj = matmul(node, deriv)
    call C3D8_get_inverse_matrix(xj, inv, det)
    !dndx = matmul(deriv, inv)
  end subroutine C3D8_get_global_deriv

  subroutine C3D8_shapefunc(local, func)
    implicit none
    real(8) ::  local(3), func(8)

    func(1) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0-local(3))
    func(2) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0-local(3))
    func(3) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0-local(3))
    func(4) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0-local(3))
    func(5) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0+local(3))
    func(6) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0+local(3))
    func(7) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0+local(3))
    func(8) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0+local(3))
  end subroutine C3D8_shapefunc

  subroutine C3D8_shapefunc_deriv(local, func)
    implicit none
    real(8) :: local(3), func(8,3)

    func(1,1) = -0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
    func(2,1) =  0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
    func(3,1) =  0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
    func(4,1) = -0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
    func(5,1) = -0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
    func(6,1) =  0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
    func(7,1) =  0.125d0*(1.0d0+local(2))*(1.0d0+local(3))
    func(8,1) = -0.125d0*(1.0d0+local(2))*(1.0d0+local(3))

    func(1,2) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
    func(2,2) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
    func(3,2) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
    func(4,2) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
    func(5,2) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(3))
    func(6,2) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
    func(7,2) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
    func(8,2) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(3))

    func(1,3) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(2,3) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(3,3) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(4,3) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
    func(5,3) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(6,3) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(7,3) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(8,3) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
  end subroutine C3D8_shapefunc_deriv

  subroutine C3D8_get_inverse_matrix(xj, inv, det)
    implicit none
    real(8) :: xj(3,3), inv(3,3), det, detinv

    det = xj(1,1) * xj(2,2) * xj(3,3) &
        + xj(2,1) * xj(3,2) * xj(1,3) &
        + xj(3,1) * xj(1,2) * xj(2,3) &
        - xj(3,1) * xj(2,2) * xj(1,3) &
        - xj(2,1) * xj(1,2) * xj(3,3) &
        - xj(1,1) * xj(3,2) * xj(2,3)

    if(det == 0.0d0) stop "determinant = 0.0"

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
end program main

