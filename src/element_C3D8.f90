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
      u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
      u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
      u(3,i)  = var%u(3*in  ) + var%du(3*in  )
    enddo

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x0, r, dndx, det)
      call C3D8_Bmat(u, dndx, B)
      call C3D8_Dmat(param%E, param%mu, D)
      call C3D8_Kmat(var%gauss(i,icel)%stress, dndx, D, B, wg, det, stiff)
    enddo
  end subroutine C3D8_stiff

  subroutine C3D8_Bmat(u, dndx, B)
    implicit none
    integer(kint) :: i, i1, i2, i3
    real(kdouble) :: u(3,8), B(6,24), dndx(8,3), dudx(3,3)

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

    if(isNLGeom)then
      dudx = 0.0d0
      dudx = matmul(u, dndx)
      do i = 1, 8
        i1 = 3*i-2
        i2 = 3*i-1
        i3 = 3*i
        B(1,i1) = B(1,i1) + dudx(1,1)*dndx(i,1)
        B(1,i2) = B(1,i2) + dudx(2,1)*dndx(i,1)
        B(1,i3) = B(1,i3) + dudx(3,1)*dndx(i,1)
        B(2,i1) = B(2,i1) + dudx(1,2)*dndx(i,2)
        B(2,i2) = B(2,i2) + dudx(2,2)*dndx(i,2)
        B(2,i3) = B(2,i3) + dudx(3,2)*dndx(i,2)
        B(3,i1) = B(3,i1) + dudx(1,3)*dndx(i,3)
        B(3,i2) = B(3,i2) + dudx(2,3)*dndx(i,3)
        B(3,i3) = B(3,i3) + dudx(3,3)*dndx(i,3)
        B(4,i1) = B(4,i1) + dudx(1,2)*dndx(i,1) + dudx(1,1)*dndx(i,2)
        B(4,i2) = B(4,i2) + dudx(2,2)*dndx(i,1) + dudx(2,1)*dndx(i,2)
        B(4,i3) = B(4,i3) + dudx(3,2)*dndx(i,1) + dudx(3,1)*dndx(i,2)
        B(5,i1) = B(5,i1) + dudx(1,2)*dndx(i,3) + dudx(1,3)*dndx(i,2)
        B(5,i2) = B(5,i2) + dudx(2,2)*dndx(i,3) + dudx(2,3)*dndx(i,2)
        B(5,i3) = B(5,i3) + dudx(3,2)*dndx(i,3) + dudx(3,3)*dndx(i,2)
        B(6,i1) = B(6,i1) + dudx(1,3)*dndx(i,1) + dudx(1,1)*dndx(i,3)
        B(6,i2) = B(6,i2) + dudx(2,3)*dndx(i,1) + dudx(2,1)*dndx(i,3)
        B(6,i3) = B(6,i3) + dudx(3,3)*dndx(i,1) + dudx(3,1)*dndx(i,3)
      enddo
    endif
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

  subroutine C3D8_Kmat(stress, dndx, D, B, wg, det, stiff)
    implicit none
    integer(kint) :: i, j, k
    real(kdouble) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det
    real(kdouble) :: stress(6), S(9,9), BN(9,24), SBN(9,24), dndx(8,3)

    DB = matmul(D, B)

    do i = 1, 24
      do j = 1, 24
        do k = 1, 6
          stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
        enddo
      enddo
    enddo

    if(isNLGeom)then
      BN = 0.0d0
      do j = 1, 8
        BN(1, 3*j-2) = dndx(j, 1)
        BN(2, 3*j-1) = dndx(j, 1)
        BN(3, 3*j  ) = dndx(j, 1)
        BN(4, 3*j-2) = dndx(j, 2)
        BN(5, 3*j-1) = dndx(j, 2)
        BN(6, 3*j  ) = dndx(j, 2)
        BN(7, 3*j-2) = dndx(j, 3)
        BN(8, 3*j-1) = dndx(j, 3)
        BN(9, 3*j  ) = dndx(j, 3)
      enddo

      S = 0.0d0
      do j = 1, 3
        S(j  , j  ) = stress(1)
        S(j  , j+3) = stress(4)
        S(j  , j+6) = stress(6)
        S(j+3, j  ) = stress(4)
        S(j+3, j+3) = stress(2)
        S(j+3, j+6) = stress(5)
        S(j+6, j  ) = stress(6)
        S(j+6, j+3) = stress(5)
        S(j+6, j+6) = stress(3)
      enddo

      SBN = matmul(S, BN)
      do i = 1, 24
        do j = 1, 24
          do k = 1, 9
            stiff(j,i) = stiff(j,i) + BN(k,j)*SBN(k,i)*wg*det
          enddo
        enddo
      enddo
    endif
  end subroutine C3D8_Kmat

  subroutine C3D8_update(mesh, var, param, icel, q)
    implicit none
    type(meshdef) :: mesh
    type(vardef) :: var
    type(paramdef) :: param
    integer(kint) :: i, in, icel
    real(kdouble) :: x0(3,8), u(3,8), r(3), dndx(8,3), xj(3,3), D(6,6), B(6,24)
    real(kdouble) :: strain(6), stress(6), q(24), det
    real(kdouble) :: Yield0, mises, plstrain(1), PPStrain(6), ystress(6)

    q = 0.0d0

    do i = 1, 8
      in = mesh%elem(i,icel)
      x0(1,i) = mesh%node(1,in)
      x0(2,i) = mesh%node(2,in)
      x0(3,i) = mesh%node(3,in)
      u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
      u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
      u(3,i)  = var%u(3*in  ) + var%du(3*in  )
    enddo

    do i = 1, 8
      call monolis_C3D8_integral_point(i, r)
      call monolis_C3D8_get_global_deriv(x0, r, dndx, det)
      call C3D8_Bmat(u, dndx, B)
      call C3D8_Dmat(param%E, param%mu, D)
      xj = matmul(u, dndx)
      strain(1) = xj(1,1)
      strain(2) = xj(2,2)
      strain(3) = xj(3,3)
      strain(4) =(xj(1,2) + xj(2,1))
      strain(5) =(xj(2,3) + xj(3,2))
      strain(6) =(xj(3,1) + xj(1,3))

      if(isNLGeom)then
        strain(1) = strain(1) + 0.5d0*dot_product(xj(:, 1), xj(:, 1))
        strain(2) = strain(2) + 0.5d0*dot_product(xj(:, 2), xj(:, 2))
        strain(3) = strain(3) + 0.5d0*dot_product(xj(:, 3), xj(:, 3))
        strain(4) = strain(4) + (xj(1,1)*xj(1,2) + xj(2,1)*xj(2,2) + xj(3,1)*xj(3,2))
        strain(5) = strain(5) + (xj(1,2)*xj(1,3) + xj(2,2)*xj(2,3) + xj(3,2)*xj(3,3))
        strain(6) = strain(6) + (xj(1,1)*xj(1,3) + xj(2,1)*xj(2,3) + xj(3,1)*xj(3,3))
      endif

      stress = matmul(D, strain)

!      call get_mises(stress(1:6), mises)
!      plstrain = var%gauss(i,icel)%equival

!      Yield0 = 200.d0
!      if(mises > Yield0) then
!        !call BackwardEuler( param, var, mises, stress, Yield0, plstrain, PPStrain, ystress )
!        !var%gauss(i,icel)%strain = plstrain
!        var%gauss(i,icel)%stress = ystress
!        var%gauss(i,icel)%equival = plstrain
!        q = q + matmul(ystress, B) * det
!      else
        var%gauss(i,icel)%strain = strain
        var%gauss(i,icel)%stress = stress
        q = q + matmul(stress, B)*det
!      endif
    enddo
  end subroutine C3D8_update

  subroutine C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)
    implicit none
    type(vardef) :: var
    integer(kint) :: i, j, k, icel
    real(kdouble) :: inv(8,8)
    real(kdouble) :: nstrain(8,6), nstress(8,6)
    real(kdouble) :: estrain(6), estress(6)

    nstrain  = 0.0d0
    nstress  = 0.0d0
    estrain  = 0.0d0
    estress  = 0.0d0

    do i = 1, 8
      do j = 1, 8
        do k = 1, 6
          nstrain(i,k) = nstrain(i,k) + inv(i,j) * var%gauss(j,icel)%strain(k)
          nstress(i,k) = nstress(i,k) + inv(i,j) * var%gauss(j,icel)%stress(k)
        enddo
      enddo
    enddo

    do i = 1, 8
      do j = 1, 6
        estrain(j) = estrain(j) + var%gauss(i,icel)%strain(j)
        estress(j) = estress(j) + var%gauss(i,icel)%stress(j)
      enddo
    enddo
    estrain = estrain/8.0d0
    estress = estress/8.0d0
  end subroutine C3D8_get_nodal_values

  subroutine BackwardEuler(mises, stress, param, var, YIeld0, plstrain, PPStrain, ystress)
    type(paramdef) :: param
    type(vardef) :: var
    real(kdouble) :: stress(6), ystress(6)   !< trial->real stress
    !real(=kreal)   :: plstrain    !< plastic strain till current substep
    !integer, intent(inout)           :: istat       !< plastic state
    real(kdouble), parameter :: tol = 1.d-3
    integer, parameter :: MAXITER = 15
    real(kdouble) :: dlambda, f, mises
    integer(kint) :: i, j
    real(kdouble) :: youngs, poisson, plstrain, dum
    real(kdouble) :: J1, J2, J3, H, dd, G, K, devia(6)
    real(kdouble) :: PPStrain(6)

    youngs = param%E
    poisson = param%mu

    f = 0.0d0
    dum = YIeld0
    J1 = (stress(1) + stress(2) + stress(3))/3.0d0
    devia(1:3) = stress(1:3) - J1
    devia(4:6) = stress(4:6)

    if(youngs == 0.0d0) stop "YOUNG's ratio==0"
    G = youngs/(2.0d0*(1.0d0 + poisson))
    K = youngs/(3.0d0*(1.0d0 - 2.0d0*poisson))
    dlambda = 0.0d0

    !> Mises or. Isotropic
    do i = 1, MAXITER
      !> H should be given in a tabulated way
      ! H= calHardenCoeff( matl, pstrain+dlambda )
      H = 200.0d0
      f = mises - 3.0d0*G*dlambda - dum
      dd = 3.d0*G + H
      dlambda = dlambda + f/dd
      !dum = calCurrYield(plstrain + dlambda, YIeld0)

      if(dabs(f) < tol*tol) exit
    enddo

    !> Update equivalent plastic strain and stress
    plstrain = plstrain + dlambda

    !> Update strain components
    PPStrain(1:3) = PPStrain(1:3) + 1.5d0*dlambda*devia(1:3)/mises
    PPStrain(4:6) = PPStrain(4:6) + 3.0d0*dlambda*devia(4:6)/mises

    !>Update stress
    devia(:) = dum/mises*devia(:)
    ystress(1:3) = devia(1:3) + J1
    ystress(4:6) = devia(4:6)
  end subroutine BackwardEuler

  !> This function calcualtes current yield stress
  function calCurrYield(plstrain, YIeld0)
    real(kdouble), intent(in) :: plstrain !< plastic strain
    real(kdouble), intent(in) :: YIeld0   !< plastic strain
    real(kdouble) :: calCurrYield

    calCurrYield = - 1.0d0
    !> Linear hardening,
    calCurrYield = YIeld0 + 200.0d0*plstrain
  end function
end module mod_soild_c3d8
