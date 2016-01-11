module physics
  use constants

  implicit none

  private

  public :: initial_speeds, compute_force, compute_initial_variables, compute_energy, integrate
contains

  subroutine initial_speeds (r, v)
    real(kind = xp), dimension(:, :), intent(in)  :: r
    real(kind = xp), dimension(:, :), intent(out) :: v

    real(kind = xp) :: omega

    omega = 0.5_xp / pi
    v(:, 1) =   omega * r(:, 2)
    v(:, 2) = - omega * r(:, 1)
    v(:, 3) = 0._xp

  end subroutine initial_speeds
  
  subroutine compute_initial_variables()
    ! We take 1/20th of the initial caracteristic distance
    epsilon2 = npoints**(2._xp/3._xp) / 400._xp
  end subroutine compute_initial_variables
  
  subroutine compute_force (m, r, a)
    real(kind=xp), dimension(:), intent(in)     :: m
    real(kind=xp), dimension(:, :), intent(in)  :: r

    real(kind=xp), dimension(:, :), intent(out) :: a

    real(kind=xp), dimension(3) :: vec, tmp

    integer :: i, j

    a = 0._xp
    do i = 1, npoints
       do j = i+1, npoints
          vec = r(i, :) - r(j, :)
          tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp
          a(i, :) = a(i, :) - tmp*m(j)*vec
          a(j, :) = a(j, :) + tmp*m(i)*vec
       end do
    end do

  end subroutine compute_force

  subroutine compute_energy (m, r, v, Ec, Ep, E)
    real(kind=xp), dimension(:), intent(in)    :: m
    real(kind=xp), dimension(:, :), intent(in) :: r, v

    real(kind=xp), intent(out) :: Ec, Ep, E

    integer :: i, j
    Ec = 0._xp
    Ep = 0._xp
    do i = 1, npoints
       Ec = Ec + 0.5_xp * m(i) * norm2(v(i,:))**2
       do j = i+1, npoints
          Ep = Ep - G*m(j)*m(i)/sqrt(norm2(r(i, :) - r(j, :))**2 + epsilon2)
       end do

    end do

    E = Ec + Ep
  end subroutine compute_energy

  subroutine integrate(f, df, dt)
    implicit none

    real(xp), dimension(:, :), intent(inout) :: f, df
    real(xp),                  intent(in)    :: dt

    f = f + df * dt

  end subroutine integrate

end module physics
