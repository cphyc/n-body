module physics
use constants

implicit none

private

public :: initial_speeds, compute_force, compute_force_omp, compute_initial_variables, compute_energy, integrate
contains

   subroutine initial_speeds (r, v)
      implicit none

      real(kind = xp), intent(in)  :: r(:,:)
      real(kind = xp), intent(out) :: v(:,:)

      real(kind = xp) :: omega

      omega = 0.5_xp / pi
      v(:, 1) =   omega * r(:, 2)
      v(:, 2) = - omega * r(:, 1)
      v(:, 3) = 0._xp

   end subroutine initial_speeds

   subroutine compute_initial_variables()
      implicit none

      ! We take 1/20th of the initial caracteristic distance
      epsilon2 = npoints**(2._xp/3._xp) / 400._xp

   end subroutine compute_initial_variables

   subroutine compute_force (m, r, a)
      implicit none

      real(kind=xp), intent(in)  :: m(:)
      real(kind=xp), intent(in)  :: r(:,:)

      real(kind=xp), intent(out) :: a(:,:)

      real(kind=xp), dimension(3) :: vec, tmp
      integer :: i, j

      a = 0._xp

      !$OMP PARALLEL PRIVATE(j, vec, tmp)
      !$OMP DO SCHEDULE(DYNAMIC, 100)
      do i = 1, npoints

         do j = i+1, npoints

            vec = r(i, :) - r(j, :)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            !$OMP CRITICAL
            a(i, :) = a(i, :) - tmp*m(j)*vec
            a(j, :) = a(j, :) + tmp*m(i)*vec
            !$OMP END CRITICAL

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force


   subroutine compute_force_omp (m, r, a)
     implicit none

     real(kind=xp), dimension(:),    intent(in)  :: m
     real(kind=xp), dimension(:, :), intent(in)  :: r

     real(kind=xp), dimension(:, :), intent(out) :: a

     real(kind=xp), dimension(3) :: vec, tmp
     integer :: i, j

     a = 0._xp

     !$OMP PARALLEL PRIVATE(j, vec, tmp)
     !$OMP DO SCHEDULE(GUIDED)
     do i = 1, npoints

        do j = 1, npoints

           vec = r(i, :) - r(j, :)
           tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

           a(i, :) = a(i, :) - tmp*m(j)*vec

        end do

     end do
     !$OMP END DO
     !$OMP END PARALLEL

   end subroutine compute_force_omp

   subroutine compute_energy (m, r, v, Ec, Ep, E)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:), v(:,:)

      real(kind=xp), intent(out) :: Ec, Ep, E

      integer :: i, j

      !$OMP PARALLEL REDUCTION(+:Ep) REDUCTION(+:Ec) PRIVATE(j)
      !$OMP DO SCHEDULE(DYNAMIC)
      do i = 1, npoints

         Ec = Ec + 0.5_xp * m(i) * norm2(v(i,:))**2

         do j = i+1, npoints

            Ep = Ep - G * m(j) * m(i) / sqrt(norm2(r(i, :) - r(j, :))**2 + epsilon2)

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      E = Ec + Ep

   end subroutine compute_energy

   subroutine integrate(f, df, dt)
      implicit none

      real(xp), intent(in) :: dt
      real(xp), intent(inout) :: f(:,:), df(:,:)

      f = f + df * dt

   end subroutine integrate

end module physics
