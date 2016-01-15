module physics
use constants

implicit none

private

public :: initial_speeds, integrate, integrate_omp, &
          compute_force, compute_force_omp, compute_force_omp_diag, &
          compute_energy, compute_energy_omp, compute_energy_omp_diag

contains

   subroutine initial_speeds (r, v)
      implicit none

      real(kind = xp), intent(in)  :: r(:,:)
      real(kind = xp), intent(out) :: v(:,:)

      v(1, :) =   r(2, :)
      v(2, :) = - r(1, :)
      v(3, :) = 0._xp

   end subroutine initial_speeds

   subroutine integrate(f, df, dt)
      implicit none

      real(xp), intent(in) :: dt
      real(xp), intent(inout) :: f(:,:), df(:,:)

      f = f + df * dt

   end subroutine integrate

   subroutine integrate_omp(f, df, dt)
      implicit none

      real(xp), intent(in) :: dt
      real(xp), intent(inout) :: f(:,:), df(:,:)

      !$OMP PARALLEL
      !$OMP WORKSHARE
      f = f + df * dt
      !$OMP END WORKSHARE
      !$OMP END PARALLEL

   end subroutine integrate_omp

   subroutine compute_force (m, r1, r2, length, a)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r1(:,:), r2(:,:)
      integer                   :: length

      real(kind=xp), intent(out) :: a(:,:)

      real(kind=xp), dimension(3) :: vec, tmp
      integer :: i, j

      a = 0._xp

      do i = 1, length

         do j = i+1, length

            vec = r1(:, i) - r2(:, j)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            a(:, i) = a(:, i) - tmp*m(j)*vec
            a(:, j) = a(:, j) + tmp*m(i)*vec

         end do

      end do

   end subroutine compute_force

   subroutine compute_force_omp (m, r1, r2, length, a)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r1(:, :), r2(:, :)
      integer                   :: length

      real(kind=xp), intent(out) :: a(:, :)

      real(kind=xp), dimension(3) :: vec, tmp
      integer :: i, j

      a = 0._xp

      !$OMP PARALLEL PRIVATE(j, vec, tmp)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = 1, length

         do j = 1, length

            if (i /= j) then
               vec = r1(:, i) - r2(:, j)
               tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

               a(:, i) = a(:, i) - tmp*m(j)*vec
            end if

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force_omp

   subroutine compute_force_omp_diag (m, r1, r2, length, a)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r1(:,:), r2(:,:)
      integer                   :: length

      real(kind=xp), intent(out) :: a(:,:)

      real(kind=xp), dimension(3) :: vec, tmp
      integer :: i, j, k

      a = 0._xp

      !$OMP PARALLEL PRIVATE(j, k, vec, tmp) REDUCTION(+:a)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = 1, length/2

         do j = 1, i-1

            vec = r1(:, i) - r2(:, j)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            a(:, i) = a(:, i) - tmp*m(j)*vec
            a(:, j) = a(:, j) + tmp*m(i)*vec

         end do

         k = length + 1 - i

         do j = 1, k-1

            vec = r1(:, k) - r2(:, j)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            a(:, k) = a(:, k) - tmp*m(j)*vec
            a(:, j) = a(:, j) + tmp*m(k)*vec

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force_omp_diag

   subroutine compute_energy (m, r, v, Ec, Ep, E)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:), v(:,:)

      real(kind=xp), intent(out) :: Ec, Ep, E

      integer :: i, j

      Ec = 0._xp
      Ep = 0._xp

      do i = 1, npoints

         Ec = Ec + 0.5_xp * m(i) * norm2(v(:,i))**2

         do j = i+1, npoints

            Ep = Ep - G * m(j) * m(i) / sqrt(norm2(r(:, i) - r(:, j))**2 + epsilon2)

         end do

      end do

      E = Ec + Ep

   end subroutine compute_energy

   subroutine compute_energy_omp (m, r, v, Ec, Ep, E)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:), v(:,:)

      real(kind=xp), intent(out) :: Ec, Ep, E

      integer :: i, j

      Ec = 0._xp
      Ep = 0._xp

      !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = 1, npoints

         Ec = Ec + 0.5_xp * m(i) * norm2(v(:,i))**2

         do j = 1, npoints

            if (i /= j) then
               Ep = Ep - 0.5_xp * G * m(j) * m(i) / sqrt(norm2(r(:, i) - r(:, j))**2 + epsilon2)
            end if

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      E = Ec + Ep

   end subroutine compute_energy_omp

   subroutine compute_energy_omp_diag (m, r, v, Ec, Ep, E)
     implicit none

     real(kind=xp), intent(in) :: m(:)
     real(kind=xp), intent(in) :: r(:,:), v(:,:)

     real(kind=xp), intent(out) :: Ec, Ep, E

     integer :: i, j, k

     Ec = 0._xp
     Ep = 0._xp

     !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j,k)
     !$OMP DO SCHEDULE(RUNTIME)
     do i = 1, npoints/2

        Ec = Ec + 0.5_xp * m(i) * norm2(v(:,i))**2

        do j = 1, i-1

           Ep = Ep -  G * m(j) * m(i) / sqrt(norm2(r(:, i) - r(:, j))**2 + epsilon2)

        end do

        k = npoints + 1 - i
        Ec = Ec + 0.5_xp * m(k) * norm2(v(:,k))**2

        do j = 1, k-1

           Ep = Ep -  G * m(j) * m(k) / sqrt(norm2(r(:, k) - r(:, j))**2 + epsilon2)

        end do

     end do
     !$OMP END DO
     !$OMP END PARALLEL

     E = Ec + Ep

   end subroutine compute_energy_omp_diag

end module physics
