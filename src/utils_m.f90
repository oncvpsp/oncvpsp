module utils_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: linspace, nonuniform_trapezoid, fermi_dirac, &
      get_pseudo_linear_mesh_parameters
contains

subroutine linspace(first, last, step, n, array)
   implicit none
   ! Input parameters
   real(dp), intent(in) :: first
   real(dp), intent(in) :: last
   real(dp), intent(in) :: step
   ! Output parameters
   integer, intent(out) :: n
   real(dp), allocatable, intent(out) :: array(:)
   ! Local variables
   integer :: i

   n = nint((last - first) / step) + 1
   allocate(array(n))
   do i = 1, n
      array(i) = first + real(i - 1, dp) * step
   end do
end subroutine linspace

function nonuniform_trapezoid(mmax, rr, f) result(integral)
   implicit none
   ! Input variables
   integer, intent(in) :: mmax
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: f(mmax)

   ! Output variable
   real(dp) :: integral

   ! Local variables
   integer :: ii

   integral = 0.0_dp
   do ii = 1, mmax - 1
      integral = integral + 0.5_dp * (f(ii) + f(ii + 1)) * (rr(ii + 1) - rr(ii))
   end do
   return
end function nonuniform_trapezoid

function fermi_dirac(x, mu, sigma) result(y)
   implicit none
   ! Input variables
   real(dp), intent(in) :: x
   real(dp), intent(in), optional :: mu
   real(dp), intent(in), optional :: sigma

   ! Output variable
   real(dp) :: y

   ! Local variables
   real(dp) :: mu_loc, sigma_loc

   if (present(mu)) then
      mu_loc = mu
   else
      mu_loc = 0.0_dp
   end if

   if (present(sigma)) then
      sigma_loc = sigma
   else
      sigma_loc = 1.0_dp
   end if

   y = 1.0_dp / (1.0_dp + exp((x - mu_loc) / sigma_loc))
   return
end function fermi_dirac

subroutine get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, n1, n2, n3, n4)
   implicit none
   ! Input variables
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl

   ! Output variables
   integer, intent(out) :: n1, n2, n3, n4

   ! Local variables
   integer :: l1
   real(dp) :: al

   al = 0.01_dp * log(rr(101) / rr(1))
   n1 = int(log(drl / rr(1)) / al + 1.0_dp)
   n2 = int(log(real(nrl, dp) * drl / rr(1)) / al + 1.0_dp)
   n3 = 0
   do l1 = 1, lmax + 1
      n3 = max(n3, irc(l1) - 1)
   end do
   n4 = min(n2, int(log((rr(n3) + 1.0_dp) / rr(1)) / al))
   n3 = int(log(1.1_dp * rr(n3) / rr(1)) / al + 1.0d0)
end subroutine get_pseudo_linear_mesh_parameters

end module utils_m
