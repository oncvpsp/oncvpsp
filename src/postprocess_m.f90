!
! Copyright (c) 1989-2025 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

module postprocess_m
   use, intrinsic :: iso_fortran_env, only: stderr => output_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use output_text_m, only: get_pseudo_linear_mesh_parameters
   implicit none
   private
   public :: get_wavefunctions, get_bound_wavefunctions, get_scattering_wavefunctions
   public :: run_test_configurations
contains

subroutine get_wavefunctions(zz, srel, mmax, rr, vfull, lloc, vp, lmax, &
                             irc, nproj, drl, nrl, mxprj, epa, npa, vkb, evkb, &
                             uu_ae, up_ae, uu_ps, up_ps, mch_ae, mch_ps, e_ae, e_ps, &
                             sign_ae, sign_ps, is_scattering)
   implicit none
   ! Input variables
   !> Atomic number
   real(dp), intent(in) :: zz
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> All-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> Angular momentum of local potential
   integer, intent(in) :: lloc
   !> Semi-local potentials
   real(dp), intent(in) :: vp(mmax,5)
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Number of projectors for each angular momentum
   integer, intent(in) :: nproj(6)
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Guess all-electron state energies
   real(dp), intent(in) :: epa(mxprj,6)
   !> Principal quantum numbers
   integer, intent(in) :: npa(mxprj,6)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax,mxprj,4)
   !> Coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj,4)

   ! Output variables
   !> All-electron wavefunctions
   real(dp), intent(out) :: uu_ae(mmax,mxprj,4)
   !> Pseudo wavefunctions
   real(dp), intent(out) :: uu_ps(mmax,mxprj,4)
   !> All-electron wavefunction radial derivatives
   real(dp), intent(out) :: up_ae(mmax,mxprj,4)
   !> Pseudo wavefunction radial derivatives
   real(dp), intent(out) :: up_ps(mmax,mxprj,4)
   !> Matching points for all-electron wavefunctions
   integer, intent(out) :: mch_ae(mxprj,4)
   !> Matching points for pseudo wavefunctions
   integer, intent(out) :: mch_ps(mxprj,4)
   !> All-electron state energies
   real(dp), intent(out) :: e_ae(mxprj,4)
   !> Pseudo state energies
   real(dp), intent(out) :: e_ps(mxprj,4)
   !> Signs of all-electron wavefunctions at matching points
   real(dp), intent(out) :: sign_ae(mxprj,4)
   !> Signs of pseudo wavefunctions at matching points
   real(dp), intent(out) :: sign_ps(mxprj,4)
   !> .true. for scattering states, .false. for bound states
   logical, intent(out) :: is_scattering(mxprj,4)

   ! Local variables
   integer :: l1, ll, npr, iproj
   integer :: n1, n2, n3, n4

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, n1, n2, n3, n4)

   do l1 = 1, lmax + 1
      ll = l1 - 1
      if (l1 == lloc + 1) then
         npr = 0
      else
         npr = nproj(l1)
      end if
      do iproj = 1, nproj(l1)
         if (epa(iproj, l1) < 0.0_dp) then
            is_scattering(iproj, l1) = .false.
            call get_bound_wavefunctions(epa(iproj, l1), npa(iproj, l1), ll, npr, iproj, zz, srel, mmax, rr, &
                                         vfull, vp(:, lloc + 1), mxprj, vkb(:, :, l1), evkb(:, l1), &
                                         uu_ae(:, iproj, l1), up_ae(:, iproj, l1), &
                                         uu_ps(:, iproj, l1), up_ps(:, iproj, l1), &
                                         mch_ae(iproj, l1), mch_ps(iproj, l1), &
                                         e_ae(iproj, l1), e_ps(iproj, l1), &
                                         sign_ae(iproj, l1), sign_ps(iproj, l1))
         else
            is_scattering(iproj, l1) = .true.
            call get_scattering_wavefunctions(epa(iproj, l1), npa(iproj, l1), ll, nproj(l1), iproj, zz, srel, mmax, rr, &
                                              vfull, vp(:, lloc + 1), mxprj, vkb(:, :, l1), evkb(:, l1), n2, &
                                              uu_ae(:, iproj, l1), up_ae(:, iproj, l1), &
                                              uu_ps(:, iproj, l1), up_ps(:, iproj, l1), &
                                              e_ae(iproj, l1), e_ps(iproj, l1), &
                                              sign_ae(iproj, l1), sign_ps(iproj, l1))
            mch_ae(iproj, l1) = n2
            mch_ps(iproj, l1) = n2
         end if
      end do
   end do
end subroutine get_wavefunctions

subroutine get_bound_wavefunctions(e_in, nn, ll, nproj, iproj, zz, srel, mmax, rr, vfull, vloc, mxprj, vkb, evkb, &
                                   uu_ae, up_ae, uu_ps, up_ps, mch_ae, mch_ps, e_ae, e_ps, sign_ae, sign_ps)
   implicit none
   ! Input variables
   !> guess all-electron state energy
   real(dp), intent(in) :: e_in
   !> principal quantum number
   integer, intent(in) :: nn
   !> angular momentum
   integer, intent(in) :: ll
   !> number of projectors for this l
   integer, intent(in) :: nproj
   !> index of projector at this l
   integer, intent(in) :: iproj
   !> atomic number
   real(dp), intent(in) :: zz
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel
   !> size of log radial mesh
   integer, intent(in) :: mmax
   !> logarithmic radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> local potential
   real(dp), intent(in) :: vloc(mmax)
   !> maximum number of projectors
   integer, intent(in) :: mxprj
   !> VKB projectors for this l
   real(dp), intent(in) :: vkb(mmax, mxprj)
   !> coefficients of VKB projectors for this l
   real(dp), intent(in) :: evkb(mxprj)

   ! Output variables
   real(dp), intent(out) :: uu_ae(mmax)
   real(dp), intent(out) :: up_ae(mmax)
   real(dp), intent(out) :: uu_ps(mmax)
   real(dp), intent(out) :: up_ps(mmax)
   integer, intent(out) :: mch_ae
   integer, intent(out) :: mch_ps
   real(dp), intent(out) :: e_ps, e_ae
   real(dp), intent(out) :: sign_ae, sign_ps

   ! Local variables
   integer :: ierr
   real(dp) :: e_max, e_min

   e_ae = e_in
   call lschfb(nn, ll, ierr, e_ae, rr, vfull, uu_ae, up_ae, zz, mmax, mch_ae, srel)
   if (ierr /= 0) then
      write(stderr,'(/a,3i4)') 'run_plot: lschfb convergence ERROR n,l,ierr=', nn, ll, ierr
      stop
   end if

   e_max = 0.9_dp * e_ae
   e_min = 1.1_dp * e_ae
   e_ps = e_ae
   call lschvkbb(ll + iproj, ll, nproj, ierr, e_ps, e_min, e_max, &
                 rr, vloc, vkb, evkb, uu_ps, up_ps, mmax, mch_ps)
   if (ierr /= 0) then
      write(stderr,'(/a,3i4)') 'run_plot: lschvkbb convergence ERROR l,iproj,ierr=', ll, iproj, ierr
      stop
   end if

   sign_ae = 1.0_dp
   sign_ps = 1.0_dp
end subroutine get_bound_wavefunctions

subroutine get_scattering_wavefunctions(e_in, n, l, nproj, iproj, z, srel, mmax, rr, &
                                        vfull, vloc, mxprj, vkb, evkb, mch, &
                                        uu_ae, up_ae, uu_ps, up_ps, e_ae, e_ps, sign_ae, sign_ps)
   implicit none
   ! Input variables
   !> guess all-electron state energy
   real(dp), intent(in) :: e_in
   !> principal quantum number
   integer, intent(in) :: n
   !> angular momentum
   integer, intent(in) :: l
   !> number of projectors for this l
   integer, intent(in) :: nproj
   !> index of projector at this l
   integer, intent(in) :: iproj
   !> atomic number
   real(dp), intent(in) :: z
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel
   !> size of log radial mesh
   integer, intent(in) :: mmax
   !> logarithmic radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> local potential
   real(dp), intent(in) :: vloc(mmax)
   !> maximum number of projectors
   integer, intent(in) :: mxprj
   !> VKB projectors for this l
   real(dp), intent(in) :: vkb(mmax, mxprj)
   !> coefficients of VKB projectors for this l
   real(dp), intent(in) :: evkb(mxprj)
   !> matching point
   integer, intent(in) :: mch

   ! Output variables
   real(dp), intent(out) :: uu_ae(mmax)
   real(dp), intent(out) :: up_ae(mmax)
   real(dp), intent(out) :: uu_ps(mmax)
   real(dp), intent(out) :: up_ps(mmax)
   real(dp), intent(out) :: e_ps, e_ae
   real(dp), intent(out) :: sign_ae, sign_ps

   ! Local variables
   integer :: ierr, n_effective

   e_ae = e_in
   call lschfs(n_effective, l, ierr, e_ae, rr, vfull, uu_ae, up_ae, z, mmax, mch, srel)
   if (ierr /= 0) then
      write(stderr, '(/a,3i4)') 'run_plot: lschfs convergence ERROR n,l,ierr=', n, l, ierr
      stop
   end if

   e_ps = e_ae
   call lschvkbs(l, nproj, e_ps, rr, vloc, vkb, evkb, uu_ps, up_ps, mmax, mch)
   if (ierr /= 0) then
      write(stderr,'(/a,3i4)') 'run_plot: lschvkbs convergence ERROR n,l,ierr=', n, l, ierr
      stop
   end if

   sign_ae = sign(1.0_dp, uu_ae(mch))
   sign_ps = sign(1.0_dp, uu_ps(mch))
end subroutine get_scattering_wavefunctions

subroutine run_test_configurations(ncnf,nacnf,lacnf,facnf,nc,nvcnf,rho,rhomod,rr,zz, &
                                   rcmax,mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
                                   lloc,vkb,evkb,srel,nvt,nat,lat,fat,eat,eatp,eaetst,etsttot)
   implicit none
   ! Input variables
   !> Number of test configurations
   integer, intent(in) :: ncnf
   !> size of log grid
   integer, intent(in) :: mmax
   !> dimension of number of projectors
   integer, intent(in) :: mxprj
   !> exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> number of core states
   integer, intent(in) :: nc
   !> index-1 of local potential
   integer, intent(in) :: lloc
   !> principal quantum number array, all configurations
   integer, intent(in) :: nacnf(30,5)
   !> angular-momenta array, all config.
   integer, intent(in) :: lacnf(30,5)
   !> number of valence states, all config.
   integer, intent(in) :: nvcnf(5)
   !> number of VKB projectors to use for each l
   integer, intent(in) :: nproj(5)
   !> reference configuration total energy
   real(dp), intent(in) :: etot
   !> pseudoatom total energy
   real(dp), intent(in) :: epstot
   !> maximum core radius for psp
   real(dp), intent(in) :: rcmax
   !> atomic number
   real(dp), intent(in) :: zz
   !> occupation number array, all configurations
   real(dp), intent(in) :: facnf(30,5)
   !> reference configurattion eigenvaluess
   real(dp), intent(in) :: ea(30)
   !> valence pseudo-charge of reference configuration
   real(dp), intent(in) :: rho(mmax)
   !> log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> unscreened semi-local pseudopotentials (plus different vloc if lloc==4)
   real(dp), intent(in) :: vpuns(mmax,5)
   !> Vanderbilt-Kleinman-Bylander projectors
   real(dp), intent(in) :: vkb(mmax,mxprj,4)
   !> VKB projector coefficients
   real(dp), intent(in) :: evkb(mxprj,4)
   !> model core charge
   real(dp), intent(in) :: rhomod(mmax,5)
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !> Output variables
   integer, intent(out) :: nvt(5)
   !> principal quantum number array for test configuration
   integer, intent(out) :: nat(30, 5)
   !> angular momentum array for test configuration
   integer, intent(out) :: lat(30, 5)
   !> occupation number array for test configuration
   real(dp), intent(out) :: fat(30, 3, 5)
   !> all-electron energy of test state
   real(dp), intent(out) :: eat(30, 3, 5)
   !> pseudo energy of test state
   real(dp), intent(out) :: eatp(30, 5)
   !> test configuration all-electron total energy
   real(dp), intent(out) :: eaetst(5)
   !> test configuration pseudoatom total energy
   real(dp), intent(out) :: etsttot(5)

   !> Local variables
   integer :: jj, ierr
   real(dp) :: rhot(mmax)

   nvt(:) = 0
   nat(:, :) = 0
   lat(:, :) = 0
   fat(:, :, :) = 0.0_dp
   eat(:, :, :) = 0.0_dp
   eatp(:, :) = 0.0_dp
   eaetst(:) = 0.0_dp
   etsttot(:) = 0.0_dp

   do jj = 1, ncnf + 1
      ! charge density is initialized to that of reference configuration
      rhot(:) = rho(:)
      call run_config(jj,nacnf,lacnf,facnf,nc,nvcnf,rhot,rhomod,rr,zz, &
                      rcmax,mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
                      lloc,vkb,evkb,srel,nvt(jj),nat(:,jj),lat(:,jj),fat(:,:,jj),eat(:,:,jj),eatp(:,jj),eaetst(jj), &
                      etsttot(jj))
   end do

end subroutine run_test_configurations

end module postprocess_m
