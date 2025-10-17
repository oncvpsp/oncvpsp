!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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

module modcore4_m
   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use utils_m, only: nonuniform_trapezoid, fermi_dirac
   use modcore3_m, only: modcore3, get_modcore3_match
   implicit none
   private
   public :: modcore4

   interface
      function objective_f(iexc, zion, mmax, nc, nv, la, rr, &
                           rhoae, rhotae, rhops, rhotps, rhocore, &
                           amp, scale) result(oo)
         use, intrinsic :: iso_fortran_env, only: dp => real64
         integer, intent(in) :: iexc
         real(dp), intent(in) :: zion
         integer, intent(in) :: mmax
         integer, intent(in) :: nc
         integer, intent(in) :: nv
         integer, intent(in) :: la(30)
         real(dp), intent(in) :: rr(mmax)
         real(dp), intent(in) :: rhoae(mmax, nv)
         real(dp), intent(in) :: rhotae(mmax)
         real(dp), intent(in) :: rhops(mmax, nv)
         real(dp), intent(in) :: rhotps(mmax)
         real(dp), intent(in) :: rhocore(mmax)
         real(dp), intent(in) :: amp
         real(dp), intent(in) :: scale
         real(dp) :: oo
      end function objective_f
   end interface

contains

subroutine modcore4(mmax,rr,nc,nv,la,zion,irps,iexc, &
                    rhops,rhotps,rhocore,rhoae,rhotae, &
                    objective_name, &
                    n_amp, n_scale, amp_params, scale_params, &
                    grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, &
                    nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, nm_iter, &
                    rhomod)
   implicit none
   ! Input parameters
   integer, intent(in) :: mmax
   integer, intent(in) :: nc
   integer, intent(in) :: nv
   integer, intent(in) :: la(nc + nv)
   real(dp), intent(in) :: zion
   integer, intent(in) :: irps
   integer, intent(in) :: iexc
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhops(mmax, nv)
   real(dp), intent(in) :: rhotps(mmax)
   real(dp), intent(in) :: rhocore(mmax)
   real(dp), intent(in) :: rhoae(mmax, nc + nv)
   real(dp), intent(in) :: rhotae(mmax)
   character(len=*), intent(in) :: objective_name
   integer, intent(in) :: n_amp
   integer, intent(in) :: n_scale
   real(dp), intent(in) :: amp_params(n_amp)
   real(dp), intent(in) :: scale_params(n_scale)

   ! Output parameters
   real(dp), intent(out) :: grid_objective(n_amp, n_scale)
   real(dp), intent(out) :: grid_opt_amp_param
   real(dp), intent(out) :: grid_opt_scale_param
   real(dp), intent(out) :: grid_opt_objective
   real(dp), intent(out) :: nm_opt_amp_param
   real(dp), intent(out) :: nm_opt_scale_param
   real(dp), intent(out) :: nm_opt_objective
   integer, intent(out) :: nm_iter
   real(dp), intent(out) :: rhomod(mmax, 5)

   ! Local variables
   integer :: ircc
   real(dp) :: rmatch
   real(dp) :: rhocmatch
   real(dp) :: nm_atol
   ! procedure(objective_f), pointer :: objective

   call get_modcore3_match(mmax, rr, rhocore, rhotps, ircc, rmatch, rhocmatch)

   select case(objective_name)
    case ('d2exc_rmse')
      ! objective => d2exc_rmse_objective
      nm_atol = 1.0e-4_dp * d2exc_rmse_objective(iexc, zion, mmax, nc, nv, la, rr, &
                                                 rhoae, rhotae, rhops, rhotps, rhocore, &
                                                 0.0_dp, 1.0_dp)
      call teter_grid_search(mmax, nv, nc, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             d2exc_rmse_objective, &
                             n_amp, amp_params, &
                             n_scale, scale_params, &
                             grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, rhomod(:, 1))

      call teter_nelder_mead(grid_opt_amp_param, grid_opt_scale_param, rhocmatch, rmatch, nm_atol, &
                             mmax, nc, nv, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             d2exc_rmse_objective, &
                             nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, rhomod(:, 1), nm_iter)
    case ('d2exc_rmse_iminus')
      ! objective => d2exc_rmse_iminus_objective
      ! By definition, I- == 0 if rhomod = 0, so using a relative tolerance w.r.t. the objective
      ! evaluated without a model core charge is not meaningful.
      ! The objective as a whole is on a log scale, so a small absolute tolerance corresponds to
      ! a very small tolerance
      nm_atol = 1.0e-4_dp
      call teter_grid_search(mmax, nv, nc, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             d2exc_rmse_iminus_objective, &
                             n_amp, amp_params, &
                             n_scale, scale_params, &
                             grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, rhomod(:, 1))

      call teter_nelder_mead(grid_opt_amp_param, grid_opt_scale_param, rhocmatch, rmatch, nm_atol, &
                             mmax, nc, nv, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             d2exc_rmse_iminus_objective, &
                             nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, rhomod(:, 1), nm_iter)
    case default
      error stop 'modcore4: ERROR unknown objective_name = ' // trim(objective_name)
   end select

   ! call teter_grid_search(mmax, nv, nc, zion, iexc, la, rr, &
   !                        rhotae, rhoae, rhotps, rhops, rhocore, &
   !                        objective, &
   !                        n_amp, amp_params, &
   !                        n_scale, scale_params, &
   !                        grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, rhomod(:, 1))

   ! call teter_nelder_mead(grid_opt_amp_param, grid_opt_scale_param, rhocmatch, rmatch, nm_atol, &
   !                        mmax, nc, nv, zion, iexc, la, rr, &
   !                        rhotae, rhoae, rhotps, rhops, rhocore, &
   !                        objective, &
   !                        nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, rhomod(:, 1), nm_iter)

   call modcore3(3,rhops,rhotps,rhocore,rhoae,rhotae,rhomod, &
                 nm_opt_amp_param,nm_opt_scale_param,irps,mmax,rr,nc,nv,la,zion,iexc)
end subroutine modcore4

subroutine teter_grid_search(mmax, nv, nc, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             objective, &
                             n_teter_amp, teter_amp_params, &
                             n_teter_scale, teter_scale_params, &
                             objective_grid, opt_amp_param, opt_scale_param, opt_objective, rhomod)
   implicit none
   ! Input variables
   !> Size of log radial mesh
   integer, intent(in) :: mmax
   !> Number of valence states
   integer, intent(in) :: nv
   !> Number of core states
   integer, intent(in) :: nc
   !> Pseudo atomic charge
   real(dp), intent(in) :: zion
   !> Exchange-correlation functional index
   integer, intent(in) :: iexc
   !> Angular momenta
   integer, intent(in) :: la(nv+nc)
   !> Log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> All-electron total charge
   real(dp), intent(in) :: rhotae(mmax)
   !> All-electron state-by-state charge
   real(dp), intent(in) :: rhoae(mmax,nv)
   !> Pseudo total charge
   real(dp), intent(in) :: rhotps(mmax)
   !> Pseudo state-by-state charge
   real(dp), intent(in) :: rhops(mmax,nv)
   !> True all-electron core charge
   real(dp), intent(in) :: rhocore(mmax)
   ! !> Objective function to minimize
   ! procedure(objective_f) :: objective
   !> Number of Teter amplitude parameters to scan
   integer, intent(in) :: n_teter_amp
   !> Teter amplitude parameters to scan
   real(dp), intent(in) :: teter_amp_params(n_teter_amp)
   !> Number of Teter scale parameters to scan
   integer, intent(in) :: n_teter_scale
   !> Teter scale parameters to scans
   real(dp), intent(in) :: teter_scale_params(n_teter_scale)

   ! Output variables
   !> d2Exc RMSE grid for all parameter combinations
   real(dp), intent(out) :: objective_grid(n_teter_amp, n_teter_scale)
   !> Optimal Teter amplitude parameter
   real(dp), intent(out) :: opt_amp_param
   !> Optimal Teter scale parameter
   real(dp), intent(out) :: opt_scale_param
   !> Minimum d2Exc RMSE found
   real(dp), intent(out) :: opt_objective
   !> Model core charge for optimal parameters
   real(dp), intent(out) :: rhomod(mmax)

   ! Local variables
   !> Current amplitude parameter
   real(dp) :: amp_param
   !> Current scale parameter
   real(dp) :: scale_param
   !> d2Exc RMSE for current parameters
   real(dp) :: d2exc_rmse
   !> Teter function input
   real(dp) :: xx
   !> Teter function output
   real(dp) :: yy
   !> Loop indices
   integer :: ii, jj, kk

   interface
      function objective(o_iexc, o_zion, o_mmax, o_nc, o_nv, o_la, o_rr, &
                         o_rhoae, o_rhotae, o_rhops, o_rhotps, o_rhocore, &
                         o_amp, o_scale) result(oo)
         use, intrinsic :: iso_fortran_env, only: dp => real64
         integer, intent(in) :: o_iexc
         real(dp), intent(in) :: o_zion
         integer, intent(in) :: o_mmax
         integer, intent(in) :: o_nc
         integer, intent(in) :: o_nv
         integer, intent(in) :: o_la(30)
         real(dp), intent(in) :: o_rr(o_mmax)
         real(dp), intent(in) :: o_rhoae(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotae(o_mmax)
         real(dp), intent(in) :: o_rhops(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotps(o_mmax)
         real(dp), intent(in) :: o_rhocore(o_mmax)
         real(dp), intent(in) :: o_amp
         real(dp), intent(in) :: o_scale
         real(dp) :: oo
      end function objective
   end interface

   opt_objective = huge(0.0_dp)
   do kk = 1, n_teter_scale
      scale_param = teter_scale_params(kk)
      do jj = 1, n_teter_amp
         amp_param = teter_amp_params(jj)
         d2exc_rmse = objective(iexc, zion, mmax, nc, nv, la, rr, &
                                rhoae, rhotae, rhops, rhotps, rhocore, &
                                amp_param, scale_param)
         objective_grid(jj, kk) = d2exc_rmse

         if (d2exc_rmse < opt_objective) then
            opt_amp_param = amp_param
            opt_scale_param = scale_param
            opt_objective = d2exc_rmse
         end if

      end do
   end do
end subroutine teter_grid_search

subroutine teter_nelder_mead(init_amp_param, init_scale_param, rhocmatch, rmatch, atol, &
                             mmax, nc, nv, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             objective, &
                             opt_amp_param, opt_scale_param, opt_objective, rhomod, iter)
   implicit none
   ! Input variables
   !> Initial guess for amplitude parameter
   real(dp), intent(in) :: init_amp_param
   !> Initial guess for scale parameter
   real(dp), intent(in) :: init_scale_param
   !> Core charge at crossover point
   real(dp), intent(in) :: rhocmatch
   !> Crossover radius
   real(dp), intent(in) :: rmatch
   !> Absolute convergence tolerance
   real(dp), intent(in) :: atol
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> Ion charge
   real(dp), intent(in) :: zion
   !> Exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> Angular momenta
   integer, intent(in) :: la(30)
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Total all-electron charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> State-by-state all-electron charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> Total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> State-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> True all-electron core charge
   real(dp), intent(in) :: rhocore(mmax)
   ! !> Objective function to minimize
   ! procedure(objective_f) :: objective

   ! Output variables
   !> Optimized amplitude parameter
   real(dp), intent(out) :: opt_amp_param
   !> Optimized scale parameter
   real(dp), intent(out) :: opt_scale_param
   !> Optimized objective function value
   real(dp), intent(out) :: opt_objective
   !> Model core charge density
   real(dp), intent(out) :: rhomod(mmax)
   !> Number of iterations taken
   integer, intent(out) :: iter


   ! Local variables
   !> Teter function output
   real(dp) :: gg
   !> Teter function argument
   real(dp) :: yy
   !> Loop indices
   integer :: ii, jj, kk
   ! 2-dimensional Nelder-Mead variables
   real(dp), parameter :: ALPHA=1.0_dp
   real(dp), parameter :: GAMMA=2.0_dp
   real(dp), parameter :: RHO=-0.5_dp
   real(dp), parameter :: SIGMA=0.5_dp
   integer, parameter :: MAX_ITER=101
   real(dp) :: xx(2,3),ff(3),xt(2),ft,x0(2),xr(2),fr,xe(2),fe,xc(2),fc
   real(dp) :: d2ref,fta(10),d2exc_rmse_min

   interface
      function objective(o_iexc, o_zion, o_mmax, o_nc, o_nv, o_la, o_rr, &
                         o_rhoae, o_rhotae, o_rhops, o_rhotps, o_rhocore, &
                         o_amp, o_scale) result(oo)
         use, intrinsic :: iso_fortran_env, only: dp => real64
         integer, intent(in) :: o_iexc
         real(dp), intent(in) :: o_zion
         integer, intent(in) :: o_mmax
         integer, intent(in) :: o_nc
         integer, intent(in) :: o_nv
         integer, intent(in) :: o_la(30)
         real(dp), intent(in) :: o_rr(o_mmax)
         real(dp), intent(in) :: o_rhoae(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotae(o_mmax)
         real(dp), intent(in) :: o_rhops(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotps(o_mmax)
         real(dp), intent(in) :: o_rhocore(o_mmax)
         real(dp), intent(in) :: o_amp
         real(dp), intent(in) :: o_scale
         real(dp) :: oo
      end function objective
   end interface

   ! Initial Nelder-Mead simplex from coarse-search minimum
   xx(1,1)=init_amp_param  ! fcfact * rhocmatch
   xx(2,1)=init_scale_param  ! rcfact * rmatch
   xx(1,2)=xx(1,1)+0.25_dp*rhocmatch
   xx(2,2)=xx(2,1)
   xx(1,3)=xx(1,1)
   xx(2,3)=xx(2,1)+0.05_dp*rmatch
   xx(1,1)=xx(1,1)-0.125_dp*rhocmatch
   xx(2,1)=xx(2,1)-0.025_dp*rmatch

   gg = 0.0_dp
   yy = 0.0_dp

   ! Fill function values for initial simplex
   do kk=1,3
      xt(:)=xx(:,kk)
      ff(kk) = objective(iexc, zion, mmax, nc, nv, la, rr, &
                         rhoae, rhotae, rhops, rhotps, rhocore, &
                         xt(1), xt(2))
   end do
   ! Nelder-Mead iteration loop
   ! write(stdout,'(/a)') 'Nelder-Mead iteration'
   do kk = 1, MAX_ITER
      ! (1) Order (dumb bubble sort)
      do jj=1,3
         if(ff(3)<ff(2)) then
            ft=ff(2) ; xt(:)=xx(:,2)
            ff(2)=ff(3) ; xx(:,2)=xx(:,3)
            ff(3)=ft ; xx(:,3)=xt(:)
         end if
         if(ff(2)<ff(1)) then
            ft=ff(1) ; xt(:)=xx(:,1)
            ff(1)=ff(2) ; xx(:,1)=xx(:,2)
            ff(2)=ft ; xx(:,2)=xt(:)
         end if
      end do
      ! stopping criterion
      if(abs(ff(3)-ff(1)) < atol) then
         ! write(stdout,'(a,i4,a)') ' converged in',kk-1,' steps'
         ! write(stdout,'(a)') 'amplitude prefactor, scale prefactor'
         ! write(stdout,'(2f10.4)') xx(1,1)/rhocmatch,xx(2,1)/rmatch
         iter = kk - 1
         exit
      else if(kk==MAX_ITER) then
         ! write(stdout,'(a,i4,a)') ' WARNING: not fully converged in 100 steps'
         ! write(stdout,'(a)') 'amplitude prefactor, scale prefactor'
         ! write(stdout,'(2f10.4)') xx(1,1)/rhocmatch,xx(2,1)/rmatch
         iter = kk - 1
         exit
      end if
      ! (2) Centroid
      x0(:)=0.5_dp*(xx(:,1)+xx(:,2))
      ! (3) Reflection
      xr(:)=x0(:)+ALPHA*(x0(:)-xx(:,3))
      fr = objective(iexc, zion, mmax, nc, nv, la, rr, &
                     rhoae, rhotae, rhops, rhotps, rhocore, &
                     xr(1), xr(2))
      if(ff(1)<= fr .and. fr<ff(2)) then
         ff(3)=fr ; xx(:,3)=xr(:)
         cycle  !kk loop
      end if
      ! (4) Expansion
      if(fr<ff(1)) then
         xe(:)=x0(:)+GAMMA*(x0(:)-xx(:,3))
         fe = objective(iexc, zion, mmax, nc, nv, la, rr, &
                        rhoae, rhotae, rhops, rhotps, rhocore, &
                        xe(1), xe(2))
         if(fe<fr) then
            ff(3)=fe ; xx(:,3)=xe(:)
            cycle  !kk
         else
            ff(3)=fr ; xx(:,3)=xr(:)
            cycle  !kk
         end if
      end if
      ! (5) Contraction
      xc(:)=x0(:)+RHO*(x0(:)-xx(:,3))
      fc = objective(iexc, zion, mmax, nc, nv, la, rr, &
                     rhoae, rhotae, rhops, rhotps, rhocore, &
                     xc(1), xc(2))
      if(fc<ff(3)) then
         ff(3)=fc ; xx(:,3)=xc(:)
         cycle  !kk
      end if
      ! (6) Reduction
      do jj=2,3
         xx(:,jj)=xx(:,1)+SIGMA*(xx(:,jj)-xx(:,1))
         ff(jj) = objective(iexc, zion, mmax, nc, nv, la, rr, &
                            rhoae, rhotae, rhops, rhotps, rhocore, &
                            xx(1,jj), xx(2,jj))
      end do  !jj
   end do  !kk
   opt_amp_param = xx(1,1)
   opt_scale_param = xx(2,1)
   opt_objective = ff(1)
end subroutine teter_nelder_mead

function d2exc_rmse_objective(iexc, zion, mmax, nc, nv, la, rr, &
                              rhoae, rhotae, rhops, rhotps, rhocore, &
                              amp, scale) result(d2exc_rmse)
   implicit none
   ! Input variables
   integer, intent(in) :: iexc
   real(dp), intent(in) :: zion
   integer, intent(in) :: mmax
   integer, intent(in) :: nc
   integer, intent(in) :: nv
   integer, intent(in) :: la(30)
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhoae(mmax, nv)
   real(dp), intent(in) :: rhotae(mmax)
   real(dp), intent(in) :: rhops(mmax, nv)
   real(dp), intent(in) :: rhotps(mmax)
   real(dp), intent(in) :: rhocore(mmax)
   real(dp), intent(in) :: amp
   real(dp), intent(in) :: scale

   ! Output variable
   real(dp) :: d2exc_rmse

   ! Local variables
   real(dp):: rhomod(mmax)
   real(dp) :: d2exc_ae(nv, nv)
   real(dp) :: d2exc_ps(nv, nv)
   real(dp) :: d2exc_dummy(nv, nv)
   real(dp) :: d2exc_rmse_dummy
   real(dp) :: xx, yy
   integer :: ii

   d2exc_ps(:, :) = 0.0_dp
   d2exc_ae(:, :) = 0.0_dp

   do ii = 1, mmax
      xx = rr(ii) / scale
      call gg1cc(yy, xx)
      rhomod(ii) = amp * yy
   end do

   ! Compute d2Exc for all-electron atom
   call der2exc(rhotae, rhocore, rhoae, rr, d2exc_ae, d2exc_dummy, d2exc_rmse_dummy, &
                zion, iexc, nc, nv, la, mmax, mmax)
   ! Compute d2Exc for pseudo atom with model core charge
   ! and RMSE w.r.t. all-electron d2Exc
   call der2exc(rhotps, rhomod, rhops, rr, d2exc_ps, d2exc_ae, d2exc_rmse, &
                zion, iexc, nc, nv, la, mmax, mmax)
   return
end function d2exc_rmse_objective

function d2exc_rmse_iminus_objective(iexc, zion, mmax, nc, nv, la, rr, &
                                     rhoae, rhotae, rhops, rhotps, rhocore, &
                                     amp, scale) result(objective)
   implicit none
   ! Constants
   real(dp), parameter :: MU = -3.0_dp
   real(dp), parameter :: SIGMA = 1.0_dp
   real(dp), parameter :: IMINUS_SHIFT = 1.0e-1_dp

   ! Input variables
   integer, intent(in) :: iexc
   real(dp), intent(in) :: zion
   integer, intent(in) :: mmax
   integer, intent(in) :: nc
   integer, intent(in) :: nv
   integer, intent(in) :: la(30)
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhoae(mmax, nv)
   real(dp), intent(in) :: rhotae(mmax)
   real(dp), intent(in) :: rhops(mmax, nv)
   real(dp), intent(in) :: rhotps(mmax)
   real(dp), intent(in) :: rhocore(mmax)
   real(dp), intent(in) :: amp
   real(dp), intent(in) :: scale

   ! Output variable
   real(dp) :: objective

   ! Local variables
   real(dp) :: rhomod(mmax)
   real(dp) :: d2exc_ae(nv, nv)
   real(dp) :: d2exc_ps(nv, nv)
   real(dp) :: d2exc_dummy(nv, nv)
   real(dp) :: d2exc_rmse_dummy
   real(dp) :: d2exc_rmse
   real(dp) :: iminus
   real(dp) :: xx, yy
   integer :: ii

   d2exc_ps(:, :) = 0.0_dp
   d2exc_ae(:, :) = 0.0_dp

   do ii = 1, mmax
      xx = rr(ii) / scale
      call gg1cc(yy, xx)
      rhomod(ii) = amp * yy
   end do

   ! Compute d2Exc for all-electron atom
   call der2exc(rhotae, rhocore, rhoae, rr, d2exc_ae, d2exc_dummy, d2exc_rmse_dummy, &
                zion, iexc, nc, nv, la, mmax, mmax)
   ! Compute d2Exc for pseudo atom with model core charge
   ! and RMSE w.r.t. all-electron d2Exc
   call der2exc(rhotps, rhomod, rhops, rr, d2exc_ps, d2exc_ae, d2exc_rmse, &
                zion, iexc, nc, nv, la, mmax, mmax)

   iminus = get_iminus(mmax, rr, rhocore, rhomod)

   objective = log10(d2exc_rmse) * fermi_dirac(log10(iminus + IMINUS_SHIFT), mu=MU, sigma=SIGMA)
   return
end function d2exc_rmse_iminus_objective

function get_iminus(mmax, rr, rhocore, rhomod) result(iminus)
   implicit none
   ! Input variables
   integer, intent(in) :: mmax
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhocore(mmax)
   real(dp), intent(in) :: rhomod(mmax)

   ! Output variable
   real(dp) :: iminus

   ! Local variables
   integer :: ii
   real(dp) :: integrand(mmax)

   do ii = 1, mmax
      integrand(ii) = max(0.0_dp, rhomod(ii), -  rhocore(ii))
   end do
   iminus = nonuniform_trapezoid(mmax, rr, integrand)
   return
end function get_iminus

end module modcore4_m
