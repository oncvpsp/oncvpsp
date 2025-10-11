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

module modcore3_m
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: modcore3, get_modcore3_match, &
      teter_grid_search, teter_nelder_mead, &
      d2exc_rmse_objective, d2exc_iminus_rmse_objective
contains

!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore3(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
                    teter_amp,teter_scale,irps,mmax,rr,nc,nv,la,zion,iexc)

   ! icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
   ! rhops  state-by-state pseudocharge density
   ! rhotps  total pseudocharge density
   ! rhoc  core-charge density
   ! rhoae  state-by-state all-electron valence charge density
   ! rhotae  total all-electron valence charge density
   ! rhomod  model core density and 4 derivatives
   ! fcfact  prefactor for model amplitude (multiplies crossover value)
   ! rcfact  prefactor for model scale (multiplies crossover radius)
   ! irps  rr index of maximum rc
   ! mmax  dimension of log grid
   ! rr log radial grid
   ! nc  number of core states
   ! nv  number of valence states
   ! la  angular-momenta
   ! zion  ion charge
   ! iexc  exchange-correlation function to be used
   implicit none

   ! Input variables
   integer :: icmod,nv,nc,iexc,irps,mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
   real(dp) :: rhotps(mmax),rhoc(mmax),rr(mmax)
   real(dp) :: zion,teter_amp,teter_scale
   logical :: srel

   ! Output variables
   real(dp) :: rhomod(mmax,5)

   ! Local variables
   real(dp) :: al,eeel,eexc
   real(dp) :: d2exc_rmse_rhom,rmatch,rhocmatch,r0,rcross
   real(dp) :: gg,tt,yy
   real(dp) :: drint,rtst,rint(20),fint(20) !ad-hoc smoothing variables
   integer :: ii,ierr,ircc,ircross,irmod,iter,jj,kk
   integer :: iint !ad-hoc smoothing variables
   real(dp) :: xx(2,3)

   call get_modcore3_match(mmax, rr, rhoc, rhotps, ircc, rmatch, rhocmatch)

   gg=0.d0
   yy=0.d0

   xx(1,1)=teter_amp
   xx(2,1)=teter_scale

   r0=1.5d0*xx(2,1)
   do ii=mmax,1,-1
      if(rr(ii)<r0) then
         call gg1cc(gg,yy)
         if( xx(1,1)*gg<rhoc(ii)) then
            rcross=rr(ii)
            ircross=ii
            exit
         end if
      end if
   end do

   !  blend the Teter function tail into the all-electron rhoc
   !  first two derivatives are filled in analytically before the blend starts
   rhomod(:,:)=0.0d0
   do ii=1,mmax
      yy=rr(ii)/xx(2,1)
      call gg1cc(gg,yy)
      tt=(rr(ii)-r0-2.0d0*rr(ircross))/(r0-rr(ircross))

      rhomod(ii,1)=xx(1,1)*gg
      if(ii<ircross) then
         call gp1cc(gg,yy)
         rhomod(ii,2)=xx(1,1)*gg/xx(2,1)
         call gpp1cc(gg,yy)
         rhomod(ii,3)=xx(1,1)*gg/xx(2,1)**2
      end if
   end do

   ! 7-point numerical first derivatives applied successively
   ! skip non-blended section for 1st and 2nd derivatives

   al = 0.01d0 * dlog(rr(101) / rr(1))

   do jj=2,3
      do ii=ircross-6,mmax-3
         rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
            &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
            &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
            &     /(60.d0*al*rr(ii))
      end do
   end do

   ! ad-hoc treatment of numerical noise near origin
   ! set up a mesh on which 2nd derivative will have a stable
   ! polynomial representation
   ! assumes dpnint remains 7th order
   drint=0.02d0*rr(ircross)
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,3)
            fint(5-jj)=rhomod(ii,3)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,3),iint-1)

   jj=4
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
         &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
         &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
         &     /(60.d0*al*rr(ii))
   end do

   ! set up a mesh on which 3rd derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,4)
            fint(5-jj)=-rhomod(ii,4)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,4),iint-1)

   jj=5
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
         &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
         &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
         &     /(60.d0*al*rr(ii))
   end do

   ! set up a mesh on which 4th derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,5)
            fint(5-jj)=rhomod(ii,5)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,5),iint-1)

   return
end subroutine modcore3

subroutine get_modcore3_match(mmax, rr, rhoc, rhotps, ircc, rmatch, rhocmatch)
   implicit none
   ! Input variables
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> True all-electron core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)

   ! Output variables
   !> Index of crossover point
   integer, intent(out) :: ircc
   !> Crossover radius
   real(dp), intent(out) :: rmatch
   !> Core charge at crossover point
   real(dp), intent(out) :: rhocmatch

   ! Local variables
   integer :: ii

   !  find valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax,1,-1
      if(rhoc(ii) .gt. rhotps(ii)) then
         ircc=ii
         rmatch=rr(ircc)
         rhocmatch=rhoc(ircc)
         exit
      end if
   end do

   if(ircc .eq. 0) then
      write(stdout,'(/a)') 'modcore3: ERROR ircc (core-valence charge crossover) &
         &        not found'
      stop
   end if

end subroutine get_modcore3_match

subroutine teter_grid_search(mmax, nv, nc, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             objective, &
                             n_teter_amp, teter_amp_params, &
                             n_teter_scale, teter_scale_params, &
                             d2excae, d2excps, &
                             d2exc_rmse_grid, opt_amp_param, opt_scale_param, d2exc_rmse_min, rhomod)
   implicit none

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
         integer, intent(in) :: o_la(o_nc + o_nv)
         real(dp), intent(in) :: o_rr(o_mmax)
         real(dp), intent(in) :: o_rhoae(o_mmax, o_nc + o_nv)
         real(dp), intent(in) :: o_rhotae(o_mmax)
         real(dp), intent(in) :: o_rhops(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotps(o_mmax)
         real(dp), intent(in) :: o_rhocore(o_mmax)
         real(dp), intent(in) :: o_amp
         real(dp), intent(in) :: o_scale
         real(dp) :: oo
      end function objective
   end interface

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
   !> Number of Teter amplitude parameters to scan
   integer, intent(in) :: n_teter_amp
   !> Teter amplitude parameters to scan
   real(dp), intent(in) :: teter_amp_params(n_teter_amp)
   !> Number of Teter scale parameters to scan
   integer, intent(in) :: n_teter_scale
   !> Teter scale parameters to scans
   real(dp), intent(in) :: teter_scale_params(n_teter_scale)
   !> d2Exc for all-electron atom
   real(dp), intent(in) :: d2excae(nv, nv)

   ! Output variables
   !> d2Exc RMSE grid for all parameter combinations
   real(dp), intent(out) :: d2exc_rmse_grid(n_teter_amp, n_teter_scale)
   !> Optimal Teter amplitude parameter
   real(dp), intent(out) :: opt_amp_param
   !> Optimal Teter scale parameter
   real(dp), intent(out) :: opt_scale_param
   !> Minimum d2Exc RMSE found
   real(dp), intent(out) :: d2exc_rmse_min
   !> Model core charge for optimal parameters
   real(dp), intent(out) :: rhomod(mmax)

   ! Local variables
   !> d2Exc atom with model core charge built with current parameters
   real(dp) :: d2excps(nv, nv)
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

   d2exc_rmse_min = huge(0.0_dp)
   do kk = 1, n_teter_scale
      scale_param = teter_scale_params(kk)
      do jj = 1, n_teter_amp
         amp_param = teter_amp_params(jj)
         d2exc_rmse = objective(iexc, zion, mmax, nc, nv, la, rr, &
                                rhoae, rhotae, rhops, rhotps, rhocore, &
                                amp_param, scale_param)
         d2exc_rmse_grid(jj, kk) = d2exc_rmse

         if (d2exc_rmse < d2exc_rmse_min) then
            opt_amp_param = amp_param
            opt_scale_param = scale_param
            d2exc_rmse_min = d2exc_rmse
         end if

      end do
   end do
end subroutine teter_grid_search

subroutine teter_nelder_mead(init_amp_param, init_scale_param, rhocmatch, rmatch, &
                             mmax, nc, nv, zion, iexc, la, rr, &
                             rhotae, rhoae, rhotps, rhops, rhocore, &
                             objective, &
                             d2excae, d2exc_rmse_no_rhom, &
                             opt_amp_param, opt_scale_param, rhomod, iter)
   implicit none

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
         integer, intent(in) :: o_la(o_nc + o_nv)
         real(dp), intent(in) :: o_rr(o_mmax)
         real(dp), intent(in) :: o_rhoae(o_mmax, o_nc + o_nv)
         real(dp), intent(in) :: o_rhotae(o_mmax)
         real(dp), intent(in) :: o_rhops(o_mmax, o_nv)
         real(dp), intent(in) :: o_rhotps(o_mmax)
         real(dp), intent(in) :: o_rhocore(o_mmax)
         real(dp), intent(in) :: o_amp
         real(dp), intent(in) :: o_scale
         real(dp) :: oo
      end function objective
   end interface

   ! Input variables
   !> Initial guess for amplitude parameter
   real(dp), intent(in) :: init_amp_param
   !> Initial guess for scale parameter
   real(dp), intent(in) :: init_scale_param
   !> Core charge at crossover point
   real(dp), intent(in) :: rhocmatch
   !> Crossover radius
   real(dp), intent(in) :: rmatch
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
   !> d2Exc for all-electron atom
   real(dp), intent(in) :: d2excae(nv, nv)
   !> d2Exc RMSE without model core charge
   real(dp), intent(in) :: d2exc_rmse_no_rhom

   ! Output variables
   !> Optimized amplitude parameter
   real(dp), intent(out) :: opt_amp_param
   !> Optimized scale parameter
   real(dp), intent(out) :: opt_scale_param
   !> Model core charge density
   real(dp), intent(out) :: rhomod(mmax)
   !> Number of iterations taken
   integer, intent(out) :: iter


   ! Local variables
   !> Teter function output
   real(dp) :: gg
   !> Teter function argument
   real(dp) :: yy
   !> Starting radius for crossing search
   real(dp) :: r0
   !> Crossover radius found in search
   real(dp) :: rcross
   !>
   real(dp) :: d2excps(nv, nv)
   !>
   real(dp) :: d2exc_rmse_rhom
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

   ! Initial Nelder-Mead simplex from coarse-search minimum
   xx(1,1)=init_amp_param
   xx(2,1)=init_scale_param
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
      if(ff(3)-ff(1)<1.0d-4*d2exc_rmse_no_rhom) then
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
   integer, intent(in) :: la(nc + nv)
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhoae(mmax, nc + nv)
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

function d2exc_iminus_rmse_objective(iexc, zion, mmax, nc, nv, la, rr, &
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
   integer, intent(in) :: la(nc + nv)
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhoae(mmax, nc + nv)
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
end function d2exc_iminus_rmse_objective

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

end module modcore3_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
