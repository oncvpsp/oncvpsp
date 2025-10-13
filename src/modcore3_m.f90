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
   public :: modcore3, get_modcore3_match
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
   real(dp) :: d2exc_rmse_rhom,r0,rcross
   real(dp) :: gg,tt,yy
   real(dp) :: drint,rtst,rint(20),fint(20) !ad-hoc smoothing variables
   integer :: ii,ierr,ircc,ircross,irmod,iter,jj,kk
   integer :: iint !ad-hoc smoothing variables
   real(dp) :: xx(2,3)

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

end module modcore3_m
