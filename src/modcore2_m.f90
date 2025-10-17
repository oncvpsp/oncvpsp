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

module modcore2_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: modcore2, get_modcore2_match
contains

!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore2(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
                    fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

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
   real(dp) :: zion,fcfact,rcfact
   logical :: srel

   ! Output variables
   real(dp) :: rhomod(mmax,5)
   real(dp) :: a0
   real(dp) :: b0

   ! convergence criterion
   real(dp), parameter :: EPS=1.0d-7

   ! Local variables
   real(dp) :: al,eeel,eexc
   real(dp) :: xx,yy,dy
   real(dp) :: x0max,x0min,r0,tt,ymatch,ytrial
   real(dp) :: fmatch(5)
   real(dp) :: drint,rtst,rint(20),fint(20) !ad-hoc smoothing variables
   integer :: ii,ierr,ircc,iter,jj,kk
   integer :: iint !ad-hoc smoothing variables

   ! find scaled valence pseudocharge - core charge crossover
   call get_modcore2_match(mmax, rr, rhoc, rhotps, fcfact, ircc, a0, b0)

   rhomod(:,:)=0.0d0
   do ii=1,mmax
      xx=b0*rr(ii)
      if(xx<3.0d0) then
         call gg1cc(yy,xx)
         rhomod(ii,1)=a0*yy
         call gp1cc(yy,xx)
         rhomod(ii,2)=a0*yy*b0
         call gpp1cc(yy,xx)
         rhomod(ii,3)=a0*yy*b0**2
      end if
   end do

   ! 7-point numerical first derivatives applied successively
   ! do jj=4,5
   !  do ii=3* j j-8,mmax-3
   !    if(rhomod(ii,1)== 0 .0d0) exit
   !     rhomod(ii,jj)=(-rhomod(ii-3 , jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
   ! &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
   ! &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
   ! &     /(60.d0*al*rr(ii))
   !  end do
   ! end do

   ! ad-hoc treatment of numerical noise near origin
   ! set up a mesh on which 2nd derivative will have   a stable
   ! polynomial representation
   ! assumes dpnint remains 7t h  order
   drint=0.05d0*rr(ircc)
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircc
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
      do ii=1,ircc
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
      do ii=1,ircc
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
end subroutine modcore2

subroutine get_modcore2_match(mmax,rr,rhoc,rhotps,fcfact,ircc,a0,b0)
   implicit none
   ! Constants
   real(dp), parameter :: EPS=1.0d-7

   ! Input variables
   integer, intent(in) :: mmax
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rhoc(mmax)
   real(dp), intent(in) :: rhotps(mmax)
   real(dp), intent(in) :: fcfact

   ! Output variables
   integer, intent(out) :: ircc
   real(dp), intent(out) :: a0
   real(dp), intent(out) :: b0

   ! Local variables
   integer :: ii,jj
   real(dp) :: rhomod(mmax,5)
   real(dp) :: al,xx,yy,dy
   real(dp) :: x0max,x0min,tt,ymatch,ytrial
   real(dp) :: fmatch(5)

   ircc = 0
   do ii = mmax,1,-1
      if(rhoc(ii) .gt. fcfact*rhotps(ii)) then
         ircc=ii
         exit
      end if
   end do

   if(ircc .eq. 0) then
      write(6,'(/a)') 'modcore2: ERROR ircc (core-valence charge crossover) &
         &        not found'
      stop
   end if

   al = 0.01d0 * dlog(rr(101) / rr(1))

   ! core charge density for Louie-Froyen-Cohen correction
   rhomod(:,1) = rhoc(:)
   xx=rr(ircc)
   fmatch(1)=rhoc(ircc)

   ! core charge derivatives
   ! 7-point numerical first derivative
   jj=2
   ii=ircc
   rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
      &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
      &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
      &     /(60.d0*al*rr(ii))
   fmatch(jj)=rhomod(ircc,jj)

   ! Fit Teter function to value and slope at ircc
   ymatch=rr(ircc)*fmatch(2)/fmatch(1)

   ! interval-halving search for dimensionless match point
   x0max=1.48d0
   x0min=0.0d0
   do jj=1,50
      xx=0.5d0*(x0max+x0min)
      call gg1cc(yy,xx)
      call gp1cc(dy,xx)
      ytrial=xx*dy/yy
      if(abs(ytrial-ymatch)<EPS) exit
      if(ytrial<ymatch) then
         x0max=xx
      else
         x0min=xx
      end if
   end do
   b0=xx/rr(ircc)
   a0=fmatch(1)/yy
end subroutine get_modcore2_match

end module modcore2_m
