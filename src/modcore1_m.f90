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

module modcore1_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: modcore1, get_modcore1_match
contains

!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore1(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
   &                   fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc,iter)

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

   integer, parameter :: dp=kind(1.0d0)

   ! Input variables
   integer :: icmod,nv,nc,iexc,irps,mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
   real(dp) :: rhotps(mmax),rhoc(mmax),rr(mmax)
   real(dp) :: zion,fcfact,rcfact
   logical :: srel

   ! Output variables
   real(dp) :: rhomod(mmax,5)
   integer :: iter

   ! convergence criterion
   real(dp), parameter :: eps=1.0d-7

   ! Local variables
   real(dp) :: a0,al,et,yy,gg,a0min,a0max,dermax,psum,sf,eeel,eexc
   real(dp) :: d2mdiff,rmatch,rhocmatch
   real(dp) :: emin,emax,sls,ss,hh,rgg,xx
   real(dp) :: aco(5),polym(5,5),work(5,5),constm(5,5),xpow(9),fmatch(5)
   integer :: ii,ierr,ircc,irmod,jj,kk,ll,l1,mch
   integer :: ipvt(5),nodes(4)

   call get_modcore1_match(irps,fcfact,mmax,rhoc,rhotps,ircc,irmod)

   al = 0.01d0 * dlog(rr(101) / rr(1))

   !  core charge density for Louie-Froyen-Cohen correction
   rhomod(:,1)= rhoc(:)

   xx=rr(ircc)
   fmatch(1)=rhoc(ircc)

   ! core charge derivatives
   do jj=2,5
      ! 7-point numerical first derivatives applied successively
      do ii=ircc-25+3*jj,mmax-3
         rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
                        -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
                        -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
            /(60.d0*al*rr(ii))
      end do
      fmatch(jj)=rhomod(ircc,jj)
   end do

   ! constants for polynomial matrix
   do jj=1,5
      constm(1,jj)=1.0d0
   end do
   do jj=1,5
      constm(2,jj)=jj+2
   end do
   do ii=3,5
      do jj=1,5
         constm(ii,jj)=constm(ii-1,jj)*(constm(2,jj)-ii+2)
      end do
   end do

   ! powers of xx=rr(ircc)
   xpow(1)=0.0d0
   xpow(2)=1.0d0
   do ii=3,9
      xpow(ii)=xx*xpow(ii-1)
   end do

   ! polynomial matrix
   do jj=1,5
      do ii=1,5
         polym(ii,jj)=constm(ii,jj)*xpow(jj-ii+5)
      end do
   end do


   ! begin iteration for monotonic polynomial match
   a0=fmatch(1)
   a0min=a0
   a0max=0.0d0

   do iter=1,50

      work(:,:)=polym(:,:)

      ! fill target value vector
      aco(1)=fmatch(1)-a0
      do jj=2,5
         aco(jj)=fmatch(jj)
      end do

      ! solve linear equations for coefficients

      call dgesv(5,1,work,5,ipvt,aco,5,kk)
      if(kk/= 0) then
         if(kk>0) write(6,'(a,i4)') &
            &      'modcore:ERROR stop - singular polym matrix',kk
         if(kk<0) write(6,'(a,i4)') &
            &      'modcore:ERROR stop - dgesv input error',kk
         stop
      end if

      ! find maximum derivative
      ! note that xpow end up properly re-written for continued iteation
      dermax=-1.0d20
      do kk=1,ircc
         do ii=3,9
            xpow(ii)=rr(kk)*xpow(ii-1)
         end do
         ii=2
         psum=0.0d0
         do jj=1,5
            psum=psum+constm(ii,jj)*xpow(jj-ii+5)*aco(jj)
         end do
         dermax=dmax1(dermax,psum)
      end do

      ! test maximum derivative and adjust a0
      ! interval-halving search after achieving monotonicity to get
      ! barely monotonic model
      if(dermax .gt. 0.0d0) then
         a0min=a0
      else
         a0max=a0
      end if

      ! success when maximum derivative bracketed just above zero
      if(dermax>0.0d0 .and. dermax<eps*dabs(fmatch(2))) exit

      if(a0max==0.0d0) then
         a0=2.5d0*a0
      else
         a0=0.5d0*(a0max+a0min)
      end if

   end do  !iterr



   ! fill in model and derivatives with fitted polynomial
   do kk=1,ircc-1
      do ii=3,9
         xpow(ii)=rr(kk)*xpow(ii-1)
      end do
      do ii=1,5
         if(ii .eq. 1) then
            psum=a0
         else
            psum=0.0d0
         end if
         do jj=1,5
            psum=psum+constm(ii,jj)*xpow(jj-ii+5)*aco(jj)
         end do
         rhomod(kk,ii)=psum
      end do
   end do
end subroutine modcore1

subroutine get_modcore1_match(irps,fcfact,mmax,rhoc,rhotps,ircc,irmod)
   integer, intent(in) :: irps
   real(dp), intent(in) :: fcfact
   integer, intent(in) :: mmax
   real(dp), intent(in) :: rhoc(mmax)
   real(dp), intent(in) :: rhotps(mmax)

   integer, intent(out) :: ircc
   integer, intent(out) :: irmod

   integer :: ii
   ircc = 0

   do ii = mmax,1,-1
      if(rhoc(ii) .gt. fcfact*rhotps(ii)) then
         ircc=ii
         exit
      end if
   end do
   if(ircc .eq. 0) then
      write(6,'(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
         &        not found'
      stop
   end if

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod=max(irps,ircc)
end subroutine get_modcore1_match

end module modcore1_m
