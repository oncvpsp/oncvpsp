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
subroutine vploc(rr,vv,vp,dvloc0,irc,mmax,lpopt)

!polynomial extrapolation of all-electron potential to rr=0

!rr log radial mesh
!vv all-electron  potential
!vp pseudopotentials, vp(:,5) is local to be created by this routine
!dvloc0  amplitude at rr==0 of smooth "null space" function to be added
!irc  matching radius for polynomial potential
!mmax  dimension of rr
!lpopt  1-5, determines polynomial to be used
!  1) match 2 derivatives, r**0,2,4
!  2) match 2 derivatives, r**0,4,6
!  3) match 3 derivatives, r**0,4,5,6
!  4) match 3 derivatives, r**0,4,6,8
!  5) match 3 derivatives, r**0,2,4,6

   implicit none
   integer, parameter :: dp=kind(1.0d0)

! Input variables
   integer :: irc,mmax,lpopt
   real(dp) :: dvloc0
   real(dp) :: rr(mmax),vv(mmax)

! Output variables
   real(dp) :: vp(mmax,5)

!Local variables
   integer :: ii
   real(dp) :: aco,al,bco,cco,dco,d3vv,x

   real(dp), allocatable :: dvv(:),d2vv(:)

   allocate(dvv(mmax),d2vv(mmax))

   al = 0.01d0 * dlog(rr(101) / rr(1))

! derivatives of vv

   do ii = irc - 4, irc + 4
      dvv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
      &         -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
   end do


   do ii = irc - 2, irc + 2
      d2vv(ii)=(2.d0*dvv(ii-2)-16.d0*dvv(ii-1)+16.d0*dvv(ii+1) &
      &         -2.d0*dvv(ii+2))/(24.d0*al*rr(ii))
   end do

   ii = irc
   d3vv=(2.d0*d2vv(ii-2)-16.d0*d2vv(ii-1)+16.d0*d2vv(ii+1) &
   &       -2.d0*d2vv(ii+2))/(24.d0*al*rr(ii))


   if(lpopt == 1) then

      bco=(3*dvv(ii)-d2vv(ii)*rr(ii))/(4*rr(ii))
      cco=( -dvv(ii)+d2vv(ii)*rr(ii))/(8*rr(ii)**3)
      aco=vv(ii)-bco*rr(ii)**2-cco*rr(ii)**4

   else if (lpopt == 2) then

      bco=( 5*dvv(ii)-d2vv(ii)*rr(ii))/(8*rr(ii)**3)
      cco=(-3*dvv(ii)+d2vv(ii)*rr(ii))/(12*rr(ii)**5)
      aco=vv(ii)-bco*rr(ii)**4-cco*rr(ii)**6

   else if(lpopt == 3) then

      bco=( 20*dvv(ii)-8*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(8*rr(ii)**3)
      cco=(-15*dvv(ii)+7*d2vv(ii)*rr(ii)-d3vv*rr(ii)**2)/(5*rr(ii)**4)
      dco=( 12*dvv(ii)-6*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(12*rr(ii)**5)
      aco = vv(ii)-bco*rr(ii)**4-cco*rr(ii)**5-dco*rr(ii)**6

   else if(lpopt == 4) then

      bco=(35*dvv(ii)-11*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(32*rr(ii)**3)
      cco=(-21*dvv(ii)+9*d2vv(ii)*rr(ii)-d3vv*rr(ii)**2)/(24*rr(ii)**5)
      dco=(15*dvv(ii)-7*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(64*rr(ii)**7)
      aco = vv(ii)-bco*rr(ii)**4-cco*rr(ii)**6-dco*rr(ii)**8

   else if(lpopt == 5) then

      bco=(15*dvv(ii)-7*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(16*rr(ii))
      cco=(-5*dvv(ii)+5*d2vv(ii)*rr(ii)-d3vv*rr(ii)**2)/(16*rr(ii)**3)
      dco=( 3*dvv(ii)-3*d2vv(ii)*rr(ii)+d3vv*rr(ii)**2)/(48*rr(ii)**5)
      aco = vv(ii)-bco*rr(ii)**2-cco*rr(ii)**4-dco*rr(ii)**6

   end if

! create polynomial potential inside irc


   if(lpopt == 1) then
      do ii=1,irc
         vp(ii,5)=aco+bco*rr(ii)**2+cco*rr(ii)**4
      end do
   else if(lpopt == 2) then
      do ii=1,irc
         vp(ii,5)=aco+bco*rr(ii)**4+cco*rr(ii)**6
      end do
   else if(lpopt == 3) then
      do ii=1,irc
         vp(ii,5)=aco+bco*rr(ii)**4+cco*rr(ii)**5+dco*rr(ii)**6
      end do
   else if(lpopt == 4) then
      do ii=1,irc
         vp(ii,5)=aco+bco*rr(ii)**4+cco*rr(ii)**6+dco*rr(ii)**8
      end do
   else if(lpopt == 5) then
      do ii=1,irc
         vp(ii,5)=aco+bco*rr(ii)**2+cco*rr(ii)**4+dco*rr(ii)**6
      end do
   end if

! add compatible supplementary function to change value at rr==0

   do ii = 1, irc
      x=rr(ii)/rr(irc)
      if(lpopt == 1) then
         vp(ii,5)=vp(ii,5) + dvloc0*(1.0d0 - x**2)**3
      else if(lpopt == 2) then
         vp(ii,5)=vp(ii,5) + dvloc0*(1.0d0 - x**4)**3
      else if(lpopt == 3 .or. lpopt == 4) then
         vp(ii,5)=vp(ii,5) + dvloc0*(1.0d0 - x**4)**4
      else if(lpopt == 5) then
         vp(ii,5)=vp(ii,5) + dvloc0*(1.0d0 - x**4)**4
      end if
   end do

   do ii=irc+1,mmax
      vp(ii,5)=vv(ii)
   end do

   deallocate(dvv,d2vv)

   return
end subroutine vploc
