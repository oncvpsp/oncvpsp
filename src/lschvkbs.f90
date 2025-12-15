!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
 subroutine lschvkbs(ll,ivkb,ee,rr,vloc,vkb,evkb,uu,up,mmax,mch)

! Normalized scattering state for pseudopotential, fully non-local unless 
! ivkb=0, which should only happen if ll=lloc

!ll  angular-momentum quantum number
!ivkb  = 0, 1 or 2 VKB proectors to be used
!ee  scattering state energy
!rr  log radial mesh
!vloc  local pseudopotential
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch  index of radius to which uu is computed

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: mmax,mch
 real(dp) :: rr(mmax),vloc(mmax),vkb(mmax,2),evkb(2)
 real(dp) :: ee
 integer :: ivkb,ll
 
!Output variables
 real(dp) :: uu(mmax),up(mmax)


!Local variables
 real(dp) :: amesh,al
 real(dp) :: cn,ro,sn
 integer :: ii,node


 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0

 call vkboutwf(ll,ivkb,ee,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)

!perform normalization sum

 ro=rr(1)/dsqrt(amesh)
 sn=ro**(2*ll+3)/(2*ll+3)

 do ii=1,mch-3
   sn=sn+al*rr(ii)*uu(ii)**2
 end do

 sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 &
&           + 28.0d0*rr(mch-1)*uu(mch-1)**2 &
&          +  9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0

!normalize u

 cn=1.0d0/dsqrt(sn)

 do ii=1,mch
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mch+1,mmax
   uu(ii)=0.0d0
 end do

 return
 end  subroutine lschvkbs
