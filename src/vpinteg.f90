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
 subroutine vpinteg(gg,hh,nn,mm,ss,rr)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

! integrals that go into construction of Vanderbilt separable pseudopotential

! product of functions gg*hh goes like rr**mm at rr -> 0
! integral on usual log mesh from 1 to nn

!Input variables
 real(dp) :: gg(nn),hh(nn),rr(nn)
 integer :: nn,mm

!Output variable
 real(dp) :: ss

!Local variables
 real(dp) :: r0,amesh,al
 integer :: ii

 al = 0.01d0 * dlog(rr(101)/rr(1))
 amesh = exp(al)

 r0=rr(1)/dsqrt(amesh)
 ss=r0**(mm+1)*(gg(1)*hh(1)/rr(1)**mm)/dfloat(mm+1)

 do ii = 4, nn - 3
   ss =  ss + al*gg(ii)*hh(ii)*rr(ii)
 end do

 ss=ss + al*(23.d0*rr(nn-2)*gg(nn-2)*hh(nn-2) &
&        + 28.d0*rr(nn-1)*gg(nn-1)*hh(nn-1) &
&        +  9.d0*rr(nn  )*gg(nn  )*hh(nn  ))/24.d0


 return
 end subroutine vpinteg
