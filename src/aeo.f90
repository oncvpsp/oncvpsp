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
! Adams extrapolation and interpolation formulas for
! outward and inward integration, Abramowitz and
! Stegun, p. 896

 function aeo(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aeo,yy(*)
 integer :: jj

 aeo=(4.16666666667d-2)*(55.0d0*yy(jj)-59.0d0*yy(jj-1) &
& +37.0d0*yy(jj-2)-9.0d0*yy(jj-3))
 return
 end function aeo

 function aio(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aio,yy(*)
 integer :: jj

 aio=(4.16666666667d-2)*(9.0d0*yy(jj+1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj-1)+yy(jj-2))
 return
 end function aio

 function aei(yy,jj)

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aei,yy(*)
 integer :: jj

 aei=-(4.16666666667d-2)*(55.0d0*yy(jj)-59.0d0*yy(jj+1) &
& +37.0d0*yy(jj+2)-9.0d0*yy(jj+3))
 return
 end function aei

 function aii(yy,jj)

 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: aii,yy(*)
 integer :: jj

 aii=-(4.16666666667d-2)*(9.0d0*yy(jj-1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj+1)+yy(jj+2))
 return
 end function aii

