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
 subroutine excpwca(rho,vxc,exc,mmax)

!calculates Perdew-Wang-Ceperly_Alder exchange-correlation potential 
! and energy density

! Ceperley - Alder exchange-correlation potential and energy
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)

!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!mmax  dimension of log radial grid

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0

!Input variables
 integer :: mmax
 real(dp) :: rho(mmax)

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii
 real(dp) :: conrs
 real(dp) :: rs,rh,sqrs,rsl,den

 conrs = (3.d0/(4.d0*pi))**thrd


 do ii=1,mmax
   rh=rho(ii)/pi4
   if(rh .le. 0.0d0) then
     vxc(ii)=0.d0
     exc(ii)=0.d0
   else if(rh .lt. 0.23873241d0) then
     rs=0.62035049d0*rh**(-0.3333333333333333d0)
     sqrs=dsqrt(rs)
     den=1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
     exc(ii)=-0.4582d0/rs - 0.1423d0/den
!          exc(i)=-0.4581652932831428d0/rs - 0.1423d0/den
     vxc(ii)=exc(ii) - rs*(0.15273333d0/rs**2 &
&      + (0.02497128d0/sqrs + 0.01581427d0)/den**2)
   else
     rs=0.62035049d0*rh**(-0.3333333333333333d0)
     rsl=dlog(rs)
     exc(ii)=-0.4582d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs &
&      + 0.002d0*rs*rsl
     vxc(ii)=exc(ii) - rs*(0.15273333d0/rs**2 &
&      + 0.01036667d0/rs - 0.003866667d0 &
&      + 0.00066667d0*(1.0d0 + rsl))
   end if
 end do

 return
 end subroutine excpwca
