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
 subroutine exchdl(rho,vxc,exc,mmax)

!calculates Hedin-Lundquist exchange-correlation potential and energy
!density

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0

!Inpupt variables
 integer :: mmax
 real(dp) :: rho(mmax)

!Output variables
 real(dp) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii
 real(dp) :: conrs
 real(dp) :: rs,ecp,aln,xx,rh

 conrs = (3.d0/(4.d0*pi))**thrd

! Hedin-Lundqvist correlation

 do ii=1,mmax
     rh=rho(ii)*pi4i
     rs=conrs/rh**thrd
     xx=rs/21.0d0
     aln=dlog(1.0d0 + 1.0d0/xx)
     ecp = aln+(xx**3*aln-xx*xx)+xx/2-1.0d0/3.0d0
     exc(ii)=-0.458175d0/rs - 0.0225d0*ecp
     vxc(ii)=-0.6109d0/rs - 0.0225d0*aln
 end do

 return
 end subroutine exchdl
