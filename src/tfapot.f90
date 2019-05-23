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
! tfapot
 function tfapot(rr,zz)

! generalized Thomas-Fermi atomic potential

!...to an article of N. H. March ( "The Thomas-Fermi Approximation in 
! Quantum Mechanics", Adv. Phys. 6, 1 (1957)). He has the formula,
! but it is not the result of his work. The original publication is: 
!     R. Latter, Phys. Rev. 99, 510 (1955).
! He says that it''s an analytic fit to an improved calculation of the 
! potential distribution of a Thomas-Fermi atom without exchange first 
! performed by Miranda (C. Miranda, Mem. Acc. Italia 5, 285 (1934)).
!                                 Alexander Seidl, TU Munich

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr,zz
 
!Output function
 real(dp) :: tfapot

!Local variables
 real(dp) :: bb,tt,xx,xs

 bb=(0.69395656d0/zz)**.33333333d0
 xx=rr/bb
 xs=sqrt(xx)

 tt=zz/(1.0d0+xs*(0.02747d0 - xx*(0.1486d0 - 0.007298d0*xx)) &
&   + xx*(1.243d0 + xx*(0.2302d0 + 0.006944d0*xx)))

 if(tt<1.0d0) tt=1.0d0
 tfapot=-tt/rr

 return
 end function tfapot


