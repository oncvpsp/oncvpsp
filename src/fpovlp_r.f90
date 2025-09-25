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
! calculates all-electron wave function overlap intergrals

subroutine fpovlp_r(gg,hh,nn,ll,kap,zz,ss,rr)

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!gg Dirac wave function one
!hh Dirac wave function two
!nn outer integration limit
!ll angular momentum
!kap  =l, -(l+1) for j=l -/+ 1/2
!zz nuclear charge
!ss overlap output
!rr log radial mesh

! product of all-elecctron relativistic wave functions gg*hh
! goes like rr**(2*gamma) as rr -> 0
! integral on usual log mesh from rr=0 to rr(nn)

!Input variables
   real(dp) :: zz
   real(dp) :: gg(nn,2),hh(nn,2),rr(nn)
   integer :: nn,ll,kap

!Output variable
   real(dp) :: ss

!Local variables
   real(dp) :: r0,amesh,al,gam,cc,cci
   integer :: ii

   al = 0.01d0 * dlog(rr(101)/rr(1))
   amesh = exp(al)

   cc=137.036d0
   cci=1.0d0/cc

   gam=sqrt(kap**2-(zz*cci)**2)

   r0=rr(1)/dsqrt(amesh)
   ss=r0**(2.0d0*gam+1.0d0)/(2.d0*gam+1.0d0)

! original version used large and small components, now deprecated (see docs)

! ss=ss*(gg(1,1)*hh(1,1)+gg(1,2)*hh(1,2))/rr(1)**(2.d0*gam)

! do ii = 1, nn - 3
!   ss =  ss + al*(gg(ii,1)*hh(ii,1)+gg(ii,2)*hh(ii,2))*rr(ii)
! end do

! ss=ss+al*(23.d0*rr(nn-2)*(gg(nn-2,1)*hh(nn-2,1)+gg(nn-2,2)*hh(nn-2,2)) &
!&        + 28.d0*rr(nn-1)*(gg(nn-1,1)*hh(nn-1,1)+gg(nn-1,2)*hh(nn-1,2)) &
!&        +  9.d0*rr(nn  )*(gg(nn  ,1)*hh(nn  ,1)+gg(nn  ,2)*hh(nn  ,2)))/24.d0

! current version based on large-component only

   ss=ss*(gg(1,1)*hh(1,1))/rr(1)**(2.d0*gam)

   do ii = 1, nn - 3
      ss =  ss + al*(gg(ii,1)*hh(ii,1))*rr(ii)
   end do

   ss=ss+al*(23.d0*rr(nn-2)*(gg(nn-2,1)*hh(nn-2,1)) &
   &        + 28.d0*rr(nn-1)*(gg(nn-1,1)*hh(nn-1,1)) &
   &        +  9.d0*rr(nn  )*(gg(nn  ,1)*hh(nn  ,1)))/24.d0

   return
end subroutine fpovlp_r
