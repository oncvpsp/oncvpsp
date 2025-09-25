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
subroutine renorm_r(uu,rr,kap,zz,mmax,cnorm)

! renormalize Dirac wave function so that large component is normalized

!uu  Dirac wave function
!rr  log radial grid
!ll  angular momentum
!kap  Dirac kappa
!zz  atomic number
!mmax  size of radial grid
!cnorm  renormalization coefficient

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   real(dp) :: rr(mmax)
   real(dp) :: zz
   integer :: kap,mmax

!Output variable
   real(dp) :: cnorm

!Input/Output variable
   real(dp) :: uu(mmax,2)

!Local variables
   real(dp) :: sn
   real(dp) :: al,amesh
   real(dp) :: r0,gam,cc,cci
   integer :: ii,nin

   al = 0.01d0 * dlog(rr(101)/rr(1))
   amesh = exp(al)

   cc=137.036d0
   cci=1.0d0/cc

   gam=sqrt(kap**2-(zz*cci)**2)

   do ii=mmax,1,-1
      if(dabs(uu(ii,1))>0.0d0) then
         nin=ii
         exit
      end if
   end do

   r0=rr(1)/dsqrt(amesh)
   sn=r0**(2.0d0*gam+1.0d0)/(2.d0*gam+1.0d0)
   sn=sn*(uu(1,1)**2)/rr(1)**(2.d0*gam)

   do ii=1,nin-3
      sn=sn+al*rr(ii)*uu(ii,1)**2
   end do

   sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2,1)**2 &
   &          + 28.0d0*rr(nin-1)*uu(nin-1,1)**2 &
   &          +  9.0d0*rr(nin  )*uu(nin  ,1)**2)/24.0d0


   cnorm=sqrt(1.0d0/sn)

   uu(:,:)=cnorm*uu(:,:)

   return
end subroutine renorm_r
