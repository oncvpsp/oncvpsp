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
! calculates spherical Bessel functions using recursion algorithm

subroutine sbf8(nm,xx,sb_out)

!nm  maximum angular momentum wanted plus 1
!xx  argument of spherical bessel function
!sb_out  output of sbf_l(xx) for l=0,...,nm-1

 implicit none
 integer, parameter :: dp=kind(1.0d0)

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi2=2.0_dp*pi

!Input variables
  integer :: nm
  real(dp) :: xx

!Output variables
  real(dp) :: sb_out(nm)

!Local variables
 integer :: nlim,nn
 real(dp) :: cc,fn,sn,ss,xi,xn,xs,xmod
 real(dp),allocatable :: sb(:)


 if(xx<= 1.0d-36) then
!  zero argument section
   sb_out(:)=0.0d0
   sb_out(1)=1.0d0
 else if(xx<1.0d-3) then
!  small argument section
   xn=1.0d0
   xs=0.5d0*xx**2
   do nn=1,nm
     sb_out(nn)=xn*(1.0d0 - xs*(1.0d0 - xs/(4*nn+6))/(2*nn+1))
     xn=xx*xn/(2*nn+1)
   end do
!  trigonometric section (small l values)
 else if(nm==1) then
   xmod=mod(xx,pi2)
   ss=sin(xmod)
   sb_out(1)=ss/xx 
 else if(nm==2) then
   xmod=mod(xx,pi2)
   ss=sin(xmod)
   cc=cos(xmod)
   sb_out(1)=ss/xx 
   sb_out(2)=(ss-xx*cc)/xx**2
 else if(nm==3) then
   xmod=mod(xx,pi2)
   ss=sin(xmod)
   cc=cos(xmod)
   sb_out(1)=ss/xx 
   sb_out(3)=((3.0d0-xx**2)*ss - 3.0d0*xx*cc)/xx**3
 else if(nm==4) then
   xmod=mod(xx,pi2)
   ss=sin(xmod)
   cc=cos(xmod)
   sb_out(1)=ss/xx 
   sb_out(3)=((3.0d0-xx**2)*ss - 3.0d0*xx*cc)/xx**3
   sb_out(4)=(15.d0*ss-15.d0*xx*cc - 6.d0*xx**2*ss + xx**3*cc )/xx**4
 else
!  recursion method (accurate for large l, slow for large arguments)
   if(xx<1.0d0) then
     nlim=nm+int(15.0d0*xx)+1
   else
     nlim=nm+int(1.36d0*xx)+15
   end if
   allocate(sb(nlim+1))
   nn=nlim
   xi=1.0d0/xx
   sb(nn+1)=0.0d0
   sb(nn)=1.0d-18
   sn=dble(2*nn-1)*1.0d-36
   do nn=nlim-1,1,-1
     sb(nn)=dble(2*nn+1)*xi*sb(nn+1) - sb(nn+2)
   end do
   do nn=1,nlim-1
     sn=sn + dble(2*nn-1)*sb(nn)*sb(nn)
   end do
   fn=1.0d0/sqrt(sn)
   sb_out(:)=fn*sb(1:nm)
   deallocate(sb)
 end if
 return
end subroutine sbf8
