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
 subroutine fphsft(ll,epsh2,depsh,pshf,rr,vv,zz,mmax,mch,npsh,srel)

! computes full potential scattering log derivatives
! returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
! and easier to compare in plots than the log derivatives themselves
! Pauli-type scalar-relativistic calculation

!ll  angular momentum
!epsh2  upper limit of energy scan
!depsh  increment of scan
!pshf  log derivatives "angles", as above
!rr  radial log grid
!vv  all-electron potential
!zz  atomic number
!mmax  dimension of rr, etc.
!mch  index of radius for log der test
!npsh  number of energy points in scan
!srel .true. for scalar-relativistic, .false. for non-relativistic

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi2=2.0d0*pi

!Input variables
 integer :: ll,mmax,npsh,mch
 real(dp) :: rr(mmax),vv(mmax)
 real(dp) :: depsh,epsh2,zz
 logical :: srel

!Output variables
 real(dp) :: pshf(npsh)

!Local variables
 real(dp) :: al,epsh,phi,phip,pshoff
 integer :: ii,ierr

 real(dp), allocatable :: uu(:),up(:)



 allocate(uu(mmax),up(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))

 pshoff = 0.0d0

 do ii = 1,npsh
   epsh = epsh2-(ii-1)*depsh

   call lschfs(ll,ierr,epsh,rr,vv,uu,up,zz,mmax,mch,srel)

   phi = uu(mch)/rr(mch)
   phip = (up(mch)-al*uu(mch))/(al*rr(mch)**2)
   pshf(ii) = atan2(rr(mch)*phip,phi) + pshoff
   if((ii>1) .and. (pshf(ii)<pshf(ii-1))) then
     pshoff = pshoff+pi2
     pshf(ii) = pshf(ii)+pi2
   end if
 end do

 deallocate(uu,up)
 return
 end subroutine fphsft
