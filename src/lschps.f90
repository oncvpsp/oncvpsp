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
subroutine lschps(ll,ierr,ee,rr,vv,uu,up,mmax,mch)

! outward integration of Srcroedinger equation for semi-local pseudopotential
! on a logarithmic mesh

!nn  principal quantum number (not used actually)
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!rr  log radial mesh
!vv  semi-local atomic pseudopotential
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: mmax,mch
   real(dp) :: rr(mmax),vv(mmax)
   integer :: ll

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   real(dp) :: ee
   integer :: ierr


!Local variables
   real(dp) :: amesh,al
   real(dp) :: aeo, aio, als, cn
   real(dp) :: ro
   real(dp) :: sls, sn, uout, upout
   integer :: ii, it

   real(dp), allocatable :: upp(:),cf(:)

   allocate(upp(mmax),cf(mmax))


   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = dexp(al)

   ierr = 0

   sls=ll*(ll+1)

! null arrays to remove leftover garbage

   uu(:)=0.0d0
   up(:)=0.0d0
   upp(:)=0.0d0

   als=al**2

! coefficient array for u in differential eq.
   do ii=1,mmax
      cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do

! start wavefunction with series

   do ii=1,4
      uu(ii)=rr(ii)**(ll+1)
      up(ii)=al*(ll+1)*rr(ii)**(ll+1)
      upp(ii)=al*up(ii)+cf(ii)*uu(ii)
   end do

! outward integration using predictor once, corrector
! twice

   do ii=4,mch-1
      uu(ii+1)=uu(ii)+aeo(up,ii)
      up(ii+1)=up(ii)+aeo(upp,ii)
      do it=1,2
         upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
         up(ii+1)=up(ii)+aio(upp,ii)
         uu(ii+1)=uu(ii)+aio(up,ii)
      end do
   end do

   uout=uu(mch)
   upout=up(mch)

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
   uout=cn*uout
   upout=cn*upout

   do ii=1,mch
      up(ii)=cn*up(ii)
      uu(ii)=cn*uu(ii)
   end do
   do ii=mch+1,mmax
      uu(ii)=0.0d0
   end do

   deallocate(upp,cf)

   return
end subroutine lschps
