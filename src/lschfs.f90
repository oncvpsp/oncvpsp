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
subroutine lschfs(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,srel)

! integrates radial Pauli-type scalar-relativistic equation
! on a logarithmic mesh
! modified routine to be used in finding norm-conserving
! pseudopotential

!nn  effective principal quantum number based on nodes inside mch (output)
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!rr  log radial mesh
!vv  local atomic potential
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations
!srel .true. for scalar-relativistic, .false. for non-relativistic

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: mmax
   real(dp) :: rr(mmax),vv(mmax)
   real(dp) :: zz
   integer :: ll,mch
   logical :: srel

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   real(dp) :: ee
   integer :: ierr,nn


!Local variables
   real(dp) :: amesh,al
   real(dp) :: aeo, aio, als, cn
   real(dp) :: fss, tfss, gamma, ro
   real(dp) :: sls, sn, uout, upout
   integer :: ii, it,node

   real(dp), allocatable :: upp(:),cf(:),dv(:),fr(:),frp(:)

   allocate(upp(mmax),cf(mmax),dv(mmax),fr(mmax),frp(mmax))


   al = 0.01d0 * log(rr(101) / rr(1))
   amesh = exp(al)

   ierr = 0

! relativistic - non-relativistic switch
   if(srel) then
      fss=(1.0d0/137.036d0)**2
   else
      fss=1.0d-12
   end if

   if(ll==0) gamma=sqrt(1.0d0-fss*zz**2)
   if(ll>0) gamma=(ll*sqrt(ll**2-fss*zz**2) &
   & +(ll+1)*sqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

   sls=ll*(ll+1)

! null arrays to remove leftover garbage

   uu(:)=0.0d0
   up(:)=0.0d0
   upp(:)=0.0d0

   node=0

   als=al**2

! coefficient array for u in differential eq.
   do ii=1,mmax
      cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do

! calculate dv/dr for darwin correction
   dv(:)=0.0d0

   dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
   &       -6.d0*vv(5))/(24.d0*al*rr(1))
   dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
   &       +2.d0*vv(5))/(24.d0*al*rr(2))

   do ii=3,mmax-2
      dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
      &         -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
   end do

!  relativistic coefficient arrays for u (fr) and up (frp).
   do ii=1,mmax
      tfss=fss
      fr(ii)=als*(rr(ii)**2)*(-tfss*(vv(ii)-ee)**2 + 0.5d0*tfss*dv(ii)/ &
      &   (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
      frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
   end do

! start wavefunction with series

   do ii=1,4
      uu(ii)=rr(ii)**gamma
      up(ii)=al*gamma*rr(ii)**gamma
      upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
   end do

! outward integration using predictor once, corrector
! twice

   do ii=4,mch-1
      uu(ii+1)=uu(ii)+aeo(up,ii)
      up(ii+1)=up(ii)+aeo(upp,ii)
      do it=1,2
         upp(ii+1)=(al+frp(ii+1))*up(ii+1)+(cf(ii+1)+fr(ii+1))*uu(ii+1)
         up(ii+1)=up(ii)+aio(upp,ii)
         uu(ii+1)=uu(ii)+aio(up,ii)
      end do
      if(uu(ii+1)*uu(ii) <= 0.0d0) node=node+1
   end do

   uout=uu(mch)
   upout=up(mch)

!perform normalization sum

   ro=rr(1)/dsqrt(amesh)
   sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)

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

!calculate effective principal quantum number as if this were a bound
!state with a barrier at mch
   nn=node+ll+1

   deallocate(upp,cf,dv,fr,frp)

   return
end subroutine lschfs
