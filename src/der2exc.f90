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
! computes the 2nd derivative of the contribution to the exchange-
! correlation energy from the region where the valence pseudo wave functions
! differ from the all-electron wave functions
!
subroutine der2exc(rhotot,rhoc,rho,rr,d2exc,d2ref,d2mdiff, &
&                   zion,iexc,nc,nv,la,ircut,mmax)

! rhotot  total valence charge, all-electron of pseudo
! rhoc  core charge, all-electron or model
! rho  valence state-by-state charge (one-electron)
! rr  log radial grid
! d2exc  Exc 2nd-derivative matrix
! d2ref  reference matrix
! d2mdiff root-mean-squared differrenc between d2exc and d2ref
! zion  ion potential
! iexc  exchange-correlation type
! nc  number of core states
! nv  number of valence states
! la  array of l values for all-electron atom
! ircut  maximum-radius for which all-electron and pseudo charges differ
! mmax  dimensiion of log grid

   implicit none
   integer, parameter :: dp=kind(1.0d0)

! Input variables
   real(dp) :: rhotot(mmax),rhoc(mmax),rho(mmax,nv),rr(mmax)
   real(dp) :: d2ref(nv,nv)
   real(dp) :: zion
   integer :: la(nv+nc)
   integer :: iexc,nc,nv,ircut,mmax

! Output variables
   real(dp) :: d2exc(nv,nv)
   real(dp) :: d2mdiff

! Local variables
   real(dp) :: hh,eeel,eexc,ss
   real(dp), allocatable :: vo(:),vxct(:),rhot(:),dvxc(:,:)
   integer :: jj,kk,l1

   allocate(vo(mmax),vxct(mmax),rhot(mmax),dvxc(mmax,nv))

   hh=0.15d0

   do kk=1,nv
      dvxc(:,kk)=0.0d0

      rhot(:)=rhotot(:)-2.d0*hh*rho(:,kk)
      call vout(1,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc)
      dvxc(:,kk)=dvxc(:,kk)+( 2.0d0/(24.0d0*hh))*vxct(:)

      rhot(:)=rhotot(:)     -hh*rho(:,kk)
      call vout(1,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc)
      dvxc(:,kk)=dvxc(:,kk)+(-16.0d0/(24.0d0*hh))*vxct(:)

      rhot(:)=rhotot(:)     +hh*rho(:,kk)
      call vout(1,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc)
      dvxc(:,kk)=dvxc(:,kk)+( 16.0d0/(24.0d0*hh))*vxct(:)

      rhot(:)=rhotot(:)-2.d0*hh*rho(:,kk)
      call vout(1,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc)
      dvxc(:,kk)=dvxc(:,kk)+(-2.0d0/(24.0d0*hh))*vxct(:)

   end do  !kk

! compute Exc 2nd-derivative wrt occupation numbers matrix

   do kk=1,nv
      l1=la(nc+kk)+1
      rhot(:)=rho(:,kk)*rr(:)**2
      do jj=1,nv
         call vpinteg(rhot,dvxc(1,jj),ircut,2*l1,d2exc(kk,jj),rr)
      end do  !jj
   end do  !kk

   ss=0.0d0
   do kk=1,nv
      do jj=1,nv
         ss=ss+(d2exc(jj,kk)-d2ref(jj,kk))**2
      end do
   end do
   d2mdiff=sqrt(ss/nv**2)

   deallocate(vo,vxct,rhot,dvxc)
end subroutine der2exc
