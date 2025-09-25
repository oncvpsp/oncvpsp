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
subroutine lschpse(nn,ll,ierr,ee,uld,rr,vv,uu,up,mmax,mch)

! integrates radial Schroedinger equation for pseudopotential
! on a logarithmic mesh finding energy at which desired log derivative
! uld is matched at point mch

!nn  principal quantum number
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!uld  log derivative to be matched
!rr  log radial mesh
!vv  semi-local pseudopotential
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: mmax,mch
   real(dp) :: rr(mmax),vv(mmax)
   real(dp) :: uld
   integer :: nn,ll

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   real(dp) :: ee
   integer :: ierr

!Local variables
   real(dp) :: aeo, aio, als, cn
   real(dp) :: de, emax, emin
   real(dp) :: eps, ro
   real(dp) :: amesh,al
   real(dp) :: sls, sn, uout, upin, upout
   integer :: ii,  it, nin, nint, node

   real(dp), allocatable :: upp(:),cf(:)
   allocate(upp(mmax),cf(mmax))

   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = dexp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
   eps=1.0d-10
   ierr = 60

   sls=ll*(ll+1)

   emin=0.0d0
   do ii=1,mmax
      emin=dmin1(emin,rr(ii)+0.5d0*sls/rr(ii)**2)
   end do
   emax=ee +  10.0d0
   emin=emin - 10.0d0

! null arrays to remove leftover garbage
   uu(:)=0.0d0
   up(:)=0.0d0
   upp(:)=0.0d0

   als=al**2

! return point for bound state convergence
   do nint=1,60

! coefficient array for u in differential eq.
      do ii=1,mmax
         cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
      end do

      nin=mch

! start wavefunction with series

      do ii=1,4
         uu(ii)=rr(ii)**(ll+1)
         up(ii)=al*(ll+1)*rr(ii)**(ll+1)
         upp(ii)=al*up(ii)+cf(ii)*uu(ii)
      end do

! outward integration using predictor once, corrector
! twice
      node=0

      do ii=4,mch-1
         uu(ii+1)=uu(ii)+aeo(up,ii)
         up(ii+1)=up(ii)+aeo(upp,ii)
         do it=1,2
            upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
            up(ii+1)=up(ii)+aio(upp,ii)
            uu(ii+1)=uu(ii)+aio(up,ii)
         end do
         if(uu(ii+1)*uu(ii) <= 0.0d0) node=node+1
      end do

      uout=uu(mch)
      upout=up(mch)


      if(node-nn+ll+1==0) then

         upin=uld*uout

! perform normalization sum

         ro=rr(1)/dsqrt(amesh)
         sn=ro**(2*ll+3)/dfloat(2*ll+3)

         do ii=1,nin-3
            sn=sn+al*rr(ii)*uu(ii)**2
         end do

         sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
         &              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
         &              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0

! normalize u

         cn=1.0d0/dsqrt(sn)
         uout=cn*uout
         upout=cn*upout
         upin=cn*upin

         do ii=1,nin
            up(ii)=cn*up(ii)
            uu(ii)=cn*uu(ii)
         end do
         do ii=nin+1,mmax
            uu(ii)=0.0d0
         end do

! perturbation theory for energy shift

         de=0.5d0*uout*(upout-upin)/(al*rr(mch))

! convergence test and possible exit

         if(dabs(de)<dmax1(dabs(ee),0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if(de>0.0d0) then
            emin=ee
         else
            emax=ee
         end if
         ee=ee+de
         if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)

      else if(node-nn+ll+1<0) then
! too few nodes
         emin=ee
         ee=0.5d0*(emin+emax)

      else
! too many nodes
         emax=ee
         ee=0.5d0*(emin+emax)
      end if

   end do

   deallocate(upp,cf)
   return

end subroutine lschpse
