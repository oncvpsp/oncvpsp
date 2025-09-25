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
subroutine lschvkbbe(nn,ll,nvkb,ierr,ee,uld,emin,emax, &
&   rr,vloc,vkb,evkb,uu,up,mmax,mch)

! integrates radial schroedinger equation for pseudopotential with
! Vanderbilt-Kleinman-Bylander non-local projectors finding energy at
! which desired log derivative uld is matched at point mch


!nn  principal quantum number
!ll  angular-momentum quantum number
!nvkb  = number of VKB projectors to be used
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!uld  bound-state energy, input guess and output calculated value
!emin  externally generaated estimate of lower bound for ee
!emax  externally generaated estimate of upper bound for ee
!rr  log radial mesh
!vloc  local part of psp
!vkb  VKB projectors
!evkb coefficients of VKB projectors
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input Variables
   real(dp) :: emin,emax,uld
   real(dp) :: rr(mmax),vloc(mmax),vkb(mmax,nvkb),evkb(nvkb)
   integer :: nn,ll,nvkb,mmax

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   real(dp) :: ee  !in/out, really - needs starting guess
   integer :: ierr,mch


!Local variables
   real(dp) :: cn
   real(dp) :: de
   real(dp) :: eps,ro
   real(dp) :: sls,sn,uout,upin,upout
   real(dp) :: amesh,al,als
   integer :: ii,nin,nint,node

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


   if(ee>emax) ee=1.25d0*emax
   if(ee<emin) ee=0.75d0*emin
   if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
   uu(:)=0.0d0
   up(:)=0.0d0
   upp(:)=0.0d0

   als=al**2

! return point for bound state convergence
   do nint=1,60

! coefficient array for u in differential eq.
      do ii=1,mmax
         cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ee)*rr(ii)**2
      end do

      nin=mch

! outward integration
      call vkboutwf(ll,nvkb,ee,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)

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

end subroutine lschvkbbe
