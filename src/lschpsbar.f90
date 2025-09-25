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
subroutine lschpsbar(nn,ll,ierr,ee,emin,emax,rr,vv,uu,up,mmax,mbar,tht)

! Finds bound states of a semi-local pseudopotential

!nn  principal quantum number
!ll  angular-momentum quantum number
!ierr  non-zero return if error/ input for box boundary condition (0 or 1)
!ee  bound-state energy, input guess and output calculated value (in/our)
!emin  lower bound, potential minimum if ==0.0
!emax  upper bound
!rr  log radial mesh
!vv  local psp
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mbar  mesh point for infinite barrier
!tht  phast-shift angle for boundary condition, units of pi

   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp


!Input variables
   real(dp) :: rr(mmax),vv(mmax)
   real(dp) :: emin,emax,tht
   integer :: nn,ll,mmax,mbar

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   real(dp) :: ee
   integer :: ierr

!Local Variables

   real(dp) :: aei,aeo,aii,aio  !functions in aeo.f90
   real(dp) :: de,uumax
   real(dp) :: eps,ro,sc,stht,ctht
   real(dp) :: sls,sn,cn,uout,upin,upout,xkap,xkap2,xx
   real(dp) :: amesh,al,als
   integer :: ii,imin,it,mch,mchb,nint,node,nin

   real(dp), allocatable :: upp(:),cf(:)
   allocate(upp(mmax),cf(mmax))

   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = dexp(al)
   stht=sin(pi*tht)
   ctht=cos(pi*tht)

   nin=0
   mch=0
   mchb=0
! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
   eps=1.0d-8
   ierr = 60

   sls=ll*(ll+1)

   if(emin==0.0d0) then
      emin=emax
      do ii=1,mbar
         if(emin>vv(ii)+0.5d0*sls/rr(ii)**2) then
            emin=vv(ii)+0.5d0*sls/rr(ii)**2
            imin=ii
         end if
      end do
   end if
! avoid error in classical turning point search
   imin=max(imin,2)

   if(ee<emin) ee=dmin1(emin+5.0d0, 0.5d0*(emin+emax))
   node=0
   de=0.0d0
   mch=0

! null arrays to remove leftover garbage
   uu(:)=0.0d0
   up(:)=0.0d0
   upp(:)=0.0d0
   cf(:)=0.0d0
!maxset=.false.

   als=al**2
! return point for bound state convergence
   do nint=1,60

      uu(:)=0.0d0 ; up(:)=0.0d0 ; upp(:)=0.0d0

      do ii=1,mmax
         cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
      end do

! find classical turning point for matching
      mch=mbar-1
      if(cf(mbar)>0.0d0) then
         do ii=mbar-1,imin,-1
            if(cf(ii)<=0.d0) then
               mch=ii
               exit
            end if
         end do
      end if

      de=0.0d0
! start wavefunction with series

      do ii=1,4
         uu(ii)=rr(ii)**(ll+1)
         up(ii)=al*(ll+1)*rr(ii)**(ll+1)
         upp(ii)=al*up(ii)+cf(ii)*uu(ii)
      end do

! outward integration using predictor once, corrector
! twice
      node=0
      uumax=0.0d0
      do ii=4,mch-1
         uu(ii+1)=uu(ii)+aeo(up,ii)
         up(ii+1)=up(ii)+aeo(upp,ii)
         do it=1,2
            upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
            up(ii+1)=up(ii)+aio(upp,ii)
            uu(ii+1)=uu(ii)+aio(up,ii)
         end do
         uumax=dmax1(uumax,abs(uu(ii+1)))
      end do

      do ii=4,mch-1
         if(abs(uu(ii+1))>eps*uumax .or. abs(uu(ii))>eps*uumax) then
            if(uu(ii+1)*uu(ii) < 0.0d0) node=node+1
         end if
      end do

      mchb=mch
      uout=uu(mchb)
      upout=up(mchb)

      de=0.0d0
      if(node-nn+ll+1==0) then

! start inward integration at 10*classical turning
! point with simple exponential

         nin=int(mch+2.3d0/al)
         if(nin<mbar-4) then
            if(nin+4>mmax) nin=mmax-4
            xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))

            do ii=nin,nin+4
               uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
               up(ii)=-rr(ii)*al*xkap*uu(ii)
               upp(ii)=al*up(ii)+cf(ii)*uu(ii)
            end do

         else

            nin=mbar

!   start inward integration at barrier
!   point with sinh or sin

            xkap2=sls/rr(mbar)**2 + 2.0d0*(vv(mbar)-ee)
            xkap=dsqrt(dabs(xkap2))


            do ii=mbar,mbar+4
               xx=xkap*(rr(ii)-rr(mbar))

               if(xkap2>0) then
                  uu(ii)=rr(ii)*(ctht*sinh(-xx)/xkap + stht*cosh(-xx))

                  up(ii)=al*rr(ii)**2*(-ctht*cosh(-xx) - xkap*stht*sinh(-xx)) &
                  &                 +al*uu(ii)
               else
                  uu(ii)=rr(ii)*(ctht*sin(-xx)/xkap + stht*cos(-xx))

                  up(ii)=al*rr(ii)**2*(-ctht*cos(-xx) + xkap*stht*sin(-xx)) &
                  &                 +al*uu(ii)
               end if
               upp(ii)=al*up(ii)+cf(mbar)*uu(ii)
            end do

         end if  !nin<mbar

! integrate inward

         do ii=nin,mchb+1,-1
            uu(ii-1)=uu(ii)+aei(up,ii)
            up(ii-1)=up(ii)+aei(upp,ii)
            do it=1,2
               upp(ii-1)=al*up(ii-1)+cf(ii-1)*uu(ii-1)
               up(ii-1)=up(ii)+aii(upp,ii)
               uu(ii-1)=uu(ii)+aii(up,ii)
            end do
         end do

         do ii=nin+1,mmax
            uu(ii)=0.0d0
            up(ii)=0.0d0
            upp(ii)=0.0d0
         end do

! scale outside wf for continuity

         sc=uout/uu(mchb)

         do ii=mchb,nin
            up(ii)=sc*up(ii)
            uu(ii)=sc*uu(ii)
         end do

         upin=up(mchb)

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

         de=0.5d0*uout*(upout-upin)/(al*rr(mchb))

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

   end do  !nint

   deallocate(upp,cf)
   return

end subroutine lschpsbar
