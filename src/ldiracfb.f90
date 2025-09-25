!
! Copyright (c) 1989-2012 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
subroutine ldiracfb(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

! Finds relativistic bound states of an al-electron potential

!nn  principal quantum number
!ll  angular-momentum quantum number
!kap =l, -(l+1) for j=l -/+ 1/2
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!rr  log radial mesh
!zz  atomic number
!vv  local psp
!uu(mmax,jj)  output radial wave functions (*rr) jj=1,2 for large, small
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   real(dp) :: rr(mmax),vv(mmax)
   real(dp) :: zz
   integer :: nn,ll,kap,mmax

!Output variables
   real(dp) :: uu(mmax,2),up(mmax,2)
   real(dp) :: ee
   integer :: ierr,mch

!Local Variables

   real(dp), allocatable :: gu(:),fu(:),gup(:),fup(:),cf(:)

   real(dp) :: aei,aeo,aii,aio  !functions in aeo.f90
   real(dp) :: cc,cci,gam,cof
   real(dp) :: de,emax,emin
   real(dp) :: eps,ro,sc
   real(dp) :: sls,sn,cn,uout,upin,upout,xkap
   real(dp) :: amesh,al
   integer :: ii,kk,nint,node,nin

   cc=137.036d0
   cci=1.0d0/cc

   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = exp(al)

! convergence factor for solution of Dirac eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
   eps=1.0d-10
   ierr = 60

! check arguments
   if(ll>nn-1) then
      write(6,'(/a,i4,a,i4)') 'ldiracfb: ERROR ll =',ll,' > nn =',nn
      ierr=1
      return
   end if
   if(kap/=ll .and. kap/=-(ll+1)) then
      write(6,'(/a,i4,a,i4)') 'ldiracfb: ERROR kap =',kap,' ll =',ll
      ierr=2
      return
   end if
   if(zz<1.0d0) then
      write(6,'(/a,f12.8)') 'ldiracfb: ERROR zz =',zz
      ierr=3
      return
   end if

   sls=ll*(ll+1)

! usual limits from non-relativistic Schroedinger eq. should be OK
   emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
   emin=emax
   do ii=1,mmax
      emin=dmin1(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
   end do
   emin=dmax1(emin,-zz**2/nn**2)
   if(ee>emax) ee=0.5d0*(emax+emin)
   if(ee<emin) ee=0.5d0*(emax+emin)
   if(ee>emax) ee=0.5d0*(emax+emin)

   allocate(gu(mmax),fu(mmax),gup(mmax),fup(mmax),cf(mmax))

! null arrays
   gu(:)=0.0d0; fu(:)=0.0d0; gup(:)=0.0d0; fup(:)=0.0d0; cf(:)=0.0d0

   gam=sqrt(kap**2-(zz*cci)**2)

   cof=(kap+gam)*(cc/zz)

! return point for bound state convergence
   do nint=1,60

! coefficient array for u in Schroedinger eq.
      do ii=1,mmax
         cf(ii)=sls + 2.0d0*(vv(ii)-ee)*rr(ii)**2
      end do

! find classical turning point for matching
      mch=0
      do ii=mmax,2,-1
         if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
            mch=ii
            exit
         end if
      end do

      if(mch==0) then
         ierr=-1
         return
      end if

! start wavefunctions with series
      do ii=1,5
         gu(ii)=rr(ii)**gam
         fu(ii)=cof*gu(ii)

         fup(ii)= al*rr(ii)*(kap*fu(ii)/rr(ii) &
         &            + cci*(ee + vv(ii))*gu(ii))

         gup(ii)= al*rr(ii)*(-kap*gu(ii)/rr(ii) &
         &            + cci*(2.0d0*cc**2 - ee - vv(ii))*fu(ii))
      end do

! outward integration using predictor once, corrector
! twice
      node=0
      do ii=5,mch-1

         fu(ii+1)=fu(ii)+aeo(fup,ii)
         gu(ii+1)=gu(ii)+aeo(gup,ii)

         do kk=1,2
            fup(ii+1)= al*rr(ii+1)*(kap*fu(ii+1)/rr(ii+1) &
            &             - cci*(ee - vv(ii+1))*gu(ii+1))

            gup(ii+1)= al*rr(ii+1)*(-kap*gu(ii+1)/rr(ii+1) &
            &             + cci*(2.0d0*cc**2 + ee - vv(ii+1))*fu(ii+1))

            fu(ii+1)=fu(ii)+aio(fup,ii)
            gu(ii+1)=gu(ii)+aio(gup,ii)
         end do
         if(gu(ii+1)*gu(ii) <= 0.0d0) node=node+1
      end do

      uout=gu(mch)
      upout=gup(mch)

      if(node-nn+ll+1==0) then

! start inward integration at 10*classical turning
! point with simple exponential

         nin=int(mch+2.3d0/al)
         if(nin+4>mmax) nin=mmax-4
         xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))

         do ii=nin,nin+4
            gu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
            fu(ii)=0.5d0*cci*(-xkap+kap/rr(ii))*gu(ii)

            fup(ii)= al*rr(ii)*(kap*fu(ii)/rr(ii) &
            &               + cci*(ee + vv(ii))*gu(ii))
            gup(ii)= al*rr(ii)*(-kap*gu(ii)/rr(ii) &
            &               + cci*(2.0d0*cc**2 - ee - vv(ii))*fu(ii))
         end do

! integrate inward
         do ii=nin,mch+1,-1

            fu(ii-1)=fu(ii)+aei(fup,ii)
            gu(ii-1)=gu(ii)+aei(gup,ii)

            do kk=1,2
               fup(ii-1)= al*rr(ii-1)*(kap*fu(ii-1)/rr(ii-1) &
               &               - cci*(ee - vv(ii-1))*gu(ii-1))

               gup(ii-1)= al*rr(ii-1)*(-kap*gu(ii-1)/rr(ii-1) &
               &               + cci*(2.0d0*cc**2 + ee - vv(ii-1))*fu(ii-1))

               fu(ii-1)=fu(ii)+aii(fup,ii)
               gu(ii-1)=gu(ii)+aii(gup,ii)
            end do
         end do

! scale outside wf for continuity

         sc=uout/gu(mch)

         do ii=mch,nin
            gup(ii)=sc*gup(ii)
            gu(ii)=sc*gu(ii)
            fup(ii)=sc*gup(ii)
            fu(ii)=sc*fu(ii)
         end do

         upin=gup(mch)

! perform normalization sum

         ro=rr(1)/dsqrt(amesh)
         sn=((1.0d0+cof**2)*ro**(2.0d0*gam+1.0d0))/(2.0d0*gam+1.0d0)

         do ii=1,nin-3
            sn=sn+al*rr(ii)*(gu(ii)**2 + fu(ii)**2)
         end do

         sn=sn + al*(23.0d0*rr(nin-2)*(gu(nin-2)**2 + fu(nin-2)**2) &
         &              + 28.0d0*rr(nin-1)*(gu(nin-1)**2 + fu(nin-1)**2) &
         &              +  9.0d0*rr(nin  )*(gu(nin  )**2 + fu(nin  )**2))/24.0d0

! normalize u

         cn=1.0d0/dsqrt(sn)
         uout=cn*uout
         upout=cn*upout
         upin=cn*upin

         do ii=1,nin
            gup(ii)=cn*gup(ii)
            gu(ii)=cn*gu(ii)
            fup(ii)=cn*fup(ii)
            fu(ii)=cn*fu(ii)
         end do
         do ii=nin+1,mmax
            gu(ii)=0.0d0
            fu(ii)=0.0d0
         end do

! perturbation theory for energy shift based on large component
! continuity of small component should follow

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

! copy local arrays for output
   uu(:,1)=gu(:)
   up(:,1)=gup(:)
   uu(:,2)=fu(:)
   up(:,2)=fup(:)

!fix sign to be positive at rr->oo
   if(uu(mch,1)<0.0d0) then
      uu(:,:)=-uu(:,:)
      up(:,:)=-up(:,:)
   end if

   deallocate(gu,fu,gup,fup)
   deallocate(cf)
   return

end subroutine ldiracfb
