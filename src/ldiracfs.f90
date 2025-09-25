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
subroutine ldiracfs(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

! Finds relativistic scattering states of an al-electron potential

!nn  effective principal quantum number based on nodes inside mch (output)
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
   integer :: ll,kap,mmax

!Output variables
   real(dp) :: uu(mmax,2),up(mmax,2)
   real(dp) :: ee
   integer :: ierr,mch,nn

!Local Variables

   real(dp), allocatable :: gu(:),fu(:),gup(:),fup(:),cf(:)

   real(dp) :: aeo,aio  !functions in aeo.f90
   real(dp) :: cc,cci,gam,cof
   real(dp) :: eps,ro
   real(dp) :: sn,cn,uout,upout
   real(dp) :: amesh,al
   integer :: ii,kk,node

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
   if(kap/=ll .and. kap/=-(ll+1)) then
      write(6,'(/a,i4,a,i4)') 'ldiracfs: ERROR kap =',kap,' ll =',ll
      ierr=2
      return
   end if
   if(zz<1.0d0) then
      write(6,'(/a,f12.8)') 'ldiracfs: ERROR zz =',zz
      ierr=3
      return
   end if

   allocate(gu(mmax),fu(mmax),gup(mmax),fup(mmax),cf(mmax))

! null arrays
   gu(:)=0.0d0; fu(:)=0.0d0; gup(:)=0.0d0; fup(:)=0.0d0; cf(:)=0.0d0

   node=0
   ierr=0

   gam=sqrt(kap**2-(zz*cci)**2)

   cof=(kap+gam)*(cc/zz)

! start wavefunctions with series
   do ii=1,5
      gu(ii)=rr(ii)**gam
      fu(ii)=cof*gu(ii)

      fup(ii)= al*rr(ii)*(kap*fu(ii)/rr(ii) &
                          + cci*(ee + vv(ii))*gu(ii))

      gup(ii)= al*rr(ii)*(-kap*gu(ii)/rr(ii) &
                          + cci*(2.0d0*cc**2 - ee - vv(ii))*fu(ii))
   end do

! outward integration using predictor once, corrector
! twice
   do ii=5,mch-1

      fu(ii+1)=fu(ii)+aeo(fup,ii)
      gu(ii+1)=gu(ii)+aeo(gup,ii)

      do kk=1,2
         fup(ii+1)= al*rr(ii+1)*(kap*fu(ii+1)/rr(ii+1) &
                                 - cci*(ee - vv(ii+1))*gu(ii+1))

         gup(ii+1)= al*rr(ii+1)*(-kap*gu(ii+1)/rr(ii+1) &
                                 + cci*(2.0d0*cc**2 + ee - vv(ii+1))*fu(ii+1))

         fu(ii+1)=fu(ii)+aio(fup,ii)
         gu(ii+1)=gu(ii)+aio(gup,ii)
      end do
      if(gu(ii+1)*gu(ii) <= 0.0d0) node=node+1
   end do

   uout=gu(mch)
   upout=gup(mch)


! perform normalization sum

   ro=rr(1)/dsqrt(amesh)
   sn=((1.0d0+cof**2)*ro**(2.0d0*gam+1.0d0))/(2.0d0*gam+1.0d0)

   do ii=1,mch-3
      sn=sn+al*rr(ii)*(gu(ii)**2 + fu(ii)**2)
   end do

   sn=sn + al*(23.0d0*rr(mch-2)*(gu(mch-2)**2 + fu(mch-2)**2) &
               + 28.0d0*rr(mch-1)*(gu(mch-1)**2 + fu(mch-1)**2) &
   &              +  9.0d0*rr(mch  )*(gu(mch  )**2 + fu(mch  )**2))/24.0d0

! normalize u

   cn=1.0d0/dsqrt(sn)
   uout=cn*uout
   upout=cn*upout

   do ii=1,mch
      gup(ii)=cn*gup(ii)
      gu(ii)=cn*gu(ii)
      fup(ii)=cn*fup(ii)
      fu(ii)=cn*fu(ii)
   end do
   do ii=mch+1,mmax
      gu(ii)=0.0d0
      fu(ii)=0.0d0
   end do

! copy local arrays for output
   uu(:,1)=gu(:)
   up(:,1)=gup(:)
   uu(:,2)=fu(:)
   up(:,2)=fup(:)

!calculate effective principal quantum number as if this were a bound
!state with a barrier at mch
   nn=node+ll+1

   deallocate(gu,fu,gup,fup)
   deallocate(cf)
   return

end subroutine ldiracfs
