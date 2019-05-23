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
 subroutine lschfb(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,srel)

!Finds bound states of an all-electron atomic potential using
!Pauli-type  scalar-relativistic Schroedinger equation

!nn  principal quantum number
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
 real(dp) :: rr(mmax),vv(mmax)
 real(dp) :: zz
 integer :: nn,ll
 integer :: mmax
 logical :: srel

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr,mch

!Local variables

 real(dp) :: aei,aeo,aii,aio,als !functions in aeo.f90
 real(dp) :: de,emax,emin
 real(dp) :: eps,fss,tfss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,xx
 integer :: ii,it,nint,node,nin

 real(dp), allocatable :: upp(:),cf(:),dv(:),fr(:),frp(:)
 allocate(upp(mmax),cf(mmax),dv(mmax),fr(mmax),frp(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-10
 ierr = 60

! check arguments
 if(ll>nn-1) then
  write(6,'(/a,i4,a,i4)') 'lschfb: ERROR ll =',ll,' > nn =',nn
  ierr=1
  return
 end if

! relativistic - non-relativistic switch
 if(srel) then
  fss=(1.0d0/137.036d0)**2
 else
  fss=1.0d-12
 end if

 if(ll==0) gamma=dsqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*dsqrt(ll**2-fss*zz**2) + &
& (ll+1)*dsqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 sls=ll*(ll+1)

 emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
 emin=emax
 do ii=1,mmax
   if(ll==0) then
! ad-hock step to eliminate probelms from absurd emin for ll=0
     emin=dmin1(emin,vv(ii)+0.25d0/rr(ii)**2)
   else
     emin=dmin1(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
   end if
 end do
 if(ee==0.0d0) ee=0.5d0*(emax+emin)
 if(ee<emin) ee=0.5d0*(emax+emin)
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
     cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do
  
! calculate dv/dr for darwin correction
   dv(:)=0.0d0

   dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&         -6.d0*vv(5))/(24.d0*al*rr(1))
   dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&         +2.d0*vv(5))/(24.d0*al*rr(2))
  
   do ii=3,mmax-2
     dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&          -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
   end do
  
!  relativistic coefficient arrays for u (fr) and up (frp).
   do ii=1,mmax
     tfss=fss
     fr(ii)=als*(rr(ii)**2)*(-tfss*(vv(ii)-ee)**2 + 0.5d0*tfss*dv(ii)/ &
&     (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
     frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
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
  
! start wavefunction with series
  
   do ii=1,4
     uu(ii)=rr(ii)**gamma
     up(ii)=al*gamma*rr(ii)**gamma
     upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
   end do
  
! outward integration using predictor once, corrector
! twice
   node=0
  
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=(al+frp(ii+1))*up(ii+1)+(cf(ii+1)+fr(ii+1))*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
     if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
   end do
  
   uout=uu(mch)
   upout=up(mch)
  
  
   if(node-nn+ll+1==0) then
  
! start inward integration at 10*classical turning
! point with simple exponential
  
     nin=mch+2.3d0/al
     if(nin+4>mmax) nin=mmax-4
     xkap=dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))
  
     do ii=nin,nin+4
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
     end do
  
! integrate inward
  
     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=(al+frp(ii-1))*up(ii-1)+(cf(ii-1)+fr(ii-1))*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do
  
! scale outside wf for continuity
  
     sc=uout/uu(mch)
  
     do ii=mch,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
  
     upin=up(mch)
  
! perform normalization sum
  
     ro=rr(1)/dsqrt(amesh)
     sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
  
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

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 deallocate(upp,cf,dv,fr,frp)
 return

 end subroutine lschfb
