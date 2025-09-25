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
! Creates monotonic polynomial model core charge matching all-electron
! core charge and 4 derivatives at "crossover" radius.
! Polynomial is 8th-order with no linear term.

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

subroutine modcore2(rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   fcfact,mmax,rr,nc,nv,la,zion,iexc)

!icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
!rhops  state-by-state pseudocharge density
!rhotps  total pseudocharge density
!rhoc  core-charge density
!rhoae  state-by-state all-electron valence charge density
!rhotae  total all-electron valence charge density
!rhomod  model core density and 4 derivatives
!fcfact  prefactor for model amplitude (multiplies crossover value)
!rcfact  prefactor for model scale (multiplies crossover radius)
!irps  rr index of maximum rc
!mmax  dimension of log grid
!rr log radial grid
!nc  number of core states
!nv  number of valence states
!la  angular-momenta
!zion  ion charge
!iexc  exchange-correlation function to be used

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: nv,nc,iexc,mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
   real(dp) :: rhotps(mmax),rhoc(mmax),rr(mmax)
   real(dp) :: zion,fcfact

!Output variables
   real(dp) :: rhomod(mmax,5)

!convergence criterion
   real(dp), parameter :: eps=1.0d-7

!Local variables
   real(dp) :: al
   real(dp) :: d2mdiff,rmatch,rhocmatch
   real(dp) :: xx,yy,dy
   real(dp) :: x0max,x0min,a0,b0,ymatch,ytrial
   real(dp) :: fmatch(5)
   real(dp) :: drint,rtst,rint(20),fint(20)  !ad-hoc smoothing variables
   real(dp), allocatable :: vxcae(:),vxcpsp(:),vo(:),d2excae(:,:),d2excps(:,:)
   real(dp), allocatable :: dvxcae(:,:),dvxcps(:,:),vxct(:)
   integer :: ii,ircc,irmod,jj,kk
   integer :: iint  !ad-hoc smoothing variables

   allocate(vxcae(mmax),vxcpsp(mmax),vo(mmax))
   allocate(dvxcae(mmax,nv),dvxcps(mmax,nv),vxct(mmax))
   allocate(d2excae(nv,nv),d2excps(nv,nv))

   d2excae(:,:)=0.0d0
   d2excps(:,:)=0.0d0

! find valence pseudocharge - core charge crossover
! this is needed for compatability with icmod=3 option
   ircc = 0
   do ii = mmax,1,-1
      if(rhoc(ii) > rhotps(ii)) then
         ircc=ii
         rmatch=rr(ircc)
         rhocmatch=rhoc(ircc)
         exit
      end if
   end do

! find scaled valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax,1,-1
      if(rhoc(ii) > fcfact*rhotps(ii)) then
         ircc=ii
         exit
      end if
   end do

   if(ircc == 0) then
      write(6,'(/a)') 'modcore2: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   end if

!set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod=mmax

   write(6,'(/a/a)') 'Model core correction analysis',&
   &                  '  based on d2Exc/dn_idn_j'


! get derivatives of all-electron xc energy

   call der2exc(rhotae,rhoc,rhoae,rr,d2excae,d2excps,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excae - all-electron derivatives'
   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excae(kk,jj),jj=1,nv)
   end do


! set model charge to zero
   rhomod(:,:)=0.0d0

! compute d2excps with no core correction
   rhomod(:,:)=0.0d0

   call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
   end do
   write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

   al = 0.01d0 * dlog(rr(101) / rr(1))

! core charge density for Louie-Froyen-Cohen correction

   rhomod(:,1)= rhoc(:)

   xx=rr(ircc)
   fmatch(1)=rhoc(ircc)

! core charge derivatives

!7-point numerical first derivative
   jj=2; ii=ircc
   rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
   &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
   &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
   &     /(60.d0*al*rr(ii))
   fmatch(jj)=rhomod(ircc,jj)

! Fit Teter function to value and slope at ircc

   ymatch=rr(ircc)*fmatch(2)/fmatch(1)

! interval-halving search for dimensionless match point

   x0max=1.48d0
   x0min=0.0d0

   do jj=1,50

      xx=0.5d0*(x0max+x0min)

      call gg1cc(yy,xx)
      call gp1cc(dy,xx)
      ytrial=xx*dy/yy

      if(abs(ytrial-ymatch)<eps) exit

      if(ytrial<ymatch) then
         x0max=xx
      else
         x0min=xx
      end if
   end do

   b0=xx/rr(ircc)
   a0=fmatch(1)/yy

   write(6,'(/a)') 'amplitude prefactor, scale prefactor'
   write(6,'(2f10.4)') a0/rhocmatch,1.0d0/(b0*rmatch)


   rhomod(:,:)=0.0d0
   do ii=1,mmax
      xx=b0*rr(ii)
      if(xx<3.0d0) then
         call gg1cc(yy,xx)
         rhomod(ii,1)=a0*yy
         call gp1cc(yy,xx)
         rhomod(ii,2)=a0*yy*b0
         call gpp1cc(yy,xx)
         rhomod(ii,3)=a0*yy*b0**2
      end if
   end do

! test model
   call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excps - pseudofunction derivatives with core correction'

   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
   end do
   write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

!7-point numerical first derivatives applied successively

! do jj=4,5
!  do ii=3*jj-8,mmax-3
!    if(rhomod(ii,1)==0.0d0) exit
!     rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
!&     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
!&     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
!&     /(60.d0*al*rr(ii))
!  end do
! end do

!ad-hoc treatment of numerical noise near origin
!set up a mesh on which 2nd derivative will have a stable
!polynomial representation
!assumes dpnint remains 7th order
   drint=0.05d0*rr(ircc)
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircc
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,3)
            fint(5-jj)=rhomod(ii,3)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,3),iint-1)

   jj=4
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
      &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
      &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
      &     /(60.d0*al*rr(ii))
   end do

!set up a mesh on which 3rd derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircc
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,4)
            fint(5-jj)=-rhomod(ii,4)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,4),iint-1)

   jj=5
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
      &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
      &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
      &     /(60.d0*al*rr(ii))
   end do

!set up a mesh on which 4th derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircc
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,5)
            fint(5-jj)=rhomod(ii,5)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,5),iint-1)


   deallocate(vxcae,vxcpsp,vo)
   deallocate(dvxcae,dvxcps)
   deallocate(d2excae,d2excps)

   return
end subroutine modcore2
