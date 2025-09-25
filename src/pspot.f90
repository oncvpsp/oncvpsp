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
subroutine pspot(ipr,ll,rr,irc,mmax,al,nbas,qroot,eig,uu,pswfopt_sb, &
&           psopt,vae,vpsp,vkb,ekin_num)

!calculates optimized pseudopotential from coefficients of optimized
!pseudo wave function and all-electron wave function

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!INPUT
!ipr 1 for first projector, 2 for second which usually has a node
!ll  angular momentum
!rr  log radial mesh
!irc  index of rc
!mmax  last point of rr
!al  log of mesh spacing factor
!nbas  number of sbfs
!qroot  q's for set of sbf
!eig  eigenvalue of all-electron wave function
!uu  rr*all-electron wave function
!pswfopt_sb  sbf coefficients for optimized pseudo wave function
!vae  all-electron potential

!OUTPUT
!psopt  rr * pseudo wave function on log grid
!vpsp  pseudopotential
!vkb  (epsilon-T)*phi component of chi projector
!ekin_num  total radial kinetic energy

!Input variables
   integer :: ipr,ll,irc,mmax,nbas
   real(dp) :: rr(mmax),qroot(nbas),uu(mmax),pswfopt_sb(nbas),vae(mmax)
   real(dp) :: al,eig

!Output variables
   real(dp) :: psopt(mmax)
   real(dp) :: vpsp(mmax),vkb(mmax)
   real(dp) :: ekin_num

!Local variables
   real(dp), allocatable :: work(:,:)
   real(dp) :: sbfder(5)
   real(dp) :: amesh,tt,ddt,ske,ro
   integer :: ii,jj,ll1
   integer :: ixp(4)

   allocate(work(mmax,5))
   work(:,:)=0.0d0


   ll1=ll+1
   ixp(1)=2
   ixp(2)=4
   ixp(3)=6
   ixp(4)=8

   do ii=1,mmax
      work(ii,1)=uu(ii)
   end do

!numerical derivatives of all-electron wave function
!note that we are taking derivatives of uu here, not uu/rr as for sbf matching
   do jj=2,3

! 5-point numerical first derivatives applied successively
!  do ii=1+2*jj,mmax-2*jj
!     work(ii,jj)=(2.d0*work(ii-2,jj-1)-16.d0*work(ii-1,jj-1)&
!&     +16.d0*work(ii+1,jj-1)-2.d0*work(ii+2,jj-1))&
!&     /(24.d0*al*rr(ii))

! 7-point numerical first derivatives applied successively
      do ii=1+3*jj,mmax-3*jj
         work(ii,jj)=(-work(ii-3,jj-1)+ 9.d0*work(ii-2,jj-1)&
         &     -45.d0*work(ii-1,jj-1)+45.d0*work(ii+1,jj-1)&
         &     -9.d0*work(ii+2,jj-1)+work(ii+3,jj-1))&
         &     /(60.d0*al*rr(ii))
      end do
   end do

!fill [0, rc] with rr * analytic pswf and 2nd derivatives
   do ii=1,irc
      tt=0.0d0
      ddt=0.0d0
      do jj=1,nbas
         call sbf_rc_der(ll,qroot(jj),rr(ii),sbfder)
         tt=tt+pswfopt_sb(jj)*rr(ii)*sbfder(1)
         ddt=ddt+pswfopt_sb(jj)*(rr(ii)*sbfder(3)+2.0d0*sbfder(2))
      end do
      work(ii,1)=tt
      work(ii,3)=ddt
   end do

   do ii=1,mmax
      psopt(ii)=work(ii,1)
   end do

!integration for radial kinetic energy

   amesh = exp(al)
   ro=rr(5)/sqrt(amesh)

!kinetic operator on psopt including centrifugal term
   work(:,4)=0.5d0*(ll*ll1*work(:,1)/rr(:)**2-work(:,3))


!pseudopotential and node test
   do ii=1,irc
      if(abs(work(ii,1))>0.0d0) then
         if(ipr<=1 .and. work(ii,1)*work(ii+1,1)<0.0d0) then
            write(6,'(a)') ' ERROR pspot:  first pseudo wave function has node, &
            &         program will stop'
            write(6,'(a)') ' ERROR pspot: try changing psp parameters for this l'
            stop
         end if
         vpsp(ii)=-work(ii,4)/work(ii,1) + eig
         vkb(ii)=-work(ii,4) + eig*work(ii,1)
      else
         vpsp(ii)=0.0d0
      end if
   end do

   do ii=irc+1,mmax
      vpsp(ii)=vae(ii)
      vkb(ii)=vae(ii)*work(ii,1)
   end do

!kinetic energy integrand
   work(:,5)=work(:,1)*work(:,4)

   ro=rr(5)/sqrt(amesh)
   ske=((work(5,5))/rr(5)**ixp(ll1)) * ro**ixp(ll1)/ixp(ll1)

   do ii=5,mmax-9
      ske=ske+al*rr(ii)*work(ii,5)
   end do

   ekin_num=ske

   deallocate(work)
   return
end subroutine pspot
