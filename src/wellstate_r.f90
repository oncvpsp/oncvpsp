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
subroutine wellstate_r(nnin,ll,kap,irc,ep,rr,vfull,vwell, &
&                      uu,up,zz,mmax,mch)

!creates quantum well to confine positive-energy state, and calculates
!the resulting all-electron wave function

!nn  principal quantum number of well state
!ll  angular momentum
!kap =l, -(l+1) for j=l -/+ 1/2
!irc  index of core radius
!ep  target energy for well state (>0)
!rr  log radial mesh
!vfull  all-electron potential
!vwell  potential with binding well for ikap=2 calculation
!uu  all-electron well-bound wave function
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log radial mesh
!mch matching mesh point for inward-outward integrations
!srel .true. for scalar-relativistic, .false. for non-relativistic

   implicit none

   integer, parameter :: dp=kind(1.0d0)

!Input variables
   real(dp) :: rr(mmax),vfull(mmax)
   real(dp) :: ep,zz
   integer :: nnin,ll,kap,irc,mmax  !(nnin is actually in/out)

!Output variables
   real(dp) :: uu(mmax,2),up(mmax,2),vwell(mmax)
   integer :: mch

!Local variables
   real(dp) :: al,cwell,et,xx,rwell,rwmax,rwmin,rwscale,umax

   real(dp),parameter :: eps=1.0d-8
   integer :: ii,ierr,iumax,l1,itrwell,nn,nnloop
   logical :: convg

   l1=ll+1
   al = 0.01d0 * dlog(rr(101) / rr(1))

!check for bound state with specified ll and nnin

   uu(:,:)=0.0d0 ; up(:,:)=0.0d0
   et=-0.1d0
   call ldiracfb(nnin,ll,kap,ierr,et,rr,zz,vfull,uu,up,mmax,mch)

!if bound state is found, check its localization
   if(ierr==0) then
      umax=0.0d0
      do ii=mmax,1,-1
         if(dabs(uu(ii,1))>=umax) then
            umax=dabs(uu(ii,1))
         else
            iumax=ii+1
            exit
         end if
      end do

!if bound state is localized compared to rc, use it and its energy
      if(rr(iumax)<0.75d0*rr(irc)) then
         ep=et
         write(6,'(/a,i2,a,i2/a,f12.8/a)') &
         &          'WARNING wellstate: localized bound state found for n=', &
         &           nnin,' l=', ll,'WARNING this state with energy',ep, &
         &          'WARNING will be used for the first projector'

!make sure vwell exists for ikap=2 calculation in this case
         vwell(:)=vfull(:)

         return

!if moderately delocalized use very low-energy scattering state
      else if(rr(iumax)<1.5d0*rr(irc)) then

         ep=0.05d0
         write(6,'(/a,i2,a,i2/a,f12.8/a)') &
         &          'WARNING wellstate: moderately localized bound state found for n=', &
         &          nnin,' l= ',ll,'WARNING scattering state with energy',ep, &
         &          'WARNING will be used for the first projector'
      end if
   end if

!if specified ep is retained and is negative
   if(ep<0.0d0) then
      ep=0.25d0
      write(6,'(/a,i2,a,i2/a,f12.8/a)') &
      &         'WARNING wellstate: negative energy specified for n=', &
      &         nnin,' l= ',ll, &
      &         'WARNING scattering state with energy',ep, &
      &         'WARNING will be used for the first projector'
   end if


!scale factor to decrease well width until trial energy exceeds target
   rwscale=0.8d0

!asymptotic potential value to give target bound state a "binding" energy
!of 0.5 Ha
   cwell=ep+0.5d0

!loop which will increment number of nodes before well wall gets too steep
   do nnloop=0,10
      nn=nnin+nnloop

      rwell=8.0d0*rr(irc)
      rwmax=0.0d0
      rwmin=0.0d0

      convg=.false.

      do itrwell=1,100


!create well potential
         do ii=1,mmax
            vwell(ii)=vfull(ii)
         end do

!   start outside rc to keep numerical derivatives at rc accurate
         do ii=irc+6,mmax
            xx=(rr(ii)-rr(irc+5))/rwell
            vwell(ii)=vwell(ii)+cwell*xx**3/(1.0d0+xx**3)
         end do

!find bound state in well
         et=ep
         call ldiracfb(nn,ll,kap,ierr,et,rr,zz,vwell,uu,up,mmax,mch)

!     if(ierr < 0) then
!       write(6,'(a,i3,a,i3,a,i3)') &
!&            'wellstate_r: ERROR ldiracfb no classical turning point, &
!&             n=',nn,' l=',ll,' kap=',kap
!       stop
!     end if

         if(abs(et-ep)<eps) then
            ep=et
            convg=.true.
            exit
         end if

!Interval-halving search after proper rwell has been bracketed
         if(rwmin>0.0d0 .and. rwmax>0.0d0) then
            if(et>ep) then
               rwmin=rwell
            else
               rwmax=rwell
            end if
            rwell=0.5d0*(rwmax+rwmin)
            cycle
         end if

         if(et>ep) then
            rwmin=rwell
            if(rwmax>0.0d0) then
               rwell=0.5d0*(rwmax+rwmin)
            else
               rwell=rwell/rwscale
            end if
            cycle
         end if

         if(et<ep) then
            rwmax=rwell
            if(rwmin>0.0d0) then
               rwell=0.5d0*(rwmax+rwmin)
            else
               rwell=rwell*rwscale
               if(rwell<dmax1(1.0d0,0.5d0*rr(irc))) exit
            end if
            cycle
         end if

      end do  !itrwell

      if(convg) exit
   end do  !nnloop

   if(.not. convg) then
      write(6,'(a,a,i3,a,i3,a,i3a,f8.4)') &
      &   'ERROR wellstate_r: well potential iteration failed', &
      &   ' to converge, n=',nn,' l=',ll,' kap=', kap,' ep=',ep
      stop
   end if

   if(cwell>0.0d0) then
      write(6,'(/a,i3,a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
      &          '   Wellstate for l =',ll,'  n =',nn,' kap=',kap, &
      &          '     eigenvalue = ',et, &
      &          '     asymptotic potential = ',cwell, &
      &          '     half-point radius =',rr(irc)+rwell
   end if

!carry over node count

   nnin=nn

   return
end subroutine wellstate_r
