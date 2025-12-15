!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

 subroutine vout(mode,rho,rhoc,vo,vxc,zion,eeel,eexc, &
&                rr,mmax,iexc)

!mode  1=> add rhoc to rho, 0 => don't
!rho  total charge density or valence/pseudovalence charge density
!rhoc  core charge density, full or model
!vo  output total potential
!vc  output exchange-correlation potential
!zion  charge of screened ion
!eeel  electron-electron interaction energy
!eexc  exchange-correlation correction to eigenvalue sum for total energy
!rr  log radial mesh
!output electrostatic and exchange-correlation potential

 implicit none

 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0_dp*pi

!Input vaiables
 integer :: iexc,mmax,mode
 real(dp) :: rho(mmax),rhoc(mmax),rr(mmax)
 real(dp) :: zion

!Output variables
 real(dp) :: vo(mmax),vxc(mmax)
 real(dp) :: eeel,eexc

!Local variables
 integer ii
 real(dp) :: al,tv

!Function
 real(dp) :: aii

 real(dp), allocatable :: rvp(:),rv(:),excca(:),difxc(:)
 real(dp), allocatable :: exca(:),vxcd(:),rhot(:)
 real(dp), allocatable :: foutv(:),foutx(:),foutfc(:)

 allocate(rvp(mmax),rv(mmax),excca(mmax),difxc(mmax))
 allocate(exca(mmax),vxcd(mmax),rhot(mmax))
 allocate(foutv(mmax),foutx(mmax),foutfc(mmax))


 al = 0.01d0 * dlog(rr(101) / rr(1))

 foutfc(:)=0.0d0

! integration for electrostatic potential
 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

 rv(mmax)=zion
 rv(mmax-1)=zion
 rv(mmax-2)=zion

 do ii=mmax-2,2,-1
   rv(ii-1)=rv(ii)+aii(rvp,ii)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=mmax-2,2,-1
   tv=tv+aii(rvp,ii)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vo(ii)=rv(ii)/rr(ii)
 end do

! electron-electron interaction for total energy
 do ii=1,mmax
  foutv(ii)=rho(ii)*vo(ii)*rr(ii)**3
 end do

 eeel=(9.0d0*foutv(1) + 28.0d0*foutv(2) &
&   + 23.0d0*foutv(3))/24.0d0
 do ii=4,mmax
   eeel=eeel + foutv(ii)
 end do
 eeel=al*eeel + foutv(1)/3.0d0

 if(mode .eq. 0) then
   do ii = 1, mmax
     rhot(ii) = rho(ii)
   end do
 else if(mode .eq. 1) then
   do ii = 1, mmax
     rhot(ii) = rho(ii) + rhoc(ii)
   end do
 else
   write(6,'(/a,i4)') 'vout: bad input mode =',mode
   stop
 end if

! exchange-correlation potential added

 if(iexc .eq. 1) then
   call excwig(rhot,vxc,exca,mmax)
 else if(iexc .eq. 2) then
   call exchdl(rhot,vxc,exca,mmax)
 else if(iexc .eq. 3) then
   call excpwca(rhot,vxc,exca,mmax)
 else if(iexc .eq. 4) then
   call excggc(rhot,vxc,exca,rr,mmax)
 else
   write(6,'(/a,i4)') 'vout: bad input iexc =',iexc
   stop
 end if

 do ii = 1, mmax
   difxc(ii)=exca(ii) - vxc(ii)
   vo(ii)=vo(ii) + vxc(ii)
 end do

! exchange-correlation correction for total energy

 do ii=1,mmax
  foutx(ii)=rho(ii)*difxc(ii)*rr(ii)**3
 end do
 eexc=(9.0d0*foutx(1) + 28.0d0*foutx(2) &
&   + 23.0d0*foutx(3))/24.0d0
 do ii=4,mmax
   eexc=eexc + foutx(ii)
 end do

 if(mode .eq. 1 .and. rhoc(1) .ne. 0.0d0) then
   if(iexc .eq. 1) then
     call excwig(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 2) then
     call exchdl(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 3) then
     call excpwca(rhoc,vxcd,excca,mmax)
   else if(iexc .eq. 4) then
     call excggc(rhoc,vxcd,excca,rr,mmax)
   else
     write(6,'(/a,i4)') 'vout: bad input iexc =',iexc
     stop
   end if
   do ii=1,mmax
    foutfc(ii)=rhoc(ii)*(exca(ii) - excca(ii))*rr(ii)**3
   end do
   eexc=eexc + (9.0d0*foutfc(1) + 28.0d0*foutfc(2) &
&     + 23.0d0*foutfc(3))/24.0d0
   do ii=4,mmax
     eexc=eexc + foutfc(ii)
   end do
 end if

 eexc=al*eexc + foutx(1)/3.0d0


 if(mode .eq. 1) then
   eexc=eexc + foutfc(1)/3.0d0
 end if

 deallocate(rvp,rv,excca,difxc)
 deallocate(exca,vxcd,rhot)
 deallocate(foutv,foutx,foutfc)
 return

 end subroutine vout
