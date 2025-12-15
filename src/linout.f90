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
! interpolates various arrays onto linear radial mesh to create file
! for Abinit input using pspcod=8

 subroutine linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                  zz,zion,mmax,iexc,fcfact,nrl,drl,atsym)


!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local potential 
!  if linear combination is used)
!rho  valence pseudocharge
!rhomod  model core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!iexc  type of exchange-correlation
!fcfact  parameter determining model core charge crossover (0.0 if not used)
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,iexc,mmax,nrl
 integer :: nproj(6)
 real(dp) :: drl,fcfact,zz,zion
 real(dp) :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,2,4)
 real(dp) :: rhomod(mmax,5)
 real(dp):: rc(6),evkb(2,4)
 character*2 :: atsym

!Output variables - printing only

!Local variables
 integer :: ii,jj,ll,l1,ixc_abinit
 integer :: dtime(8)
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:),vpl(:,:)
 real(dp), allocatable :: vkborl(:,:,:)
 character*2 :: pspd(3)


 allocate(rhol(nrl),rl(nrl),vkbl(nrl,2,4),vpl(nrl,5),rhomodl(nrl,5))
 allocate(vkborl(mmax,2,4))

! divide Kleinman-Bylander / Vanderbilt potentials by their r -> 0 dependence

 do l1 = 1, lmax + 1
   ll = l1 - 1
   if (ll .ne. lloc) then
     vkborl(:,1,l1)=vkb(:,1,l1)/rr(:)**(ll+1)
     vkborl(:,2,l1)=vkb(:,2,l1)/rr(:)**(ll+1)
   end if
 end do
!
! interpolation of everything onto linear output mesh
! cubic extrapolation to zero from nearest 4 points of linear mesh
!
 do  ii=1,nrl
   rl(ii)=drl*dble(ii-1)
 end do
!
 do l1=1,max(lmax+1,lloc+1)
   call dp3int(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)
   vpl(1,l1)=4.0d0*vpl(2,l1)-6.0d0*vpl(3,l1)  &
&           +4.0d0*vpl(4,l1)-      vpl(5,l1)
   if(l1 .ne. lloc + 1) then

     call dp3int(rr,vkborl(1,1,l1),mmax,rl,vkbl(1,1,l1),nrl)
     call dp3int(rr,vkborl(1,2,l1),mmax,rl,vkbl(1,2,l1),nrl)

     vkbl(1,1,l1)=4.0d0*vkbl(2,1,l1)-6.0d0*vkbl(3,1,l1)  &
&                +4.0d0*vkbl(4,1,l1)-      vkbl(5,1,l1)
     vkbl(1,2,l1)=4.0d0*vkbl(2,2,l1)-6.0d0*vkbl(3,2,l1)  &
&                +4.0d0*vkbl(4,2,l1)-      vkbl(5,2,l1)
   end if
 end do
 
! rescale linear-mesh vkbl to restore actual r dependence

 do l1=1,lmax+1
   ll=l1-1
   if(ll.ne.lloc)then
     do ii = 1, nrl
       vkbl(ii,1,l1)=vkbl(ii,1,l1)*rl(ii)**(ll+1)
       vkbl(ii,2,l1)=vkbl(ii,2,l1)*rl(ii)**(ll+1)
     end do
   end if
 end do

 call dp3int(rr,rho,mmax,rl,rhol,nrl)

 rhol(1)= 4.0d0*rhol(2)-6.0d0*rhol(3) &
&        +4.0d0*rhol(4)-      rhol(5)

 do jj=1,5
   call dp3int(rr,rhomod(1,jj),mmax,rl,rhomodl(1,jj),nrl)

   rhomodl(1,jj)=4.0d0*rhomodl(2,jj)-6.0d0*rhomodl(3,jj) &
&               +4.0d0*rhomodl(4,jj)-      rhomodl(5,jj)
 end do

! Output for Abinit input using pspcod=8

 if(iexc==1) then
   ixc_abinit=4
 else if(iexc==2) then
   ixc_abinit=5
 else if(iexc==3) then
   ixc_abinit=2
 else if(iexc==4) then
   ixc_abinit=11
 end if

 call date_and_time(VALUES=dtime)
 ii=dtime(1)-2000
 if(ii<10) then
   write(pspd(1),'(a,i1)') '0',ii
 else
   write(pspd(1),'(i2)') ii
 end if
 ii=dtime(2)
 if(ii<10) then
   write(pspd(2),'(a,i1)') '0',ii
 else
   write(pspd(2),'(i2)') ii
 end if
 ii=dtime(3)
 if(ii<10) then
   write(pspd(3),'(a,i1)') '0',ii
 else
   write(pspd(3),'(i2)') ii
 end if

 write(6,'(/a)') 'Begin PSPCODE8'
 write(6,'(3a,4f6.2)') atsym,'    ONCVPSP' &
&  ,'  r_core=',(rc(l1),l1=1,lmax+1)
 write(6,'(2f12.4, 5a)') zz,zion, '      ', pspd,  &
&  '    zatom,zion,pspd'
 write(6,'(6i6, a)') 8, ixc_abinit,lmax,lloc, &
&  nrl, 0, '    pspcod,pspxc,lmax,lloc,mmax,r2well'
 write(6,'(3f12.8, a)') rl(nrl),fcfact, 0.0,  &
&  '    rchrg fchrg qchrg'
 write(6,'(5i6, a)') (nproj(l1),l1=1,5),  &
&  '    nproj'
 write(6,'(i6, a)') 0, &
&  '                 extension_switch'

 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(ll==lloc) then
     write(6,'(i4)') ll
     do ii = 1,nrl
       write(6,'(i6,1p,2d21.13)') ii,rl(ii),vpl(ii,l1)
     end do
   elseif(nproj(l1)==1) then
     write(6,'(i4,23x,1p,d21.13)') ll,evkb(1,l1)
     do ii = 1,nrl
       write(6,'(i6,1p,2d21.13)') ii,rl(ii),vkbl(ii,1,l1)
     end do
   elseif(nproj(l1)==2) then
     write(6,'(i4,23x,1p,2d21.13)') ll,evkb(1,l1), &
&      evkb(2,l1)
     do ii=1,nrl
       write(6,'(i6,1p,3d21.13)') ii,rl(ii),vkbl(ii,1,l1), &
&        vkbl(ii,2,l1)
     end do
   endif
 end do

 if(fcfact>0.0d0) then
   do ii=1,nrl
     write(6,'(i6,1p,6d21.13)') ii,rl(ii),  &
&      (rhomodl(ii,jj),jj=1,5)
   end do
 end if
  
 deallocate(rhol,rl,vkbl,vpl,rhomodl)
 deallocate(vkborl)

 return
 end subroutine linout
