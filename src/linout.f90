!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
&                  rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&                  epsh1,epsh2,depsh,rlmax,psfile)


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
!rhotae  all-electron valence charge
!rhoc  all-electron core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!mxprj  dimension of number of projectors
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol
!remaining input variables to be echoed:
!  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf
!  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact
!  epsh1,epsh2,depsh,rlmax,psfile
!psfile  should be 'psp8' or 'both'

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,iexc,mmax,mxprj,nrl,icmod
 integer :: nproj(6)
 real(dp) :: drl,fcfact,rcfact,zz,zion
 real(dp) :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,mxprj,4)
 real(dp) :: rhotae(mmax),rhoc(mmax)
 real(dp) :: rhomod(mmax,5)
 real(dp):: rc(6),evkb(mxprj,4)
 character*2 :: atsym

!additional input for psp8 output to echo input file, all as defined
! in the main progam
 integer :: na(30),la(30),ncon(6),nbas(6)
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer :: nc,nv,lpopt,ncnf
 real(dp) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp) :: dvloc0,epsh1,epsh2,depsh,rlmax
 character*4 :: psfile

!Output variables - printing only

!Local variables
 integer :: ii,iprj,jj,ll,l1,ixc_abinit
 integer :: dtime(8)
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:),vpl(:,:)
 real(dp),allocatable :: rhotael(:),rhocl(:)
 character*2 :: pspd(3)


 allocate(rhol(nrl),rl(nrl),vkbl(nrl,mxprj,4),vpl(nrl,5),rhomodl(nrl,5))
 allocate(rhotael(nrl),rhocl(nrl))

!
! interpolation of everything onto linear output mesh
 
 do  ii=1,nrl
   rl(ii)=drl*dble(ii-1)
 end do
!
 do l1=1,max(lmax+1,lloc+1)
   call dpnint(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)

! override dpnint extrapolation to zero for vpl
   vpl(1,l1)=vpuns(1,l1)

   if(l1 .ne. lloc + 1) then
     do iprj=1,nproj(l1)

       call dpnint(rr,vkb(1,iprj,l1),mmax,rl,vkbl(1,iprj,l1),nrl)

     end do
   end if
 end do

 call dpnint(rr,rho,mmax,rl,rhol,nrl)
 call dpnint(rr,rhotae,mmax,rl,rhotael,nrl)
 call dpnint(rr,rhoc,mmax,rl,rhocl,nrl)

 do jj=1,5
   call dpnint(rr,rhomod(1,jj),mmax,rl,rhomodl(1,jj),nrl)
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
 else if(iexc<0) then
   ixc_abinit=iexc
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
 write(6,'(3a,4f10.5)') atsym,'    ONCVPSP-4.0.1' &
&  ,'  r_core=',(rc(l1),l1=1,lmax+1)
 write(6,'(2f12.4, 5a)') zz,zion, '      ', pspd,  &
&  '    zatom,zion,pspd'
 write(6,'(i6,i8,i4,3i6, a)') 8, ixc_abinit,lmax,lloc, &
&  nrl, 0, '    pspcod,pspxc,lmax,lloc,mmax,r2well'
 write(6,'(3f12.8, a)') rl(nrl),fcfact, 0.0,  &
&  '    rchrg fchrg qchrg'
 write(6,'(5i6, a)') (nproj(l1),l1=1,5),  &
&  '    nproj'
 write(6,'(2i6, a)') 1,1, &
&  '           extension_switch'

! write the VKB projectors and the local potential
 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(ll==lloc) then
     write(6,'(i4)') ll
     do ii = 1,nrl
       write(6,'(i6,1p,2e21.13)') ii,rl(ii),vpl(ii,l1)
     end do
   else if(nproj(l1)>0) then
     write(6,'(i4,23x,1p,6e21.13)') ll,(evkb(jj,l1),jj=1,nproj(l1))
     do ii = 1,nrl
       if(rl(ii)<=rc(l1)+2.0d0*drl) then
         write(6,'(i6,1p,6e21.13)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
       else
         write(6,'(i6,f6.2,6f4.0)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
       end if
     end do
   endif
 end do


! write the model core charge if called for
 if(fcfact>0.0d0) then
   do ii=1,nrl
     write(6,'(i6,1p,6e21.13)') ii,rl(ii),  &
&      (rhomodl(ii,jj),jj=1,5)
   end do
 end if

! write pseudo valence, all-electron valence and all-electron core charges
 do ii=1,nrl
   write(6,'(i6,1p,4e21.13)') ii,rl(ii),rhol(ii),rhotael(ii),rhocl(ii)
 end do

 write(6,'(a)') '<INPUT>'
 write(6,'(a/a/a/a)') &
&    '#', &
&    '#ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&    '#scalar-relativistic version 4.0.1 0r/27/2018', &
&    '#'

 write(6,'(a/a/a/a)') &
&    '#While it is not required under the terms of the GNU GPL, it is',&
&    '#suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&    '#in any publication utilizing these pseudopotentials.', &
&    '#'

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',trim(psfile)
 write(6,'(a/a)') '#','#   n    l    f'
 do ii=1,nc+nv
   write(6,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc0(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc0(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact, rcfact'
 write(6,'(i5,2f10.5)') icmod,fcfact,rcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.2)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

! write termination signal
 write(6,'(a)')'</INPUT>'
 write(6,'(/a)') 'END_PSP'
  
 deallocate(rhol,rl,vkbl,vpl,rhomodl,rhotael,rhocl)

 return
 end subroutine linout
