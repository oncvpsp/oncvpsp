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
 subroutine run_plot(lmax,np,fp,ep,lloc,rc,npx,lpx,fpx,epx,nxtra, &
&                    vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,drl,nrl, &
&                    rho,rhoc,rhomod,pswf,srel,cvgplt)


! write output for plotting pseudopotentials, model core charge, and
! all-electron and pseudo wave functions

!lmax  maximum angular momentum
!np  principal quantum number for corresponding all-electron state
!fp  occupancy
!ep  bound-state or scattering state reference energies for vkb potentials
!lloc  l for local potential
!rc  core radii
!npx  principal quantum number for corresponding all-electron state
!lpx  l's for extra valence states
!fpx  occupancies for extra valence states
!epx  energies for extra valence states
!nxtra  number of extra valence states
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials
!vpuns  unscreened vp
!zz  atomic number
!mmax  size of radial grid
!drl  spacing of linear radial mesh
!nrl  number of points in radial mesh
!rho  valence pseudocharge
!rhoc core charge
!rhomod  model core charge
!pswf  pseudo wave functions for 1st and 2nd projectors
!srel .true. for scalar-relativistic, .false. for non-relativistic
!cvgplt  Energy per electron error vs. cutoff 

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,mmax,nlim,nxtra,nrl
 integer :: np(6),npx(6),lpx(6),nproj(6)
 real(dp) :: zz,drl
 real(dp) :: rr(mmax),vp(mmax,5),vpuns(mmax,5),vfull(mmax),vkb(mmax,2,4)
 real(dp) :: rho(mmax),rhoc(mmax),rhomod(mmax,5),pswf(mmax,2,4)
 real(dp):: ep(6),fp(6),rc(6),epx(6),fpx(6),evkb(2,4),cvgplt(2,4,2,4)
 logical :: srel

!Output variables - printing only

!Local variables
 integer :: ll,l1,ii,jj,ierr,mch,mchf,n1,n2,node
 real(dp) :: al,emax,emin,etest,sls,rmx,sgnae,sgnps
 real(dp), allocatable :: uu(:),u2(:),up(:)

 allocate(uu(mmax),u2(mmax),up(mmax))

 al = 0.01d0 * dlog(rr(101)/rr(1))

 write(6,'(/a)') 'DATA FOR PLOTTING'

 n1 = idint(dlog(drl/rr(1))/al+1.0d0)
 n2 = idint(dlog(dfloat(nrl)*drl/rr(1))/al+1.0d0)


 write(6,'(/2a/)') ' radii, charge,',' pseudopotentials (ll=0, 1, lmax)'

 do ii=n1,n2
   write(6,'(a,6f12.7)') '!p',rr(ii),rho(ii),(vpuns(ii,l1),l1=1,lmax+1)
 end do

 if(lloc .eq. 4) then
   do ii=n1,n2
     write(6,'(a, 2f12.7)') ' !L',rr(ii),vpuns(ii,lloc + 1)
   end do
 end if

! output for analysis of core charge and models
 rmx=0
 do ii=n1,n2
  rmx=max(rmx,5.0d0*rho(ii),rhomod(ii,1))
 end do

 write(6,'(//2a/)') ' radii, charge,',' core charge, model core charge'
 do ii=n1,n2
  if(rhoc(ii)<rmx)then
   write(6,'(a,8f12.7)') '!r',rr(ii),rho(ii),rhoc(ii), &
&   rhomod(ii,1)
  else if(rr(ii)>0.01d0) then
   write(6,'(a,8f12.7)') '!r',rr(ii),rho(ii),rmx, &
&   rhomod(ii,1)
  end if
 end do

! loop for wave function output

 do l1 = 1, lmax + 1
   ll = l1 - 1
 write(6,'(//a,i2,a,i2,a,a/)') 'n=',ll+1,',  l=',ll, &
& ', all-electron wave function', ', pseudo w-f'
   if(fp(l1) .ne. 0.0d0 .or. ep(l1)<0.0d0) then
     etest=ep(l1)
     call lschfb(np(l1),ll,ierr,etest, &
&                rr,vfull,uu,up,zz,mmax,mch,srel)
     sls=(l1-1)*l1
     emax=vp(mmax,l1)+0.5d0*sls/rr(mmax)**2
     emin=emax
     do ii=1,mmax
       emin=dmin1(emin,vp(ii,l1)+0.5d0*sls/rr(ii)**2)
     end do
     call lschvkbb(ll+1,ll,nproj(l1),ierr,etest,emin,emax, &
&                  rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                  u2,up,mmax,mch)
   else
     call lschfs(ll,ierr,ep(l1), &
&                rr,vfull,uu,up,zz,mmax,n2,srel)
     call lschvkbs(ll,nproj(l1),ep(l1),rr,vp(1,lloc+1), &
&                  vkb(1,1,l1),evkb(1,l1),u2,up,mmax,n2)
   end if
  sgnae=sign(1.0d0,uu(n2))
  sgnps=sign(1.0d0,u2(n2))
   do ii = n1, n2
     write(6,'(a,i6,3f12.6)') '&',ll,rr(ii),sgnae*uu(ii),&
&     sgnps*u2(ii)
   end do

 end do
!
! section for extra valence wavefunction
!
 if(nxtra > 0) then
   do jj = 1, nxtra
    ll = lpx(jj)
    etest = epx(jj)
    call lschfb(npx(jj),ll,ierr,etest, &
&               rr,vfull,uu,up,zz,mmax,mch,srel)
    l1 = ll + 1
    sls=(l1-1)*l1
    emax=vp(mmax,l1)+0.5d0*sls/rr(mmax)**2
    emin=emax
    do ii=1,mmax
      emin=dmin1(emin,vp(ii,l1)+0.5d0*sls/rr(ii)**2)
    end do
    call lschvkbb(ll+2,ll,nproj(l1),ierr,etest,emin,emax, &
&                 rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                 u2,up,mmax,mch)
    sgnae=sign(1.0d0,uu(n2))
    sgnps=sign(1.0d0,u2(n2))
    write(6,'(//a,i2,a,i2,a,a/)') 'n=',ll+2,',  l=',ll, &
&    ', all-electron wave function', ', pseudo w-f'
    do ii = n1, n2
      write(6,'(a,i6,3f12.6)') '&',ll+4,rr(ii), &
&      sgnae*uu(ii),sgnps*u2(ii)
    end do
   end do
 end if

 do l1=1,lmax+1
   ll=l1-1
   if(nproj(l1)==2 .and. ll/=lloc) then
     write(6,'(//a,2i2,a,i2,a/)') 'n=',ll+1,ll+2,'  l=',ll, &
&     ', projecctor pseudo wave functions, well or 2nd valence'
     do ii = n1, n2
       write(6,'(a,i6,3f12.6)') '@',ll,rr(ii),pswf(ii,1,l1),&
&            pswf(ii,2,l1)
     end do
   end if 
 end do

! convergence profile plots

 write(6,'(/a)') 'convergence profiles, (ll=0,lmax)'
 do l1=1,lmax+1
   ll=l1-1
   do jj=1,4
     if(cvgplt(1,jj,1,l1)/=0.d0) then
       write(6,'(a,i6,3f12.6)') '!C',ll,cvgplt(1,jj,1,l1),cvgplt(2,jj,1,l1)
     end if
   end do
 end do

 deallocate(uu,u2,up)
 return
 end subroutine run_plot
