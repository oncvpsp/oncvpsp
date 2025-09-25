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
subroutine run_plot(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,mxprj,drl,nrl, &
&                    rho,rhoc,rhomod,srel,cvgplt)


! write output for plotting pseudopotentials, model core charge, and
! all-electron and pseudo wave functions

!lmax  maximum angular momentum
!npa  principal quantum number for corresponding all-electron state
!epa  bound-state or scattering state reference energies for vkb potentials
!lloc  l for local potential
!irc  indices of core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials
!vpuns  unscreened vp
!zz  atomic number
!mmax  size of radial grid
!mxprj  dimension of number of projector
!drl  spacing of linear radial mesh
!nrl  number of points in radial mesh
!rho  valence pseudocharge
!rhoc core charge
!rhomod  model core charge
!srel .true. for scalar-relativistic, .false. for non-relativistic
!cvgplt  Energy per electron error vs. cutoff

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: lmax,lloc,mmax,mxprj,nrl
   integer :: npa(mxprj,6),irc(6),nproj(6)
   real(dp) :: zz,drl
   real(dp) :: rr(mmax),vp(mmax,5),vpuns(mmax,5),vfull(mmax),vkb(mmax,mxprj,4)
   real(dp) :: rho(mmax),rhoc(mmax),rhomod(mmax,5)
   real(dp):: epa(mxprj,6),evkb(mxprj,4),cvgplt(2,7,mxprj,4)
   logical :: srel

!Output variables - printing only

!Local variables
   integer :: ll,l1,ii,jj,ierr,mch,n1,n2,n3,n4,nn
   integer :: iprj,npr
   real(dp) :: al,emax,emin,etest,rmx,sgnae,sgnps
   real(dp) :: r0,dr
   real(dp), allocatable :: uu(:),u2(:),up(:)

   allocate(uu(mmax),u2(mmax),up(mmax))

   al = 0.01d0 * dlog(rr(101)/rr(1))

   write(6,'(/a)') 'DATA FOR PLOTTING'

   n1 = idint(dlog(drl/rr(1))/al+1.0d0)
   n2 = idint(dlog(dfloat(nrl)*drl/rr(1))/al+1.0d0)
   n3=0
   do l1=1,lmax+1
      n3=max(n3,irc(l1)-1)
   end do
   n4=min(n2,idint(dlog((rr(n3)+1.0d0)/rr(1))/al))
   n3=idint(dlog(1.1d0*rr(n3)/rr(1))/al+1.0d0)

   dr=dmin1(drl,rr(n4)/200)
   write(6,'(/2a/)') ' radii, charge,',' pseudopotentials (ll=0, 1, lmax)'
   r0=rr(n1)-dr
   do ii=n1,n4
      if(rr(ii)>r0) then
         write(6,'(a,6(f12.7,1x))') '!p',rr(ii),rho(ii),(vpuns(ii,l1),l1=1,lmax+1)
         r0=rr(ii)+dr
      end if
   end do

   if(lloc == 4) then
      r0=rr(n1)-dr
      do ii=n1,n4
         if(rr(ii)>r0) then
            write(6,'(a, 2(f12.7,1x))') ' !L',rr(ii),vpuns(ii,lloc + 1)
            r0=rr(ii)+dr
         end if
      end do
   end if

! output for analysis of core charge and models
   rmx=0
   do ii=n1,n2
      rmx=max(rmx,5.0d0*rho(ii),rhomod(ii,1))
   end do

   write(6,'(//2a/)') ' radii, charge,',' core charge, model core charge'
   dr=dmin1(drl,rr(n2)/200.0d0)
   r0=rr(n1)-dr
   do ii=n1,n2
      if(rr(ii)>r0) then
         if(rhoc(ii)<rmx)then
            write(6,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rhoc(ii), &
            &       rhomod(ii,1)
         else if(rr(ii)>0.01d0) then
            write(6,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rmx, &
            &       rhomod(ii,1)
         end if
         r0=rr(ii)+dr
      end if
   end do

! loop for wave function output

   do l1 = 1, lmax + 1
      ll = l1 - 1
      npr=nproj(l1)
      if(l1==lloc+1) npr=0

      do iprj=1,nproj(l1)
         if(epa(iprj,l1)<0.0d0) then
            write(6,'(//a,i2,a,i2,a,a/)') 'n=',npa(iprj,l1),',  l=',ll, &
            &           ', all-electron wave function', ', pseudo w-f'
            etest=epa(iprj,l1)
            call lschfb(npa(iprj,l1),ll,ierr,etest, &
            &                  rr,vfull,uu,up,zz,mmax,mch,srel)
            if(ierr /= 0) then
               write(6,'(/a,3i4)') 'run_plot: lschfb convergence ERROR n,l,ierr=', &
               &         npa(iprj,l1),ll,ierr
               stop
            end if

            emax=0.9d0*etest
            emin=1.1d0*etest

            call lschvkbb(ll+iprj,ll,npr,ierr,etest,emin,emax, &
            &                    rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
            &                    u2,up,mmax,mch)
            sgnae=1.0d0
            sgnps=1.0d0
         else
            write(6,'(//a,i2,a,i2,a,a/)') 'scattering, iprj=',iprj,',  l=',ll, &
            &           ', all-electron wave function', ', pseudo w-f'
            call lschfs(nn,ll,ierr,epa(iprj,l1), &
            &                  rr,vfull,uu,up,zz,mmax,n2,srel)
            call lschvkbs(ll,npr,epa(iprj,l1),rr,vp(1,lloc+1), &
            &                    vkb(1,1,l1),evkb(1,l1),u2,up,mmax,n2)
            sgnae=sign(1.0d0,uu(n2))
            sgnps=sign(1.0d0,u2(n2))
         end if

         dr=dmin1(drl,rr(n2)/200.0d0)
         r0=rr(n1)-dr
         do ii = n1, n2
            if(rr(ii)>r0) then
               write(6,'(a,i5,i1,3(f12.6,1x))') '&',iprj,ll,rr(ii),sgnae*uu(ii),&
               &          sgnps*u2(ii)
               r0=rr(ii)+dr
            end if
         end do
      end do  !iproj

! orthonormal projector plots

      write(6,'(/)')
      dr=dmin1(drl,rr(n3)/200.0d0)
      r0=rr(n1)-dr
      do ii=n1,n3
         if(rr(ii)>r0) then
            write(6,'(a,i6,6(f12.6,1x))') '!J',ll,rr(ii), &
            &       (vkb(ii,jj,l1),jj=1,nproj(l1))
            r0=rr(ii)+dr
         end if
      end do

   end do  !l1
!
! convergence profile plots

   write(6,'(/a)') 'convergence profiles, (ll=0,lmax)'
   write(6,*) 'lmax',lmax
   do l1=1,lmax+1
      ll=l1-1
      do jj=1,7
         if(cvgplt(1,jj,1,l1)/=0.d0) then
            write(6,'(a,i6,3(f12.6,1x))') '!C',ll,cvgplt(1,jj,1,l1),cvgplt(2,jj,1,l1)
         end if
      end do
   end do

   deallocate(uu,u2,up)
   return
end subroutine run_plot
