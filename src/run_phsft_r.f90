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
 subroutine run_phsft_r(lmax,lloc,nproj,ep,epsh1,epsh2,depsh,vkb,evkb, &
&                     rr,vfull,vp,zz,mmax,irphs)

! computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
! at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
! results for 1 and 2 projectors, or the semi-local pseudpotential
! when that is the local potential for some l
! the computed quantity is reminiscent of a scattering phase shift, but isn't
! This version is for fully-relativistic pseudopotentials

!lmax  maximum angular momentum
!lloc  l for local potential
!nproj  number ov VKB projectors for  each l
!ep  bound-state or scattering state reference energies for vkb potentials
!epsh1  low energy limit for "phase shift" calculation
!epsh2  high energy limit for "phase shift" calculation
!depsh  energy increment
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!zz  atomic number
!mmax  size of radial grid
!irphs  index of rr beyond which all vp==vlocal

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,mmax,irphs
 integer :: nproj(6)
 real(dp) :: epsh1,epsh2,depsh,zz
 real(dp) :: rr(mmax),vp(mmax,5,2),ep(6,2)
 real(dp) :: vfull(mmax),vkb(mmax,2,4,2),evkb(2,4,2)

!Output variables - printing only

!Local variables
 integer :: ii,ll,l1,npsh
 integer :: ikap,kap,mkap
 real(dp) :: epsh

 real(dp),allocatable :: pshf(:),pshp(:)

 npsh=(epsh2-epsh1)/depsh+1

 allocate(pshf(npsh),pshp(npsh))

! loop for phase shift calculation -- full, then local or Kleinman-
! Bylander / Vanderbilt
 
 do l1 = 1, lmax+1
   ll = l1 - 1
  if(l1==1) then
   mkap=1
  else
   mkap=2
  end if
! loop on J = ll +/- 1/2
  do ikap=1,mkap
   if(ikap==1) kap=-(ll+1)
   if(ikap==2) kap=  ll
   call fphsft_r(ll,kap,epsh2,depsh,pshf,rr,vfull,zz,mmax,irphs,npsh)
   if(ll .eq. lloc) then  
     call  vkbphsft(ll,0,epsh2,depsh,ep(l1,ikap),pshf,pshp, &
&                   rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   mmax,irphs,npsh)
   else
     call  vkbphsft(ll,nproj(l1),epsh2,depsh,ep(l1,ikap),pshf,pshp, &
&                   rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   mmax,irphs,npsh)
   end if

   write(6,'(/a,i2)') 'log derivativve data for plotting, l=',ll
   write(6,'(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=',rr(irphs)
   write(6,'(a/)') 'l, energy, all-electron, pseudopotential'
   do ii = 1, npsh
    epsh = epsh2 - depsh * dfloat(ii - 1)
    if(ikap==1) then
     write(6,'(a, i6, 3f12.6)') '! ',-ll, epsh, pshf(ii), pshp(ii)
    else
     write(6,'(a, i6, 3f12.6)') '! ',ll, epsh, pshf(ii), pshp(ii)
    end if
   end do
  end do !ikap
 end do !l1
 deallocate(pshf,pshp)
 return
 end subroutine run_phsft_r
