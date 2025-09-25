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
! self-consistent all-electron and pseudopotential atomic calculations
! compared for reference and tests atomic configurations


subroutine run_config_r(jj,nacnf,lacnf,facnf,nc,nvcnf,rhov,rhomod,rr,zz, &
&                  mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
&                  lloc,vkb,evkb)

!jj  index of current configufation
!nacnf  principal quantum number array, all configurations
!lacnf  angular-momenta array, all config.
!nc  number of core states
!nvcnf  number of valence states, all config.
!rhov valence pseudo-charge of reference configuration
!rhomod  model core charge
!rr  log radial mesh
!zz  atomic number
!rcmax  maximum core radius for psp
!mmax  size of log grid
!mxprj  dimension of number of projectors
!iexc  exchange-correlation function to be used
!ea  reference configurattion eigenvalues
!etot  reference configuration total energy
!epstot  pseudoatom total energy
!nproj  number of VKB projectors to use for each l
!vpuns  unscreened semi-local pseudopotentials (plus different vloc if lloc==4)
!lloc  index-1 of local potential
!vkb   Vanderbilt-Kleinman-Bylander projectors
!evkb VKB projector coefficients
!srel .true. for scalar-relativistic, .false. for non-relativistic

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables

   integer :: jj,mmax,mxprj,iexc,nc,lloc
   integer :: nacnf(30,5),lacnf(30,5),nvcnf(5),nproj(5)
   real(dp) :: etot,epstot,zz
   real(dp) :: facnf(30,5),ea(30,2),rhov(mmax),rr(mmax)
   real(dp) :: vpuns(mmax,5),vkb(mmax,mxprj,4,2),evkb(mxprj,4,2),rhomod(mmax,5)

!Output variables  only printing

!Local variables
   integer :: ii,it,kk,ll,l1,ierr,mch,nvt
   integer :: ikap,kap,mkap
   integer :: nat(30),lat(30),natp(30),latp(30),nav(4)
   integer ::  indxr(30),indxe(30)
   real(dp) :: et,eaetst,etsttot,fj
   real(dp) :: eat(30,2,3),fat(30,3),rpk(30,2),eatp(30,2),fatp(30,3)

   real(dp),allocatable :: rho(:),rhoc(:),rhocps(:),vi(:),vfull(:)
   real(dp),allocatable :: uu(:,:),up(:,:)

   allocate(rho(mmax),rhoc(mmax),vi(mmax),vfull(mmax),rhocps(mmax))
   allocate(uu(mmax,2),up(mmax,2))

! atom tests comparing reference state and excited configuration for
! all-electron and pseudo atoms.  Total excitation energies and excited-
! -configuration eigenvalues are compared

! For multi-projector potentials, the bound-state solver lschvkbb has poor
! stability, and requires a trial energy sufficiently close to the correct
! result for a given local potential. This appears to be due to the fact
! that intermediate stages towards the solution of the radial Schroedinger
! equation with a non-local potential may not increase their node count
! monotonically with increasing energy (the same phenonenon which gives
! rise to ghost states).

! To aviod these problems, the pseudoatom calculation is started with the
! screened local potential and valence eigenvalues of the reference
! configuration.  This is adiabatically changed in two steps.  In the first,
! the occupation numbers are changed by 2% within the potential iteration
! loop for each of the first 50 steps to those of a maximally-ionized
! intermediate configuration.

! The second step is started with the local potential of this maximally-
! ionized psuedo-atom, and eigenvalues from an all-electron
! calculation of this configuration including unoccupied states that
! are to be adiabatically filled to produce the final configuration.
! The same 2% strategy is used filling pseudo-atom states in this step.

! set up pseudo atom arrays
!arrange data arrays for adiabatic configufation switching of pseudo atom

   fat(:,:)=0.0d0
   nat(:)=0
   lat(:)=0

! set up cross indces of reference (r) and excited (e) cofigurations
   indxr(:)=0
   indxe(:)=0
   do ii=1,nc+nvcnf(1)
      do kk=1,nc+nvcnf(jj)
         if(nacnf(kk,jj)==nacnf(ii,1) .and. &
         &         lacnf(kk,jj)==lacnf(ii,1)) then
            indxr(ii)=kk
            indxe(kk)=ii
         end if
      end do
   end do

! set up merged list with reference levels first
! levels not present in one confguration or the other get zero occupancy
   eat(:,:,:)=0.0d0
   do ii=1,nc+nvcnf(1)
      nat(ii)=nacnf(ii,1)
      lat(ii)=lacnf(ii,1)
      fat(ii,1)=facnf(ii,1)
      eat(ii,:,1)=ea(ii,:)
      if(indxr(ii)>0) then
         fat(ii,2)=facnf(indxr(ii),jj)
      else
         fat(ii,2)=0.0d0
      end if
      fat(ii,3)=fat(ii,2)
   end do

   nvt=nvcnf(1)
   do kk=1,nc+nvcnf(jj)
      if(indxe(kk)==0) then
         nvt=nvt+1
         nat(nc+nvt)=nacnf(kk,jj)
         lat(nc+nvt)=lacnf(kk,jj)
         fat(nc+nvt,1)=0.0d0
         fat(nc+nvt,3)=facnf(kk,jj)
      end if
   end do

!set up array for maximally ionized intermediate state
   do kk=1,nc+nvt
      fat(ii,2)=dmin1(fat(ii,1),fat(ii,3))
   end do



   eat(:,:,2)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
   eaetst=0.d0

!all-electron atom solution for maximally-ionized state

   call relatom(nat,lat,eat(1,1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
   &            rr,vfull,zz,mmax,iexc,eaetst,ierr)
   if(ierr/=0) then
      write(6,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
      &       ' WARNING no output for configuration',jj
      if(ierr==-1) then
         write(6,'(a,i4)') &
         &         'run_config_r: WARNING no classical turning point error, iteration',it
      else
         write(6,'(a)') 'run_config_r: WARNING self-consistency failed to converge'
      end if
      deallocate(rho,rhoc,rhocps,vi,vfull)
      deallocate(uu,up)
      return
   end if

!fill in energies for empty levels of maximally ionized state

   do kk=1,nc+nvt
      if(eat(kk,1,2)==0.0d0) then
         et=0.0d0
         ll=lat(kk)
         kap=ll+1
         call ldiracfb(nat(kk),ll,kap,ierr,et, &
         &                rr,zz,vfull,uu,up,mmax,mch)
         if(ierr /= 0) then
            write(6,'(/a,4i4)') &
            &            'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
            &            nat(kk),ll,kap,ierr
            deallocate(rho,rhoc,rhocps,vi,vfull)
            deallocate(uu,up)
            return
         end if
         eat(kk,1,2)=et
         if(ll==0) then
            eat(kk,2,2)=et
         else
            kap=ll
            call ldiracfb(nat(kk),ll,kap,ierr,et, &
            &                  rr,zz,vfull,uu,up,mmax,mch)
            if(ierr>0) then
               write(6,'(/a,4i4)') &
               &              'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
               &              nat(kk),ll,kap,it
               deallocate(rho,rhoc,rhocps,vi,vfull)
               deallocate(uu,up)
               return
            end if
            eat(kk,2,2)=et
         end if
      end if
   end do

   eat(:,:,3)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
   eaetst=0.d0

!all-electron atom solution for excited state

   call relatom(nat,lat,eat(1,1,3),fat(1,3),rpk,nc,nc+nvt,it,rhoc,rho, &
   &            rr,vfull,zz,mmax,iexc,eaetst,ierr)
   if(ierr/=0) then
      write(6,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
      &       ' WARNING no output for configuration',jj
      if(ierr==-1) then
         write(6,'(a,i4)') &
         &        'run_config_r: WARNING no classical turning point error, iteration',it
      else
         write(6,'(a)') 'run_config_r: WARNING self-consistency failed to converge'
      end if
      deallocate(rho,rhoc,rhocps,vi,vfull)
      deallocate(uu,up)
      return
   end if


!first pseudoatom run from reference to maximally-ionized configuration
   do l1=1,4
      nav(l1)=l1-1
   end do

   do kk=1,nvt
      latp(kk)=lat(nc+kk)
      l1=latp(kk)+1
      nav(l1)=nav(l1)+1
      natp(kk)=nav(l1)
      eatp(kk,:)=eat(nc+kk,:,1)
      fatp(kk,1)=fat(nc+kk,1)
      fatp(kk,2)=fat(nc+kk,2)
   end do

   rho(:)=rhov(:)

   rhocps(:)=rhomod(:,1)

   call psatom_r(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
   &           rr,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
   &           vkb,evkb,ierr)

   if(ierr/=0) then
      write(6,'(a/a,i2)') 'run_config_r: WARNING for fully non-local PS atom,', &
      &       ' WARNING no output for configuration',jj
      deallocate(rho,rhoc,rhocps,vi,vfull)
      return
   end if

!second pseudoatom run from maximally-ionized to excited configuration

   do kk=1,nvt
      eatp(kk,:)=eat(nc+kk,:,2)
      fatp(kk,1)=fat(nc+kk,2)
      fatp(kk,2)=fat(nc+kk,3)
   end do

   call psatom_r(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
   &           rr,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
   &           vkb,evkb,ierr)

   if(ierr/=0) then
      write(6,'(a/a,i2)') 'run_config_r: WARNING for fully non-local PS atom,', &
      &       ' WARNING no output for configuration',jj
      deallocate(rho,rhoc,rhocps,vi,vfull)
      return
   end if


   write(6,'(/a)') '   n   l kap    f        eae           eps        diff'
   do ii=1,nc
      ll = lat(ii)
      if(ll==0) then
         mkap=1
      else
         mkap=2
      end if
! loop on J = ll +/- 1/2
      do ikap=1,mkap
         if(ikap==1) then
            kap=-(ll+1)
            fj=dble(ll+1)*fat(ii,3)/dble(2*ll+1)
         else
            kap=  ll
            fj=dble(ll)*fat(ii,3)/dble(2*ll+1)
         end if
         write(6,'(3i4,f8.4,f14.8)') nat(ii),lat(ii),kap,fj,eat(ii,ikap,3)
      end do
   end do

   do ii=1,nvt
      ll = lat(ii+nc)
      if(ll==0) then
         mkap=1
      else
         mkap=2
      end if
! loop on J = ll +/- 1/2
      do ikap=1,mkap
         if(ikap==1) then
            kap=-(ll+1)
            fj=dble(ll+1)*fat(ii+nc,3)/dble(2*ll+1)
         else
            kap=  ll
            fj=dble(ll)*fat(ii+nc,3)/dble(2*ll+1)
         end if
         write(6,'(3i4,f8.4,2f14.8,1p,e12.2)') nat(ii+nc),lat(ii+nc),kap,fj, &
         &    eat(ii+nc,ikap,3),eatp(ii,ikap),eatp(ii,ikap)-eat(ii+nc,ikap,3)
      end do
   end do


   write(6,'(/a)') '    Total energies and differences'
   write(6,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      AE_ref=',etot, &
   & '  AE_tst=',eaetst,'  dif=',eaetst-etot
   write(6,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      PS_ref=',epstot, &
   & '  PS_tst=',etsttot,'  dif=',etsttot-epstot
   write(6,'(a,1p,e10.2)') '      PSP excitation error=', &
   & eaetst-etot-etsttot+epstot

   deallocate(rho,rhoc,rhocps,vi,vfull)
   deallocate(uu,up)
   return
end subroutine run_config_r
