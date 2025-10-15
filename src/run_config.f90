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

!> self-consistent all-electron and pseudopotential atomic calculations
!> compared for reference and tests atomic configurations
!> atom tests comparing reference state and excited configuration for
!> all-electron and pseudo atoms.  Total excitation energies and excited-
!> -configuration eigenvalues are compared
!>
!> For multi-projector potentials, the bound-state solver lschvkbb has poor
!> stability, and requires a trial energy sufficiently close to the correct
!> result for a given local potential. This appears to be due to the fact
!> that intermediate stages towards the solution of the radial Schroedinger
!> equation with a non-local potential may not increase their node count
!> monotonically with increasing energy (the same phenonenon which gives
!> rise to ghost states).
!>
!> To aviod these problems, the pseudoatom calculation is started with the
!> screened local potential and valence eigenvalues of the reference
!> configuration.  This is adiabatically changed in two steps.  In the first,
!> the occupation numbers are changed by 2% within the potential iteration
!> loop for each of the first 50 steps to those of a maximally-ionized
!> intermediate configuration.
!>
!> The second step is started with the local potential of this maximally-
!> ionized psuedo-atom, and eigenvalues from an all-electron
!> calculation of this configuration including unoccupied states that
!> are to be adiabatically filled to produce the final configuration.
!> The same 2% strategy is used filling pseudo-atom states in this step.
subroutine run_config(jj,nacnf,lacnf,facnf,nc,nvcnf,rhov,rhomod,rr,zz, &
                      rcmax,mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
                      lloc,vkb,evkb,srel,nvt,nat,lat,fat,eat,eatp,eaetst,etsttot)
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   ! Input variables
   !> index of current configufation
   integer, intent(in) :: jj
   !> size of log grid
   integer, intent(in) :: mmax
   !> dimension of number of projectors
   integer, intent(in) :: mxprj
   !> exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> number of core states
   integer, intent(in) :: nc
   !> index-1 of local potential
   integer, intent(in) :: lloc
   !> principal quantum number array, all configurations
   integer, intent(in) :: nacnf(30,5)
   !> angular-momenta array, all config.
   integer, intent(in) :: lacnf(30,5)
   !> number of valence states, all config.
   integer, intent(in) :: nvcnf(5)
   !> number of VKB projectors to use for each l
   integer, intent(in) :: nproj(5)
   !> reference configuration total energy
   real(dp), intent(in) :: etot
   !> pseudoatom total energy
   real(dp), intent(in) :: epstot
   !> maximum core radius for psp
   real(dp), intent(in) :: rcmax
   !> atomic number
   real(dp), intent(in) :: zz
   !> occupation number array, all configurations
   real(dp), intent(in) :: facnf(30,5)
   !> reference configurattion eigenvaluess
   real(dp), intent(in) :: ea(30)
   !> valence pseudo-charge of reference configuration
   real(dp), intent(in) :: rhov(mmax)
   !> log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> unscreened semi-local pseudopotentials (plus different vloc if lloc==4)
   real(dp), intent(in) :: vpuns(mmax,5)
   !> Vanderbilt-Kleinman-Bylander projectors
   real(dp), intent(in) :: vkb(mmax,mxprj,4)
   !> VKB projector coefficients
   real(dp), intent(in) :: evkb(mxprj,4)
   !> model core charge
   real(dp), intent(in) :: rhomod(mmax,5)
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !> Output variables
   integer, intent(out) :: nvt
   !> principal quantum number array for test configuration
   integer, intent(out) :: nat(30)
   !> angular momentum array for test configuration
   integer, intent(out) :: lat(30)
   !> occupation number array for test configuration
   real(dp), intent(out) :: fat(30,3)
   !> all-electron energy of test state
   real(dp), intent(out) :: eat(30,3)
   !> pseudo energy of test state
   real(dp), intent(out) :: eatp(30)
   !> test configuration all-electron total energy
   real(dp), intent(out) :: eaetst
   !> test configuration pseudoatom total energy
   real(dp), intent(out) :: etsttot

   ! Local variables
   integer :: ii,it,kk,l1,ierr,mch
   integer :: natp(30),latp(30),nav(4)
   integer ::  indxr(30),indxe(30)
   real(dp) :: et
   real(dp) :: rpk(30),fatp(30,3)

   real(dp),allocatable :: rho(:),rhoc(:),rhocps(:),vi(:),vfull(:),vxc(:)
   real(dp),allocatable :: uu(:),up(:),uae(:,:),upae(:,:)

   allocate(rho(mmax),rhoc(mmax),vi(mmax),vfull(mmax),vxc(mmax),rhocps(mmax))
   allocate(uu(mmax),up(mmax),uae(mmax,30),upae(mmax,30))

   ! set up pseudo atom arrays
   ! arrange data arrays for adiabatic configufation switching of pseudo atom

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
   eat(:,1)=0.0d0
   do ii=1,nc+nvcnf(1)
      nat(ii)=nacnf(ii,1)
      lat(ii)=lacnf(ii,1)
      fat(ii,1)=facnf(ii,1)
      eat(ii,1)=ea(ii)
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

   ! set up array for maximally ionized intermediate state
   do kk=1,nc+nvt
      fat(ii,2)=dmin1(fat(ii,1),fat(ii,3))
   end do



   eat(:,2)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
   eaetst=0.d0

   ! all-electron atom solution for maximally-ionized state

   call sratom(nat,lat,eat(1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
      &        rr,vfull,vxc,zz,mmax,iexc,eaetst,ierr,srel,uae,upae)
   if(ierr/=0) then
      write(stdout,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
         &       ' WARNING no output for configuration',jj
      if(ierr==-1) then
         write(stdout,'(a,i4)') 'no classical turning point error, iteration',it
      else
         write(stdout,'(a)') 'run_config: WARNING self-consistency failed to converge'
      end if
      deallocate(rho,rhoc,rhocps,vi,vfull)
      deallocate(uu,up)
      return
   end if

   ! fill in energies for empty levels of maximally ionized state

   do kk=1,nc+nvt
      if(abs(eat(kk,2)) < tiny(0.0_dp)) then
         et=0.0d0
         call lschfb(nat(kk),lat(kk),ierr,et, &
            &                rr,vfull,uu,up,zz,mmax,mch,srel)
         if(ierr .ne. 0) then
            write(stdout,'(/a,3i4)') &
               &            'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
               &            nat(kk),lat(kk),it
            deallocate(rho,rhoc,rhocps,vi,vfull)
            deallocate(uu,up)
            return
         end if
         eat(kk,2)=et
      end if
   end do

   eat(:,3)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
   eaetst=0.d0

   ! all-electron atom solution for excited state

   call sratom(nat,lat,eat(1,3),fat(1,3),rpk,nc,nc+nvt,it,rhoc,rho, &
      &        rr,vfull,vxc,zz,mmax,iexc,eaetst,ierr,srel,uae,upae)
   if(ierr/=0) then
      write(stdout,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
         &       ' WARNING no output for configuration',jj
      if(ierr==-1) then
         write(stdout,'(a,i4)') 'no classical turning point error, iteration',it
      else
         write(stdout,'(a)') 'run_config: WARNING self-consistency failed to converge'
      end if
      deallocate(rho,rhoc,rhocps,vi,vfull)
      deallocate(uu,up)
      return
   end if


   ! first pseudoatom run from reference to maximally-ionized configuration
   do l1=1,4
      nav(l1)=l1-1
   end do

   do kk=1,nvt
      latp(kk)=lat(nc+kk)
      l1=latp(kk)+1
      nav(l1)=nav(l1)+1
      natp(kk)=nav(l1)
      eatp(kk)=eat(nc+kk,1)
      fatp(kk,1)=fat(nc+kk,1)
      fatp(kk,2)=fat(nc+kk,2)
   end do

   rho(:)=rhov(:)

   rhocps(:)=rhomod(:,1)

   call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
      &           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
      &           vkb,evkb,ierr)

   if(ierr/=0) then
      write(stdout,'(a,a/a,i2)') 'run_config: WARNING for fully non-local PS atom,', &
         &       ' stg. 1', &
         &       ' WARNING no output for configuration',jj
      deallocate(rho,rhoc,rhocps,vi,vfull)
      return
   end if

   ! second pseudoatom run from maximally-ionized to excited configuration

   do kk=1,nvt
      eatp(kk)=eat(nc+kk,2)
      fatp(kk,1)=fat(nc+kk,2)
      fatp(kk,2)=fat(nc+kk,3)
   end do
   do ii=1,nvt
   end do

   call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
      &           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
      &           vkb,evkb,ierr)

   if(ierr/=0) then
      write(stdout,'(a,a/a,i2)') 'run_config: WARNING for fully non-local PS atom,', &
         &       ' stg. 2', &
         &       ' WARNING no output for configuration',jj
      deallocate(rho,rhoc,rhocps,vi,vfull)
      return
   end if



   ! write(stdout,'(/a)') '   n   l     f        eae           eps        diff'
   ! do ii=1,nc
   !    write(stdout,'(2i4,f8.4,f14.8)') nat(ii),lat(ii),fat(ii,3),eat(ii,3)
   ! end do
   ! do ii=1,nvt
   !    if(abs(fat(ii+nc,3)) < tiny(0.0_dp)) cycle
   !    write(stdout,'(2i4,f8.4,2f14.8,1p,d12.2)') nat(ii+nc),lat(ii+nc),fat(ii+nc,3), &
   !       &  eat(ii+nc,3),eatp(ii),eatp(ii)-eat(ii+nc,3)
   ! end do

   ! write(stdout,'(/a)') '    Total energies and differences'
   ! write(stdout,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=',etot, &
   !    & '  AE_tst=',eaetst,'  dif=',eaetst-etot
   ! write(stdout,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=',epstot, &
   !    & '  PS_tst=',etsttot,'  dif=',etsttot-epstot
   ! write(stdout,'(a,1p,d10.2)') '      PSP excitation error=', &
   !    & eaetst-etot-etsttot+epstot

   deallocate(rho,rhoc,rhocps,vi,vfull,vxc)
   deallocate(uu,up,uae,upae)
   return
end subroutine run_config
