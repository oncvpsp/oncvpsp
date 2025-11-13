!
! Copyright (c) 1989-201r by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
subroutine sratom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
                  rr,vi,vxc,zz,mmax,iexc,etot,ierr,srel,u,up)
!! self-consistent scalar-relativistic all-electron atom
!! calculation using log mesh (non-relativistic when srel=.false.)
   implicit none
   integer, parameter :: dp=kind(1.0d0)

   ! Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> nc  number of core states
   integer, intent(in) :: nc
   !> ncv  number of core+valence states
   integer, intent(in) :: ncv
   !> na  principal quantum number array, dimension ncv
   integer, intent(in) :: na(ncv)
   !> la  angular-momenta
   integer, intent(in) :: la(ncv)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> fa  occupancies
   real(dp), intent(in) :: fa(ncv)
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> srel  .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   ! Output variables
   !> it  number of iterations (output)
   integer, intent(out) :: it
   !> ierr  error flag
   integer, intent(out) :: ierr
   !> etot  all-electron total energy (output)
   real(dp), intent(out) :: etot
   !> ea  eigenvalues (output)
   real(dp), intent(out) :: ea(ncv)
   !> rpk  radius of outermost peak of wave function
   real(dp), intent(out) :: rpk(ncv)
   real(dp), intent(out) :: rho(mmax)
   real(dp), intent(out) :: rhoc(mmax)
   !> vi  all-electron potential (output)
   real(dp), intent(out) :: vi(mmax)
   real(dp), intent(out) :: vxc(mmax)
   real(dp), intent(out) :: u(mmax, ncv)
   real(dp), intent(out) :: up(mmax, ncv)

   ! Local function
   real(dp) :: tfapot

   ! Local variables
   integer :: nin,mch
   real(dp) :: amesh,al
   real(dp) :: dr,eeel,eexc,et,rl,rl1,sd,sf,sn,eeig
   real(dp) :: thl,vn,zion
   integer :: ii,jj
   logical :: convg

   real(dp), allocatable :: vo(:),vi1(:),vo1(:)

   ! blend parameter for Anderson iterative potential mixing
   real(dp), parameter ::  bl=0.5d0

   allocate(vo(mmax),vi1(mmax),vo1(mmax))

   ! why all this is necessary is unclear, but it seems to be
   u(:,:)=0.d0; up(:,:)=0.d0; vo(:)=0.d0; vi1(:)=0.d0; vo1(:)=0.d0; vxc(:)=0.d0
   dr=0.d0; eeel=0.d0; eexc=0.d0; et=0.d0; rl=0.d0; rl1=0.d0
   sd=0.d0; sf=0.d0; sn=0.d0; eeig=0.d0; thl=0.d0; vn=0.d0; zion=0.d0
   nin=0; mch=0

   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = dexp(al)

   do ii=1,mmax
      vi(ii)=tfapot(rr(ii),zz)
   end do

   ! starting approximation for energies
   sf=0.0d0
   do ii=1,ncv
      sf=sf+fa(ii)
      zion=zz+1.0d0-sf
      ea(ii)=-0.5d0*(zion/na(ii))**2
      if(ea(ii)>vi(mmax)) ea(ii)=2.0d0*vi(mmax)
   end do

   ! big self  self-consietency loop

   do it=1,100
      convg=.true.

      rhoc(:) = 0.0d0
      rho(:)=0.0d0

      ! solve for bound states in turn
      eeig=0.0d0
      do ii=1,ncv

         ! skip unoccupied states
         if(fa(ii)==0.0d0) then
            ea(ii)=0.0d0
            cycle
         end if
         et=ea(ii)
         ierr = 0
         call lschfb(na(ii),la(ii),ierr,et,rr,vi,u(:,ii),up(:,ii),zz,mmax,mch,srel)
         if(ierr .ne. 0) then
            write(6,'(/a,3i4)') 'sratom123: lschfb convergence ERROR n,l,iter=', &
               &       na(ii),la(ii),it
            stop
         end if

         ! overall convergence criterion based on eps within lschfb
         if(ea(ii)/= et) convg=.false.
         ea(ii)=et

         ! accumulate charge and eigenvalues
         eeig = eeig + fa(ii) * ea(ii)
         rho(:)=rho(:) + fa(ii)*(u(:,ii)/rr(:))**2
         if(ii<=nc) then
            rhoc(:)=rhoc(:) + fa(ii)*(u(:,ii)/rr(:))**2
         end if


         ! find outermost peak of wavefunction
         do jj=mch-1,1,-1
            if(up(jj,ii)*up(jj+1,ii)<0.0d0) then
               rpk(ii)=rr(jj)
               exit
            end if
         end do

      end do

      if(ierr/=0) then
         exit
      end if


      ! output potential
      call vout(0,rho,rhoc,vo,vxc,sf-zz,eeel,eexc, &
         &            rr,mmax,iexc)

      etot =  eeig + eexc - 0.5d0*eeel

      ! generate next iteration using d. g. anderson''s
      ! method
      thl=0.0d0
      if(it>1) then
         sn=0.0d0
         sd=0.0d0
         do ii=1,mmax
            rl=vo(ii)-vi(ii)
            rl1=vo1(ii)-vi1(ii)
            dr=rl-rl1
            sn=sn + rl*dr*rr(ii)**2
            sd=sd + dr*dr*rr(ii)**2
         end do
         thl=sn/sd
      end if

      do ii=1,mmax
         vn=(1.0d0-bl)*((1.0d0-thl)*vi(ii) + thl*vi1(ii)) &
            &   + bl*((1.0d0-thl)*vo(ii) + thl*vo1(ii))
         vi1(ii)=vi(ii)
         vo1(ii)=vo(ii)
         vi(ii)=vn
      end do

      if(convg) exit

      if(it==100 .and. .not. convg) then
         write(6,'(/a)') 'sratom: WARNING failed to converge'
      end if

   end do !it

   if(.not. convg .and. ierr==0) then
      ierr=100
   end if

   ! total energy output

   ! output potential for e-e interactions

   call vout(0,rho,rhoc,vo,vxc,sf,eeel,eexc, &
      &          rr,mmax,iexc)

   etot =  eeig + eexc - 0.5d0*eeel


   deallocate(vo,vi1,vo1)
   return

end subroutine sratom
