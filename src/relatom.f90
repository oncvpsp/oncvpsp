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
 subroutine relatom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&           rr,vi,zz,mmax,iexc,etot,ierr)

! self-consistent fully-relativistic all-electron atom
! partially-occupied orbitals are weighted proportionally to 2j+1
! calculation using log mesh (non-relativistic when srel=.false.)

!na  principal quantum number array, dimension ncv
!la  angular-momenta
!ea  eigenvalues (output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nc  number of core states
!ncv  number of core+valence states
!it  number of iterations (output)
!rr  log radial mesh
!vi  all-electron potential (output)
!zz  atomic number
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  all-electron total energy (output)
!ierr  error flag
!srel  .true. for scalar-relativistic, .false. for non-relativistic

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 
 integer :: mmax,iexc,nc,ncv
 integer :: na(ncv),la(ncv)
 real(dp) :: zz
 real(dp) :: fa(ncv),rr(mmax)

!Output variables
 integer :: it,ierr
 real(dp) :: etot
 real(dp) :: ea(30,2),rpk(30,2)
 real(dp) :: rho(mmax),rhoc(mmax),vi(mmax)

!Local function
 real(dp) :: tfapot

!Local variables
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,rl,rl1,sd,sf,sn,eeig
 real(dp) :: thl,vn,zion,fj
 integer :: ii,jj,ll,ikap,kap,mkap
 logical :: convg

 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp), allocatable :: vo(:),vi1(:),vo1(:),vxc(:)

! blend parameter for Anderson iterative potential mixing
 real(dp), parameter ::  bl=0.5d0

 allocate(uu(mmax,2),up(mmax,2))
 allocate(vo(mmax),vi1(mmax),vo1(mmax),vxc(mmax))

! why all this is necessary is unclear, but it seems to be
 uu(:,:)=0.d0; up(:,:)=0.d0; vo(:)=0.d0; vi1(:)=0.d0; vo1(:)=0.d0; vxc(:)=0.d0
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
   do ikap=1,2
     ea(ii,ikap)=-0.5d0*(zion/na(ii))**2
     if(ea(ii,ikap)>vi(mmax)) ea(ii,ikap)=2.0d0*vi(mmax)
   end do
   if(la(ii)==0) ea(ii,2)=0.0d0
 end do

! big self  self-consietency loop

 do it=1,100
   convg=.true.

   rhoc(:) = 0.0d0
   rho(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do ii=1,ncv
    ll = la(ii)
    if(ll==0) then
     mkap=1
    else
     mkap=2
    end if
! loop on J = ll +/- 1/2
    do ikap=1,mkap
     if(ikap==1) then
       kap=-(ll+1)
       fj=dble(ll+1)*fa(ii)/dble(2*ll+1)
     else
       kap=  ll
       fj=dble(ll)*fa(ii)/dble(2*ll+1)
     end if
     et=ea(ii,ikap)
     ierr = 0

     call ldiracfb(na(ii),ll,kap,ierr,et, &
&                rr,zz,vi,uu,up,mmax,mch)

     if(ierr>0) then
       write(6,'(/2a,4i4)') 'relatom: ldiracfb convergence error', &
&           ' n,l,kap,iter=', &
&       na(ii),ll,kap,it
       exit
     end if

! overall convergence criterion based on eps within lschfb
     if(ea(ii,ikap)/= et) convg=.false.
     ea(ii,ikap)=et

! accumulate charge and eigenvalues
     eeig = eeig + fj*ea(ii,ikap)
     rho(:)=rho(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
     if(ii<=nc) then
      rhoc(:)=rhoc(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
     end if


! find outermost peak of wavefunction
     do jj=mch-1,1,-1
       if(up(jj,1)*up(jj+1,1)<0.0d0) then
         rpk(ii,ikap)=rr(jj)
         exit
       end if
     end do
    end do
   end do

   if(ierr/=0) then
    exit
   end if


! output potential
   call vout(0,rho,rhoc,vo,vxc,sf-zz,eeel,eexc, &
&            rr,mmax,iexc)

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
     write(6,'(/a)') 'relatom: WARNING failed to converge'
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

 deallocate(uu,up)
 deallocate(vo,vi1,vo1,vxc)
 return

 end subroutine relatom
