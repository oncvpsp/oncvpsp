!:
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
 subroutine run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
&                    vkb,evkb,nproj,rr,vp,mmax,mxprj)

! Two tests for ghosts;  The local potential is terminated by a hard wall
! barrier at a radius of mbfact*rc, presently 3*rc.  A basis set is formed
! from eigenstates this combination and used to compute a hamiltonian matrix
! whose off-diagonal terms come exclusively from the non-local pojectors.
! A cutoff is chosen to simulate a reasonably well converged plane-wave
! calculation.  The hamiltoian is diagonalized, and the results used
! as follows:
!
! Test 1) Negative eigenvalues are compared with those of the radial
! Schroedinger equation for the full-range pseudopotential.  Those 
! lying below the Schroedinger Eq. results and which presumably have
! more nodes are reported and tagged as GHOST(-).
!
! Test2) For barrier states with  positive eigenvalues, the average radius
! is computed.  Such states could be ghost resonances, and highly localized
! ones with <r>/rc<1 are reported as GHOST(+).  Real fairly localized
! resonances can occur, particularly at lower energies, and the GHOST(+)
! report is an indication that the log-derivative plot should be examined
! carefully around the reported energy to see if the all-electron plot
! has a corresponding drop of pi  (arctan(log der), of course).


!lmax  maximum angular momentum
!la  angular momenta of all-electron states
!nc  number of core electrons
!nv  number of valence electrons
!ea  all-electron bound-state energies
!lloc  l for local potential
!irc  core radii indices
!vkb  VKB projectors
!qmsbf maximum q in sbf basis for each l
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!mmax  size of radial grid
!mxprj  dimension of number of projectors

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
 integer :: nc,nv,lmax,lloc,mmax,mxprj
 integer :: la(30),irc(6),nproj(6)
 real(dp) :: rr(mmax),vp(mmax,5),vkb(mmax,mxprj,4)
 real(dp):: qmsbf(6),ea(30),evkb(mxprj,4)

!Output variables - printing only

!Local variables
 integer :: ll,l1,ii,jj,kk,ierr,mbfact,mbar,mch,nn,npgh
 integer :: iprj,nprj,nbas,nbmax,lwork,info
 integer :: nodes(4),lpgh(20)
 real(dp) :: al,amesh,basfact,ee,ecut,eps,et,emin,emax,sn,ro,tht
 real(dp) :: cpgh(20),ebcut(6),epgh(20),rpgh(20)
 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: gg0(:,:),hmat(:,:),hmev(:),work(:)
 real(dp), allocatable :: uua(:,:)

 allocate(uu(mmax),up(mmax))

 al = 0.01d0 * dlog(rr(101)/rr(1))
 amesh = dexp(al)
 eps=1.0d-4 !threshold for negative-energy ghost report
 tht=0.0d0 !sets barrier boundary conditions to zero

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors

   write(6,'(/a)') ' Testing for bound ghosts'
   write(6,'(5a)') '   n','   l','    E NL Schr. Eq', &
&        '   E Basis Diag.','   E Cutoff'

 nodes(:)=0
 npgh=0
 do l1 = 1, lmax + 1
   nprj=nproj(l1)  
   if(nprj==0) cycle

   write(6,'()')
   ll = l1 - 1

!Construct basis of local-pseidupotential plus hard-wall barrier eigenfunctions
!Barrier is set at mbfact*rc

   mbfact=3
   mbar = idint(dlog(mbfact*rr(irc(l1))/rr(1))/al+1.0d0)

! Basis cutoff energy is based on largest q in sbf optimization basis
  basfact=1.25d0
  ebcut(l1)= basfact*0.5d0*qmsbf(l1)**2

! Estimate maximum basis set size nbmax based on ll=0 sbf
! with 20% allowance for binding by local potential

   nbmax=idint(2.0d0*dsqrt(1.2d0*ebcut(l1))*rr(mbar)/pi)

   lwork=3*nbmax

   allocate(gg0(nprj,nbmax),hmat(nbmax,nbmax),hmev(nbmax),work(lwork))
   allocate(uua(mmax,nbmax))

   gg0(:,:)=0.0d0
   hmat(:,:)=0.0d0
   hmev(:)=0.0d0

   
   nbas=0
   do jj=1,nbmax
     nn=ll+jj

     emin=0.0d0
     emax=ebcut(l1)
     ee=0.0d0

     call lschpsbar(nn,ll,ierr,ee,emin,emax,rr,vp(1,lloc+1),uu,up, &
&                   mmax,mbar,tht)

     if(ierr/=0) then
         if(emax/=ebcut(l1)) then
         write(6,'(a,3i4,3f12.6)')  &
&              'run_ghosts 143: lschpsbar ERROR ierr,nn,ll,ee,emin,emax', &
&              ierr,nn,ll,ee,emin,emax
         stop
         end if
       nbas=jj-1
       exit
     end if

! barrier potential eigenvalues are contirbution to hamiltonian matrix diagonal
     hmat(jj,jj)=ee

     uua(:,jj)=uu(:)

! oerlap of eigenfunction jj and projectors
     do iprj=1,nprj
       call vpinteg(uu,vkb(1,iprj,l1),irc(l1),2*ll+2,gg0(iprj,jj),rr)
     end do

   end do !jj

   ecut=ebcut(l1)

! add non-local terms to hamiltonian matrix
   do jj=1,nbas
     do ii=1,jj
       do iprj=1,nprj
         hmat(ii,jj)=hmat(ii,jj)+evkb(iprj,l1)*gg0(iprj,ii)*gg0(iprj,jj)
       end do
       hmat(jj,ii)=hmat(ii,jj)
     end  do
   end do

! find eigenvalues of the hamiltonian matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', nbas, hmat, nbmax, hmev, work, lwork, info )
     if(info .ne. 0) then
      write(6,'(a,i4)') 'run_ghosts: hamiltonian matrix eigenvalue ERROR, &
&          info=',info
      stop
     end if
!treat valence states for this l
     nn=ll
     jj=0
     do kk=1,nv
       if(la(nc+kk)/=ll) cycle
       jj=jj+1
       nn=nn+1
       ee=ea(nc+kk)
       emax=0.9d0*ea(nc+kk)
       emin=1.1d0*ea(nc+kk)
       call lschvkbb(nn,ll,nprj,ierr,ee,emin,emax, &
&                    rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                    uu,up,mmax,mch)
       if(ierr/=0) then
         write(6,'(a,3i4,2f12.6)') 'run_ghosts: lschvkbb ERROR', &
&              nn,ll,ierr,ea(nc+kk),ee
       end if

       do ii=1,10
         et=dmin1(ee-eps,ee*(1.0d0+eps))
         if(hmev(jj)<et) then
           write(6,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&                '  WARNING - GHOST(-)'
           jj=jj+1
         else
           exit
         end if
       end do !ii

       write(6,'(2i4,2f16.6,f10.2)') nn,ll,ee,hmev(jj),ecut
     end do !kk

!if there were no valence states, look for any nodeless bound state
!if this fails, use vacuum level (0.0) to compare to basis-set states

     if(jj==0) then
       ee=-1.0d5
       do kk=1,nv
         ee=max(ee,ea(nc+kk))
       end do
       emax=0.0d0
       emin=2.0d0*ee
       call lschvkbb(ll+1,ll,nprj,ierr,ee,emin,emax, &
&                    rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                    uu,up,mmax,mch)
       if(ierr/=0) ee=0.0d0
       do jj=1,10
         et=dmin1(ee-eps,ee*(1.0d0+eps))
         if(hmev(jj)<et) then
           write(6,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&                '  WARNING - GHOST(-)'
         else if(ee<0.0d0) then
           write(6,'(2i4,2f16.6,f10.2)') ll+1,ll,ee,hmev(jj),ecut
           exit
         else
           exit
         end if
       end do
     end if

!calculate mean radius of all eigenstates to look for highly localized
!states which indicate postive-energy ghosts
     do jj=1,nbas
       uu(:)=0.0d0
       do ii=1,nbas
         uu(:)=uu(:)+hmat(ii,jj)*uua(:,ii)
       end do

       ro=rr(1)/dsqrt(amesh)
       sn=(ro*uu(1))**2/dfloat(2*ll+4)

       do ii=1,mbar-3
         sn=sn+al*(rr(ii)*uu(ii))**2
       end do

       sn=sn + al*(23.0d0*(rr(mbar-2)*uu(mbar-2))**2 &
&                + 28.0d0*(rr(mbar-1)*uu(mbar-1))**2 &
&                +  9.0d0*(rr(mbar  )*uu(mbar  ))**2)/24.0d0

!      if(sn<rr(irc(l1)) .and. hmev(jj)>0.0d0 .and. npgh<20) then
       if(sn<rr(irc(l1)) .and. hmev(jj)>ee+0.05d0 .and. npgh<20) then
         npgh=npgh+1
         lpgh(npgh)=ll
         epgh(npgh)=hmev(jj)
         rpgh(npgh)=sn/rr(irc(l1))
         cpgh(npgh)=ecut
       end if
     end do


   deallocate(gg0,hmat,hmev,work)
   deallocate(uua)
 end do !l1

   write(6,'(/a)') ' Testing for highly-localized positive-energy ghosts'
   write(6,'(5a)') '    ','   l','    <radius>/rc  ', &
&        '   E Basis Diag.','   E Cutoff'

 if(npgh>0) then
   do jj=1,npgh
           write(6,'(/4x,i4,2f16.6,f10.2,a)') lpgh(jj),rpgh(jj), &
&                epgh(jj),cpgh(jj),'  WARNING - GHOST(+)'
   end do
 end if

 deallocate(uu,up)
 return
 end subroutine run_ghosts
