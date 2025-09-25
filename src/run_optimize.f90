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
subroutine run_optimize(eig,ll,mmax,mxprj,rr,uua,qq,&
&                        irc,qcut,qmsbf,ncon_in,nbas_in,npr, &
&                        psopt,vpsp,vkb,vae,cvgplt)

! calls routines to generate optimized pseudo wave function and semi-local
! pseudopotential and prints diagnostic information on process and
! convergence performance

!eig  energy at which pseudo wave function is computed
!ll  angular momentum
!mmax  dimension of log grid
!mxprj  dimension of number of projectors
!rr  log radial grid
!uu  all-electron wave function for first projector
!qq  2x2 matrix of norms and overlaps of uu and uu2
!uu2  all-electron wave function for second projector
!irc  index of core radius
!qcut  q cutoff defining residual energy
!qmsbf  maximum q in sbfs for this ll
!ncon_in  number of constraints for matching psuedo and AE wave functions
!nbas_in  number of basis functions for pseudo wave function
!np_in   number of projectors = 0,1,2
!psopt  optimized pseudo wave function(s)
!vpsp  corresponding pseudopotential
!vkb  Vanderbilt-Kleinman-Bylander projectors (without local v correction)
!vae  all-electron potential
!cvgplt  energy error vs. cutoff energy for plotting

   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: Ha_eV=27.21138386d0  ! 1 Hartree, in eV

!Input variables
   integer :: ll,mmax,mxprj,irc,ncon_in,nbas_in,npr
   real(dp) :: rr(mmax)
   real(dp) :: uua(mmax,mxprj),vae(mmax),qq(mxprj,mxprj)
   real(dp) :: eig(mxprj),qcut

!Output variables
   real(dp) :: qmsbf
   real(dp) :: psopt(mmax,mxprj),vpsp(mmax),vkb(mmax,mxprj)

!Local variables
   real(dp) :: uord(6)
   real(dp) :: al,rc,ulgd
   real(dp) :: ps0norm
   real(dp) :: err,lerr,qerr
   real(dp) :: cons(6),cvgplt(2,7,mxprj)
   real(dp), allocatable :: orbasis(:,:)
   real(dp), allocatable :: orbasis_der(:,:),pswf0_sb(:),pswfnull_sb(:,:)
   real(dp), allocatable :: pswf0_or(:),pswfnull_or(:,:)
   real(dp), allocatable :: work(:)
   integer :: ii,iprj,jj,l1,nbas,ncon,nconmx

   integer :: nnull,nqout

   real(dp), allocatable :: qroot(:),sbasis(:)
   real(dp), allocatable :: eresid0(:),eresiddot(:,:)
   real(dp), allocatable :: eresidmat(:,:,:),pswfresid_sb(:,:)
   real(dp), allocatable :: qout(:),eresq(:),leresq(:)

   real(dp), allocatable :: pswfopt_sb(:),pswfopt_or(:)
   real(dp) :: eresidmin
   real(dp) :: ekin_anal,ekin_num,qmax

!parameters whose default values here give well-converged results
! increasing dq and/or dr will speed up the code at some cost
! in optimization accuracy
   real(dp) :: dq,dr,dqout
   dq=0.02d0  !linear integration mesh spacing for above
   dr=0.02d0  !linear integration spacing for Fourier transforms
!dr=0.002d0  !linear integration spacing for Fourier transforms
   dqout=0.5000001d0  !linear mesh spacing for final E_resid(q) interpolation
   !increment is to avoid exact commnesurability with dq mesh

   al = 0.01d0 * dlog(rr(101) / rr(1))
   rc=rr(irc)
   l1=ll+1

   nconmx=ncon_in+npr-1

   allocate(work(mmax))

!maximum basis size
   nbas=nbas_in+npr-1
   allocate(qroot(nbas))

!calculate derivatives of all-electron wave function at r_c
!for the first projector to set sbf wave vectors for all projectors

   call wf_rc_der(rr,uua(1,1),al,rc,irc,mmax,uord)
   ulgd=uord(2)/uord(1)

   qmax=200.0d0  !very large value should always allow nbas q's to be found

!select wave vectors for spherical Bessel functions
   call qroots(ll,rc,ulgd,nbas,dq,qmax,qroot)

   qmsbf=qroot(nbas)

!write(6,'(/a)') 'qroots'
!write(6,'(6f12.6)') (qroot(ii),ii=1,nbas)

   psopt(:,:)=0.0d0

! loop over projectors
   do iprj=1,npr

      write(6,'(/a,i4)') 'Calculating optimized projector #',iprj, &
      &        ' for l=',ll

!basis size and number of constraints for this projector
      nbas=nbas_in+iprj-1
      ncon=ncon_in+iprj-1
      nnull=nbas-ncon

      allocate(sbasis(nbas))
      allocate(orbasis(nbas,nbas))
      allocate(orbasis_der(ncon,nbas))
      allocate(pswf0_sb(nbas),pswf0_or(nbas))
      allocate(pswfopt_sb(nbas),pswfopt_or(nbas))
      allocate(pswfresid_sb(nbas,nnull))
      allocate(pswfnull_sb(nbas,nnull))
      allocate(pswfnull_or(nbas,nnull))

!calculate derivatives of all-electron wave function at r_c
!to define basis set for current projector and derivative constraints

      call wf_rc_der(rr,uua(1,iprj),al,rc,irc,mmax,uord)

!q considered infinity for E_resid calculation
      qmax=dq*int(2.0d0*qroot(nbas)/dq)
      qmax=max(qmax,20.0d0)


      nqout=2 + int(0.75d0*qmax/dqout)

      allocate(qout(nqout),eresq(nqout),leresq(nqout),eresid0(nqout))
      allocate(eresiddot(nnull,nqout),eresidmat(nnull,nnull,nqout))

      qout(1)=qcut
      do ii=nqout,2,-1
         qout(ii)=(ii-2)*dqout
      end do

!calculate orthogonal basis and constraint matrix
      call sbf_basis_con(ll,rr,mmax,irc,nbas,qroot,psopt,orbasis,orbasis_der, &
      &                     iprj,mxprj,ncon,ncon_in)


!load constraint vector for value/derivative matching
      cons(:)=0.0d0
      do jj=1,ncon_in
         cons(jj)=uord(jj)
      end do
!load overlap constraints if present
      if(iprj>=2) then
         do jj=1,iprj-1
            cons(ncon_in+jj)=qq(jj,iprj)
         end do
      end if

!calculate constrained basis for residual energy minimization
!satisfying off-diagonal norm conservation
      call const_basis(nbas,ncon,cons,orbasis,orbasis_der, &
      &                   pswf0_or,pswfnull_or, &
      &                   pswf0_sb,pswfnull_sb,ps0norm)

!calculate eigenvectors, eigenvalues, and inhomogeneous terms for
!the residual energy for a set of q lower cutoffs

      write(6,'(/a,f10.6)') '    Fraction of norm inside rc',qq(iprj,iprj)

      call eresid(ll,irc,nnull,nbas,mmax,rr,qmax,qroot, &
      &                    uua(1,iprj),pswf0_sb,pswfnull_sb, nqout,qout, &
      &                    eresid0,eresiddot,eresidmat)

      write(6,'(a,f7.2,a,f7.2,a)') '    Optimizing pswf for qcut=',qout(1), &
      &   ' a_B^-1,  ecut=',0.5d0*qout(1)**2,' Ha'
      write(6,'(a,f6.2,a,f7.1,a)') '    q_infinity defining residual KE=',&
      &   qmax,'  (E_inf=',0.5d0*qmax**2,' Ha)'

!find the null-basis coefficients which minimize the eresid while
!satisfying diagonal norm conservation

      call optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
      &                eresid0,eresiddot,eresidmat,&
      &                pswfnull_sb,pswfnull_or,qq(iprj,iprj),ps0norm,eresidmin, &
      &                pswfopt_sb,pswfopt_or,ekin_anal,eresq)

      write(6,'(a,1p,d10.2,a)') '    Residual kinetic energy error=', &
      &        eresidmin,' Ha'

! find the Vanderbilt projectors and optimized wave functions

      call pspot(iprj,ll,rr,irc,mmax,al,nbas,qroot,eig(iprj),uua(1,iprj), &
      &             pswfopt_sb,psopt(1,iprj),vae,work,vkb(1,iprj),ekin_num)

!semi-local potential
      if(iprj==1) then
         do ii=1,irc
            vpsp(ii)=work(ii)
         end do
         do ii=irc+1,mmax
            vpsp(ii)=vae(ii)
         end do
      end if

      write(6,'(/a)') '    Total kinetic energy consistency test'
      write(6,'(a)') '      Fourier integration compared to d^2/dr^2 integral'
      write(6,'(a,f12.8,a,f12.8,a,f12.8)') '      Fourier',ekin_anal, &
      &    '  r-space=',ekin_num,'  ratio=',ekin_anal/ekin_num
      write(6,'(a)') '    Potential consistency test at r_c'
      write(6,'(a,f12.8,a,f12.8,a,1p,d10.2)') '    "vpsp"=',work(irc),'  vae=', &
      &    vae(irc),'  difference=',vae(irc)-work(irc)

!interpolate the convergence behavior of the optimized pseudo wave function

      write(6,'(/a)') '    Energy error per electron        Cutoff'
      write(6,'(a)') '         Ha          eV             Ha'
      leresq(:)=log(eresq(:))
      err=0.01d0
      do jj=1,7
         lerr=log(err)
         do ii=3,nqout
            if(leresq(ii-1)>lerr .and. leresq(ii)<=lerr) then
               qerr=qout(ii)-dqout*(lerr-leresq(ii))/(leresq(ii-1)-leresq(ii))
               cvgplt(1,jj,iprj)= 0.5d0*qerr**2
               cvgplt(2,jj,iprj)= err
               if(mod(jj,2)/=0) then
                  write(6,'(4x,2f12.5,f12.2)') err,err*Ha_eV,0.5d0*qerr**2
               end if
            end if
         end do
         err=sqrt(0.1d0)*err
      end do

      deallocate(sbasis)
      deallocate(orbasis)
      deallocate(orbasis_der)
      deallocate(pswf0_sb,pswf0_or)
      deallocate(pswfopt_sb,pswfopt_or)
      deallocate(pswfresid_sb)
      deallocate(pswfnull_sb)
      deallocate(pswfnull_or)
      deallocate(eresiddot,eresidmat)
      deallocate(qout,eresq,leresq,eresid0)

   end do  !iprj

   deallocate(work,qroot)

   return
end subroutine run_optimize
