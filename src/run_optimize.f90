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
 subroutine run_optimize(eig,eig2,ll,mmax,rr,uu,uu2,qq,&
&                        irc,qcut,ncon_in,nbas,npr, &
&                        psopt,vpsp,vkb,vae,cvgplt)

! calls routines to generate optimized pseudo wave function and semi-local
! pseudopotential and prints diagnostic information on process and 
! convergence performance

!eig  energy at which pseudo wave function is computed
!eig2  energy at which 2nd pseudo wave function is computed
!ll  angular momentum
!mmax  dimension of log grid
!rr  log radial grid
!uu  all-electron wave function for first projector
!qq  2x2 matrix of norms and overlaps of uu and uu2
!uu2  all-electron wave function for second projector
!irc  index of core radius
!qcut  q cutoff defining residual energy
!ncon_in  number of constraints for matching psuedo and AE wave functions
!nbas  number of basis functions fo pseudo wave function
!npr   number of projectors = 0,1,2
!psopt  optimized pseudo wave function(s)
!vpsp  corresponding pseudopotential
!vkb  Vanderbilt-Kleinman-Bylander projectors (without local v correction)
!vae  all-electron potential
!cvgplt  energy error vs. cutoff energy for plotting

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV

!Input variables
 integer :: ll,mmax,irc,ncon_in,nbas,npr
 real(dp) :: rr(mmax),sbf(mmax)
 real(dp) :: uu(mmax),uu2(mmax),vae(mmax),qq(2,2)
 real(dp) :: eig,eig2,qcut

!Output variables
 real(dp) :: psopt(mmax,2),vpsp(mmax),vkb(mmax,2)

!Local variables
 real(dp) :: uord(6)
 real(dp) :: al,rc,ulgd,tt
 real(dp) :: sn,ps0norm,amesh,ro
 real(dp) :: err,lerr,qerr,ehaerr,eeverr
 real(dp) :: qroot(11)
 real(dp) :: sbasis(11)
 real(dp) :: cons(6),cvgplt(2,4,2)
 real(dp), allocatable :: orbasis(:,:)
 real(dp), allocatable :: orbasis_der(:,:),pswf0_sb(:),pswfnull_sb(:,:)
 real(dp), allocatable :: pswf0_or(:),pswfnull_or(:,:)
 real(dp), allocatable :: work(:)
 integer :: ii,jj,l1,lmax,ibas,ncon
 logical :: found_root

 integer :: nnull,nqout

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
 dqout=0.5000001d0  !linear mesh spacing for final E_resid(q) interpolation
                    !increment is to avoid exact commnesurability with dq mesh

 psopt(:,:)=0.0d0

 write(6,'(//a,i4)') 'Calculating first optimized projector for l=',ll

 l1=ll+1
 ncon=ncon_in

 nnull=nbas-ncon

 allocate(pswfresid_sb(nbas,nnull))
 allocate(orbasis(nbas,nbas))
 allocate(orbasis_der(6,nbas),pswf0_sb(nbas),pswfnull_sb(nbas,nnull))
 allocate(pswf0_or(nbas),pswfnull_or(nbas,nnull))
 allocate(pswfopt_sb(nbas),pswfopt_or(nbas))
 allocate(work(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))

 rc=rr(irc)

!calculate derivatives of all-electron wave function at r_c
 call wf_rc_der(rr,uu,al,rc,irc,mmax,uord)

 ulgd=uord(2)/uord(1)

 qmax=200.0d0  !very large value should always allow nbas q's to be found

!select wave vectors for spherical Bessel functions
 call qroots(ll,rc,ulgd,nbas,dq,qmax,qroot)

!q considered infinity for E_resid calculation
 qmax=dq*int(2.0d0*qroot(nbas)/dq) 
 qmax=max(qmax,20.0d0)

 nqout=2 + 0.75d0*qmax/dqout
 allocate(qout(nqout),eresq(nqout),leresq(nqout))
 allocate(eresid0(nqout),eresiddot(nnull,nqout),eresidmat(nnull,nnull,nqout))

 qout(1)=qcut
 do ii=nqout,2,-1
  qout(ii)=(ii-2)*dqout
 end do

!write(6,'(/a)') 'qroots'
!write(6,'(6f12.6)') (qroot(ii),ii=1,nbas)

!calculate orthogonal basis
 call sbf_basis(ll,rr,mmax,irc,nbas,qroot,sbasis,orbasis,orbasis_der)

!load constraint vector
 cons(:)=0.0d0
 do jj=1,ncon
  cons(jj)=uord(jj)
 end do

!calculate constrained basis for residual energy minimization
 call const_basis(nbas,ncon,cons,orbasis,orbasis_der, &
&                 pswf0_or,pswfnull_or, &
&                 pswf0_sb,pswfnull_sb,ps0norm)

!calculate eigenvectors, eigenvalues, and inhomogeneous terms for
!the residual energy for a set of q lower cutoffs

 write(6,'(/a,f10.6)') '    Fraction of norm inside rc',qq(1,1)

 call eresid(ll,irc,nnull,nbas,mmax,rr,dr,dq,qmax,qroot, &
&                  uu,pswf0_sb,pswfnull_sb, nqout,qout, &
&                  eresid0,eresiddot,eresidmat)

 write(6,'(a,f5.2,a,f5.2,a)') '    Optimizing pswf for qcut=',qout(1), &
& ' a_B^-1,  ecut=',0.5d0*qout(1)**2,' Ha'
 write(6,'(a,f6.2,a,f7.1,a)') '    q_infinity defining residual KE=',&
& qmax,'  (E_inf=',0.5d0*qmax**2,' Ha)'

!calculate eigenvectors, eigenvalues, and inhomogeneous terms for
!the residual energy for a set of q lower cutoffs

 call optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
&              eresid0,eresiddot,eresidmat,&
&              pswfnull_sb,pswfnull_or,qq(1,1),ps0norm,eresidmin,pswfopt_sb, &
&              pswfopt_or,ekin_anal,eresq)

 write(6,'(a,1p,d10.2,a)') '    Residual kinetic energy error=', &
&      eresidmin,' Ha'

!find the semi-local pseudopotential / Kleinman-Bylander projector

 call pspot(1,ll,rr,irc,mmax,al,nbas,qroot,eig,uu,pswfopt_sb, &
&           psopt(1,1),vae,vpsp,vkb(1,1),ekin_num)

 write(6,'(/a)') '    Total kinetic energy consistency test'
 write(6,'(a)') '      Fourier integration compared to d^2/dr^2 integral'
 write(6,'(a,f12.8,a,f12.8,a,f12.8)') '      Fourier',ekin_anal, &
&  '  r-space=',ekin_num,'  ratio=',ekin_anal/ekin_num
 write(6,'(a)') '    Potential consistency test at r_c'
 write(6,'(a,f12.8,a,f12.8,a,1p,d10.2)') '      vpsp=',vpsp(irc),'  vae=', &
&  vae(irc),'  difference=',vae(irc)-vpsp(irc)
 do ii=irc+1,mmax
   vpsp(ii)=vae(ii)
 end do

!interpolate the convergence behavior of the optimized pseudo wave function

 write(6,'(/a)') '    Energy error per electron        Cutoff'
 write(6,'(a)') '         Ha          eV             Ha'
 leresq(:)=log(eresq(:))
 err=0.01d0
 do jj=1,4
  lerr=log(err)
  do ii=3,nqout
   if(leresq(ii-1)>lerr .and. leresq(ii)<=lerr) then
    qerr=qout(ii)-dqout*(lerr-leresq(ii))/(leresq(ii-1)-leresq(ii))
    cvgplt(1,jj,1)= 0.5d0*qerr**2
    cvgplt(2,jj,1)= err
    write(6,'(4x,2f12.5,f12.2)') err,err*Ha_eV,0.5d0*qerr**2
   end if
  end do
  err=0.1d0*err
 end do

 if(npr==0 .or. npr==1) then
  deallocate(eresid0,eresiddot,eresidmat)
  deallocate(qout)
  deallocate(orbasis)
  deallocate(orbasis_der,pswf0_sb,pswfnull_sb)
  deallocate(pswf0_or,pswfnull_or)
  deallocate(pswfopt_sb,pswfopt_or)
  return
 end if

! begin section for optimized second projector wavefunction

 write(6,'(/a,i4)') 'Calculating second optimized projector for l=',ll

 call wf_rc_der(rr,uu2,al,rc,irc,mmax,uord)

! load constraint vector for value/derivative matching
 cons(:)=0.0d0
 do jj=1,ncon
  cons(jj)=uord(jj)
 end do

! last constraint is off-diagonal overlap
 ncon=ncon+1
 nnull=nbas-ncon
 cons(ncon)=qq(1,2)

! last row of or-basis derivative matrix is loaded with optimized
! first-projector coefficients in or basis
 orbasis_der(ncon,:)=pswfopt_or(:)

! resize arrays
 deallocate(pswfresid_sb)
 deallocate(pswfnull_sb)
 deallocate(pswfnull_or)
 deallocate(eresiddot,eresidmat)

 allocate(pswfresid_sb(nbas,nnull))
 allocate(pswfnull_sb(nbas,nnull))
 allocate(pswfnull_or(nbas,nnull))
 allocate(eresiddot(nnull,nqout),eresidmat(nnull,nnull,nqout))

!calculate constrained basis for residual energy minimization
!satisfying off-diagonal norm conservation
 call const_basis(nbas,ncon,cons,orbasis,orbasis_der, &
&                 pswf0_or,pswfnull_or, &
&                 pswf0_sb,pswfnull_sb,ps0norm)

!calculate eigenvectors, eigenvalues, and inhomogeneous terms for
!the residual energy for a set of q lower cutoffs

 write(6,'(/a,f10.6)') '    Fraction of norm inside rc',qq(2,2)

 call eresid(ll,irc,nnull,nbas,mmax,rr,dr,dq,qmax,qroot, &
&                  uu2,pswf0_sb,pswfnull_sb, nqout,qout, &
&                  eresid0,eresiddot,eresidmat)

 write(6,'(a,f5.2,a,f5.2,a)') '    Optimizing pswf for qcut=',qout(1), &
& ' a_B^-1,  ecut=',0.5d0*qout(1)**2,' Ha'
 write(6,'(a,f6.2,a,f7.1,a)') '    q_infinity defining residual KE=',&
& qmax,'  (E_inf=',0.5d0*qmax**2,' Ha)'

!find the null-basis coefficients which minimize the eresid while
!satisfying diagonal norm conservation

 call optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
&              eresid0,eresiddot,eresidmat,&
&              pswfnull_sb,pswfnull_or,qq(2,2),ps0norm,eresidmin,pswfopt_sb, &
&              pswfopt_or,ekin_anal,eresq)

 write(6,'(a,1p,d10.2,a)') '    Residual kinetic energy error=', &
&      eresidmin,' Ha'

! find the Vanderbilt second projector

 call pspot(2,ll,rr,irc,mmax,al,nbas,qroot,eig2,uu2,pswfopt_sb, &
&           psopt(1,2),vae,work,vkb(1,2),ekin_num)

 write(6,'(/a)') '    Total kinetic energy consistency test'
 write(6,'(a)') '      Fourier integration compared to d^2/dr^2 integral'
 write(6,'(a,f12.8,a,f12.8,a,f12.8)') '      Fourier',ekin_anal, &
&  '  r-space=',ekin_num,'  ratio=',ekin_anal/ekin_num
 write(6,'(a)') '    Potential consistency test at r_c'
 write(6,'(a,f12.8,a,f12.8,a,1p,d10.2)') '    "vpsp"=',work(irc),'  vae=', &
&  vae(irc),'  difference=',vae(irc)-work(irc)
 do ii=irc+1,mmax
   vpsp(ii)=vae(ii)
 end do

!interpolate the convergence behavior of the optimized pseudo wave function

 write(6,'(/a)') '    Energy error per electron        Cutoff'
 write(6,'(a)') '         Ha          eV             Ha'
 leresq(:)=log(eresq(:))
 err=0.01d0
 do jj=1,4
  lerr=log(err)
  do ii=3,nqout
   if(leresq(ii-1)>lerr .and. leresq(ii)<=lerr) then
    qerr=qout(ii)-dqout*(lerr-leresq(ii))/(leresq(ii-1)-leresq(ii))
    cvgplt(1,jj,2)= 0.5d0*qerr**2
    cvgplt(2,jj,2)= err
    write(6,'(4x,2f12.5,f12.2)') err,err*Ha_eV,0.5d0*qerr**2
   end if
  end do
  err=0.1d0*err
 end do

 deallocate(eresid0,eresiddot,eresidmat)
 deallocate(qout)
 deallocate(orbasis)
 deallocate(orbasis_der,pswf0_sb,pswfnull_sb)
 deallocate(pswfopt_sb,pswfopt_or)
 deallocate(work)

 return
 end subroutine run_optimize
