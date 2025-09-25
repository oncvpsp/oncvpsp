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
! calculates convergence-optimized pseudo-wave-function coefficients
! in orthonormal and spherical Bessel function bases

subroutine optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
&                    eresid0,eresiddot,eresidmat,&
&                    pswfnull_sb,pswfnull_or,uunorm,ps0norm,eresidmin,&
&                    pswfopt_sb,pswfopt_or,ekin_anal,eresq)

!nnull  number of unconstrained basis vectors for residual minimization
!nbas  number of sbf basis functions
!pswf0_sb(nbas)  sbf basis  coefficients for constrained part of pseudo
!                    wave function
!pswf0_or(nbas)  or basis  coefficients
!nqout  number of q values in [0,qmax] for which results are to be saved
!qout(nqout)  q values in [0,qmax] for which residuals are to be calculated
!eresid0(nqout) set of <pswf0| E_resid |pswf0> matrix elements
!eresiddot(nnull,nqout) set of <pswfnull| E_resid |pswf0> matrix elements
!eresidmat(nnull,nnull,nqout) set of <pswfnull | E_resid |pswfnull'>
!pswfnull_sb(nbas,nnull)  sbf coefficients of null-space eigenfunctions
!pswfnull_sb(nbas,nnull)  or basis coefficients of null-space eigenfunctions
!uunorm  all-electron charge inside rc
!ps0norm  ps0 charge inside rc
!pswfopt_sb optimized pseudowavefunction sbf coefficients
!orbasis_it  transformation matrix from sb to or basis
!pswfopt_or optimized pseudowavefunction or basis coefficients
!ekin_anal  total pswf kinetic energy calculated analytically from E_resid
!eresq(nqout)  E_resid for optimized pswf as a function of cutoff

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: nnull,nbas,nqout
   real(dp) :: pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
   real(dp) :: pswf0_or(nbas),pswfnull_or(nbas,nnull)
   real(dp) :: eresid0(nqout),eresiddot(nnull,nqout),eresidmat(nnull,nnull,nqout)
   real(dp) :: qout(nqout),eresq(nqout)
   real(dp) :: uunorm,ps0norm

!Output variables
   real(dp) :: eresidmin,ekin_anal
   real(dp) :: pswfopt_sb(nbas),pswfopt_or(nbas)

!Local variables
   real(dp) :: yy,tt,emin,rtsgn
   real(dp) :: x1,x1min,x1max
   real(dp), allocatable :: work(:),wmat(:,:),wev(:),wvec(:)
   real(dp), allocatable :: eresidevec_nu(:,:)
   real(dp), allocatable :: pswfopt(:),pswfopt_nu(:)
   real(dp), allocatable :: eresiddot_ev(:),eresideval(:)
   real(dp), parameter :: eps=1.d-10
   integer :: ii,iter,jj,kk,info
   integer, parameter :: niter=100
   logical :: converged

!find eigenvalues and eigenvectors of saved E_resid matrix for qout(1)
!which is meant to be the cutoff used for optimization

   allocate(work(5*nnull),wmat(nnull,nnull),wev(nnull),wvec(nnull))
   allocate(eresidevec_nu(nnull,nnull),eresideval(nnull))
   allocate(eresiddot_ev(nnull))

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   wmat(:,:)=eresidmat(:,:,1)

   call dsyev( 'V', 'U', nnull, wmat, nnull, wev, work, 5*nnull, info )
   if(info /= 0) then
      write(6,'(a,i4)') 'optimize: residual-energy eigenvalue ERROR, info=',info
      stop
   end if

!convert <ps0| E_resid |psnull> vector from the null-space representation
!to the eigenvector representation

   wvec(:)=eresiddot(:,1)
   eresiddot_ev(:)=0.0d0
   do jj=1,nnull
      do kk=1,nnull
         eresiddot_ev(jj)=eresiddot_ev(jj)+wvec(kk)*wmat(kk,jj)
      end do
   end do

!write(6,'(/a,1p,4d14.4)') 'eresideval',wev(:)
!write(6,'(/a,1p,4d14.4)') 'eresiddot_ev',eresiddot_ev(:)
!write(6,'(/a,1p,4d14.4)') 'eresid0',eresid0(1)

   eresideval(:)=wev(:)
   eresidevec_nu(:,:)=wmat(:,:)

   deallocate(work,wmat,wev,wvec)

!find coefficients of E_resid eigenvectors which must be added to pswf0
!to minimize E_resid

!note that E_resid is an exteremely simple function in this representation,
!consisting of a constant plus as sum of linear and quadratic terms in each of
!these coefficients

   allocate(pswfopt(nnull),pswfopt_nu(nnull))
   pswfopt(:)=0.d0
   pswfopt_nu(:)=0.d0

!write(6,'(/a,1p,3d16.6)') 'initial emin, final, diff',emin0,emin,emin0-emin

! Interval-halving search to minimize E*r

! This is all based on Eqs. (14) - (16) in my paper as correccted in the
! erratum.  It is obvious from Eq.(14) that to minimize E^r, the signs of x_i
! must be the opposite of those of f_i.  Eq.(16) is derived by treating
! x_1 as depenent upon x_2...x_N-M through Eq.(15).  We can then turn
! around and treat those x's as functions all dependent on x_1.  Since the
! signs of x_1 and f_1 must be opposite at the minimum, and since e_1
! is the smallest eigenvalue, the denominator of Eq.(16) is always positive,
! and the magnitudes of the x_i are all monotonically increasing functions
! of x_1.  The limits of |x_1| are zero and D_norm (yy here), so a simple
! interval-halving search varying |x_1| to satisfy the norm constraint
! (Eq.(15)) should be robust.

! Translations between the notation in the paper and here are as follows:
!   pswfopt(ii) = x_i
!   eresiddot_ev(ii) = f_i
!   eresideval(ii) = e_i
!   D_norm == yy
!   E^r_00 = eresid0(1)
!   N-M = nnull


   rtsgn=-sign(1.0d0, eresiddot_ev(1))
   if(ps0norm>uunorm) then
      write(6,'(/a)') 'optimize: ERROR ps0norm > uunorm, code will stop'
      stop
   end if
   yy=sqrt((uunorm-ps0norm))

   x1min=0.0d0
   x1max=yy
   converged=.false.

   do iter=1,niter

      x1=0.5d0*(x1max+x1min)

      tt=x1**2-yy**2
      do ii=2,nnull
         pswfopt(ii)=-eresiddot_ev(ii)/(eresideval(ii)-eresideval(1) &
         &               +abs(eresiddot_ev(1))/x1)
         tt=tt+pswfopt(ii)**2
      end do

      if(x1max-x1min<eps) then
         pswfopt(1)=rtsgn*x1
         converged=.true.
         exit
      end if

      if(tt<0.0d0) then
         x1min=x1
      else
         x1max=x1
      end if
   end do

   if(converged .eqv. .false.) then
      write(6,'(/a)') 'optimize: ERROR interval-halving search not converged'
      stop
   end if

   emin=eresid0(1)
   do ii=1,nnull
      emin=emin+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
      &            +eresideval(ii)*pswfopt(ii)**2)
   end do
   eresidmin=emin

   ekin_anal=0.0d0
!find cutoff-dependence of E_resid for optimized pswf
!must transform pswfopt to the null basis representation to be
!consistent with the input dot product and
   allocate(work(nnull))
   do ii=1,nnull
      pswfopt_nu(ii)=0.0d0
      do jj=1,nnull
         pswfopt_nu(ii)=pswfopt_nu(ii)+eresidevec_nu(ii,jj)*pswfopt(jj)
      end do
   end do

   do kk=1,nqout
      emin=eresid0(kk)
      work(:)=0.0d0
      do ii=1,nnull
         emin=emin+2.0d0*eresiddot(ii,kk)*pswfopt_nu(ii)
         work(:)=work(:)+eresidmat(:,ii,kk)*pswfopt_nu(ii)
      end do
      do ii=1,nnull
         emin=emin+pswfopt_nu(ii)*work(ii)
      end do
      eresq(kk)=emin
      if(qout(kk)==0.0d0) ekin_anal=emin
   end do

!calculate optimized pswf in sbf representation
   pswfopt_sb(:)=pswf0_sb(:)
   pswfopt_or(:)=pswf0_or(:)
   do ii=1,nnull
      pswfopt_sb(:)=pswfopt_sb(:)+pswfopt_nu(ii)*pswfnull_sb(:,ii)
      pswfopt_or(:)=pswfopt_or(:)+pswfopt_nu(ii)*pswfnull_or(:,ii)
   end do

   deallocate(pswfopt,pswfopt_nu)
   return
end subroutine optimize
