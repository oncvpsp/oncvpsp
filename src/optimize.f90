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
 real(dp) :: orbasis_it(nbas,nbas)
 real(dp) :: qout(nqout),eresq(nqout)
 real(dp) :: uunorm,ps0norm

!Output variables
 real(dp) :: eresidmin,ekin_anal
 real(dp) :: pswfopt_sb(nbas),pswfopt_or(nbas)
 real(dp) :: eres(nqout)

!Local variables
 real(dp) :: dpso,yy,tt,emin,rtsgn,sn
 real(dp) :: eropt(2)
 real(dp), allocatable :: work(:),wmat(:,:),wev(:),wvec(:)
 real(dp), allocatable :: eresidevec_nu(:,:)
 real(dp), allocatable :: pswfopt(:),pswfomin(:),pswfopt_nu(:)
 real(dp), allocatable :: eresiddot_ev(:),eresideval(:)
 real(dp), parameter :: eps=1.d-6
 integer :: ii,iter,jj,kk,ll,nn,info
 integer, parameter :: nopt=20,niter=5


!find eigenvalues and eigenvectors of saved E_resid matrix for qout(1)
!which is meant to be the cutoff used for optimization

 allocate(work(5*nnull),wmat(nnull,nnull),wev(nnull),wvec(nnull))
 allocate(eresidevec_nu(nnull,nnull),eresideval(nnull))
 allocate(eresiddot_ev(nnull))

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

 wmat(:,:)=eresidmat(:,:,1)

 call dsyev( 'V', 'U', nnull, wmat, nnull, wev, work, 5*nnull, info )
 if(info .ne. 0) then
  write(6,'(a,i4)') 'optimize: redisual-energy eigenvalue error, info=',info
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

 allocate(pswfopt(nnull),pswfomin(nnull),pswfopt_nu(nnull))

!Search a coarse grid of values allowed by the norm constraint for
!coefficients 2,...,nnull, calculating coefficient 1 from the norm
!constraint.  1 corresponds to the smallest eigenvalue, so it is
!appropriate to make it "dependent."
 emin=1.0d5
 yy=sqrt((uunorm-ps0norm))
 dpso=yy/nopt
 if(nnull==2) then
  do kk=-nopt,nopt
   pswfopt(2)=dpso*kk
   pswfopt(1)=sqrt(yy**2-pswfopt(2)**2)
   eropt(1)=eresid0(1)
   do ii=1,nnull
    eropt(1)=eropt(1)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                +eresideval(ii)*pswfopt(ii)**2)
   end do
   if(eropt(1)<emin) then
    rtsgn=1.0d0
    emin=eropt(1)
    pswfomin(:)=pswfopt(:)
   end if
   pswfopt(1)=-sqrt(yy**2-pswfopt(2)**2)
   eropt(2)=eresid0(1)
   do ii=1,nnull
    eropt(2)=eropt(2)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                +eresideval(ii)*pswfopt(ii)**2)
   end do
   if(eropt(2)<emin) then
    rtsgn=-1.0d0
    emin=eropt(2)
    pswfomin(:)=pswfopt(:)
   end if
  end do
  pswfopt(:)=pswfomin(:)
 else if(nnull==3) then
  do kk=-nopt,nopt
   pswfopt(2)=dpso*kk
   do jj=-nopt,nopt
    pswfopt(3)=dpso*jj
    tt=yy**2-pswfopt(2)**2-pswfopt(3)**2
    if(tt>0.0d0) then
     pswfopt(1)=sqrt(tt)
     eropt(1)=eresid0(1)
     do ii=1,nnull
      eropt(1)=eropt(1)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                  +eresideval(ii)*pswfopt(ii)**2)
     end do
     if(eropt(1)<emin) then
      rtsgn=1.0d0
      emin=eropt(1)
      pswfomin(:)=pswfopt(:)
     end if
     pswfopt(1)=-sqrt(tt)
     eropt(2)=eresid0(1)
      do ii=1,nnull
      eropt(2)=eropt(2)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                  +eresideval(ii)*pswfopt(ii)**2)
      end do
      if(eropt(2)<emin) then
       rtsgn=-1.0d0
       emin=eropt(2)
       pswfomin(:)=pswfopt(:)
      end if
    end if
   end do
  end do
  pswfopt(:)=pswfomin(:)
 else if(nnull==4) then
  do kk=-nopt,nopt
   pswfopt(2)=dpso*kk
   do jj=-nopt,nopt
    pswfopt(3)=dpso*jj
    do nn=-nopt,nopt
     pswfopt(4)=dpso*nn
     tt=yy**2-pswfopt(2)**2-pswfopt(3)**2-pswfopt(4)**2
     if(tt>0.0d0) then
      pswfopt(1)=sqrt(tt)
      eropt(1)=eresid0(1)
      do ii=1,nnull
       eropt(1)=eropt(1)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                   +eresideval(ii)*pswfopt(ii)**2)
      end do
      if(eropt(1)<emin) then
       rtsgn=1.0d0
       emin=eropt(1)
       pswfomin(:)=pswfopt(:)
      end if
      pswfopt(1)=-sqrt(tt)
      eropt(2)=eresid0(1)
       do ii=1,nnull
       eropt(2)=eropt(2)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                   +eresideval(ii)*pswfopt(ii)**2)
       end do
       if(eropt(2)<emin) then
        rtsgn=-1.0d0
        emin=eropt(2)
        pswfomin(:)=pswfopt(:)
       end if
     end if
    end do
   end do
  end do
  pswfopt(:)=pswfomin(:)
 else if(nnull==5) then
  do kk=-nopt,nopt
   pswfopt(2)=dpso*kk
   do jj=-nopt,nopt
    pswfopt(3)=dpso*jj
    do ll=-nopt,nopt
     pswfopt(4)=dpso*ll
     do nn=-nopt,nopt
      pswfopt(5)=dpso*nn
      tt=yy**2-pswfopt(2)**2-pswfopt(3)**2-pswfopt(4)**2-pswfopt(5)**2
      if(tt>0.0d0) then
       pswfopt(1)=sqrt(tt)
       eropt(1)=eresid0(1)
       do ii=1,nnull
        eropt(1)=eropt(1)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                    +eresideval(ii)*pswfopt(ii)**2)
       end do
       if(eropt(1)<emin) then
        rtsgn=1.0d0
        emin=eropt(1)
        pswfomin(:)=pswfopt(:)
       end if
       pswfopt(1)=-sqrt(tt)
       eropt(2)=eresid0(1)
        do ii=1,nnull
        eropt(2)=eropt(2)+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&                    +eresideval(ii)*pswfopt(ii)**2)
        end do
        if(eropt(2)<emin) then
         rtsgn=-1.0d0
         emin=eropt(2)
         pswfomin(:)=pswfopt(:)
        end if
      end if
     end do
    end do
   end do
  end do
  pswfopt(:)=pswfomin(:)
 else
  write(6,'(/a,a)') 'optimize: only set up  for E_resid optimization on', &
&  ' <=5  variables'
  write(6,'(a,i3)') 'nnull=',nnull
  stop
 end if

!Iterative refinement exploiting the typical small value of eigenvalue 1.
!Treating coefficient 1 as a constant, dE_resid/d(coeff.) = 0 for 2-nnull
!is solved analytically.  Coeff. 1 is then updated from the norm condition.
!This converges very rapidly.  It is necessary that the proper sign for
!the square root is found first from the coarse search.

!write(6,'(/a,d18.8)') 'initial sweep result emin=',emin

 do iter=1,niter
  tt=yy**2
  do ii=2,nnull
   pswfopt(ii)=-eresiddot_ev(ii)/(eresideval(ii)+eresideval(1) &
&             -0.5d0*rtsgn*eresiddot_ev(1)/pswfopt(1))
   tt=tt-pswfopt(ii)**2 
   if(tt<0.0d0) then 
    write(6,'(/a)') 'optimize: norm constraint violated (line 260)'
    stop
   end if
   pswfopt(1)=rtsgn*sqrt(tt)
  end do

  emin=eresid0(1)
  do ii=1,nnull
   emin=emin+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&             +eresideval(ii)*pswfopt(ii)**2)
  end do
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

 deallocate(pswfopt,pswfomin,pswfopt_nu)
 return
 end subroutine optimize

