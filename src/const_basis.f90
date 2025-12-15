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
!calculates basis for residual energy minimization

 subroutine const_basis(nbas,ncon,cons,orbasis,orbasis_der,&
&                       pswf0_or,pswfnull_or, &
&                       pswf0_sb,pswfnull_sb,ps0norm)

!ncon  number of constraints
!nbas  number of orthogonalized sbf basis functions
!cons  constraints (value at rc and derivatives up to 4th)
!orbasis  coefficients of sbfs in orthogonal basis
!orbasis_der  values and derivatives of orthogonalized sbf basis 
!              functions at rc (constraint matrix)
!pswf0_sb  linear combination of sbfs that matches cons at rc
!pswfnull_sb  nbas-ncon linear combinations of sbfs satisfying value and ncon-1
!           derivatives at rc == 0 (null space of constraint matrix)
!          These form an orthonormal set, and are all orthogonal to
!          the pswf0 combination (which is not normalized)
!ps0norm  charge contained in ps0 inside rc

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!input variables
 integer :: ncon,nbas
 real(dp) :: orbasis(nbas,nbas),orbasis_der(6,nbas),cons(6)

!Output variables
 real(dp) :: pswf0_or(nbas),pswfnull_or(nbas,nbas-ncon)
 real(dp) :: pswf0_sb(nbas),pswfnull_sb(nbas,nbas-ncon)
 real(dp) :: ps0norm

!Local variables
 real(dp) :: usvd(6,6),ss(6),work(36),con_tst(6)
 real(dp), allocatable :: amat(:,:),vt(:,:)
 real(dp) :: sn
 real(dp), parameter :: eps=1.d-8
 integer :: ii,jj,kk,info


!DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!     $                   WORK, LWORK, INFO )
 allocate(amat(6,nbas),vt(nbas,nbas))

!save orbasis_der from overwrite
 amat(:,:)=orbasis_der(:,:) 

!singular value decomposition of constraint matrix
 call dgesvd('A','A',ncon,nbas,amat,6,ss,usvd,6,vt,nbas, &
&             work,36,info)

 if(info /= 0) then
  write(6,'(/a,i4)') 'const_basis: dgesvd info = ',info
  stop
 end if

!write(6,'(/a,6d12.4)') 'singular values',(ss(ii),ii=1,ncon)

!write(6,'(/a)') 'usvd'
!do jj=1,ncon
! write(6,'(6f12.6)') (usvd(jj,ii),ii=1,ncon)
!end do

!write(6,'(/a)') 'vt'
!do jj=1,nbas
! write(6,'(11f12.6)') (vt(jj,ii),ii=1,nbas)
!end do

!operate on constraint vector with usvd^transpose
 work(:)=0.0d0
 do ii=1,ncon
  work(ii)=0.0d0
  do jj=1,ncon
   work(ii)=work(ii)+usvd(jj,ii)*cons(jj)
  end do
 end do

!scale work with inverse transpost of singular value "matrix" ss
 do ii=1,ncon
  if(abs(ss(ii))>eps) then
   work(ii)=work(ii)/ss(ii)
  end if
 end do

!form coefficient vector for basic matching pseudo wave function
!null vectors (with ss=0) can be added to improve kinetic residual
!convergence without changing match at rc
 do ii=1,nbas
  pswf0_or(ii)=0.0d0
  do jj=1,ncon
   pswf0_or(ii)=pswf0_or(ii)+work(jj)*vt(jj,ii)
  end do
 end do

!copy coefficient vectors spanning null space of constraint matrix
!note that vt (transpose) rows become columns of pswfnull
!write(6,'(/)')
 do jj=1,nbas-ncon
  sn=0.0d0
  do ii=1,nbas
   pswfnull_or(ii,jj)=vt(ncon+jj,ii)
   sn=sn+pswfnull_or(ii,jj)**2
  end do
 end do


!now, test result
 do ii=1,ncon
  con_tst(ii)=0.0d0
  do jj=1,nbas
   con_tst(ii)=con_tst(ii)+orbasis_der(ii,jj)*pswf0_or(jj)
  end do
 end do

 write(6,'(/a)') '    Constraint consistency test'
 write(6,'(a)') '      pswf0 val/der aewf val/der     diff'
 do ii=1,ncon
  write(6,'(4x,2f14.8,1p,d12.2)') con_tst(ii),cons(ii),con_tst(ii)-cons(ii)
 end do

!calculate charge contained in ps0 inside rc
 ps0norm=0.0d0
 do ii=1,nbas
  ps0norm=ps0norm+pswf0_or(ii)**2
 end do

!convert ps0 and psnull in orthogonal representation to sbf represention
 
 pswf0_sb(:)=0.0d0
 pswfnull_sb(:,:)=0.0d0
 do jj=1,nbas
  do kk=1,nbas
   pswf0_sb(jj)     =pswf0_sb(jj)     +orbasis(jj,kk)*pswf0_or(kk)
   do ii=1,nbas-ncon
    pswfnull_sb(jj,ii)=pswfnull_sb(jj,ii)+orbasis(jj,kk)*pswfnull_or(kk,ii)
   end do
  end do
 end do

 deallocate(amat,vt)
 return
 end subroutine const_basis
