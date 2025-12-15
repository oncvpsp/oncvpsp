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
 subroutine sbf_basis(ll,rr,mmax,irc,nbas,qroot,sbasis,orbasis,orbasis_der)

!orthonormalize basis functions and derivatives at rc and find all-electron
!charge inside rc.

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!INPUT
!ll  angujlar momentum
!rr  log radial mesh
!mmax  number of points in log radial mesh
!irc  index rr such that rr(irc)=rc
!nbas  number of basis functions
!qroot  q values for j_l(q*r)

!OUTPUT
!sbasis  normalization coefficients for j_l basis set
!orbasis  matrix for orthonormal basis (rows j_l, columns basis vectors)
!orbasis_der  values and derivatives of orthonormal basis set at rc,

!Arguments
 integer :: ll,mmax,irc,nbas
 real(dp) :: rr(mmax),qroot(nbas)
 real(dp) :: sbasis(nbas),orbasis(nbas,nbas)
 real(dp) :: orbasis_der(6,nbas)

!Local variables
 integer :: ii,jj,ll1,ibas,info
 real(dp) :: al,amesh,ro,rc,sn,xx,tt
 real(dp) :: sb_out(10),sbfder(5)
 real(dp), allocatable :: sbfar(:,:),sev(:),work(:),sovlp(:,:)
 real(dp), allocatable :: sovlp_save(:,:)

 ll1=ll+1
 rc=rr(irc)
 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 allocate(sbfar(irc,nbas),sev(nbas),work(5*nbas),sovlp(nbas,nbas))
 allocate(sovlp_save(nbas,nbas))

 do ibas=1,nbas
  do jj=1,irc
   xx=qroot(ibas)*rr(jj)
   call sbf8(ll1,xx,sb_out)
   sbfar(jj,ibas)=rr(jj)*sb_out(ll1)
  end do
 end do

!perform sbf orthonormalization sum for overlap matrix

 sovlp(:,:)=0.0d0
 ro=rr(1)/sqrt(amesh)
 do ibas=1,nbas
  do jj=ibas,nbas

   sn=(sbfar(1,ibas)/rr(1)**ll)*(sbfar(1,jj)/rr(1)**ll)&
&     *ro**(2*ll+3)/dfloat(2*ll+3)
 
   do ii=1,irc-3
     sn=sn+al*rr(ii)*sbfar(ii,ibas)*sbfar(ii,jj)
   end do
 
   sn=sn + al*(23.0d0*rr(irc-2)*sbfar(irc-2,ibas)*sbfar(irc-2,jj)&
&            + 28.0d0*rr(irc-1)*sbfar(irc-1,ibas)*sbfar(irc-1,jj)&
&            +  9.0d0*rr(irc  )*sbfar(irc  ,ibas)*sbfar(irc  ,jj))/24.0d0
 
   sovlp(ibas,jj)=sn
   sovlp(jj,ibas)=sn
  end do
 end do

!normalization for j_l basis
 do ibas=1,nbas
  sbasis(ibas)=1.0d0/sqrt(sovlp(ibas,ibas))
 end do

 sovlp_save(:,:)=sovlp(:,:)
!find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

 call dsyev( 'V', 'U', nbas, sovlp, nbas, sev, work, 5*nbas, info )
 if(info .ne. 0) then
  write(6,'(a,i4)') 'sbf_basis: overlap matrix eigenvalue error, info=',info
  stop
 end if

!write(6,'(/a)') 'sbf overlap matrix eigenvalues'
!write(6,'(1p,6d12.3)') (sev(ibas),ibas=1,nbas)

!scale eigenvectors to form orthonormal basis coefficients for sbf's
!note that we are reversing order so that the leading eigenvector is the
!most linearly independent
 do ibas=1,nbas
  if(sev(ibas)>0.0d0) then
   tt=1.0d0/sqrt(sev(ibas))
  else
   write(6,'(a,f12.6)') 'sbfbasis: negative eigenvalue of overlap matrix'
   stop
  end if
  do jj=1,nbas
   orbasis(jj,nbas-ibas+1)=tt*sovlp(jj,ibas)
  end do
 end do

!find rc derivatives of basis funtction
 orbasis_der(:,:)=0.0d0
 do jj=1,nbas !sbf loop
  call sbf_rc_der(ll,qroot(jj),rc,sbfder)
  do ibas=1,nbas !orbasis loop
   do ii=1,5
    orbasis_der(ii,ibas)=orbasis_der(ii,ibas)+orbasis(jj,ibas)*sbfder(ii)
   end do
  end do
 end do

 deallocate(sbfar,sev,work,sovlp)

 return
 end subroutine sbf_basis
