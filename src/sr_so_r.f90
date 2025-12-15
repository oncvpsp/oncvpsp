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
 subroutine sr_so_r(lmax,irc,nproj,rr,mmax,evkb,vkb, &
&                       vsr,esr,vso,eso)

! reformulates non-local potentials based on j = l +/- 1/2 to scalar-
! relativistic and L dot S projectors
! uses relationship <L dot S> = (J^2 - L^2 - S^2)/2
! so L dot S = +/- l/2 for j = l +/- 1/2

!lmax  maximum angular momentum
!irc  core radii indices
!nproj  number of projectors for each l
!rr  log radial grid
!mmax  size of radial grid
!vkb  vkb projectors
!evkb  coefficients of BKB projectors
!vsr  normalized scalar projectors
!esr  energy  coefficients of vscal
!vso  normalized spin-orbig projectors
!esol  energy  coefficients of vso

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,lpopt,mmax
 integer :: irc(6),nproj(6)
 real(dp) :: rr(mmax),vkb(mmax,2,4,2),evkb(2,4,2)

!Output variables
 real(dp) :: vsr(mmax,4,4),esr(4,4),vso(mmax,4,4),eso(4,4)

!Local variables
 integer :: ii,jj,kk,ierr,ik1,ik2,ip1,ip2,ll,l1,info,nn
 real(dp) :: sn,tt,qq12
 real(dp) :: sovl(4,4),sovlev(4),ascl(4,4),aso(4,4)
 real(dp) :: asclst(4,4),wsclst(4),asost(4,4),wsost(4)
 real(dp) :: asclt(4,4),asot(4,4)
 real(dp) :: sphalf(4,4),smhalf(4,4)
 real(dp) :: fscl(2),fso(2),work(20)
 real(dp), allocatable :: vkbt(:,:),vkbst(:,:)
 logical :: sorted

 allocate(vkbt(mmax,4),vkbst(mmax,4))

 do l1=1,lmax+1
  ll=l1-1

  if(ll==0) then
   vsr(:,:,l1)=0.0d0
   vso(:,:,l1)=0.0d0
   esr(:,l1)=0.0d0
   eso(:,l1)=0.0d0
   if(nproj(l1)>=1) then
    do ii=1,nproj(l1)
     vsr(:,ii,l1)=vkb(:,ii,l1,1)
     esr(ii,l1)=evkb(ii,l1,1)
    end do
   end if
   cycle
  end if

  nn=2*nproj(l1)

  fscl(1)=(ll+1)/dble(2*ll+1)
  fscl(2)=ll/dble(2*ll+1)
  fso(1)=2/dble(2*ll+1)
  fso(2)=-2/dble(2*ll+1)

! construct overlap matrix and diagonal energy matrices

   sovl(:,:)=0.0d0 
   ascl(:,:)=0.0d0 
   aso(:,:)=0.0d0
   vkbt(:,:)=0.0d0

   do ik1=1,2
    do ip1=1,nproj(l1)
     ii=ip1+(ik1-1)*nproj(l1)
     
     ascl(ii,ii)=fscl(ik1)*evkb(ip1,l1,ik1)
     aso(ii,ii)=fso(ik1)*evkb(ip1,l1,ik1)

     vkbt(:,ii)=vkb(:,ip1,l1,ik1)

     do ik2=1,2
      do ip2=1,nproj(l1)
       jj=ip2+(ik2-1)*nproj(l1)

       call vpinteg(vkb(1,ip1,l1,ik1),vkb(1,ip2,l1,ik2),irc(l1),2*l1, &
&                   sovl(ii,jj),rr)

      end do
     end do
    end do
   end do

   call dsyev( 'V', 'U', nn, sovl, 4, sovlev, work, 20, info )

! construct S^(-1/2) AND s^(1/2)

   do jj=1,nn
    tt=sqrt(sovlev(jj))
    do ii=1,nn
     sphalf(ii,jj)=tt*sovl(ii,jj)
     smhalf(ii,jj)=sovl(ii,jj)/tt
    end do
   end do

! take linear combinations to form orthonormal basis functions

   vkbst(:,:)=0.0d0

   do jj=1,nn
    do ii=1,nn
     vkbst(:,jj)=vkbst(:,jj) + smhalf(ii,jj)*vkbt(:,ii)
    end do
   end do

! construct A^(-1)* = S^(1/2)^T A^(-1) S^(1/2)

     asclt(:,:)=0.0d0
     asclst(:,:)=0.0d0
     asot(:,:)=0.0d0
     asost(:,:)=0.0d0

     do ii=1,nn
      do jj=1,nn
       do kk=1,nn
        asclt(ii,jj)=asclt(ii,jj)+ascl(ii,kk)*sphalf(kk,jj)
        asot(ii,jj) =asot(ii,jj) + aso(ii,kk)*sphalf(kk,jj)
       end do
      end do
     end do

     do ii=1,nn
      do jj=1,nn
       do kk=1,nn
        asclst(ii,jj)=asclst(ii,jj)+asclt(kk,jj)*sphalf(kk,ii)
        asost(ii,jj) =asost(ii,jj) + asot(kk,jj)*sphalf(kk,ii)
       end do
      end do
     end do

! find eigenvalues and eigenvectors of the A* matrices

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', nn, asclst, 4, wsclst, work, 20, info )

     call dsyev( 'V', 'U', nn,  asost, 4,  wsost, work, 20, info )


     if(info .ne. 0) then
      write(6,'(a,i4)') 'sr_so_r: A* matrix eigenvalue error, info=',info
      stop
     end if

! take linear combinations to form orthonormal projectors

     vsr(:,:,l1)=0.0d0
     vso(:,:,l1)=0.0d0
     esr(:,l1)=0.0d0
     eso(:,l1)=0.0d0

     do ii=1,nn
      esr(ii,l1)=wsclst(ii)
      eso(ii,l1)=  wsost(ii)
      do jj=1,nn
       vsr(:,ii,l1)=vsr(:,ii,l1) + asclst(jj,ii)*vkbst(:,jj)
       vso(:,ii,l1)= vso(:,ii,l1) +   asost(jj,ii)*vkbst(:,jj)
      end do
     end do

! bubble-sort on coefficient magnitudes for scalar and then s-o
! (Yes, I know bubble-sort is the least-efficient sorting algorithm.)

     do ii=1,100
      sorted=.true.
      do jj=2,nn
       if(abs(esr(jj-1,l1))<abs(esr(jj,l1))) then
        tt=esr(jj,l1)
        vkbt(:,1)=vsr(:,jj,l1)
        esr(jj,l1)=esr(jj-1,l1)
        vsr(:,jj,l1)=vsr(:,jj-1,l1)
        esr(jj-1,l1)=tt
        vsr(:,jj-1,l1)=vkbt(:,1)
        sorted=.false.
       end if
      end do
      if(sorted) exit
     end do
  
     do ii=1,100
      sorted=.true.
      do jj=2,nn
       if(abs(eso(jj-1,l1))<abs(eso(jj,l1))) then
        tt=eso(jj,l1)
        vkbt(:,1)=vso(:,jj,l1)
        eso(jj,l1)=eso(jj-1,l1)
        vso(:,jj,l1)=vso(:,jj-1,l1)
        eso(jj-1,l1)=tt
        vso(:,jj-1,l1)=vkbt(:,1)
        sorted=.false.
       end if
      end do
      if(sorted) exit
     end do

     write(6,'(/a,i2/1p,4e12.4)') &
&         ' Orthonormal scalar projector coefficients, l = ',ll, &
&        (esr(jj,l1),jj=1,nn)
     write(6,'(/a,i2/1p,4e12.4)') &
&         ' Orthonormal spin-orbit projector coefficients, l = ',ll, &
&        (eso(jj,l1),jj=1,nn)

 end do ! l1

 deallocate(vkbt,vkbst)
 return
 end subroutine sr_so_r
