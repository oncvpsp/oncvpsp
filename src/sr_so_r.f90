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
subroutine sr_so_r(lmax,irc,nproj,rr,mmax,mxprj,evkb,vkb, &
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
!mmax  dimension of log grid
!mxprj  dimension of number of projectors
!vkb  vkb projectors
!evkb  coefficients of BKB projectors
!vsr  normalized scalar projectors
!esr  energy  coefficients of vscal
!vso  normalized spin-orbig projectors
!esol  energy  coefficients of vso

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: lmax,mmax,mxprj
   integer :: irc(6),nproj(6)
   real(dp) :: rr(mmax),vkb(mmax,mxprj,4,2),evkb(mxprj,4,2)

!Output variables
   real(dp) :: vsr(mmax,2*mxprj,4),esr(2*mxprj,4)
   real(dp) :: vso(mmax,2*mxprj,4),eso(2*mxprj,4)

!Local variables
   integer :: ii,jj,kk,ik1,ik2,ip1,ip2,ipk,ll,l1,info,nn
   real(dp) :: apk,tt
   real(dp) :: sovl(2*mxprj,2*mxprj),sovlev(2*mxprj),ascl(2*mxprj,2*mxprj),aso(2*mxprj,2*mxprj)
   real(dp) :: asclst(2*mxprj,2*mxprj),wsclst(2*mxprj),asost(2*mxprj,2*mxprj),wsost(2*mxprj)
   real(dp) :: asclt(2*mxprj,2*mxprj),asot(2*mxprj,2*mxprj)
   real(dp) :: sphalf(2*mxprj,2*mxprj),smhalf(2*mxprj,2*mxprj)
   real(dp) :: fscl(mxprj),fso(mxprj),work(10*mxprj)
   real(dp), allocatable :: vkbt(:,:),vkbst(:,:)
   logical :: sorted

   allocate(vkbt(mmax,2*mxprj),vkbst(mmax,2*mxprj))

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

      call dsyev( 'V', 'U', nn, sovl, 2*mxprj, sovlev, work, 10*mxprj, info )

      if(info /= 0) then
         write(6,'(a,i4)') 'sr_so_r: S matrix eigenvalue ERROR, info=',info
         stop
      end if

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

      call dsyev( 'V', 'U', nn, asclst, 2*mxprj, wsclst, work, 10*mxprj, info )

      if(info /= 0) then
         write(6,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
         stop
      end if

      call dsyev( 'V', 'U', nn,  asost, 2*mxprj,  wsost, work, 10*mxprj, info )


      if(info /= 0) then
         write(6,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
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

      write(6,'(/a,i2)') &
      &         ' Orthonormal scalar projector coefficients, l = ',ll
      write(6,'(1p,6e12.4)') (esr(jj,l1),jj=1,nn)
      write(6,'(/a,i2)') &
      &         ' Orthonormal spin-orbit projector coefficients, l = ',ll
      write(6,'(1p,6e12.4)') (eso(jj,l1),jj=1,nn)

! Set sign of projectors (physically irrelevant) so that they are positive
! at their peak (needed for compaisons apparently)

      do jj=1,nn
         apk=0.0d0
         do ii=1,mmax
            if(abs(vso(ii,jj,l1))>apk) then
               apk=abs(vso(ii,jj,l1))
               ipk=ii
            end if
         end do
         if(vso(ipk,jj,l1)<0.0d0) then
            vso(:,jj,l1)=-vso(:,jj,l1)
         end if
         apk=0.0d0
         do ii=1,mmax
            if(abs(vsr(ii,jj,l1))>apk) then
               apk=abs(vsr(ii,jj,l1))
               ipk=ii
            end if
         end do
         if(vsr(ipk,jj,l1)<0.0d0) then
            vsr(:,jj,l1)=-vsr(:,jj,l1)
         end if
      end do

   end do  ! l1

   deallocate(vkbt,vkbst)
   return
end subroutine sr_so_r
