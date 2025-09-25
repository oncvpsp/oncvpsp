
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
subroutine run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj,pswf, &
&                   vfull,vp,evkb,vkb,nlim,vr)

! computes Vanderbilt / Kleinman-Bylander non-local potentials

!lmax  maximum angular momentum
!lloc  l for local potential (lloc==4 => use linear combination)
!lpopt  choice of polynomial for lloc==4
!dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!irc  core radii indices
!nproj  number of projectors for each l
!rr  log radial grid
!mmax  size of radial grid
!mxprj  dimension of number of projectors
!pswf pseudo wave functions
!vfull  all-electron potential
!vp  semi-local pseudopotentials for first projectors
!     (vp(:,5) is local potential if lloc=4)
!vkb  semi-local "potentials"*pswf (input); VKB projectors (output)
!evkb  coefficients of BKB projectors (output)
!nlim  index of maximum rc including that of vlocal (output)
!vr  effective scalar-relativisic "potential" calculated in vrel

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: lmax,lloc,lpopt,lwork,mmax,mxprj
   integer :: irc(6),nproj(6)
   real(dp) :: rr(mmax),pswf(mmax,mxprj,5),vfull(mmax),vp(mmax,5)
   real(dp) :: dvloc0

!Input/Output variables
   real(dp) :: vkb(mmax,mxprj,4),evkb(mxprj,4)
   real(dp) :: vr(mmax,mxprj,6)
   integer :: nlim

!Local variables
   integer :: ii,ipk,jj,kk,l1,info,np
   real(dp) :: apk,sn,tt
   real(dp) :: xx,ff
   real(dp) :: bb(mxprj,mxprj),bbev(mxprj),bbi(mxprj,mxprj)
   real(dp) :: sovl(mxprj,mxprj),sovlev(mxprj),smhalf(mxprj,mxprj)
   real(dp) :: sphalf(mxprj,mxprj)
   real(dp) :: bbist(mxprj,mxprj),bbistev(mxprj),bbit(mxprj,mxprj)
   real(dp), allocatable :: vloc(:),vkbt(:,:),vkbst(:,:),work(:)
   logical :: sorted

   evkb(:,:)=0.0d0

   lwork=3*mxprj
   allocate(vloc(mmax),vkbt(mmax,mxprj),vkbst(mmax,mxprj),work(lwork))


! calculate local potential

   if(lloc==4) then
      call vploc(rr,vfull,vp,dvloc0,irc(5),mmax,lpopt)
   end if

   vloc(:)=vp(:,lloc+1)

! Vanderbilt / Kleinman-Bylander projector construction
   nlim=irc(lloc+1)
   do l1=1,lmax+1
      nlim=max(nlim,irc(l1))
      if(l1==lloc+1) cycle
      do jj=1,nproj(l1)
!incorporates scalar-relativistic "potential" in last 5% of [0:rc] to smoothly
!force projectors to zero
         do ii=1,irc(l1)
            xx=20.0d0*((rr(ii)/rr(irc(l1)))-1.0d0)
            if(xx>-1.0d0) then
               ff=(1.0d0-xx**2)**2
            else
               ff=0.0d0
            end if
            vkb(ii,jj,l1)=vkb(ii,jj,l1)-(vloc(ii)+ff*vr(ii,jj,l1))*pswf(ii,jj,l1)
         end do
         do ii=irc(l1)+1,mmax
            vkb(ii,jj,l1)=0.0d0
         end do
      end do
   end do

! Vanderbilt B-matrix construction

   do l1=1,lmax+1
      if(l1==lloc+1) cycle
      np=nproj(l1)
      bb(:,:)=0.0d0
      call vpinteg(pswf(1,1,l1),vkb(1,1,l1),irc(l1),2*l1,bb(1,1),rr)
      evkb(1,l1)=1.0d0/bb(1,1)
! this is Kleinman-Bylander result if nproj==1

      if(nproj(l1)>=2) then

         do jj=1,np
            do ii=1,np
               call vpinteg(pswf(1,ii,l1),vkb(1,jj,l1),irc(l1),2*l1,bb(ii,jj),rr)
            end do
         end do

! symmetrize exactly
         write(6,'(/a,i4)') 'B matrix Hermiticity error, ll=',l1-1
         do jj=2,np
            do ii=1,jj-1
               write(6,'(2i4,1p,e14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
               tt=0.5d0*(bb(ii,jj)+bb(jj,ii))
               bb(ii,jj)=tt
               bb(jj,ii)=tt
            end do
         end do

! find eigenvalues and eigenvectors of the B matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev( 'V', 'U', np, bb, mxprj, bbev, work, lwork, info )
         if(info /= 0) then
            write(6,'(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=',info
            stop
         end if

! take linear combinations to form diagonal projectors

         do jj=1,np
            vkbt(:,jj)=vkb(:,jj,l1)
         end do

         do jj=1,np
            vkb(:,jj,l1)=0.0d0
            evkb(jj,l1)=1.0d0/bbev(jj)
            do ii=1,np
               vkb(:,jj,l1)=vkb(:,jj,l1)+bb(ii,jj)*vkbt(:,ii)
            end do
         end do

      end if  !nproj>=2

! normalize projectors
      do jj=1,np
         call vpinteg(vkb(1,jj,l1),vkb(1,jj,l1),irc(l1),2*l1,sn,rr)
         if(vkb(irc(l1)-10,jj,l1)>=0.0d0) then
            vkb(:,jj,l1)=vkb(:,jj,l1)/sqrt(sn)
         else
            vkb(:,jj,l1)=-vkb(:,jj,l1)/sqrt(sn)
         end if
         evkb(jj,l1)=sn*evkb(jj,l1)
      end do


! calculate overlap
      if(nproj(l1)>=2) then


         bbi(:,:)=0.0d0
         do jj=1,np
            sovl(jj,jj)=1.0d0
            bbi(jj,jj)=evkb(jj,l1)
            do ii=jj+1,np
               call vpinteg(vkb(1,ii,l1),vkb(1,jj,l1),irc(l1),2*l1,sn,rr)
               sovl(ii,jj)=sn
               sovl(jj,ii)=sn
            end do
         end do
         call vpinteg(vkb(1,1,l1),vkb(1,2,l1),irc(l1),2*l1,sn,rr)

! find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev( 'V', 'U', np, sovl, mxprj, sovlev, work, lwork, info )
         if(info /= 0) then
            write(6,'(a,i4)') 'run_vkb: sovl matrix ERROR error, info=',info
            stop
         end if

!   write(6,'(/1p,2e12.4)') (sovlev(ii),ii=1,2)

! construct S^(-1/2) AND s^(1/2)

         do jj=1,np
            tt=sqrt(sovlev(jj))
            do ii=1,np
               sphalf(ii,jj)=tt*sovl(ii,jj)
               smhalf(ii,jj)=sovl(ii,jj)/tt
            end do
         end do

! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

         bbit(:,:)=0.0d0
         bbist(:,:)=0.0d0

         do ii=1,np
            do jj=1,np
               do kk=1,np
                  bbit(ii,jj)=bbit(ii,jj)+bbi(ii,kk)*sphalf(kk,jj)
               end do
            end do
         end do

         do ii=1,np
            do jj=1,np
               do kk=1,np
                  bbist(ii,jj)=bbist(ii,jj)+bbit(kk,jj)*sphalf(kk,ii)
               end do
            end do
         end do

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev( 'V', 'U', np , bbist, mxprj, bbistev, work, lwork, info )
         if(info /= 0) then
            write(6,'(a,i4)') 'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=',info
            stop
         end if


! take linear combinations to form orthonormal basis functions

         vkbt(:,:)=vkb(:,:,l1)
         vkbst(:,:)=0.0d0

         do jj=1,np
            do ii=1,np
               vkbst(:,jj)=vkbst(:,jj) + smhalf(ii,jj)*vkbt(:,ii)
            end do
         end do


! take linear combinations to form orthonormal projectors

         vkb(:,:,L1)=0.0d0

         do ii=1,np
            evkb(ii,l1)=bbistev(ii)
            do jj=1,np
               vkb(:,ii,l1)=vkb(:,ii,l1) + bbist(jj,ii)*vkbst(:,jj)
            end do
         end do

! re-order if necessary so first projector has largest magnitude coefficient
! which is consistent with relativistic sr_so_r.f90 output

         do ii=1,100
            sorted=.true.
            do jj=2,np

               if(abs(evkb(jj,l1))>abs(evkb(jj-1,l1))) then
                  tt=evkb(jj-1,l1)
                  vkbt(:,1)=vkb(:,jj-1,l1)
                  evkb(jj-1,l1)=evkb(jj,l1)
                  vkb(:,jj-1,l1)=vkb(:,jj,l1)
                  evkb(jj,l1)=tt
                  vkb(:,jj,l1)=vkbt(:,1)
                  sorted=.false.
               end if
            end do
            if(sorted) exit
         end do

         write(6,'(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
         &        (evkb(jj,l1),jj=1,np)

! Set sign of projectors (physically irrelevant) so that they are positive
! at their peak (needed for compaisons apparently)

         do jj=1,np
            apk=0.0d0
            do ii=1,mmax
               if(abs(vkb(ii,jj,l1))>apk) then
                  apk=abs(vkb(ii,jj,l1))
                  ipk=ii
               end if
            end do
            if(vkb(ipk,jj,l1)<0.0d0) then
               vkb(:,jj,l1)=-vkb(:,jj,l1)
            end if
         end do

      end if  !nproj(l1)>=2

      if(nproj(l1)<mxprj) then
         do jj=nproj(l1)+1, mxprj
            evkb(jj,l1)=0.0d0
            vkb(:,jj,l1)=0.0d0
         end do
      end if

   end do  ! l1

   deallocate(vloc,vkbt,vkbst,work)

   return
end subroutine run_vkb
