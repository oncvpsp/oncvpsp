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
subroutine vkboutwf(ll,nvkb,ep,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)

! computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions

!ll  angular momentum
!nvkb  switch for 1 or 2 projedtors
!ep  energy at which wave function is to be calculated
!vkb  Vanderbilt-Kleinman-Bylander projectors for this l
!evkb  projector coefficients
!rr  log radial mesh
!vloc  local pseudopotential
!uu  wave function
!up  1st derivative of uu
!node  count of number of nodes from 0 to rr(mch)
!mmax  dimension of log mesh
!mch  index of radius to which wave function is to be integrated

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   real(dp) :: rr(mmax),vloc(mmax),vkb(mmax,nvkb),evkb(nvkb)
   real(dp) :: ep
   integer :: nvkb,ll,mmax,mch

!Output variables
   real(dp) :: uu(mmax),up(mmax)
   integer :: node

!Local variables
   real(dp) :: rc,rn
   real(dp), allocatable ::  phi(:,:),phip(:,:)
   real(dp), allocatable ::  gg0(:),gg(:,:)
   integer, allocatable :: ipiv(:)

   integer :: ii,jj,ierr,info

   uu(:)=0.0d0
   up(:)=0.0d0

! homogeneous solution
   call lschps(ll,ierr,ep,rr,vloc,uu,up,mmax,mch)

   rc=0.0d0
   if(nvkb/=0) then
! find cutoff radius for projectors
      do ii=mmax,1,-1
         if(abs(vkb(ii,1))>0.0d0) then
            rc=rr(ii)
            exit
         end if
      end do

      allocate(phi(mmax,nvkb),phip(mmax,nvkb))
      allocate(gg(nvkb,nvkb),gg0(nvkb))
      allocate(ipiv(nvkb))

      phi(:,:)=0.0d0
      phip(:,:)=0.0d0
      gg(:,:)=0.0d0
      gg0(:)=0.0d0

! inhomogeneous solutions
      do ii=1,nvkb
         call lschkb(ll,ierr,ep,vkb(1,ii),rr,vloc,phi(1,ii),phip(1,ii),mmax,mch)
      end do


! projector matrix elements and coefficient matrix
      do jj=1,nvkb
         call vpinteg(uu,vkb(1,jj),mch,2*ll+2,gg0(jj),rr)
         gg0(jj)=evkb(jj)*gg0(jj)
         do ii=1,nvkb
            call vpinteg(phi(1,ii),vkb(1,jj),mch,2*ll+2,gg(jj,ii),rr)
            gg(jj,ii)=-evkb(jj)*gg(jj,ii)
         end do
         gg(jj,jj)=1.0d0+gg(jj,jj)
      end do

! solve linear equations for coefficients

!    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      call dgesv(nvkb, 1, gg, nvkb, ipiv, gg0, nvkb, info)
      if(info/=0) then
         write(6,'(/a,i4)') 'vkbout: dgesv ERROR, stopping info =',info
         stop
      end if

! output wave functions
      do jj=1,nvkb
         uu(:)=uu(:)+gg0(jj)*phi(:,jj)
         up(:)=up(:)+gg0(jj)*phip(:,jj)
      end do

      deallocate(phi,phip)
      deallocate(gg,gg0)
      deallocate(ipiv)
   end if

! lower cutoff for node counting to avoid small-r noise
   rn=dmax1(0.10d0*rc, 0.05d0)

   node=0
   do ii=6,mch
! note historic evolution!
!  if(rr(ii)>0.5d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.25d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.10d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.07d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
      if(rr(ii)>rn .and. uu(ii-1)*uu(ii)<0.0d0) then
         node=node+1
      end if
   end do

   return
end subroutine vkboutwf
