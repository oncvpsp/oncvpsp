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
subroutine vkbphsft(ll,ivkb,epsh2,depsh,ep,pshf,pshp, &
& rr,vloc,vkb,evkb,mmax,mch,npsh)

! computes Vanderfilt / Kleinman-Bylander scattering log derivatives
! (or semi-local if ivkb==0)
! returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
! and easier to compare in plots than the log derivatives themselves
! Pauli-type scalar-relativistic calculation

!ll  angular momentum
!ivkb  number of projectors
!epsh2  upper limit of energy scan
!depsh  increment of scan
!ep  reference energy for psp creation (bound or scattering)
!pshf  all-electron log derivatives "angles", as above (input)
!pshp  pseudopotential log derivatives "angles", as above (output)
!rr  radial log grid
!vloc  local part of psp
!vkb  VKB projectors
!mmax  dimension of rr, etc.
!mch  index of radius for log der test
!npsh  number of energy points in scan

   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi2=2.0d0*pi
   real(dp), parameter :: eps=1.0d-8

!Input variables
   integer :: ll,ivkb,mmax,npsh,mch
   real(dp) :: rr(mmax),vloc(mmax),vkb(mmax,*),evkb(*)
   real(dp) :: pshf(npsh)
   real(dp) :: depsh,epsh2,ep

!Output variables
   real(dp) :: pshp(npsh)

!Local variables
   real(dp) :: al,dnpi,epsh,phi,phip,pshoff
   real(dp) :: dmin,dmax,emin,emax,et,psmin,psmax,pst,shift,shift2
   integer :: ii,jj,node
   logical :: jump

   real(dp), allocatable :: uu(:),up(:)



   allocate(uu(mmax),up(mmax))

   al = 0.01d0 * dlog(rr(101) / rr(1))

   pshoff = 0.0d0

   do ii = 1,npsh
      epsh = epsh2-(ii-1)*depsh

      call vkboutwf(ll,ivkb,epsh,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)

      phi = uu(mch)/rr(mch)
      phip = (up(mch)-al*uu(mch))/(al*rr(mch)**2)
      pshp(ii) = atan2(rr(mch)*phip,phi) + pshoff

! Shift for continuity and to avoid false jumps suggesting ghost states

! -2 pi jump is harmless cycling from pi to -pi
      if(ii>1) then
         if(pshp(ii)<pshp(ii-1)-4.0d0) then
            pshoff = pshoff+pi2
            pshp(ii) = pshp(ii)+pi2
            ! +/- pi jump can be actual bound semi-core state, ghost state, or spurious
            ! jump from sudden change of sign of both uu and up with no real change
            ! of shape.
         else if(abs(pshp(ii)-pshp(ii-1))>2.0d0) then
            shift=sign(pi,pshp(ii)-pshp(ii-1))

            ! interval-halving search to determine if this is a discontinuous pi jump
            ! or varies continuously on some scale, indicating a real or ghost
            ! bound state / resonance
            jump=.true.
            emin=epsh
            emax=epsh+depsh
            psmin=pshp(ii)
            psmax=pshp(ii-1)
!       if((psmax-psmin)>pi) psmin=psmin+pi2
            if((psmax-psmin)>2.9d0) psmin=psmin+pi2

            do jj=1,25

               if(.not. jump) cycle
               et=0.5d0*(emin+emax)
               shift2=0.0d0

               call vkboutwf(ll,ivkb,et,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)
               phi = uu(mch)/rr(mch)
               phip = (up(mch)-al*uu(mch))/(al*rr(mch)**2)
               pst = atan2(rr(mch)*phip,phi) + pshoff
               if((psmax-pst)>2.9d0) then
                  pst=pst+pi2
                  shift2=pi2
               end if

               dmin=abs(pst-psmin)
               dmax=abs(psmax-pst)

               if(dmin>dmax) then
                  emax=et
                  psmax=pst
               else
                  emin=et
                  psmin=pst
               end if

               ! this is the test to see if an interval of change less than pi had been
               ! reached
               if(dmax<0.5d0 .and. dmin<0.5d0) jump=.false.

            end do

            ! if this is a spurious abrupt jump, restore +/- pi
            if(jump) then
               pshoff = pshoff-shift
               pshp(ii) = pshp(ii)-shift
               ! if this is a real continuous transition, check for 2 pi issue and fix
            else if(pshp(ii)<pshp(ii-1)-2.5d0) then
               pshoff = pshoff+pi2
               pshp(ii) = pshp(ii)+pi2
            end if

         end if
      end if

! calculate shift to align with all-electron results
      if(abs(ep-epsh)-0.5d0*depsh<eps) then
         dnpi = pi*(nint((pshf(ii)-pshp(ii))/pi))
      end if

   end do

   pshp(:)=pshp(:)+dnpi

   deallocate(uu,up)
   return
end subroutine vkbphsft
