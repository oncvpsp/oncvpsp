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

!> computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
!> at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
!> results for 1 and 2 projectors, or the semi-local pseudpotential
!> when that is the local potential for some l
!> the computed quantity is reminiscent of a scattering phase shift, but isn't
 subroutine run_phsft(lmax,lloc,nproj,epa,epsh1,epsh2,depsh,rxpsh,npsh,vkb,evkb, &
&                     rr,vfull,vp,zz,mmax,mxprj,irc,srel, &
&                     irpsh,rpsh,epsh,pshf,pshp)
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 implicit none

!Input variables
 !> maximum angular momentum
 integer, intent(in) :: lmax
 !> l for local potential
 integer, intent(in) :: lloc
 !> size of radial grid
 integer, intent(in) :: mmax
 !> dimension of number of projectors
 integer, intent(in) :: mxprj
 !> number ov V / KB projectors for  each l
 integer, intent(in) :: nproj(6)
 integer, intent(in) :: irc(6)
 !> low energy limit for "phase shift" calculation
 real(dp), intent(in) :: epsh1
 !> high energy limit for "phase shift" calculation
 real(dp), intent(in) :: epsh2
 !> energy increment
 real(dp), intent(in) :: depsh
 !> radius for phase shift calculation
 real(dp), intent(in) :: rxpsh
 !> Size of phase shift energies
 integer, intent(in) :: npsh
 !> atomic number
 real(dp), intent(in) :: zz
 !> log radial grid
 real(dp), intent(in) :: rr(mmax)
 !> semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
 real(dp), intent(in) :: vp(mmax,5)
 !> bound-state or scattering state reference energies for vkb potentials
 real(dp), intent(in) :: epa(mxprj,6)
 !> all-electron potential
 real(dp), intent(in) :: vfull(mmax)
 !> VKB projectors
 real(dp), intent(in) :: vkb(mmax,mxprj,4)
 !> coefficients of VKB projectors
 real(dp), intent(in) :: evkb(mxprj,4)
 !> .true. for scalar-relativistic, .false. for non-relativistic
 logical, intent(in) :: srel

!Output variables
 !> Log radial grid mesh index at which log derivative is calculated
 integer, intent(out) :: irpsh(4)
 !> Radius at which log derivative is calculated
 real(dp), intent(out) :: rpsh(4)
 !> Phase shift (log derivative) energies
 real(dp), intent(out) :: epsh(npsh)
 !> All-electron phase shifts (log derivatives)
 real(dp), intent(out) :: pshf(npsh,4)
 !> Pseudo phase shifts (log derivatives)
 real(dp), intent(out) :: pshp(npsh,4)

!Local variables
 integer :: ii,ll,l1

! loop for phase shift calculation -- full, then local or Kleinman-
! Bylander / Vanderbilt

 do ii = 1, npsh
    epsh(ii) = epsh2 - real(ii - 1, dp) * depsh
 end do

 do l1 = 1, 4
   ll = l1 - 1
   !> index of rr beyond which all vp==vlocal
   if(ll<=lmax) then
     irpsh(l1) = irc(l1) + 2
   else
     irpsh(l1) = irc(lloc + 1)
   end if
   rpsh(l1) = rr(irpsh(l1))

   if (rxpsh > 0.0_dp) then
     if (rxpsh < rpsh(l1)) then
        write(stderr, '(a,f12.5,a,i1,a,f12.5)') &
          'ERROR: run_phsft: rxpsh=', rxpsh, ' < rc(', l1 - 1, ')=', rpsh(l1)
        stop
     end if
     rpsh(l1) = rxpsh
     do ii = 1, mmax
       if (rr(ii) .ge. rpsh(l1)) then
         irpsh(l1) = ii
         exit
       end if
     end do
   end if

   call fphsft(ll,epsh2,depsh,pshf(:,l1),rr,vfull,zz,mmax,irpsh(l1),npsh,srel)
   if(ll .eq. lloc) then
     call  vkbphsft(ll,0,epsh2,depsh,epa(1,l1),pshf(:,l1),pshp(:,l1), &
&                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                   mmax,irpsh(l1),npsh)
   else
     call  vkbphsft(ll,nproj(l1),epsh2,depsh,epa(1,l1),pshf(:,l1),pshp(:,l1), &
&                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                   mmax,irpsh(l1),npsh)
   end if
 end do
 return
 end subroutine run_phsft
