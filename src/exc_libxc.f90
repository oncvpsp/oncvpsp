!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University, and M. Verstratte, University of Liege
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

subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax)
!
!  interface to libxc library for potential and energy density calculation
!
   use functionals_m
   implicit none

   integer, parameter :: dp=kind(1.0d0)

   real(dp), intent(in) :: al
   real(dp), intent(in) :: rho(mmax),rr(mmax)
   real(dp), intent(out) :: vxc(mmax),exc(mmax)
   integer, intent(in) :: mmax, iexc
   integer :: ii

! local
   integer ::  id(2), ifunc
   integer, save :: functprinted = 0
   type(xc_functl_t) :: functls(2)
   integer :: nspin, irel, deriv_method
   real(dp),allocatable :: ip(:)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

   real(dp), allocatable :: v(:), e(:)
   real(dp), allocatable :: dpr(:), dppr(:), dlap(:)
   real(dp), allocatable :: tau_dummy(:), vtau_dummy(:)
   real(dp), allocatable :: rho_loc(:)

   allocate(v(mmax), e(mmax))
   allocate(dpr(mmax), dppr(mmax), dlap(mmax))
   allocate(tau_dummy(mmax), vtau_dummy(mmax))
   allocate(rho_loc(mmax))
!
!conf = (3.d0*pi**2)**thrd
!conrs = (3.d0/(4.d0*pi))**thrd

   rho_loc = rho(1:mmax) /4.0d0 / pi

   call derivs(mmax, rho_loc, al, rr, dpr, dppr, dlap)

!
! determine xc to be used from input code (negative by convention)

   id(1) = abs(iexc)/1000
   id(2) = abs(iexc) - id(1)*1000
   nspin = 1
   irel = 0  !XC_NON_RELATIVISTIC
!deriv_method = XC_DERIV_ANALYTICAL
   deriv_method = XC_DERIV_NUMERICAL
   ! call xc_functl_init(functls(1),nspin)
   ! call xc_functl_init(functls(2),nspin)
   call xc_functl_init_functl(functls(1),id(1),nspin,deriv_method)
   call xc_functl_init_functl(functls(2),id(2),nspin,deriv_method)

   if (functprinted==0) then
      call xc_functl_write_info(functls(1), 6)
      call xc_functl_write_info(functls(2), 6)
      write (*,*)
      functprinted = 1
   end if

! no use of vtau kin en density optional argument
! ip is an ionization potential parameter for certain xc functionals, not supported here. Set to 0.5 Ha
   allocate (ip(nspin))
   ip = 0.5_dp
   tau_dummy = 0.0d0
   vtau_dummy = 0.0d0
   vxc(1:mmax) = 0.0d0
   exc(1:mmax) = 0.0d0
   do ifunc = 1, 2
      v = 0.0d0
      e = 0.0d0
      call xc_functl_get_vxc(functls(ifunc), mmax, al, rr, rho_loc, dpr, dlap, &
      &        tau_dummy, ip, v, e, vtau_dummy)
      vxc(1:mmax) = vxc(1:mmax) + v(1:mmax)
      exc(1:mmax) = exc(1:mmax) + e(1:mmax)
   end do

! go to convention for Don Hamman code - I do not understand:
! in principle libxc outputs energy per particle,
! and here we divide by rho again to get an LDA exc
! identical to the other routine...
   do ii=1,mmax
      exc(ii) = exc(ii) / max(rho_loc(ii),1.d-20)
   end do

   ! finish
   call xc_functl_end(functls(1))
   call xc_functl_end(functls(2))
   deallocate (ip)
   deallocate(v, e)
   deallocate(dpr, dppr, dlap)
   deallocate(tau_dummy, vtau_dummy)
   deallocate(rho_loc)

end subroutine exc_libxc


subroutine derivs(mmax, rho, al, rr, dpr, dppr, dlap)
! NB: presumes incoming rho is properly scaled without 4pi or r**2 factors.

   implicit none
   integer, parameter :: dp=kind(1.0d0)

   integer, intent(in) :: mmax
   real(dp), intent(in) :: al
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: rho(mmax)
   real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)

! local vars
   integer :: i
   real(dp) :: dpn(mmax), dppn(mmax)
   real(dp) :: c11,c12,c13,c14,c15
   real(dp) :: c21,c22,c23,c24,c25

   c11 =   2.0d0 / 24.0d0
   c12 = -16.0d0 / 24.0d0
   c13 =   0.0d0 / 24.0d0
   c14 =  16.0d0 / 24.0d0
   c15 =  -2.0d0 / 24.0d0
!
   c21 =   -1.0d0 / 12.0d0
   c22 =   16.0d0 / 12.0d0
   c23 =  -30.0d0 / 12.0d0
   c24 =   16.0d0 / 12.0d0
   c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!
   i=1
   dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2) &
   &         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
   dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
   &         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)
   i=2
   dpn(i) = -25.d0/12.d0*rho(i) +4.d0*rho(i+1) -3.d0*rho(i+2)  &
   &         +4.d0/3.d0*rho(i+3) -1.d0/4.d0*rho(i+4)
   dppn(i) = 15.d0/4.d0*rho(i) -77.d0/6.d0*rho(i+1) +107.d0/6.d0*rho(i+2) &
   &         -13.d0*rho(i+3) +61.d0/12.d0*rho(i+4) -5.d0/6.d0*rho(i+5)

   do i = 3, mmax - 2
      dpn(i) =  c11*rho(i-2) + c12*rho(i-1) + c14*rho(i+1) + c15*rho(i+2)
      dppn(i) = c21*rho(i-2) + c22*rho(i-1) + c23*rho(i)   + c24*rho(i+1) &
      &           +c25*rho(i+2)
   end do

   i=mmax-1
   dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
   &         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
   dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
   &          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)
   i=mmax
   dpn(i) = +25.d0/12.d0*rho(i) -4.d0*rho(i-1) +3.d0*rho(i-2) &
   &         -4.d0/3.d0*rho(i-3) +1.d0/4.d0*rho(i-4)
   dppn(i) = -15.d0/4.d0*rho(i) +77.d0/6.d0*rho(i-1) -107.d0/6.d0*rho(i-2) &
   &          +13.d0*rho(i-3) -61.d0/12.d0*rho(i-4) +5.d0/6.d0*rho(i-5)

!
! r derivatives of d
!
   do i = 1, mmax
      dpr(i) = dpn(i) / (al * rr(i))
      dppr(i) = (dppn(i) - al * dpn(i)) / (al * rr(i))**2
      dlap(i) = (dppn(i) + al * dpn(i)) / (al * rr(i))**2
   end do

end subroutine derivs
