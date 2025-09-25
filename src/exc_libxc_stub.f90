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
!  dummy substitute for interface to libxc library to allow successful
!  build without this library
!
   implicit none

   integer, parameter :: dp=kind(1.0d0)

   real(dp), intent(in) :: al
   real(dp), intent(in) :: rho(mmax),rr(mmax)
   real(dp), intent(out) :: vxc(mmax),exc(mmax)
   integer, intent(in) :: mmax, iexc

   write(6,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc,' requires libxc'
   write(6,'(a)') 'The present oncvpsp executable was built without libxc.'
   write(6,'(a)') 'Program will stop.'

   stop
   return
end subroutine exc_libxc
