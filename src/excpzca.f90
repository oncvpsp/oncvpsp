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
subroutine excpzca(rho,vxc,exc,mmax)

!calculates Perdew-Zunger-Ceperly_Alder exchange-correlation potential
! and energy density J. P. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)

!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!mmax  dimension of log radial grid

   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi4=4.0d0*pi
   real(dp), parameter :: pi4i=1.0d0/pi4
   real(dp), parameter :: thrd=1.0d0/3.0d0

!Input variables
   integer :: mmax
   real(dp) :: rho(mmax)

!Output variables
   real(dp) :: vxc(mmax),exc(mmax)

!Local variables
   integer :: ii
   real(dp) :: rs,rh,sqrs,rsl,den

   real(dp), parameter :: conrs = (3.d0/(4.d0*pi))**thrd
   real(dp), parameter :: convx = (1.5d0/pi)**(2.0d0/3.0d0)
   real(dp), parameter ::conex = 0.75d0*convx


!rs>1, Ceperley-Alder fit to QMC
   real(dp), parameter :: gam = -0.1423d0
   real(dp), parameter :: bt1 =  1.0529d0
   real(dp), parameter :: bt2 =  0.3334d0

!rs<1, Perdew-Zunger leading terms of RPA high-density expansion
!AA and BB known from Gell-Mann & Bruckner
!PZ adjust CC and DD for contiunity of exc and vxc at rs=1
   real(dp), parameter :: AA =   0.0311d0
   real(dp), parameter :: BB =  -0.0480d0
   real(dp), parameter :: CC =   0.0020d0
   real(dp), parameter :: DD =  -0.0116d0

!PZ modified coeffients from libxc
!real(dp), parameter :: CC =   0.0020191519406228d0
!real(dp), parameter :: DD =  -0.0116320663789130d0

   do ii=1,mmax
      if(rho(ii)>1.0d-20) then
         rh=rho(ii)/pi4
         rs=conrs*rh**(-thrd)
         if(rs > 1.0d0) then
            sqrs=dsqrt(rs)
            den=1.0d0 + bt1*sqrs +bt2*rs
            exc(ii)=-conex/rs + gam/den
            vxc(ii)=-convx/rs + gam*(1.0d0 + (3.5d0*bt1*sqrs &
            &            + 4.0d0*bt2*rs)*thrd)/den**2
         else
            rsl=dlog(rs)
            exc(ii)=-conex/rs + AA*rsl + BB + CC*rs*rsl + DD*rs
            vxc(ii)=-convx/rs + AA*rsl + (BB-thrd*AA) + 2.0d0*thrd*CC*rs*rsl &
            &            + thrd*(2.0d0*DD-CC)*rs
         end if
      else
         vxc(ii)=0.0d0 ; exc(ii)=0.0d0
      end if
   end do

   return
end subroutine excpzca
