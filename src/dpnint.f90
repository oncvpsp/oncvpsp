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
subroutine dpnint(xx, yy, nn, tt, ss, mm)

! local polynomial interpolation of data yy on nn points xx
! giving values ss on mm points tt
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output mm interpolated values ss on points tt

   implicit none

   integer, parameter :: dp=kind(1.0d0)

!Input variables
   real(dp) ::  xx(*),yy(*),tt(*)
   integer :: nn,mm

!Output variables
   real(dp) :: ss(*)

!Local variables
   real(dp) :: sum,term,zz
   integer :: ii,imin,imax,iprod,iy,istart,jj,kk

! set order of polynomial
   integer, parameter :: npoly=7

   if(nn<npoly+1) then
      write(6,'(/a,i6,a,i4)') 'dpnint: interpolation ERROR, n=', &
      &       nn,'< npoly=',npoly
      stop
   end if

! note: output point 1 is skipped in this version because of special
! properties of pp data (ie., log grid)
! this point receives special treatment (see below)

   ss(1:mm)=0.0d0
   imin = 1
   do jj = 2, mm
      if(tt(jj)<xx(1)) then
         write(6,'(/a)') 'dpnint: interpolation ERROR - out of range'
         stop
      end if
      if(tt(jj)>xx(nn)) then
         write(6,'(/a)') 'dpnint: interpolation ERROR - out of range'
         stop
      end if

! interval halving search for xx(ii) points bracketing tt(jj)
      if(jj>2) then
         if(tt(jj)<tt(jj-1)) imin=1
      end if
      imin = 1
      imax = nn
      do kk = 1, nn
         ii = (imin + imax) / 2
         if(tt(jj)>xx(ii)) then
            imin = ii
         else
            imax = ii
         end if
         if(imax - imin == 1) then
            exit
         end if
      end do


      zz=tt(jj)

      if(mod(npoly,2)==1) then
         istart=imin-int(real(npoly,dp)/2)
      else if(zz-xx(imin) < xx(imax)-zz) then
         istart=imin-int(real(npoly,dp)/2)
      else
         istart=imax-int(real(npoly,dp)/2)
      end if

      istart = min(istart, nn - npoly)
      istart = max(istart, 1)

      sum=0.0d0
      do iy=istart,istart+npoly
         if(yy(iy)==0.0d0) cycle
         term=yy(iy)
         do iprod=istart, istart+npoly
            if(iprod==iy) cycle
            term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
         end do
         sum=sum+term
      end do
      ss(jj)=sum

! special treatment for the origin
! Do order npoly EXTRAPOLATION to the origin using the points inerpolated
! on the linear grid rather than the log grid, since this represents an
! extrapolation of 1 grid point.  Exrapolation from the innermost points
! of the log grid would represent a huge extrapolation, and round-off
! error would not be acceptable If the fitted function is a polynomial
! of order npoly or less, this is exact.

      if(tt(1)==0.0d0) then

         istart=2
         zz=0.0d0

         sum=0.0d0
         do iy=istart,istart+npoly
            if(ss(iy)==0.0d0) cycle
            term=ss(iy)
            do iprod=istart, istart+npoly
               if(iprod==iy) cycle
               term=term*(zz-tt(iprod))/(tt(iy)-tt(iprod))
            end do
            sum=sum+term
         end do
         ss(1)=sum

      end if

   end do
   return
end subroutine dpnint
