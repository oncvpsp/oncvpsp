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
 subroutine dp3int(xx, yy, nn, tt, ss, mm)

! local cubic polynomial interpolation of data yy on nn points xx
! xx must be ordered in ascending order
! output mm interpolated values ss on points tt

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) ::  xx(*),yy(*),tt(*),ss(*)
 integer nn,mm

!Local variables
 real(dp) ::  d41,d31,d21,d42,d32,d43,t4,t3,t2,t1
 integer ii,imin,imax,jj,kk

! basis functions

 if(nn<4) then
   write(6,'(/a,i6)') 'dp3int: interpolation error, n=',nn
   stop
 end if

! note: output point 1 is skipped in this version because of special
! properties of pp data

 imin = 1
 do jj = 2, mm
   if(tt(jj)<xx(1)) then
     write(6,'(/a)') 'dp3int: interpolation error - out of range'
     stop
   end if
   if(tt(jj)>xx(nn)) then
     write(6,'(/a)') 'dp3int: interpolation error - out of range'
     stop
   end if
   if(jj>2 .and. tt(jj)<tt(jj - 1)) then
     imin = 1
   end if
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(tt(jj)>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do
   ii = imin

   ii = ii - 1
   ii = min(ii, nn - 3)
   ii = max(ii, 1)

   d41 = xx(ii + 3) - xx(ii)
   d31 = xx(ii + 2) - xx(ii)
   d21 = xx(ii + 1) - xx(ii)
   d42 = xx(ii + 3) - xx(ii + 1)
   d32 = xx(ii + 2) - xx(ii + 1)
   d43 = xx(ii + 3) - xx(ii + 2)

   t4 = tt(jj) - xx(ii + 3)
   t3 = tt(jj) - xx(ii + 2)
   t2 = tt(jj) - xx(ii + 1)
   t1 = tt(jj) - xx(ii)

   ss(jj) =  - yy(ii    ) * (t2 * t3 * t4) / (d41 * d31 * d21) &
&            + yy(ii + 1) * (t1 * t3 * t4) / (d21 * d32 * d42) &
&            - yy(ii + 2) * (t1 * t2 * t4) / (d31 * d32 * d43) &
&            + yy(ii + 3) * (t1 * t2 * t3) / (d41 * d42 * d43)

 end do
 return
 end subroutine dp3int
