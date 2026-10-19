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
 pure function smootherstep(x)
 !**********************************************
 ! 5th order "smootherstep" function: 0 for x <= 0, 1 for x >= 1, smooth in between
  implicit none
  integer, parameter :: dp=kind(1.0d0)

  real(dp), intent(in) :: x
  real(dp) :: smootherstep
  if (x <= 0.0_dp) then
   smootherstep = 0.0_dp
  else if (x >= 1.0_dp) then
   smootherstep = 1.0_dp
  else
   smootherstep = x**3 * (6.0_dp*x**2 - 15.0_dp*x + 10.0_dp)
  end if
 end function smootherstep
