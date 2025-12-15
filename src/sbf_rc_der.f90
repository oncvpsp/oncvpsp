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
 subroutine sbf_rc_der(llin,qq,rc,sbfder)

!calculate spherical Bessel function values and first four derivatives

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!INPUT
! llin  angular momentum l
! qq  wave vector
! rc  core radius
!
!OUTPUT
! sbfder  values and first 4 derivatives of j_l(qq*rc)

!Arguments
 integer :: llin
 real(dp) :: qq,rc,xx,sbfder(5)
!
!local variables
 real(dp) :: sb_out(10),sbfad(5,10)
 integer :: ii,ll

 if(llin>5) then
  write(6,'(a,i4,a)') 'sbf_rc_der: argument error, llin = ',llin,&
&  ' >5'
  stop
 end if

! calculate sbf for needed l values
 xx=qq*rc
 call sbf8(llin+5,xx,sb_out)

! derivatives based on derivative formula applied recursively
!
! (1/z d/dz)^m [z^-n j_n(z)]= (-1)^m z^(-n-m) j_(n+m)(z)
!

! store values
 do ll=llin,llin+5
  sbfad(1,ll+1)=sb_out(ll+1)
 end do

! first derivatives
 do ll=llin,llin+3
  sbfad(2,ll+1)=ll*sbfad(1,ll+1)/xx &
&                - sbfad(1,ll+2)
 end do
 
! second derivatives
 do ll=llin,llin+2
  sbfad(3,ll+1)=-ll*sbfad(1,ll+1)/xx**2 &
&              + ll*sbfad(2,ll+1)/xx&
&                  -sbfad(2,ll+2)
 end do

! third derivatives
 do ll=llin,llin+1
  sbfad(4,ll+1)= 2*ll*sbfad(1,ll+1)/xx**3 &
&               -2*ll*sbfad(2,ll+1)/xx**2 &
&                 +ll*sbfad(3,ll+1)/xx &
&                    -sbfad(3,ll+2)
 end do

! fourth derivatives
 do ll=llin,llin
  sbfad(5,ll+1)=-6*ll*sbfad(1,ll+1)/xx**4 &
&               +6*ll*sbfad(2,ll+1)/xx**3 &
&               -3*ll*sbfad(3,ll+1)/xx**2 &
&                 +ll*sbfad(4,ll+1)/xx &
&                    -sbfad(4,ll+2)
 end do

 do ii=1,5
  sbfder(ii)=sbfad(ii,llin+1)*qq**(ii-1)
 end do

 return
 end subroutine sbf_rc_der
