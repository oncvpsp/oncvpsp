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
subroutine check_data(atsym,zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl, &
&                      fa,facnf, &
&                      rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                      ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

!check input data for various violations of size or sign
!see annotated input document for all definitions
!*cnf variables are arrays of the basis atomic configuration variables
!to be used for the test configurations

   implicit none
   integer, parameter :: dp=kind(1.0d0)

! Input variables
   character(len=2) :: atsym
   character(len=4) :: psfile

   real(dp) :: zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl
   real(dp) :: fa(30),facnf(30,5),rc(6),ep(6),qcut(6),debl(5)

   integer :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf
   integer :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
   integer :: ncon(6),nbas(6),nproj(6)

! Local variables
   integer :: ierr,ii,jj,l1

   logical :: occ

   real(dp) :: sf,rcmax

   ierr=0

   if(trim(atsym)=='') then
      write(6,'(a)') 'test_data: must have non-blank atsym'
      ierr=ierr+1
   end if

   if(trim(psfile)/='psp8' .and. trim(psfile)/='upf' &
   &   .and. trim(psfile)/='both') then
      write(6,'(a)') 'test_data: psfile must == psp8 or upf or both'
      ierr=ierr+1
   end if

   if(zz<=0.0d0) then
      write(6,'(a)') 'test_data: must have positive atomic number z'
      ierr=ierr+1
   end if

   if(nc+nv<=0 .or. nc+nv>30) then
      write(6,'(a)') 'test_data: must have 1 <= nc+nv <=30'
      ierr=ierr+1
   end if

   if(iexc==0 .or. iexc>4) then
      write(6,'(a)') 'test_data: must have iexc <=4, iexc /= 0.'
      write(6,'(a)') '  negative values for libxc'
      ierr=ierr+1
   end if

   sf=0.0d0
   do ii=1,nc+nv
      sf=sf+fa(ii)
      if(la(ii)<0 .or. la(ii)>3) then
         write(6,'(a,a,i4)') 'test_data: must have 0<= l <=3,', &
         &   ' reference configuration line',ii
         ierr=ierr+1
      end if
      if(na(ii)<=la(ii)) then
         write(6,'(a,i4)') 'test_data: l > n, reference configuration line',ii
         ierr=ierr+1
      end if
      if(fa(ii)<=0.0d0) then
         write(6,'(a,i4)') 'test_data: f <= 0.0, reference configuration line',ii
         ierr=ierr+1
      end if
   end do

   if(sf>zz) then
      write(6,'(a)') 'test_data: reference configuration is negative ion'
      ierr=ierr+1
   end if

   if(lmax<0 .or. lmax>3) then
      write(6,'(a)') 'test_data: 0 <= lmax <= 3'
      ierr=ierr+1
   end if

   do l1=1,lmax+1
      if(rc(l1)<=0.0d0) then
         write(6,'(a,i4)') 'test_data: must have rc>0.0, l=',l1-1
         ierr=ierr+1
      end if
      if(ep(l1)<=0.0d0) then
         occ=.false.
         do ii=nc+1,nc+nv
            if(la(ii)==l1-1) then
               occ=.true.
            end if
         end do
!   if(.not. occ) then
!     write(6,'(a,a,i4)') 'test_data: must have ep>0.0 for scattering state', &
!&     ' l=',l1-1
!     ierr=ierr+1
!   end if
      end if

      if(ncon(l1)<3 .or. ncon(l1)>5) then
         write(6,'(a,i4)') 'test_data: must have 3 <= ncon <= 5, l=',l1-1
         ierr=ierr+1
      end if
      if(nproj(l1)==2) then
         if(nbas(l1)<ncon(l1)+3 .or. nbas(l1)>ncon(l1)+5) then
            write(6,'(a,i4/a)') 'test_data: must have ncon+3 <= nbas <= ncon+5, l=',&
            &    l1-1,' for nproj=2'
            ierr=ierr+1
         end if
      else
         if(nbas(l1)<ncon(l1)+2 .or. nbas(l1)>ncon(l1)+5) then
            write(6,'(a,i4)') 'test_data: must have ncon+2 <= nbas <= ncon+5, l=',l1-1
            ierr=ierr+1
         end if
      end if
      if(qcut(l1)<=0.0d0) then
         write(6,'(a,i4)') 'test_data: must have qcut>0.0, l=',l1-1
         ierr=ierr+1
      end if
   end do

   if(lloc==4) then
      if(lpopt<1 .or. lpopt>5) then
         write(6,'(a)') 'test_data: must have 1 <= lpopt <= 5'
         ierr=ierr+1
      end if
      if(rc(5)<=0.0d0) then
         write(6,'(a)') 'test_data: for lloc==4, must have rc > 0.0'
         ierr=ierr+1
      end if
   else if(lloc<0 .or. lloc>lmax) then
      write(6,'(a)') 'test_data: must have 0 <= lloc <= lmax or lloc == 4'
      ierr=ierr+1
   end if

   do l1=1,lmax+1
      if(rc(l1)<rc(lloc+1)) then
         write(6,'(a,i2,a)') 'test_data: rc < rc(lloc) for l =',l1-1,&
         &         '.  Not allowed.'
         ierr=ierr+1
      end if
   end do

   do l1=1,lmax+1
      if(nproj(l1)==0 .and. (lloc==4 .or. lloc==l1-1)) then
         cycle
      else if(nproj(l1)<1 .or. nproj(l1)>5) then
         write(6,'(a,i4,a)') 'test_data: must have nproj in [1,5] (0 OK lloc=4 or', &
         &        l1-1,')'
         ierr=ierr+1
      end if
      if(debl(l1)<0.0d0) then
         write(6,'(a,i4)') 'test_data: must have debl>0.0)',l1-1
         ierr=ierr+1
      end if
   end do

   if(icmod<0 .or. icmod>4) then
      write(6,'(a)') 'test_data: must have 0<= icmod <=4'
      ierr=ierr+1
   end if

   if((icmod==1 .or. icmod==3) .and. fcfact<=0.0d0) then
      write(6,'(a)') 'test_data: must have fcfact>0.0 for icmod= 1 or 3'
      ierr=ierr+1
   end if

   if((icmod==3) .and. rcfact<=0.0d0) then
      write(6,'(a)') 'test_data: must have rcfact>0.0 for icmod= 3'
      ierr=ierr+1
   end if

   if(epsh2<epsh1) then
      write(6,'(a)') 'test_data: must have epsh2 > epsh1'
      ierr=ierr+1
   end if
   if(depsh<0.0d0) then
      write(6,'(a)') 'test_data: must have depsh>0.0'
      ierr=ierr+1
   end if

   rcmax=0.0d0
   do l1=1,lmax+1
      rcmax=max(rcmax,rc(l1))
   end do
   if(lloc==4) then
      rcmax=max(rcmax,rc(5))
   end if

   if(rlmax<rcmax) then
      write(6,'(a)') 'test_data: must have rlmax>rcmax'
      ierr=ierr+1
   end if

   if(drl<=0.0d0) then
      write(6,'(a)') 'test_data: must have drl > 0.0'
      ierr=ierr+1
   end if

   if(ncnf<0 .or. ncnf>4)then
      write(6,'(a)') 'test_data: must have 0 <= ncnf <= 4, skipping config. checks'
      ierr=ierr+1
   else
      do jj=2,ncnf+1
         sf=0.0d0
         do ii=1,nc+nvcnf(jj)
            sf=sf+facnf(ii,jj)
            if(lacnf(ii,jj)<0 .or. lacnf(ii,jj)>3) then
               write(6,'(a,a,i4,a,i4)') 'test_data: must have 0<= l <=3,', &
               &     ' test configuration',jj-1, ' line',ii-nc
               ierr=ierr+1
            end if
            if(nacnf(ii,jj)<=lacnf(ii,jj)) then
               write(6,'(a,a,i4,a,i4)') 'test_data: l > n,', &
               &     ' test configuration',jj-1, ' line',ii-nc
               ierr=ierr+1
            end if
            if(facnf(ii,jj)<0.0d0) then
               write(6,'(a,a,i4,a,i4)') 'test_data: f < 0.0,', &
               &     ' test configuration',jj-1, ' line',ii-nc
               ierr=ierr+1
            end if
         end do

         if(sf>zz) then
            write(6,'(a,i4)') 'test_data: negative ion, test configuration',jj-1
            ierr=ierr+1
         end if
      end do
   end if

   if(ierr>0) then
      write(6,'(a,i4,a)') 'ERROR: test_data found',ierr,' ERROR; stopping'
      stop
   end if

   return
end subroutine check_data
