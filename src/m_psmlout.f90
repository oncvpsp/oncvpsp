!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! Copyright (c) 2015 by Alberto Garcia, ICMAB-CSIC for xmlf90-wxml calls
! in support of the PSML format
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
! Support for PSML file creation
!
module m_psmlout

   public :: psmlout, psmlout_r
   private

   integer, parameter :: dp = selected_real_kind(10,100)
   character(len=1), dimension(0:4) :: lsymb = ['s','p','d','f','g']

CONTAINS

subroutine psmlout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                  irct, srel, &
&                  zz,zion,mmax,iexc,icmod,drl,atsym, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,ep,qcut,debl,facnf,dvloc0,fcfact, &
&                  epsh1,epsh2,depsh,rlmax,psfile)

!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local potential
!  if linear combination is used)
!rho  valence pseudocharge
!rhomod  model core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy
!psfile  should be 'upf'

   ! Alberto Garcia, February 1, 2015

   use xmlf90_wxml     ! To write XML files
   use m_libxc_list  ! For ease of libxc handling

   implicit none

   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
   integer :: lmax,lloc,iexc,mmax,icmod
   integer :: nproj(6)
   integer :: irct  ! index of point at which rho_core is matched
   logical :: srel  ! whether it is scalar-relativistic or not
   real(dp) :: drl,fcfact,zz,zion
   real(dp), target :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,2,4)
   real(dp), target :: rhomod(mmax,5)
   real(dp):: rc(6),evkb(2,4)
   character(len=2) :: atsym

!additional input for upf output to echo input file, all as defined
! in the main progam
   integer :: na(30),la(30),ncon(6),nbas(6)
   integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
   integer :: nc,nv,lpopt,ncnf
   real(dp) :: fa(30),ep(6),qcut(6),debl(6),facnf(30,5)
   real(dp) :: dvloc0,epsh1,epsh2,depsh,rlmax
   character(len=4) :: psfile

!Local variables
   integer :: dtime(8)

!-------
!psml stuff

   type(xmlf_t) :: xf
   type(libxc_t) :: libxc_id(2)

   integer   :: i, j, npts
   integer   :: ii, jj, l1
   real(dp)  :: rcore
   real(dp)  :: total_valence_charge
   real(dp), pointer :: r(:), chval(:), chcore(:)
   real(dp), pointer :: vps(:), vlocal(:)
   real(dp), allocatable :: r0(:), f0(:)

   character(len=100)    :: line
   character(len=4)      :: polattrib, coreattrib
   character(len=2) :: nameat
   character(len=40):: psflavor

   character(len=1), dimension(0:4) :: lsymb = ['s','p','d','f','g']

   character(len=30) :: xcfuntype, xcfunparam
   integer          :: ncore, nval, ncp, norbs, npots

   integer, allocatable  :: n(:), l(:)
   integer, allocatable  :: nn(:), ll(:)
   real(dp), allocatable :: f(:), ff(:)
   real(dp), allocatable :: fdown(:), fup(:)
   integer :: nr
   character(len=10)     :: datestr

   logical :: polarized, there_is_core, found
   integer :: lun, stat

!---

   call date_and_time(VALUES=dtime)
   write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)

   !-------------- beginning of file echoing
   call get_unit(lun)
   open(unit=lun,file="_tmp_input",form="formatted",&
        status="unknown",position="rewind",action="write")

   write(lun,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write(lun,'(a)') '# atsym  z    nc    nv    iexc   psfile'
   write(lun,'(a,a,f6.2,3i6,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
   &      '      ',trim(psfile)
   write(lun,'(a/a)') '#','#   n    l    f        energy (Ha)'
   do ii=1,nc+nv
      write(lun,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
   end do

   write(lun,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
   write(lun,'(i5)')  lmax
   write(lun,'(a/a)') '#','#   l,   rc,     ep,   ncon, nbas, qcut'
   do l1=1,lmax+1
      write(lun,'(i5,2f8.2,2i5,f8.2)') l1-1,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
   end do

   write(lun,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
   &      '   dvloc0'
   write(lun,'(2i5,f10.2,a,f8.2)') lloc,lpopt,rc(5),'   ',dvloc0

   write(lun,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
   &      '# l, nproj, debl'
   do l1=1,lmax+1
      write(lun,'(2i5,f10.4)') l1-1,nproj(l1),debl(l1)
   end do

   write(lun,'(a/a/a)') '#','# MODEL CORE CHARGE', &
   &      '# icmod, fcfact'
   write(lun,'(i5,f8.2,2i5,f8.2)') icmod,fcfact

   write(lun,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
   &      '# epsh1, epsh2, depsh'
   write(lun,'(3f8.2)') epsh1,epsh2,depsh

   write(lun,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
   write(lun,'(2f8.2)') rlmax,drl

   write(lun,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
   write(lun,'(i5)') ncnf
   write(lun,'(a/a)') '# nvcnf','#   n    l    f'
   do jj=2,ncnf+1
      write(lun,'(i5)') nvcnf(jj)
      do ii=nc+1,nc+nvcnf(jj)
         write(lun,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
      end do
      write(lun,'(a)') '#'
   end do
   close(lun)

   call xml_OpenFile("ONCVPSP.psml",xf, indent=.false.)

   call xml_AddXMLDeclaration(xf,"UTF-8")

   call xml_NewElement(xf,"psml")
   call my_add_attribute(xf,"version","0.8")
   call my_add_attribute(xf,"energy_unit","hartree")
   call my_add_attribute(xf,"length_unit","bohr")

   call xml_NewElement(xf,"provenance")
   call my_add_attribute(xf,"creator","ONCVPSP-2.1.2+psml")
   call my_add_attribute(xf,"date",datestr)
   call xml_NewElement(xf,"input-file")
   call my_add_attribute(xf,"name","oncvpsp-input")
!
   call get_unit(lun)
   open(unit=lun,file="_tmp_input",form="formatted",&
        status="old",position="rewind",action="read")

   rewind(lun)
   do
      read(lun,fmt="(a)",iostat=stat) line
      if (stat /= 0) exit
      call xml_AddPcData(xf,trim(line),line_feed=.true.)
   enddo
   close(lun)

   call xml_EndElement(xf,"input-file")
   !
   call xml_EndElement(xf,"provenance")

!  nameat = symbol(nint(zz))
   nameat = atsym
   ncore  = nc
   nval   = nv

   norbs = ncore + nval
   allocate (n(norbs), l(norbs), f(norbs))

   total_valence_charge = 0.0_dp
   ncp = ncore + 1
   do i = 1, norbs
      n(i) = na(i)
      l(i) = la(i)
      f(i) = fa(i)
      if (i > ncore) then
         total_valence_charge =   total_valence_charge + f(i)
      endif
   enddo
   lmax = lmax

   npots = lmax + 1
   allocate (ll(npots), nn(npots), ff(npots))
   do i = 1, npots
      ll(i) = i - 1
      found = .false.
      ! look for the appropriate shell in the valence
      do j = ncp, norbs
         if (l(j) == ll(i)) then
            found = .true.
            nn(i) = n(j)
            ff(i) = f(j)
            exit
         endif
      enddo
      if (.not. found) then
         ! generate the appropriate effective n
         nn(i) = ll(i) + 1
         do j = 1, ncore
            if (l(j) == ll(i)) then
               nn(i) = nn(i) + 1
            endif
            ff(i) = 0.0_dp
         enddo
      endif
   enddo

   psflavor ="Hamann's oncvpsp"

   polarized = .false.
   polattrib = "no"
   there_is_core = (icmod >= 1)
   if (there_is_core) then
      coreattrib = "yes"
   else
      coreattrib = "no"
   endif

   ! XC name handling

   select case(iexc)

    case(1)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Wigner'
      libxc_id = [ XC_LDA_X, XC_LDA_C_WIGNER ]
    case(2)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Hedin-Lundqvist'
      libxc_id = [ XC_LDA_X, XC_LDA_C_HL ]
    case(3)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Ceperley-Alder PZ'
      libxc_id = [ XC_LDA_X, XC_LDA_C_PZ ]

    case(4)
      xcfuntype    = 'GGA'
      xcfunparam   = 'Perdew-Burke-Ernzerhof'
      libxc_id = [ XC_GGA_X_PBE, XC_GGA_C_PBE ]

    case  default

      xcfuntype    = '-----'
      xcfunparam   = '-----'
      libxc_id = [  XC_NOT_IMPL, XC_NOT_IMPL ]   ! ???

   end select
   !
   !
   call xml_NewElement(xf,"header")
   call my_add_attribute(xf,"atomic-label",nameat)
   call my_add_attribute(xf,"atomic-number",str(zz))
   call my_add_attribute(xf,"z-pseudo",str(zion))
   call my_add_attribute(xf,"flavor",psflavor)
   if (srel) then
      call my_add_attribute(xf,"relativity","scalar")
   else
      call my_add_attribute(xf,"relativity","no")
   endif
   call my_add_attribute(xf,"polarized",polattrib)
   call my_add_attribute(xf,"core-corrections",coreattrib)

   call xml_NewElement(xf,"exchange-correlation")
   call xml_NewElement(xf,"annotation")
   call my_add_attribute(xf,"oncvpsp-xc-code",str(iexc))
   call my_add_attribute(xf,"oncvpsp-xc-type",trim(xcfuntype))
   call my_add_attribute(xf,"oncvpsp-xc-authors",trim(xcfunparam))
   call xml_EndElement(xf,"annotation")

   call xml_NewElement(xf,"libxc-info")
   call my_add_attribute(xf,"number-of-functionals","2")
   do i = 1, 2
      call xml_NewElement(xf,"functional")
      call my_add_attribute(xf,"name",trim(libxc_id(i)%name))
      call my_add_attribute(xf,"type",trim(libxc_id(i)%xc_kind%str))
      call my_add_attribute(xf,"id",str(libxc_id(i)%code))
      call xml_EndElement(xf,"functional")
   enddo
   call xml_EndElement(xf,"libxc-info")
   call xml_EndElement(xf,"exchange-correlation")
   !

   !
   call do_configuration()
   call xml_EndElement(xf,"header")

!AG: save
   if(lloc==4) then
      ! fitted local potential
   else
      ! 'l_local="',lloc,'"'
   end if
!AG -- decide how to handle the case of Vlocal as one of the sl pots.


!AG: decide whether to use a single mesh (Hamann's own) or have the
! projectors use another one (linear, shorter)

   npts = mmax

   r => rr(:)
   chval => rho(:)
   chcore => rhomod(:,1)

   nr = npts + 1
   allocate(r0(nr))
   call add_zero_r(r,r,r0)
   r0(1) = 0.0_dp

   call xml_NewElement(xf,"grid")
   call my_add_attribute(xf,"npts",str(nr))

   call xml_NewElement(xf,"annotation")
   call my_add_attribute(xf,"type","oncvpsp plus r=0")
   call my_add_attribute(xf,"oncvpsp-npts",str(npts))
   call my_add_attribute(xf,"oncvpsp-r1",trim(str(r(1))))
   call my_add_attribute(xf,"oncvpsp-factor",trim(str(r(2)/r(1))))
   call xml_EndElement(xf,"annotation")

   call xml_NewElement(xf,"grid-data")
   call xml_AddArray(xf,r0(1:nr))
   call xml_EndElement(xf,"grid-data")

   call xml_EndElement(xf,"grid")

   call xml_NewElement(xf,"semilocal-potentials")
   if (srel) then
      call my_add_attribute(xf,"set","scalar_relativistic")
   else
      call my_add_attribute(xf,"set","non_relativistic")
   endif
   !
   allocate(f0(nr))
   vpsd: do i = 1, npots
      vps => vpuns(:,i)
      call add_zero_r(vps(1:npts),r,f0)
      call write_psml_item(xf, class="slps", &
                           n=nn(i), l=ll(i), &
                           rc=rc(i),flavor=psflavor,&
                           f=f0)
   enddo vpsd
   call xml_EndElement(xf,"semilocal-potentials")
!
!--------
   call xml_NewElement(xf,"valence-charge")
   call my_add_attribute(xf,"total-charge",  &
                         str(total_valence_charge))
   call xml_NewElement(xf,"radfunc")

   call xml_NewElement(xf,"data")
   call add_zero_r(chval(1:npts),r,f0)
!  call xml_AddArray(xf,chval(1:npts))
   call xml_AddArray(xf,f0(1:nr))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"valence-charge")


   if (there_is_core) then
      rcore = rr(irct)
      call xml_NewElement(xf,"pseudocore-charge")
      call my_add_attribute(xf,"matching-radius",str(rcore))
      call my_add_attribute(xf,"number-of-continuous-derivatives", &
                            str(4))
      call my_add_attribute(xf,"annotation",  &
                            "Monotonic 8th-order polynomial with no linear term")
      call xml_NewElement(xf,"radfunc")

      call xml_NewElement(xf,"data")
      call add_zero_r(chcore(1:npts),r,f0)
      call xml_AddArray(xf,f0(1:nr))
      call xml_EndElement(xf,"data")
      call xml_EndElement(xf,"radfunc")
      call xml_EndElement(xf,"pseudocore-charge")
      deallocate(chcore)
   endif

   call xml_NewElement(xf,"pseudopotential-operator")
   call my_add_attribute(xf,"version","0.1")
   call my_add_attribute(xf,"energy_unit","hartree")
   call my_add_attribute(xf,"length_unit","bohr")

   call xml_NewElement(xf,"provenance")
   call my_add_attribute(xf,"creator","oncvpsp 2.1.2+psml")
   call xml_EndElement(xf,"provenance")

   vlocal => vpuns(:,lloc+1)
   call xml_NewElement(xf,"local-potential")
   if (lloc > lmax) then
      call my_add_attribute(xf,"type","oncv-fit")
   else
      call my_add_attribute(xf,"type","l="//str(lloc))
   endif
   call xml_NewElement(xf,"radfunc")
   call xml_NewElement(xf,"data")
   call add_zero_r(vlocal(1:npts),r,f0)
   call xml_AddArray(xf,f0(1:nr))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"local-potential")

   call xml_NewElement(xf,"projectors")
   if (srel) then
      call my_add_attribute(xf,"set","scalar_relativistic")
   else
      call my_add_attribute(xf,"set","non_relativistic")
   endif

   do l1=1,lmax+1
      if(l1==lloc+1) cycle
      do jj=1,nproj(l1)
         call add_zero_r(vkb(:,jj,l1),r,f0)
         call write_psml_item(xf, class="proj", &
                              seq=jj, l=l1-1, &
                              ekb=evkb(jj,l1), &
                              type="oncv", f=f0)
      enddo
   enddo

   call xml_EndElement(xf,"projectors")

   call xml_EndElement(xf,"pseudopotential-operator")

   call xml_EndElement(xf,"psml")


   call xml_Close(xf)

   deallocate(f0,r0)

CONTAINS
subroutine do_configuration()

   call xml_NewElement(xf,"valence-configuration")
   call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
   do i = ncp, norbs
      if (f(i) < 1.0e-10_dp) cycle
      call xml_NewElement(xf,"shell")
      call my_add_attribute(xf,"n",str(n(i)))
      call my_add_attribute(xf,"l",lsymb(l(i)))
      call my_add_attribute(xf,"occupation",str(f(i)))
      if (polarized) then
         call my_add_attribute(xf,"occupation-down",str(fdown(i)))
         call my_add_attribute(xf,"occupation-up",str(fup(i)))
      endif
      call xml_EndElement(xf,"shell")
   enddo
   call xml_EndElement(xf,"valence-configuration")
end subroutine do_configuration

end subroutine psmlout
!
!=================================================
! Fully relativistic version
!
subroutine psmlout_r(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                  irct, &
&                  vsr,esr,vso,eso, &
&                  zz,zion,mmax,iexc,icmod,drl,atsym, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,ep,qcut,debl,facnf,dvloc0,fcfact, &
&                  epsh1,epsh2,depsh,rlmax,psfile)

!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local potential
!  if linear combination is used)
!rho  valence pseudocharge
!rhomod  model core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy
!psfile  'upf' or 'psp8'

   ! Alberto Garcia, February 1, 2015

   use xmlf90_wxml     ! To write XML files
   use m_libxc_list  ! For ease of libxc handling

   implicit none

   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
   integer :: lmax,lloc,iexc,mmax,icmod
   integer :: nproj(6)
   integer :: irct  ! index of point at which rho_core is matched
   real(dp) :: drl,fcfact,zz,zion
   real(dp), target :: rr(mmax),vpuns(mmax,5,2),rho(mmax),vkb(mmax,2,4,2)
   real(dp), target :: rhomod(mmax,5)
   real(dp):: rc(6),evkb(2,4,2)
   real(dp), target :: vsr(mmax,4,4),vso(mmax,4,4)
   real(dp) :: esr(4,4),eso(4,4)
   character(len=2) :: atsym

!additional input for upf output to echo input file, all as defined
! in the main progam
   integer :: na(30),la(30),ncon(6),nbas(6)
   integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
   integer :: nc,nv,lpopt,ncnf
   real(dp) :: fa(30),ep(6),qcut(6),debl(6),facnf(30,5)
   real(dp) :: dvloc0,epsh1,epsh2,depsh,rlmax
   character(len=4) :: psfile


!Local variables
   integer :: dtime(8)
   integer :: npr_so(4), npr_sr(4)

   type(xmlf_t) :: xf
   type(libxc_t) :: libxc_id(2)

   integer   :: i, j, npts
   integer   :: ii, jj, l1, jk
   real(dp)  :: rcore, jval
   real(dp)  :: total_valence_charge
   real(dp), pointer :: r(:), chval(:), chcore(:)
   real(dp), pointer :: vlocal(:)
   real(dp), allocatable :: r0(:), f0(:), vps(:)

   character(len=100)    :: line
   character(len=4)      :: polattrib, coreattrib
   character(len=2) :: nameat
   character(len=40):: psflavor

   character(len=1), dimension(0:4) :: lsymb = ['s','p','d','f','g']

   character(len=30) :: xcfuntype, xcfunparam
   integer          :: ncore, nval, ncp, norbs, npots

   integer, allocatable  :: n(:), l(:)
   integer, allocatable  :: nn(:), ll(:)
   real(dp), allocatable :: f(:), ff(:)
   real(dp), allocatable :: fdown(:), fup(:)
   integer :: nr
   character(len=10)     :: datestr

   logical :: polarized, there_is_core, found
   integer :: lun, stat

!---

   call date_and_time(VALUES=dtime)
   write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)

   !-------------- beginning of file echoing
   call get_unit(lun)
   open(unit=lun,file="_tmp_input",form="formatted",&
        status="unknown",position="rewind",action="write")

   write(lun,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write(lun,'(a)') '# atsym  z    nc    nv    iexc   psfile'
   write(lun,'(a,a,f6.2,3i6,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
   &      '      ',trim(psfile)
   write(lun,'(a/a)') '#','#   n    l    f        energy (Ha)'
   do ii=1,nc+nv
      write(lun,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
   end do

   write(lun,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
   write(lun,'(i5)')  lmax
   write(lun,'(a/a)') '#','#   l,   rc,     ep,   ncon, nbas, qcut'
   do l1=1,lmax+1
      write(lun,'(i5,2f8.2,2i5,f8.2)') l1-1,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
   end do

   write(lun,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
   &      '   dvloc0'
   write(lun,'(2i5,f10.2,a,f8.2)') lloc,lpopt,rc(5),'   ',dvloc0

   write(lun,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
   &      '# l, nproj, debl'
   do l1=1,lmax+1
      write(lun,'(2i5,f10.4)') l1-1,nproj(l1),debl(l1)
   end do

   write(lun,'(a/a/a)') '#','# MODEL CORE CHARGE', &
   &      '# icmod, fcfact'
   write(lun,'(i5,f8.2,2i5,f8.2)') icmod,fcfact

   write(lun,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
   &      '# epsh1, epsh2, depsh'
   write(lun,'(3f8.2)') epsh1,epsh2,depsh

   write(lun,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
   write(lun,'(2f8.2)') rlmax,drl

   write(lun,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
   write(lun,'(i5)') ncnf
   write(lun,'(a/a)') '# nvcnf','#   n    l    f'
   do jj=2,ncnf+1
      write(lun,'(i5)') nvcnf(jj)
      do ii=nc+1,nc+nvcnf(jj)
         write(lun,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
      end do
      write(lun,'(a)') '#'
   end do
   close(lun)


   call xml_OpenFile("ONCVPSP.psml",xf, indent=.false.)

   call xml_AddXMLDeclaration(xf,"UTF-8")

   call xml_NewElement(xf,"psml")
   call my_add_attribute(xf,"version","0.8")
   call my_add_attribute(xf,"energy_unit","hartree")
   call my_add_attribute(xf,"length_unit","bohr")

   call xml_NewElement(xf,"provenance")
   call my_add_attribute(xf,"creator","ONCVPSP-2.1.2+psml")
   call my_add_attribute(xf,"date",datestr)
   call xml_NewElement(xf,"input-file")
   call my_add_attribute(xf,"name","oncvpsp-input")
!
   call get_unit(lun)
   open(unit=lun,file="_tmp_input",form="formatted",&
        status="old",position="rewind",action="read")

   rewind(lun)
   do
      read(lun,fmt="(a)",iostat=stat) line
      if (stat /= 0) exit
      call xml_AddPcData(xf,trim(line),line_feed=.true.)
   enddo
   close(lun)

   call xml_EndElement(xf,"input-file")
   !
   call xml_EndElement(xf,"provenance")

!  nameat = symbol(nint(zz))
   nameat = atsym
   ncore  = nc
   nval   = nv

   norbs = ncore + nval
   allocate (n(norbs), l(norbs), f(norbs))

   total_valence_charge = 0.0_dp
   ncp = ncore + 1
   do i = 1, norbs
      n(i) = na(i)
      l(i) = la(i)
      f(i) = fa(i)
      if (i > ncore) then
         total_valence_charge =   total_valence_charge + f(i)
      endif
   enddo
   lmax = lmax

   npots = lmax + 1
   allocate (ll(npots), nn(npots), ff(npots))
   do i = 1, npots
      ll(i) = i - 1
      found = .false.
      ! look for the appropriate shell in the valence
      do j = ncp, norbs
         if (l(j) == ll(i)) then
            found = .true.
            nn(i) = n(j)
            ff(i) = f(j)
            exit
         endif
      enddo
      if (.not. found) then
         ! generate the appropriate effective n
         nn(i) = ll(i) + 1
         do j = 1, ncore
            if (l(j) == ll(i)) then
               nn(i) = nn(i) + 1
            endif
            ff(i) = 0.0_dp
         enddo
      endif
   enddo

   psflavor ="Hamann's oncvpsp"

   polarized = .false.
   polattrib = "no"
   there_is_core = (icmod >= 1)
   if (there_is_core) then
      coreattrib = "yes"
   else
      coreattrib = "no"
   endif

   ! XC name handling

   select case(iexc)

    case(1)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Wigner'
      libxc_id = [ XC_LDA_X, XC_LDA_C_WIGNER ]
    case(2)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Hedin-Lundqvist'
      libxc_id = [ XC_LDA_X, XC_LDA_C_HL ]
    case(3)
      xcfuntype    = 'LDA'
      xcfunparam   = 'Ceperley-Alder PZ'
      libxc_id = [ XC_LDA_X, XC_LDA_C_PZ ]

    case(4)
      xcfuntype    = 'GGA'
      xcfunparam   = 'Perdew-Burke-Ernzerhof'
      libxc_id = [ XC_GGA_X_PBE, XC_GGA_C_PBE ]

    case  default

      xcfuntype    = '-----'
      xcfunparam   = '-----'
      libxc_id = [  XC_NOT_IMPL, XC_NOT_IMPL ]   ! ???

   end select
   !
   !
   call xml_NewElement(xf,"header")
   call my_add_attribute(xf,"atomic-label",nameat)
   call my_add_attribute(xf,"atomic-number",str(zz))
   call my_add_attribute(xf,"z-pseudo",str(zion))
   call my_add_attribute(xf,"flavor",psflavor)
   call my_add_attribute(xf,"relativity","dirac")
   call my_add_attribute(xf,"polarized",polattrib)
   call my_add_attribute(xf,"core-corrections",coreattrib)

   call xml_NewElement(xf,"exchange-correlation")
   call xml_NewElement(xf,"annotation")
   call my_add_attribute(xf,"oncvpsp-xc-code",str(iexc))
   call my_add_attribute(xf,"oncvpsp-xc-type",trim(xcfuntype))
   call my_add_attribute(xf,"oncvpsp-xc-authors",trim(xcfunparam))
   call xml_EndElement(xf,"annotation")

   call xml_NewElement(xf,"libxc-info")
   call my_add_attribute(xf,"number-of-functionals","2")
   do i = 1, 2
      call xml_NewElement(xf,"functional")
      call my_add_attribute(xf,"name",trim(libxc_id(i)%name))
      call my_add_attribute(xf,"type",trim(libxc_id(i)%xc_kind%str))
      call my_add_attribute(xf,"id",str(libxc_id(i)%code))
      call xml_EndElement(xf,"functional")
   enddo
   call xml_EndElement(xf,"libxc-info")
   call xml_EndElement(xf,"exchange-correlation")
   !

   !
   call do_configuration()
   call xml_EndElement(xf,"header")

!AG: save
   if(lloc==4) then
      ! fitted local potential
   else
      ! 'l_local="',lloc,'"'
   end if
!AG -- decide how to handle the case of Vlocal as one of the sl pots.


!AG: decide whether to use a single mesh (Hamann's own) or have the
! projectors use another one (linear, shorter)

   npts = mmax

   r => rr(:)
   chval => rho(:)
   chcore => rhomod(:,1)
   allocate(vps(npts))

   nr = npts + 1
   allocate(r0(nr))
   call add_zero_r(r,r,r0)
   r0(1) = 0.0_dp

   call xml_NewElement(xf,"grid")
   call my_add_attribute(xf,"npts",str(nr))

   call xml_NewElement(xf,"annotation")
   call my_add_attribute(xf,"type","oncvpsp plus r=0")
   call my_add_attribute(xf,"oncvpsp-npts",str(npts))
   call my_add_attribute(xf,"oncvpsp-r1",trim(str(r(1))))
   call my_add_attribute(xf,"oncvpsp-factor",trim(str(r(2)/r(1))))
   call xml_EndElement(xf,"annotation")

   call xml_NewElement(xf,"grid-data")
   call xml_AddArray(xf,r0(1:nr))
   call xml_EndElement(xf,"grid-data")

   call xml_EndElement(xf,"grid")

   !
   ! Semilocal potentials
   !
   allocate(f0(nr))

   if (trim(psfile)=="psp8") then

      call xml_NewElement(xf,"semilocal-potentials")
      call my_add_attribute(xf,"set","scalar_relativistic")
! sr components
      do i = 1, npots
         l1 = i
         ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
         vps(:) = ((ll(i)+1)*vpuns(:,l1,1)+ ll(i)*vpuns(:,l1,2)) / dble(2*ll(i)+1)
         call add_zero_r(vps(1:npts),r,f0)
         call write_psml_item(xf, class="slps", &
                              n=nn(i), l=ll(i), &
                              rc=rc(i),flavor=psflavor,&
                              f=f0)
      enddo
      call xml_EndElement(xf,"semilocal-potentials")
!
! so components
!
      call xml_NewElement(xf,"semilocal-potentials")
      call my_add_attribute(xf,"set","spin_orbit")

      do i = 2, npots
         l1 = i
         vps(:) = 2*(vpuns(:,l1,1) - vpuns(:,l1,2)) / dble(2*ll(i)+1)
         call add_zero_r(vps(1:npts),r,f0)
         call write_psml_item(xf, class="slps", &
                              n=nn(i), l=ll(i), &
                              rc=rc(i),flavor=psflavor,&
                              f=f0)
      enddo

      call xml_EndElement(xf,"semilocal-potentials")

   else   ! upf

      call xml_NewElement(xf,"semilocal-potentials")
      call my_add_attribute(xf,"set","lj")

      do i = 1, npots
         l1 = i
         if (i == 1) then
            jval = 0.0
            ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
            vps(:) = vpuns(:,l1,1)
            call add_zero_r(vps(1:npts),r,f0)
            call write_psml_item(xf, class="slps", &
                                 n=nn(i), l=ll(i), j=jval, &
                                 rc=rc(i),flavor=psflavor,&
                                 f=f0)
         else
            do jj=1,2
               jval = ll(i) + (3-2*jj)*0.5  ! convert (1,2) to (1,-1)*1/2
               ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
               vps(:) = vpuns(:,l1,jj)
               call add_zero_r(vps(1:npts),r,f0)
               call write_psml_item(xf, class="slps", &
                                    n=nn(i), l=ll(i), j=jval, &
                                    rc=rc(i),flavor=psflavor,&
                                    f=f0)
            enddo
         endif
      enddo
      call xml_EndElement(xf,"semilocal-potentials")

   endif  ! upf vs psp8
!
!--------
   call xml_NewElement(xf,"valence-charge")
   call my_add_attribute(xf,"total-charge",  &
                         str(total_valence_charge))
   call xml_NewElement(xf,"radfunc")

   call xml_NewElement(xf,"data")
   call add_zero_r(chval(1:npts),r,f0)
!  call xml_AddArray(xf,chval(1:npts))
   call xml_AddArray(xf,f0(1:nr))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"valence-charge")


   if (there_is_core) then
      rcore = rr(irct)
      call xml_NewElement(xf,"pseudocore-charge")
      call my_add_attribute(xf,"matching-radius",str(rcore))
      call my_add_attribute(xf,"number-of-continuous-derivatives", &
                            str(4))
      call my_add_attribute(xf,"annotation",  &
                            "Monotonic 8th-order polynomial with no linear term")
      call xml_NewElement(xf,"radfunc")

      call xml_NewElement(xf,"data")
      call add_zero_r(chcore(1:npts),r,f0)
      call xml_AddArray(xf,f0(1:nr))
      call xml_EndElement(xf,"data")
      call xml_EndElement(xf,"radfunc")
      call xml_EndElement(xf,"pseudocore-charge")
      deallocate(chcore)
   endif


   call xml_NewElement(xf,"pseudopotential-operator")
   call my_add_attribute(xf,"version","0.1")
   call my_add_attribute(xf,"energy_unit","hartree")
   call my_add_attribute(xf,"length_unit","bohr")

   call xml_NewElement(xf,"provenance")
   call my_add_attribute(xf,"creator","oncvpsp 2.1.2+psml")
   call xml_EndElement(xf,"provenance")

   vlocal => vpuns(:,lloc+1,1)
   call xml_NewElement(xf,"local-potential")
   if (lloc > lmax) then
      call my_add_attribute(xf,"type","oncv-fit")
   else
      call my_add_attribute(xf,"type","l="//str(lloc))
   endif
   call xml_NewElement(xf,"radfunc")
   call xml_NewElement(xf,"data")
   call add_zero_r(vlocal(1:npts),r,f0)
   call xml_AddArray(xf,f0(1:nr))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"local-potential")

!
!     Scalar-relativistic projectors
!
   if (trim(psfile)=="psp8") then

! set up projector number for sr_so calculations based on non-zero coefficients
      npr_sr(:)=0
      npr_so(:)=0
      do l1=1,lmax+1
         do ii=1,4
            if(abs(esr(ii,l1))>0.0d0) npr_sr(l1)=npr_sr(l1)+1
            if(abs(eso(ii,l1))>0.0d0) npr_so(l1)=npr_so(l1)+1
         end do
      end do

      call xml_NewElement(xf,"projectors")
      call my_add_attribute(xf,"set","scalar_relativistic")

      do l1=1,lmax+1
         if(l1==lloc+1) cycle
         do jj=1,npr_sr(l1)
            call add_zero_r(vsr(:,jj,l1),r,f0)
            call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, &
                                 ekb=esr(jj,l1), &
                                 type="oncv", f=f0)
         enddo
      enddo
      call xml_EndElement(xf,"projectors")
!
!     Spin-orbit part
!
      call xml_NewElement(xf,"projectors")
      call my_add_attribute(xf,"set","spin_orbit")

      do l1=2,lmax+1
         if(l1==lloc+1) cycle
         do jj=1,npr_so(l1)
            call add_zero_r(vso(:,jj,l1),r,f0)
            call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, &
                                 ekb=eso(jj,l1), &
                                 type="oncv", f=f0)
         enddo
      enddo

      call xml_EndElement(xf,"projectors")

   else  ! upf

      call xml_NewElement(xf,"projectors")
      call my_add_attribute(xf,"set","lj")

      do l1=1,lmax+1
         if(l1==lloc+1) cycle

         if (l1 == 1) then
            ! l=0, only one j=0
            jval = 0.0

            do jj=1,nproj(l1)
               call add_zero_r(vkb(:,jj,l1,1),r,f0)
               call write_psml_item(xf, class="proj", &
                                    seq=jj, l=l1-1, j=jval, &
                                    ekb=evkb(jj,l1,1), &
                                    type="oncv", f=f0)
            enddo

         else  ! l1 /=1  (l/=0)

            ! two j values
            do jk=1,2
               jval = l1-1 + (3-2*jk)*0.5  ! convert (1,2) to (1,-1)*1/2
               do jj=1,nproj(l1)
                  call add_zero_r(vkb(:,jj,l1,jk),r,f0)
                  call write_psml_item(xf, class="proj", &
                                       seq=jj, l=l1-1, j=jval, &
                                       ekb=evkb(jj,l1,jk), &
                                       type="oncv", f=f0)
               enddo
            enddo

         endif  ! l1 == 1

      enddo  ! over l shells
      call xml_EndElement(xf,"projectors")
!
   endif  ! upf vs psp8

   call xml_EndElement(xf,"pseudopotential-operator")

   call xml_EndElement(xf,"psml")


   call xml_Close(xf)

   deallocate(f0,r0)

CONTAINS
subroutine do_configuration()

   call xml_NewElement(xf,"valence-configuration")
   call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
   do i = ncp, norbs
      if (f(i) < 1.0e-10_dp) cycle
      call xml_NewElement(xf,"shell")
      call my_add_attribute(xf,"n",str(n(i)))
      call my_add_attribute(xf,"l",lsymb(l(i)))
      call my_add_attribute(xf,"occupation",str(f(i)))
      if (polarized) then
         call my_add_attribute(xf,"occupation-down",str(fdown(i)))
         call my_add_attribute(xf,"occupation-up",str(fup(i)))
      endif
      call xml_EndElement(xf,"shell")
   enddo
   call xml_EndElement(xf,"valence-configuration")
end subroutine do_configuration


end subroutine psmlout_r

subroutine my_add_attribute(xf,name,value)
   use xmlf90_wxml, only: xmlf_t, xml_AddAttribute

   type(xmlf_t), intent(inout)   :: xf
   character(len=*), intent(in)  :: name
   character(len=*), intent(in)  :: value

   call xml_AddAttribute(xf,name,trim(value))
end subroutine my_add_attribute

subroutine add_zero_r(f,r,f0)
   ! Adds an r=0 element to a grid function f0, extrapolating

   double precision, intent(in)  :: f(:), r(:)
   double precision, intent(out) :: f0(:)

   integer :: i, npts
   double precision :: r2

   npts = size(f)
   if (size(f0) /= npts +1) stop "nr /= npts + 1 in add_zero_r"
   do i = 1, npts
      f0(i+1) = f(i)
   enddo
   r2 = r(1)/(r(2)-r(1))
   f0(1) = f(2) - (f(3)-f(2))*r2

end subroutine add_zero_r

subroutine get_unit(lun)

!     Get an available Fortran unit number

   integer, intent(out) ::  lun

   integer :: i
   logical :: unit_used

   do i = 10, 99
      lun = i
      inquire(lun,opened=unit_used)
      if (.not. unit_used) return
   enddo
   stop 'NO LUNS'
end subroutine get_unit

subroutine write_psml_item(xf,class, &
                           n, l, j, s, &
                           seq, &
                           rc, ekb, &
                           flavor, type, set, &
                           f)

   use xmlf90_wxml


   type(xmlf_t), intent(inout)   :: xf
   character(len=*), intent(in)  :: class

   integer, intent(in), optional  :: n
   integer, intent(in), optional  :: l
   real(dp), intent(in), optional  :: j
   real(dp), intent(in), optional  :: s

   ! for sl potentials
   real(dp), intent(in), optional  :: rc
   character(len=*), intent(in), optional  :: flavor

   ! for projectors
   integer, intent(in), optional   :: seq
   real(dp), intent(in), optional  :: ekb
   character(len=*), intent(in), optional  :: type

   character(len=*), intent(in), optional  :: set
   real(dp), intent(in), optional  :: f(:)

   call xml_NewElement(xf,trim(class))
   if (present(set))  call my_add_attribute(xf,"set",set)

   ! we might want to check input values
   if (present(n))  call my_add_attribute(xf,"n",str(n))
   if (present(l))  call my_add_attribute(xf,"l",lsymb(l))
   if (present(j))  call my_add_attribute(xf,"j", &
                                          str(j,format="(f3.1)"))
   ! spin: +0.5 or -0.5
   if (present(s))  call my_add_attribute(xf,"s", &
                                          str(s,format="(f4.1)"))

   if (present(seq)) call my_add_attribute(xf,"seq",str(seq))

   if (present(rc))  call my_add_attribute(xf,"rc",str(rc))
   if (present(ekb))  call my_add_attribute(xf,"ekb",str(ekb))

   if (present(flavor))  call my_add_attribute(xf,"flavor",flavor)
   if (present(type))  call my_add_attribute(xf,"type",type)

   call xml_NewElement(xf,"radfunc")
   call xml_NewElement(xf,"data")
   call xml_AddArray(xf,f(:))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")

   call xml_EndElement(xf,trim(class))

end subroutine write_psml_item

end module m_psmlout
