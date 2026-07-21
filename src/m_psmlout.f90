!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! Copyright (c) 2015-2017 by Alberto Garcia, ICMAB-CSIC, for PSML output
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

  implicit none
  
  public :: psmlout, psmlout_r, copy_input_file_for_psml
  public :: get_ec_hints
  private

  integer, parameter :: dp = selected_real_kind(10,100)

  
  character(len=1), dimension(0:4) :: lsymb = (/'s','p','d','f','g'/)
  
  character(len=*), parameter :: PSML_VERSION = "1.1"
  character(len=*), parameter :: PSML_NAMESPACE = "http://esl.cecam.org/PSML/ns/1.1"
  character(len=*), parameter, public :: PSML_PATCH = "psml-4.0.1-81"
  character(len=*), parameter, public :: PSML_CREATOR = "ONCVPSP-4.0.1" // "+" // PSML_PATCH
  character(len=*), parameter :: PSML_FILENAME = "ONCVPSPPSML"

  integer, parameter  :: POLY_ORDER_EXTRAPOL = 7  ! For extrapolation at r=0

  logical :: check_interp
  
  ! Support for cutoff hints -----
  type hint_t
     real(dp) :: err_level
     integer  :: cutoff_hint
  end type hint_t
  ! ------------------------------
  type(hint_t), dimension(7) :: ec_hint

  CONTAINS

 subroutine psmlout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                  irct, srel, &
&                  zz,zion,mmax,mxprj,iexc,icmod,nrl0,drl,atsym,epstot, &
&                  na,la,nc,nv, &
&                  fa,epa, &
&                  psfile,uupsa,ea)

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
!zion  at this point, total valence charge (becomes pseudoion charge)
!mmax  size of log radial grid
!mxprj  dimension of number of projectors
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl0 suggested size of linear radial grid in input file
!drl spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy
!psfile  should be 'upf' or 'psp8' or 'both'
!uupsa  pseudo-atomic orbital array
!ea: energy levels
!epa: actual reference energies for projectors   

  ! Alberto Garcia, February 1, 2015

  use xmlf90_wxml     ! To write XML files
  use m_uuid        ! To generate uuid

 implicit none

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
 integer :: lmax,lloc,iexc,mmax,mxprj,nrl0,icmod
 integer :: nproj(6)
 integer :: irct ! index of point at which rho_core is matched
 logical :: srel ! whether it is scalar-relativistic or not
 real(dp) :: drl,zz,zion,epstot
 real(dp), target :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,mxprj,4)
 real(dp), target :: rhomod(mmax,5)
 real(dp):: rc(6),evkb(mxprj,4)
 character*2 :: atsym
 real(dp) :: uupsa(mmax,nv)
 character(len=40) :: fname

!additional input for upf output to echo input file, all as defined
! in the main program
 integer :: na(30),la(30)
 integer :: nc,nv
 real(dp) :: fa(30),epa(mxprj,6),ea(30)
 character*4 :: psfile

!Local variables
 integer :: dtime(8), nrl, nrl_short, nrl_full
!-------
!psml stuff

  type(xmlf_t) :: xf

  integer   :: i, j, npts
  integer   :: ii, jj, l1
  real(dp)  :: rcore
  real(dp)  :: total_valence_charge
  real(dp), pointer :: r(:), chval(:), chcore(:)
  real(dp), pointer :: vps(:), vlocal(:)
  real(dp), allocatable :: r0(:), f0(:), div_by_r(:)
  character(len=100)    :: line
  character*4      :: polattrib, coreattrib
  character*1      :: pscode, char_dummy
  character(len=2) :: nameat
  character(len=40):: psflavor

  character(len=1), dimension(0:4) :: lsymb = (/'s','p','d','f','g'/)

  integer          :: ncore, nval, ncp, norbs, npots
  integer          :: n1
  real(dp)         :: energy_level, eref
  
  integer, allocatable  :: n(:), l(:)
  integer, allocatable  :: nn(:), ll(:)
  real(dp), allocatable :: f(:), ff(:)
  real(dp), allocatable :: fdown(:), fup(:)

  real(dp) :: rmax, delta, al
  integer, allocatable  :: isample(:)
  
  character(len=10)     :: datestr

  logical :: tdopsp, nonrel, polarized, there_is_core, found
  integer :: lun, stat

  logical :: write_wfns
  
  external :: dpnint   ! Resampler + extrapolator to r=0

 !---
  
  call get_psml_options(write_wfns,check_interp)
  
  call date_and_time(VALUES=dtime)
  write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)

  call xml_OpenFile(PSML_FILENAME,xf, indent=.false.)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"psml")
  call my_add_attribute(xf,"version",PSML_VERSION)
  call my_add_attribute(xf,"energy_unit","hartree")
  call my_add_attribute(xf,"length_unit","bohr")
  call my_add_attribute(xf,"uuid",generate_uuid(version=1))
  call my_add_attribute(xf,"xmlns",PSML_NAMESPACE)

  call xml_NewElement(xf,"provenance")
  if (srel) then
     call my_add_attribute(xf,"creator",PSML_CREATOR // " (scalar-relativistic)")
  else
     call my_add_attribute(xf,"creator",PSML_CREATOR // " (non-relativistic)")
  endif
  call my_add_attribute(xf,"date",datestr)
  call xml_NewElement(xf,"annotation")
       call my_add_attribute(xf,"action","semilocal-pseudopotential-generation")
       call my_add_attribute(xf,"action-cont","projectors-generation")
  call xml_EndElement(xf,"annotation")
  call xml_NewElement(xf,"input-file")
  call my_add_attribute(xf,"name","oncvpsp-input")
!
!    Use INP_COPY (generated in the main program)
!
  call cdata_section_from_file(xf,"INP_COPY")

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

  !
  !
  call xml_NewElement(xf,"pseudo-atom-spec")
  call my_add_attribute(xf,"atomic-label",nameat)
  call my_add_attribute(xf,"atomic-number",str(zz))
  call my_add_attribute(xf,"z-pseudo",str(zion))
  call my_add_attribute(xf,"flavor",psflavor)
  if (srel) then
     call my_add_attribute(xf,"relativity","scalar")
  else
     call my_add_attribute(xf,"relativity","no")
  endif
  call my_add_attribute(xf,"spin-dft",polattrib)
  call my_add_attribute(xf,"core-corrections",coreattrib)

  ! Extra non-standard pieces of information
  call xml_NewElement(xf,"annotation")
    call my_add_attribute(xf,"pseudo-energy",str(epstot))
    ! Cutoff hints a la pseudo-dojo
    ! Use the last three error levels (10^-5, sqrt(10)*10^-5, 10^-4 
    call my_add_attribute(xf,"cutoff_hint_low",str(ec_hint(5)%cutoff_hint))
    call my_add_attribute(xf,"cutoff_hint_normal",str(ec_hint(6)%cutoff_hint))
    call my_add_attribute(xf,"cutoff_hint_high",str(ec_hint(7)%cutoff_hint))
  call xml_EndElement(xf,"annotation")
    
  ! XC name handling
  call exchange_correlation_info(xf,iexc)
  !
  call configuration_info()
  call xml_EndElement(xf,"pseudo-atom-spec")

!AG: save
   if(lloc==4) then
       ! fitted local potential
   else
       ! 'l_local="',lloc,'"'
   end if
!AG -- decide how to handle the case of Vlocal as one of the sl pots.

 npts = mmax

 allocate(div_by_r(mmax))
 r => rr(:)
 chval => rho(:)
 chcore => rhomod(:,1)

 ! Compute log grid step factor
 ! r(i) = r(1)*exp(a(i-1))
 al = 0.01d0 * log(rr(101)/rr(1))

 ! Use a selection of the log-grid points, with a spacing of at least
 ! delta, and up to rmax.  rmax could be determined on the basis of
 ! the tails of the wavefunctions (to accomodate them and the charge
 ! density), but it is more straightforward to simply keep the maximum
 ! rmax in the original logarithmic grid
 
 rmax = rr(mmax)
 delta = 0.5_dp*drl
 call get_sampled_grid(mmax,rr,rmax,delta,nrl_full,isample,r0)
 
 ! Get the data length for the shorter range
 ! typically used in oncvpsp (drl*nrl0 in the linear grid)
 ! This will be used for everything except the valence charge and
 ! (possibly) the wavefunctions
 
 nrl_short = get_npts_in_range(r0,range=drl*nrl0)

  !! --------------------------------------
  !! Full range for grid and valence charge
 
  nrl = nrl_full
 
  call xml_NewElement(xf,"grid")
  call my_add_attribute(xf,"npts",str(nrl))

  call xml_NewElement(xf,"annotation")
  call my_add_attribute(xf,"type","sampled from oncvpsp log grid")
  call my_add_attribute(xf,"recipe", &
       "r(i:1..N) = r1*exp(a*(i-1)) + r=0; resampled")
  call my_add_attribute(xf,"recipe-cont","r1: scale; a: step")
  call my_add_attribute(xf,"scale",trim(str(rr(1))))
  call my_add_attribute(xf,"step",trim(str(al)))
  call my_add_attribute(xf,"delta",trim(str(delta)))
  call my_add_attribute(xf,"rmax",trim(str(rmax)))
  call xml_EndElement(xf,"annotation")

  call xml_NewElement(xf,"grid-data")
  call xml_AddArray(xf,r0(1:nrl))
  call xml_EndElement(xf,"grid-data")

  call xml_EndElement(xf,"grid")

  allocate(f0(nrl))
  !
  !
  call xml_NewElement(xf,"valence-charge")
  call my_add_attribute(xf,"total-charge",  &
                      str(total_valence_charge))
  call my_add_attribute(xf,"is-unscreening-charge","yes")
  call my_add_attribute(xf,"rescaled-to-z-pseudo","no")
  call xml_NewElement(xf,"radfunc")

  call resample(r,chval,npts,r0,isample,f0,nrl)
  call check_grid(r,chval,npts,r0,f0,nrl,"chval.check")
  call add_data(xf,f0(1:nrl))
  call xml_EndElement(xf,"radfunc")
  call xml_EndElement(xf,"valence-charge")

  !! --------------------------------------
  !! Every other magnitude uses a shorter section of the grid
  nrl = nrl_short
  
  if (there_is_core) then
     rcore = rr(irct)
     call xml_NewElement(xf,"pseudocore-charge")
     call my_add_attribute(xf,"matching-radius",str(rcore))
     call my_add_attribute(xf,"number-of-continuous-derivatives", &
                                    str(4))
     call xml_NewElement(xf,"annotation")
     select case (icmod)
        case (1)
           call my_add_attribute(xf,"model-charge-form",  &
                "Polynomial")
        case (2)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function fitted using value and slope")
        case (3)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function with specified parameters")
        case (4)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function optimized using XC hardness")
        end select
     call xml_EndElement(xf,"annotation")
        
     call xml_NewElement(xf,"radfunc")

     call resample(r,chcore,npts,r0,isample,f0,nrl)
     call check_grid(r,chcore,npts,r0,f0,nrl,"chcore.check")
     call add_data(xf,f0(1:nrl))

     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pseudocore-charge")
     deallocate(chcore)
  endif

  call xml_NewElement(xf,"semilocal-potentials")
  if (srel) then
     call my_add_attribute(xf,"set","scalar_relativistic")
  else
     call my_add_attribute(xf,"set","non_relativistic")
  endif
  !  
  vpsd: do i = 1, npots
     vps => vpuns(:,i)

     ! This call resamples vps onto r0, but avoids extrapolation to r=0,
     ! in case unscreening of the potentials causes unwanted oscillations
     ! This is the same procedure introduced in V4 of oncvpsp for upf
     ! and psp8 output. It is used below for the local potential also.
     
     call resample(r,vps,npts,r0,isample,f0,nrl, override_at_zero=.true.)
     write(fname,"(a,i1,a)") "slps.l=", ll(i), ".check"
     call check_grid(r,vps,npts,r0,f0,nrl,fname)
     call write_psml_item(xf, class="slps", &
                          n=nn(i), l=ll(i), &
                          rc=rc(i), eref= epa(1,ll(i)+1),&
                          f=f0(1:nrl))
  enddo vpsd
  call xml_EndElement(xf,"semilocal-potentials")
!

        vlocal => vpuns(:,lloc+1)
        call xml_NewElement(xf,"local-potential")
            if (lloc > lmax) then
               call my_add_attribute(xf,"type","oncv-fit")
            else
               call my_add_attribute(xf,"type","l="//str(lloc))
            endif
            call xml_NewElement(xf,"radfunc")
            call resample(r,vlocal,npts,r0,isample,f0,nrl, &
                          override_at_zero=.true.)
               call check_grid(r,vlocal,npts,r0,f0,nrl,"vlocal.check")
               call add_data(xf,f0(1:nrl))
            call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"local-potential")

      call xml_NewElement(xf,"nonlocal-projectors")
      if (srel) then
         call my_add_attribute(xf,"set","scalar_relativistic")
      else
         call my_add_attribute(xf,"set","non_relativistic")
      endif

      do l1=1,lmax+1
         if(l1==lloc+1) cycle
         do jj=1,nproj(l1)
            ! Reference energies for projector construction

            eref = epa(jj,l1)

            ! Store projectors without r factor
            div_by_r(:) = vkb(:,jj,l1)/r(:)
            call resample(r,div_by_r,npts,r0,isample,f0,nrl)
            ! For l>0, the value at r=0 should be exactly 0
            if (l1 > 1) then
               f0(1) = 0.0_dp
            endif
            write(fname,"(a,i1,a,i1,a)") "proj.l=", l1-1, ".seq=", jj,".check"
            call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
            call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, &
                                 ekb=evkb(jj,l1), &
                                 eref=eref, &
                                 type="oncv", f=f0(1:nrl))
         enddo
      enddo
      call xml_EndElement(xf,"nonlocal-projectors")

      if (write_wfns) then
         call xml_NewElement(xf,"pseudo-wave-functions")
         if (srel) then
            call my_add_attribute(xf,"set","scalar_relativistic")
         else
         call my_add_attribute(xf,"set","non_relativistic")
      endif
      nrl = nrl_full
      do ii = 1, nv
         l1 = la(nc+ii) + 1
         n1 = na(nc+ii)
         energy_level = ea(nc+ii)
         div_by_r(:) = uupsa(:,ii)/r(:)
         call resample(r,div_by_r,npts,r0,isample,f0,nrl)
         ! For l>0, the value at r=0 should be exactly 0
         if (l1 > 1) then
            f0(1) = 0.0_dp
         endif
         write(fname,"(a,i1,a,i1,a)") "wfn.n=", n1, ".l=", l1-1,".check"
         call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
         call write_psml_item(xf, class="pswf", &
              n=n1, l=l1-1, &
              energy_level=energy_level,&
              f=f0(1:nrl))
      enddo
      call xml_EndElement(xf,"pseudo-wave-functions")
   endif
   
  call xml_EndElement(xf,"psml")


  call xml_Close(xf)

  deallocate(f0,r0, div_by_r)

CONTAINS
  subroutine configuration_info()

  call xml_NewElement(xf,"valence-configuration")
  call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
  do i = ncp, norbs
     if (f(i) .lt. 1.0e-10_dp) cycle
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
end subroutine configuration_info

end subroutine psmlout
!
!=================================================
! Fully relativistic version
! Output depends on setting of 'psfile': [upf |  psp8 | both]
!
 subroutine psmlout_r(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,vpuns_sr,rho,rhomod, &
&                  irct, &
&                  vsr,esr,vso,eso, &
&                  zz,zion,mmax,mxprj,iexc,icmod,nrl0,drl,atsym,epstot, &
&                  na,la,nc,nv, &
&                  fa,epa, &
&                  psfile,uupsa,ea)

!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
! vsr, esr,  vso, eso  : projectors and coeffs for SR/SO split   
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials, per j (vp(:,5,:) is local
!  potential if linear combination is used)
!vpuns_sr  scalar-relativistic average of vpuns, i.e. the same array used
!  for psp8/upf output, passed straight through instead of being
!  independently recombined here from vpuns per j (algebraically
!  equivalent, but this avoids a redundant recomputation near r=0, where
!  the per-j terms nearly cancel)
!rho  valence pseudocharge
!rhomod  model core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes pseudoion charge)
!mmax  size of log radial grid
!mxprj  dimension of number of projectors
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!nrl size of linear radial grid
!drl spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy
!psfile  should be 'upf' or 'psp8' or 'both'
!uupsa pseudo-atomic orbital array
!ea  pseudo-atomic orbital energies
!epa actual reference energies for projectors
   
  use xmlf90_wxml     ! To write XML files
  use m_uuid        ! To generate uuid

 implicit none

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
 integer :: lmax,lloc,iexc,mmax,mxprj,nrl0,icmod
 integer :: nproj(6)
 integer :: irct ! index of point at which rho_core is matched
 real(dp) :: drl,zz,zion,epstot
 real(dp), target :: rr(mmax),vpuns(mmax,5,2),rho(mmax),vkb(mmax,mxprj,4,2)
 real(dp), target :: rhomod(mmax,5)
 real(dp) :: vpuns_sr(mmax,5)
 real(dp):: rc(6),evkb(mxprj,4,2)
 real(dp), target :: vsr(mmax,2*mxprj,4),vso(mmax,2*mxprj,4)
 real(dp) :: esr(2*mxprj,4),eso(2*mxprj,4)
 character*2 :: atsym

 integer :: na(30),la(30)
 integer :: nc,nv
 real(dp) :: fa(30),epa(mxprj,6,2),ea(30,2)
 character*4 :: psfile
 real(dp) :: uupsa(mmax,2,nv)


!Local variables
 integer :: dtime(8), nrl, nrl_short, nrl_full
 integer :: npr_so(4), npr_sr(4)

  type(xmlf_t) :: xf

  integer   :: i, j, npts
  integer   :: ii, jj, l1, jk
  real(dp)  :: rcore, jval
  real(dp)  :: total_valence_charge
  real(dp), pointer :: r(:), chval(:), chcore(:)
  real(dp), pointer :: vlocal(:)
  real(dp), allocatable :: r0(:), f0(:), vps(:), div_by_r(:)

  character(len=100)    :: line
  character*4      :: polattrib, coreattrib
  character*1      :: pscode, char_dummy
  character(len=2) :: nameat
  character(len=40):: psflavor

  character(len=1), dimension(0:4) :: lsymb = (/'s','p','d','f','g'/)

  integer          :: ncore, nval, ncp, norbs, npots
  integer          :: n1
  real(dp)         :: energy_level, eref

  integer, allocatable  :: n(:), l(:)
  integer, allocatable  :: nn(:), ll(:)
  real(dp), allocatable :: f(:), ff(:)
  real(dp), allocatable :: fdown(:), fup(:)
  
  real(dp) :: rmax, delta, al
  integer, allocatable  :: isample(:)
  
  character(len=10)     :: datestr

  logical :: tdopsp, nonrel, polarized, there_is_core, found
  integer :: lun, stat
  character(len=40) :: fname, relat_output_spec
  logical :: write_wfns, want_sr_so_form, want_lj_form
 !---
  
  call get_psml_options(write_wfns,check_interp,relat_output_spec)

  ! Output options
  ! Command-line arguments take precedence over 'psfile' in the input file
  
  if (relat_output_spec == "") then
     want_sr_so_form =  (trim(psfile)=="psp8" .or. trim(psfile) == "both")
     want_lj_form =     (trim(psfile)=="upf" .or. trim(psfile) == "both")
  else
     want_sr_so_form = (relat_output_spec=="sr-so" .or. relat_output_spec=="both" )
     want_lj_form =    ( relat_output_spec=="lj" .or. relat_output_spec=="both" )
  endif

 call date_and_time(VALUES=dtime)
 write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)

  call xml_OpenFile(PSML_FILENAME,xf, indent=.false.)

  call xml_AddXMLDeclaration(xf,"UTF-8")

  call xml_NewElement(xf,"psml")
  call my_add_attribute(xf,"version",PSML_VERSION)
  call my_add_attribute(xf,"energy_unit","hartree")
  call my_add_attribute(xf,"length_unit","bohr")
  call my_add_attribute(xf,"uuid",generate_uuid(version=1))
  call my_add_attribute(xf,"xmlns",PSML_NAMESPACE)

  call xml_NewElement(xf,"provenance")
  call my_add_attribute(xf,"creator",PSML_CREATOR // " (fully-relativistic)")
  call my_add_attribute(xf,"date",datestr)
  call xml_NewElement(xf,"annotation")
       call my_add_attribute(xf,"action","semilocal-pseudopotential-generation")
       call my_add_attribute(xf,"action-cont","projectors-generation")
  call xml_EndElement(xf,"annotation")
  call xml_NewElement(xf,"input-file")
  call my_add_attribute(xf,"name","oncvpsp-input")
!
!    Use INP_COPY (generated in the main program)
!                                                                                               
  call cdata_section_from_file(xf,"INP_COPY")
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

  call xml_NewElement(xf,"pseudo-atom-spec")
  call my_add_attribute(xf,"atomic-label",nameat)
  call my_add_attribute(xf,"atomic-number",str(zz))
  call my_add_attribute(xf,"z-pseudo",str(zion))
  call my_add_attribute(xf,"flavor",psflavor)
  call my_add_attribute(xf,"relativity","dirac")
  call my_add_attribute(xf,"spin-dft",polattrib)
  call my_add_attribute(xf,"core-corrections",coreattrib)

  ! Extra non-standard pieces of information
  call xml_NewElement(xf,"annotation")
    call my_add_attribute(xf,"pseudo-energy",str(epstot))
    ! Cutoff hints a la pseudo-dojo
    ! Use the last three error levels (10^-5, sqrt(10)*10^-5, 10^-4 
    call my_add_attribute(xf,"cutoff_hint_low",str(ec_hint(5)%cutoff_hint))
    call my_add_attribute(xf,"cutoff_hint_normal",str(ec_hint(6)%cutoff_hint))
    call my_add_attribute(xf,"cutoff_hint_high",str(ec_hint(7)%cutoff_hint))
  call xml_EndElement(xf,"annotation")

  ! XC name handling
  call exchange_correlation_info(xf,iexc)
  !
  call configuration_info()
  call xml_EndElement(xf,"pseudo-atom-spec")

!AG: save
   if(lloc==4) then
       ! fitted local potential
   else
       ! 'l_local="',lloc,'"'
   end if
!AG -- decide how to handle the case of Vlocal as one of the sl pots.

 npts = mmax

 allocate(div_by_r(mmax))
 r => rr(:)
 chval => rho(:)
 chcore => rhomod(:,1)
 allocate(vps(npts))

  ! Compute log grid step factor
  ! r(i) = r(1)*exp(a(i-1))
  al = 0.01d0 * log(rr(101)/rr(1))

 ! Use a selection of the log-grid points, with a spacing of at least
 ! delta, and up to rmax.  rmax could be determined on the basis of
 ! the tails of the wavefunctions (to accomodate them and the charge
 ! density), but it is more straightforward to simply keep the maximum
 ! rmax in the original logarithmic grid

 rmax = rr(mmax)
 delta = 0.5_dp*drl
 call get_sampled_grid(mmax,rr,rmax,delta,nrl_full,isample,r0)

 ! Get the data length for the shorter range
 ! typically used in oncvpsp (drl*nrl0 in the linear grid)
 ! This will be used for everything except the valence charge and
 ! (possibly) the wavefunctions

 nrl_short = get_npts_in_range(r0,range=drl*nrl0)

  !! --------------------------------------
  !! Full range for grid and valence charge
 
  nrl = nrl_full

  call xml_NewElement(xf,"grid")
  call my_add_attribute(xf,"npts",str(nrl))

  call xml_NewElement(xf,"annotation")
  call my_add_attribute(xf,"type","sampled from oncvpsp log grid")
  call my_add_attribute(xf,"recipe", &
       "r(i:1..N) = r1*exp(a*(i-1)) + r=0; resampled")
  call my_add_attribute(xf,"recipe-cont","r1: scale; a: step")
  call my_add_attribute(xf,"scale",trim(str(rr(1))))
  call my_add_attribute(xf,"step",trim(str(al)))
  call my_add_attribute(xf,"delta",trim(str(delta)))
  call my_add_attribute(xf,"rmax",trim(str(rmax)))
  call xml_EndElement(xf,"annotation")

  call xml_NewElement(xf,"grid-data")
  call xml_AddArray(xf,r0(1:nrl))
  call xml_EndElement(xf,"grid-data")

  call xml_EndElement(xf,"grid")

  allocate(f0(nrl))
  !
  call xml_NewElement(xf,"valence-charge")
  call my_add_attribute(xf,"total-charge",  &
                      str(total_valence_charge))
  call my_add_attribute(xf,"is-unscreening-charge","yes")
  call my_add_attribute(xf,"rescaled-to-z-pseudo","no")
  
  call xml_NewElement(xf,"radfunc")

  call resample(r,chval,npts,r0,isample,f0,nrl)
  call check_grid(r,chval,npts,r0,f0,nrl,"chval.check")
  call add_data(xf,f0(1:nrl))
  call xml_EndElement(xf,"radfunc")
  call xml_EndElement(xf,"valence-charge")

  !! --------------------------------------
  !! Every other magnitude uses a shorter section of the grid
  nrl = nrl_short

  if (there_is_core) then
     rcore = rr(irct)
     call xml_NewElement(xf,"pseudocore-charge")
     call my_add_attribute(xf,"matching-radius",str(rcore))
     call my_add_attribute(xf,"number-of-continuous-derivatives", &
          str(4))
     
     call xml_NewElement(xf,"annotation")

     select case (icmod)
        case (1)
           call my_add_attribute(xf,"model-charge-form",  &
                "Polynomial")
        case (2)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function fitted using value and slope")
        case (3)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function with specified parameters")
        case (4)
           call my_add_attribute(xf,"model-charge-form",  &
                "Teter function optimized using XC hardness")
        end select
     call xml_EndElement(xf,"annotation")
        
     call xml_NewElement(xf,"radfunc")

     call resample(r,chcore,npts,r0,isample,f0,nrl)
     call check_grid(r,chcore,npts,r0,f0,nrl,"chcore.check")
     call add_data(xf,f0(1:nrl))

     call xml_EndElement(xf,"radfunc")
     call xml_EndElement(xf,"pseudocore-charge")
     deallocate(chcore)
  endif

  !  
  ! Semilocal potentials
  !
  if (want_sr_so_form) then

     call xml_NewElement(xf,"semilocal-potentials")
     call my_add_attribute(xf,"set","scalar_relativistic")
! sr components
    do i = 1, npots
     l1 = i
     ! Use the sr average as computed for psp8/upf output (vpuns_sr) rather
     ! than recombining it here from the per-j vpuns: the two are
     ! algebraically equivalent, but this avoids a redundant, independently
     ! rounded recomputation of a value near r=0 where the per-j terms
     ! nearly cancel, so the PSML sr values are guaranteed to agree with
     ! psp8/upf's local-channel treatment bit for bit within a given build.
     vps(:) = vpuns_sr(:,l1)

     ! This call resamples vps onto r0, but avoids extrapolation to r=0,
     ! in case unscreening of the potentials causes unwanted oscillations
     ! This is the same procedure introduced in V4 of oncvpsp for upf
     ! and psp8 output. It is used below for other semilocal forms and
     ! for the local potential.

     call resample(r,vps,npts,r0,isample,f0,nrl, override_at_zero=.true.)
     write(fname,"(a,i1,a)") "slps-sr.l=", ll(i), ".check"
     call check_grid(r,vps,npts,r0,f0,nrl,fname)
     call write_psml_item(xf, class="slps", &
                          n=nn(i), l=ll(i), &
                          rc=rc(i), &
                          f=f0(1:nrl))
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
     call resample(r,vps,npts,r0,isample,f0,nrl, override_at_zero=.true.)
     write(fname,"(a,i1,a)") "slps-so.l=", ll(i), ".check"
     call check_grid(r,vps,npts,r0,f0,nrl,fname)
     call write_psml_item(xf, class="slps", &
                          n=nn(i), l=ll(i), &
                          rc=rc(i),&
                          f=f0(1:nrl))
  enddo

  call xml_EndElement(xf,"semilocal-potentials")

  endif

  if (want_lj_form) then

     call xml_NewElement(xf,"semilocal-potentials")
     call my_add_attribute(xf,"set","lj")

     do i = 1, npots
        l1 = ll(i) + 1
        if (l1 == 1) then
           jval = 0.5_dp
           ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=1/2 stored in index 1
           vps(:) = vpuns(:,l1,1)
           call resample(r,vps,npts,r0,isample,f0,nrl, override_at_zero=.true.)
           write(fname,"(a,i1,a,f3.1,a)") "slps-lj.l=", ll(i), ".j=",jval,".check"
           call check_grid(r,vps,npts,r0,f0,nrl,fname)
           call write_psml_item(xf, class="slps", &
                          n=nn(i), l=ll(i), j=jval, &
                          rc=rc(i),eref= epa(1,ll(i)+1,1),&
                          f=f0(1:nrl))
        else
           do jj=1,2
              jval = ll(i) + (3-2*jj)*0.5  ! convert (1,2) to (1,-1)*1/2
              ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=1/2 stored in index 1
              vps(:) = vpuns(:,l1,jj)
              call resample(r,vps,npts,r0,isample,f0,nrl, override_at_zero=.true.)
              write(fname,"(a,i1,a,f3.1,a)") "slps-lj.l=", ll(i), ".j=",jval,".check"
              call check_grid(r,vps,npts,r0,f0,nrl,fname)
              call write_psml_item(xf, class="slps", &
                          n=nn(i), l=ll(i), j=jval, &
                          rc=rc(i),eref= epa(1,ll(i)+1,jj),&
                          f=f0(1:nrl))
           enddo
        endif
     enddo
     call xml_EndElement(xf,"semilocal-potentials")

  endif ! lj form
!

        vlocal => vpuns(:,lloc+1,1)
        call xml_NewElement(xf,"local-potential")
            if (lloc > lmax) then
               call my_add_attribute(xf,"type","oncv-fit")
            else
               call my_add_attribute(xf,"type","l="//str(lloc))
            endif
            call xml_NewElement(xf,"radfunc")
            call resample(r,vlocal,npts,r0,isample,f0,nrl, &
                          override_at_zero=.true.)
               call check_grid(r,vlocal,npts,r0,f0,nrl,"vlocal.check")
               call add_data(xf,f0(1:nrl))
            call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"local-potential")

 if (want_sr_so_form) then
!
!     Scalar-relativistic projectors. These are obtained by
!     massaging the lj components, so 'eref' has no meaning.       
!
! set up projector number for sr_so calculations based on non-zero coefficients
  npr_sr(:)=0 
  npr_so(:)=0
  do l1=1,lmax+1
     if(nproj(l1)>0) then
        do ii=1,2*nproj(l1)
           if(abs(esr(ii,l1))>0.0d0) npr_sr(l1)=npr_sr(l1)+1
           if(abs(eso(ii,l1))>0.0d0) npr_so(l1)=npr_so(l1)+1
        end do
     end if
  end do

      call xml_NewElement(xf,"nonlocal-projectors")
      call my_add_attribute(xf,"set","scalar_relativistic")
 
      do l1=1,lmax+1
         if(l1==lloc+1) cycle
         do jj=1,npr_sr(l1)
            ! Store projectors without r factor
            div_by_r(:) = vsr(:,jj,l1)/r(:)
            call resample(r,div_by_r,npts,r0,isample,f0,nrl)
            ! For l>0, the value at r=0 should be exactly 0
            if (l1 > 1) then
               f0(1) = 0.0_dp
            endif
            write(fname,"(a,i1,a,i1,a)") "proj-sr.l=", l1-1, ".seq=", jj,".check"
            call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
            call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, &
                                 ekb=esr(jj,l1), &
                                 type="oncv", f=f0(1:nrl))
         enddo
      enddo
      call xml_EndElement(xf,"nonlocal-projectors")
!
!     Spin-orbit part
!
      call xml_NewElement(xf,"nonlocal-projectors")
      call my_add_attribute(xf,"set","spin_orbit")

      do l1=2,lmax+1
         if(l1==lloc+1) cycle
         do jj=1,npr_so(l1)
            ! Store projectors without r factor
            div_by_r(:) = vso(:,jj,l1)/r(:)
            call resample(r,div_by_r,npts,r0,isample,f0,nrl)
            ! For l>0, the value at r=0 should be exactly 0
            if (l1 > 1) then   ! always in this block
               f0(1) = 0.0_dp
            endif
            write(fname,"(a,i1,a,i1,a)") "proj-so.l=", l1-1, ".seq=", jj,".check"
            call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
            call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, &
                                 ekb=eso(jj,l1), &
                                 type="oncv", f=f0(1:nrl))
         enddo
      enddo

      call xml_EndElement(xf,"nonlocal-projectors")

 endif
   
 if (want_lj_form) then


      call xml_NewElement(xf,"nonlocal-projectors")
      call my_add_attribute(xf,"set","lj")

      do l1=1,lmax+1
         if(l1==lloc+1) cycle

         if (l1 == 1) then
            ! l=0, only one j=1/2
            jval = 0.5_dp

            do jj=1,nproj(l1)

               eref = epa(jj,l1,1)

               ! Store projectors without r factor
               div_by_r(:) = vkb(:,jj,l1,1)/r(:)
               call resample(r,div_by_r,npts,r0,isample,f0,nrl)
               ! For l>0, the value at r=0 should be exactly 0
               if (l1 > 1) then ! never in this block
                  f0(1) = 0.0_dp
               endif
               write(fname,"(a,i1,a,f3.1,a,i1,a)") "proj-lj.l=", l1-1, ".j=", jval, &
                                                   ".seq=", jj, ".check"
               call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
               call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, j=jval, &
                                 ekb=evkb(jj,l1,1), &
                                 eref=eref, &
                                 type="oncv", f=f0(1:nrl))
            enddo

         else  ! l1 /=1  (l/=0)

           ! two j values
           ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=1/2 stored in index 1
           do jk=1,2
              jval = l1-1 + (3-2*jk)*0.5  ! convert (1,2) to (1,-1)*1/2
              do jj=1,nproj(l1)

                 eref = epa(jj,l1,jk)

                 ! Store projectors without r factor
                 div_by_r(:) = vkb(:,jj,l1,jk)/r(:)
                 call resample(r,div_by_r,npts,r0,isample,f0,nrl)
                 ! For l>0, the value at r=0 should be exactly 0
                 if (l1 > 1) then  ! always in this block
                    f0(1) = 0.0_dp
                 endif
                 write(fname,"(a,i1,a,f3.1,a,i1,a)") "proj-lj.l=", l1-1, ".j=", jval, &
                                                   ".seq=", jj, ".check"
                 call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
                 call write_psml_item(xf, class="proj", &
                                 seq=jj, l=l1-1, j=jval, &
                                 ekb=evkb(jj,l1,jk), &
                                 eref=eref, &
                                 type="oncv", f=f0(1:nrl))
              enddo
           enddo

        endif  ! l1 == 1

      enddo ! over l shells
      call xml_EndElement(xf,"nonlocal-projectors")
!
   endif ! lj form

   if (write_wfns) then
      call xml_NewElement(xf,"pseudo-wave-functions")
      call my_add_attribute(xf,"set","lj")

      nrl = nrl_full
      
      do ii = 1, nv
         l1 = la(nc+ii) + 1
         n1 = na(nc+ii)
         if (l1 == 1) then
            ! l=0, only one j=1/2
            ! index 1
            jval = 0.5_dp
            energy_level = ea(nc+ii,1)
            div_by_r(:) = uupsa(:,1,ii)/r(:)
            call resample(r,div_by_r,npts,r0,isample,f0,nrl)
            write(fname,"(a,i1,a,i1,a,f3.1,a)") "wfn.n=", n1, ".l=", l1-1, ".j=",jval,".check"
            call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
            call write_psml_item(xf, class="pswf", &
                 n=n1, l=l1-1, j=jval,&
                 energy_level=energy_level,&
                 f=f0(1:nrl))
         else
            ! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=1/2 stored in index 1
            do jk=1,2
               jval = l1-1 + (3-2*jk)*0.5  ! convert (1,2) to (1,-1)*1/2
               energy_level = ea(nc+ii,jk)
               div_by_r(:) = uupsa(:,jk,ii)/r(:)
               call resample(r,div_by_r,npts,r0,isample,f0,nrl)
               ! For l>0, the value at r=0 should be exactly 0
               f0(1) = 0.0_dp
               write(fname,"(a,i1,a,i1,a,f3.1,a)") "wfn.n=", n1, ".l=", l1-1, ".j=",jval,".check"
               call check_grid(r,div_by_r,npts,r0,f0,nrl,fname)
               call write_psml_item(xf, class="pswf", &
                 n=n1, l=l1-1, j=jval,&
                 energy_level=energy_level,&
                 f=f0(1:nrl))
            enddo
         endif
      enddo
      call xml_EndElement(xf,"pseudo-wave-functions")
   endif

  call xml_EndElement(xf,"psml")


  call xml_Close(xf)

  deallocate(f0,r0,div_by_r)

CONTAINS
  subroutine configuration_info()

  call xml_NewElement(xf,"valence-configuration")
  call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
  do i = ncp, norbs
     if (f(i) .lt. 1.0e-10_dp) cycle
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
end subroutine configuration_info

end subroutine psmlout_r

subroutine add_data(xf,f)
  use xmlf90_wxml

  type(xmlf_t), intent(inout)   :: xf
  real(dp), intent(in)          :: f(:)

  call xml_NewElement(xf,"data")
  call my_add_attribute(xf,"npts",str(size(f)))
  call xml_AddArray(xf,f(:))
  call xml_EndElement(xf,"data")
end subroutine add_data

     subroutine my_add_attribute(xf,name,value)
       use xmlf90_wxml, only: xmlf_t, xml_AddAttribute

       type(xmlf_t), intent(inout)   :: xf
       character(len=*), intent(in)  :: name
       character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
     end subroutine my_add_attribute

    subroutine get_unit(lun)

!     Get an available Fortran unit number

      integer, intent(out) ::  lun

      integer i
      logical unit_used

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
                               rc, ekb, eref, energy_level,&
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
       real(dp), intent(in), optional  :: eref
       character(len=*), intent(in), optional  :: type

       ! for pseudo-wave-functions
       real(dp), intent(in), optional  :: energy_level

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
         if (present(eref))  call my_add_attribute(xf,"eref",str(eref))

         if (present(energy_level))  &
              call my_add_attribute(xf,"energy_level",str(energy_level))

         if (present(flavor))  call my_add_attribute(xf,"flavor",flavor)
         if (present(type))  call my_add_attribute(xf,"type",type)

         call xml_NewElement(xf,"radfunc")
           call xml_NewElement(xf,"data")
           call my_add_attribute(xf,"npts",str(size(f)))
           call xml_AddArray(xf,f(:))
           call xml_EndElement(xf,"data")
         call xml_EndElement(xf,"radfunc")

       call xml_EndElement(xf,trim(class))

     end subroutine write_psml_item

     ! Produces files with interpolation checks
     ! Not as useful for the sampled log grid, so
     ! this routine is only enabled by the '-c' option
     !
     subroutine check_grid(r1,v1,n1,rg,vg,ng,fname)
       integer, intent(in) :: n1, ng
       real(dp), intent(in), dimension(n1) :: r1, v1
       real(dp), intent(in), dimension(ng) :: rg, vg
       character(len=*), intent(in) :: fname

       logical, save  :: first_time = .true.
       integer, save  :: lun_summary
       
       integer  :: i, lun
       real(dp) :: v5, diff, maxdiff, rmax, v7
       real(dp) :: al, sum5, sum7, diff5, diff7

       if (.not. check_interp) return
       
       if (first_time) then
          ! Open cumulative report file
          call get_unit(lun_summary)
          open(unit=lun_summary,file='GRID_CHECKS',form="formatted", &
               status="unknown",action="write",position="rewind")
          write(lun_summary,"(a30,3a16,a12)") 'Filename', 'Norm-diff@5th',&
                                              'Norm-diff@7th', 'MaxDiff', 'r(MaxDiff)'
          first_time = .false.
       endif
       
       ! Compute log grid step factor
       ! r(i) = r(1)*exp(a(i-1))
       al = 0.01d0 * dlog(r1(101)/r1(1))

       call get_unit(lun)
       open(unit=lun,file=trim(fname),form="formatted", &
            status="unknown",action="write",position="rewind")
       
       maxdiff = 0.0_dp
       sum5 = 0
       sum7 = 0
       do i = 1, n1
          if (r1(i) > rg(ng)) exit
          call dpnint1(5,rg,vg,ng,r1(i),v5,.false.)
          call dpnint1(7,rg,vg,ng,r1(i),v7,.false.)
          write(lun,"(4es24.16)") r1(i), v1(i), v5, v7
          diff5 = abs(v1(i)-v5)
          diff7 = abs(v1(i)-v7)
          sum5 = sum5 + al*r1(i)*diff5
          sum7 = sum7 + al*r1(i)*diff7
          diff = diff5
          if (diff > maxdiff) then
             rmax = r1(i)
             maxdiff = diff
          endif
       end do
       ! Normalize functional distances
       sum5 = sum5 / rg(ng)
       sum7 = sum7 / rg(ng)
       close(lun)
       write(lun_summary,"(a30,3es16.4,f12.8)") &
                  trim(fname), sum5, sum7, maxdiff, rmax
     end subroutine check_grid
!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! 
! Modified by Alberto Garcia, March 2015
! This routine is included in this module with permission from D.R. Hamann.
!
 subroutine dpnint1(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 from routine
! dpnint by D.R. Hamann. 
! Changes:
!   -- A single value is returned
!   -- It can extrapolate, instead of stopping,
!      when called with an abscissa outside the
!      data range.
!   -- If the number of data points is less than
!      npoly+1, npoly is implicitly reduced, without
!      error, and without warning.
!   -- Debug interface 
!
! local polynomial interpolation of data yy on nn points xx
! giving value val on point r
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output interpolated value val on point r

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp), intent(in) :: xx(*),yy(*)
 real(dp), intent(in) :: r
 real(dp), intent(out) :: val
 integer, intent(in)   ::  nn,npoly
 logical, intent(in)   ::  debug

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,kk,iend

! interval halving search for xx(ii) points bracketing r

   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(r>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=r

!   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)
   iend = min(istart+npoly,nn)

 !  if (debug) print *, "istart, iend: ", istart, iend
   sum=0.0d0
   do iy=istart,iend
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, iend
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   val=sum

 end subroutine dpnint1

 ! Use a selection of the log-grid points, with a spacing
 ! of at least delta, and up to rmax.
 ! rmax is determined on the basis of the tails of the
 ! wavefunctions (to accomodate the charge density)
 ! The actual tolerance is given by CUTOFF_RMAX above
 
 subroutine get_sampled_grid(mmax,rr,rmax,delta,nrl,isample,r0)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: rmax
   real(dp), intent(in)  :: delta
   integer, intent(out)  :: nrl
   integer, allocatable, intent(out) :: isample(:)
   real(dp), allocatable, intent(out) :: r0(:)

   integer  :: is, j
   real(dp) :: rs

   ! First scan to get size of sampled grid
   is = 1
   rs = 0.0_dp
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-rs) < delta) cycle
      is = is + 1
      rs = rr(j)
   enddo
   
   nrl = is
   allocate(isample(nrl),r0(nrl))

   is = 1
   r0(is) = 0.0_dp
   isample(is) = 0
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-r0(is)) < delta) cycle
      is = is + 1
      r0(is) = rr(j)
      isample(is) = j
   enddo
 end subroutine get_sampled_grid

 function get_npts_in_range(r0,range) result (npts)
   real(dp), intent(in)  :: r0(:)
   real(dp), intent(in)  :: range
   integer               :: npts

   integer :: i
   npts = size(r0)
   do i = 1, size(r0)
      if (r0(i) >= range) then
         npts = i
         exit
      endif
   enddo
 end function get_npts_in_range

 ! Prepare a resampled representation of a function (defined on
 ! grid 'rr', with values in 'ff'), using a reduced grid 'r0', whose
 ! points are a subset of those in 'rr', plus r=0.
 ! The values of the new representation ('f0') are simply copied, except
 ! at r=0 (the first point in r0). Here the user can select whether
 ! to extrapolate, or to set f0 to the first value of ff.
 ! Note that the latter makes sense only if rr is sufficiently fine
 ! near r=0 (as is the case for logarithmic grids)

 subroutine resample(rr,ff,mmax,r0,isample,f0,nrl,override_at_zero)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: ff(:)
   integer, intent(in)   :: nrl
   real(dp), intent(in)  :: r0(:)
   integer, intent(in)   :: isample(:)
   real(dp), intent(out) :: f0(:)

   logical, intent(in), optional :: override_at_zero

   integer  :: is
   real(dp) :: val
   logical  :: extrapolate_to_zero
   
   do is = 2, nrl
      f0(is) = ff(isample(is))
   enddo
   ! Choice of treatments of point at r=0
   extrapolate_to_zero = .true.
   if (present(override_at_zero)) then
      extrapolate_to_zero = .not. override_at_zero
   endif
   if (extrapolate_to_zero) then
      ! Polynomial extrapolation with sampled points
      call dpnint1(POLY_ORDER_EXTRAPOL,r0(2:),f0(2:),nrl-1,0.0_dp,val,.false.)
      f0(1) = val
   else
      ! Simply set f0(r=0) = ff(r=r1)
      f0(1) = ff(1)
   endif
   !...
   ! Others
   ! ...
   
 end subroutine resample
   
subroutine copy_input_file_for_psml()
! Makes two copies of the input file: one for oncvpsp
! to read, and another to echo the input in the PSML file

      character(len=132) :: line
      integer            :: stat

      open(unit=55,file='INPUT_FILE',action='write',status='replace', &
           form='formatted')
      open(unit=66,file='INP_COPY',action='write',status='replace', &
           form='formatted')
      do
         read(5,fmt="(a)",iostat=stat) line
         if (stat .ne. 0) exit
         write(55,fmt="(a)") trim(line)
         write(66,fmt="(a)") trim(line)
      enddo
      close(5)
      close(55)
      close(66)

!  Now re-open INPUT_FILE as unit 5 for further processing
!  by oncvpsp
      
      open(unit=5,file='INPUT_FILE',action='read',status='old', &
           form='formatted',position='rewind')
      
end subroutine copy_input_file_for_psml

subroutine cdata_section_from_file(xf,filename)
  use xmlf90_wxml

  type(xmlf_t), intent(inout)   :: xf
  character(len=*), intent(in) :: filename

  integer :: stat
  character(len=512) :: line
  character(len=1)   :: nl = char(10)
  
  character(len=32000) :: buffer ! to accumulate characters
  
      open(44,file=trim(filename),form="formatted",status="old", &
           position="rewind",action="read")
      buffer = ""
      do
         read(44,fmt="(a)",iostat=stat) line
         if (stat .ne. 0) exit
         buffer = trim(buffer) // trim(line) // nl
      enddo
      close(44)
      
      call xml_AddCDATASection(xf,trim(buffer),line_feed=.true.)
!                                                                     
    end subroutine cdata_section_from_file
    
    subroutine exchange_correlation_info(xf,iexc)
      use xmlf90_wxml

      type(xmlf_t), intent(inout)   :: xf
      integer, intent(in)           :: iexc
      
      character(len=120) :: xcfuntype, names(2), libxc_string
      character(len=60)  :: types(2), libxc_type
      integer :: x_code, c_code
      integer :: i, nfuncs
      
      integer :: libxc_id(2)

      external :: libxc_info

      call xml_NewElement(xf,"exchange-correlation")

      select case(iexc)

      case(1) 
         xcfuntype    = 'LDA -- Wigner'
         nfuncs = 2
         libxc_id = (/ 1, 2 /)
         names(1) = "Slater exchange (LDA)"
         names(2) = "Wigner (LDA)"
         types(1) = "exchange"
         types(2) = "correlation" 
      case(2) 
         xcfuntype    = 'LDA -- Hedin-Lundqvist'
         nfuncs = 2
         libxc_id = (/ 1, 4 /)
         names(1) = "Slater exchange (LDA)"
         names(2) = "Hedin & Lundqvist (LDA)"
         types(1) = "exchange"
         types(2) = "correlation" 
      case(3) 
         xcfuntype    = 'LDA -- Ceperley-Alder Perdew-Zunger'
         nfuncs = 2
         libxc_id = (/ 1, 9 /)
         names(1) = "Slater exchange (LDA)"
         names(2) = "Perdew & Zunger (LDA)"
         types(1) = "exchange"
         types(2) = "correlation" 

      case(4) 
         xcfuntype    = 'GGA -- Perdew-Burke-Ernzerhof'
         nfuncs = 2
         libxc_id = (/ 101, 130 /)
         names(1) = "Perdew, Burke & Ernzerhof (GGA)"
         names(2) = "Perdew, Burke & Ernzerhof (GGA)"
         types(1) = "exchange"
         types(2) = "correlation" 

      case(:-1)      ! libxc encoding -XXXCCC
                     !             or -YYY for single exc functional
         xcfuntype    = 'libxc-interface'
         if (-iexc < 1000) then
            ! Special syntax for single functional
            ! (For example, Teter exch-corr functional: iexc=-020
            nfuncs = 1
            libxc_id(1) = -iexc
         else
            x_code = -iexc/1000
            c_code = -iexc - 1000*x_code
            nfuncs = 2 
            libxc_id = (/ x_code, c_code /)
         endif
     
      end select

     call xml_NewElement(xf,"annotation")
     call my_add_attribute(xf,"oncvpsp-xc-code",str(iexc))
     call my_add_attribute(xf,"oncvpsp-xc-type",trim(xcfuntype))
     call xml_EndElement(xf,"annotation")

     call xml_NewElement(xf,"libxc-info")
     select case(iexc)

     case(1:4) 
        call my_add_attribute(xf,"number-of-functionals",str(nfuncs))
        do i = 1, nfuncs
           call xml_NewElement(xf,"functional")
           call my_add_attribute(xf,"name",trim(names(i)))
           call my_add_attribute(xf,"type",trim(types(i)))
           call my_add_attribute(xf,"id",str(libxc_id(i)))
           call xml_EndElement(xf,"functional")
        enddo
        
     case(:-1)      ! libxc encoding -XXXCCC
                    !             or -YYY for single exc functional
        call my_add_attribute(xf,"number-of-functionals",str(nfuncs))
        do i = 1, nfuncs
           call xml_NewElement(xf,"functional")
           call libxc_info(libxc_id(i),libxc_string,libxc_type)
           call my_add_attribute(xf,"name",trim(libxc_string))
           call my_add_attribute(xf,"type",trim(libxc_type))
           call my_add_attribute(xf,"id",str(libxc_id(i)))
           call xml_EndElement(xf,"functional")
        enddo
        
     end select
  
        call xml_EndElement(xf,"libxc-info")
  
  call xml_EndElement(xf,"exchange-correlation")

end subroutine exchange_correlation_info

subroutine get_psml_options(write_wfns,check_interp,relat_output_spec)
  use m_getopts

  logical, intent(out)          :: write_wfns
  logical, intent(out)          :: check_interp
  character(len=*), intent(out), optional :: relat_output_spec
  
      character(len=200) :: opt_arg
      character(len=10)  :: opt_name, spec
      integer :: iostat, n_opts

      !     Process options

      write_wfns = .false.
      check_interp = .false.
      spec = ""
      n_opts = 0
      do
!        i: and t: (BE-ref's -i/-t input-file selectors) and h,v (BE-ref's
!        -h/-v help/version flags) are accepted here purely so this scan
!        does not treat them as unknown options; the main program's own
!        argument parser already handles them.
         call getopts('wr:ci:t:hv',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
           case ('c')
              check_interp = .true.
           case ('w')
              write_wfns = .true.
           case ('r')
              read(opt_arg,*) spec
           case ('i','t','h','v')
              ! Handled by the main program; ignored here.
           case ('?',':')
             write(0,*) "Invalid option: ", opt_arg(1:1)
             write(0,*) "Usage: oncvpsp [-c -w -r spec]"
             write(0,*) " -c  : produce interpolation .check files"
             write(0,*) " -w  : output pseudo-wavefunctions in psml file"
             write(0,*) " -r [ lj | sr-so | both ] : kind of relativistic output" 
             STOP
          end select
       enddo
       if (present(relat_output_spec)) then
          relat_output_spec = spec
       endif

end subroutine get_psml_options

subroutine get_ec_hints(cvgplt) 
! Find cutoff hints 'a la pseudo-dojo (ppgen starting point)'
  real(dp), intent(in) :: cvgplt(2,7,2,4)

  integer :: i
  
    ! Take the maxima over angular momenta; nearest integer
  do i = 1, 7
     ec_hint(i) = hint_t(cvgplt(2,i,1,1),nint(maxval(cvgplt(1,i,1,:))))
  enddo
  
end subroutine get_ec_hints

end module m_psmlout
