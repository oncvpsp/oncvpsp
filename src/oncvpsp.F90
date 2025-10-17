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
program oncvpsp
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!
   use, intrinsic :: iso_fortran_env, only: stdin => input_unit, stdout => output_unit, stderr => error_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use m_psmlout, only: psmlout
   use input_text_m, only: read_input_text
   use postprocess_m, only: get_wavefunctions, run_test_configurations
   use output_text_m, only: write_input_text, &
      write_reference_configuration_results_text, &
      write_test_configs_text, &
      write_rho_vpuns_text, write_vloc_text, &
      write_rho_rhoc_rhom_text, &
      write_wavefunctions_vkb_text, &
      write_convergence_profile_text, &
      write_phsft_text, &
      write_modcore1_text, &
      write_modcore2_text, &
      write_modcore3_text, &
      write_modcore4_text
   use modcore1_m, only: modcore1, get_modcore1_match
   use modcore2_m, only: modcore2, get_modcore2_match
   use modcore3_m, only: modcore3, get_modcore3_match
   use modcore4_m, only: modcore4
   use utils_m, only: linspace, get_pseudo_linear_mesh_parameters
#if (defined WITH_TOML)
   use input_toml_m, only: read_input_toml
#endif
#if (defined WITH_HDF5)
   use hdf5_utils_m, only: HID_T, hdf_open_file, hdf_close_file
   use output_hdf5_m, only: do_hdf5, hdf5_file_id
   use output_hdf5_m, only: write_input_hdf5, &
      write_log_mesh_hdf5, &
      write_sratom_hdf5, &
      write_optimize_inputs_hdf5, &
      write_reference_configuration_results_hdf5, &
      write_output_hdf5
#endif
   implicit none

   integer, parameter :: mxprj = 5
   real(dp), parameter :: amesh = 1.006_dp
#if RELATIVISTIC == 1
   logical, parameter :: srel = .true.
#elif RELATIVISTIC == 0
   logical, parameter :: srel = .false.
#endif
!
   integer :: ii,ierr,iexc,iexct,ios,iprint,irps,it,icmod,lpopt
   integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
   integer :: mch,mchf,mmax,n1,n2,n3,n4,nc,nlim,nlloc,nlmax,nrl
   integer :: nv,irct,ncnf
   integer :: iprj
   integer,allocatable :: npa(:,:)
!
   integer :: dtime(8),na(30),la(30),np(6)
   integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
   integer :: indxr(30),indxe(30)
   integer :: irc(6),nodes(4)
   integer :: nproj(6),npx(6),lpx(6)
   integer :: ncon(6),nbas(6)

   real(dp) :: al,csc,csc1,deblt,depsh,depsht,drl,eeel
   real(dp) :: eeig,eexc
   real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t,rxpsh
   real(dp) :: et,etest,emin,sls
   real(dp) :: fcfact,rcfact,dvloc0
   real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
   real(dp) :: sf,zz,zion,zval,etot
   real(dp) :: xdummy
!
   real(dp) :: cl(6),debl(6),ea(30),ep(6),fa(30),facnf(30,5)
   real(dp) :: fnp(6),fp(6)
   real(dp) :: qcut(6),qmsbf(6),rc(6),rc0(6)
   real(dp) :: rpk(30)
   real(dp) :: epx(6),fpx(6)
   real(dp) :: epstot
   real(dp), parameter :: eps=1.0d-8

   real(dp), allocatable :: evkb(:,:),cvgplt(:,:,:,:),qq(:,:)
   real(dp), allocatable :: rr(:)
   real(dp), allocatable :: rho(:),rhoc(:),rhot(:),rhozero(:)
   real(dp), allocatable :: uu(:),up(:)
   real(dp), allocatable :: vp(:,:),vfull(:),vkb(:,:,:),pswf(:,:,:)
   real(dp), allocatable :: vwell(:)
   real(dp), allocatable :: vpuns(:,:)
   real(dp), allocatable :: vo(:),vxc(:)
   real(dp), allocatable :: rhoae(:,:),rhops(:,:),rhotae(:)
   real(dp), allocatable :: uupsa(:,:) !pseudo-atomic orbitals array
   real(dp), allocatable :: epa(:,:),fpa(:,:)
   real(dp), allocatable :: uua(:,:),upa(:,:)
   real(dp), allocatable :: usratom(:,:), upsratom(:,:)
   real(dp), allocatable :: vr(:,:,:)

   logical :: isboundpa(mxprj)

   ! Model core charge optimization variables
   !> Model core charge density and its first n derivatives
   real(dp), allocatable :: rhomod(:, :)
   !> Radial cutoff index for d2Exc evaluation in model core optimization
   integer :: irmod
   !> 2nd derivate of all-electron Exc
   real(dp), allocatable :: d2excae(:, :)
   real(dp), allocatable :: d2exc_dummy(:, :)
   real(dp) :: d2exc_rmse_dummy
   !> 2nd derivate of pseudo Exc without model core
   real(dp), allocatable :: d2excps_no_rhom(:, :)
   !> RMSE of 2nd derivate difference without model core vs AE
   real(dp) :: d2exc_rmse_no_rhom
   !> 2nd derivate difference with model core
   real(dp), allocatable :: d2excps_rhom(:, :)
   !> RMSE of 2nd derivate difference with model core vs AE
   real(dp) :: d2exc_rmse_rhom
   !> Number of iterations in modcore1
   integer :: modcore1_iter
   !> Index of core charge - valence pseudocharge crossover in modcore1
   integer :: modcore1_ircc
   !> Index of core charge - scaled valence pseudocharge crossover in modcore2
   integer :: modcore2_ircc
   !> Dimensionless factor a in modcore2
   real(dp) :: modcore2_a0
   !> Dimensionless factor b in modcore2
   real(dp) :: modcore2_b0
   !> Log radial mesh index of core charge - valence pseudocharge crossover in modcore3
   integer :: modcore3_ircc
   !> Radial coordinate of core charge - valence pseudocharge crossover in modcore3
   real(dp) :: modcore3_rmatch
   !> Value of core charge at crossover in modcore3
   real(dp) :: modcore3_rhocmatch
   !> Whether teter_amp* and teter_scale* are relative to modcore3 rmatch/rhocmatch
   logical :: teter_relative
   !> Teter amplitude parameter for icmod=4
   real(dp) :: teter_amp
   !> Size of Teter amplitude search grid
   integer :: n_teter_amp
   !> Minimum amplitude prefactor for coarse grid search
   real(dp) :: teter_amp_min
   !> Maximum amplitude prefactor for coarse grid search
   real(dp) :: teter_amp_max
   !> Amplitude prefactor step size for coarse grid
   real(dp) :: teter_amp_step
   !> Teter amplitude prefactor search grid
   real(dp), allocatable :: teter_amp_prefacs(:)
   !> Teter amplitude parameter search grid
   real(dp), allocatable :: teter_amp_params(:)
   !> Teter scale parameter for icmod=3
   real(dp) :: teter_scale
   !> Size of Teter scale search grid
   integer :: n_teter_scale
   !> Minimum scale prefactor for coarse grid search
   real(dp) :: teter_scale_min
   !> Maximum scale prefactor for coarse grid search
   real(dp) :: teter_scale_max
   !> Scale prefactor step size for coarse grid
   real(dp) :: teter_scale_step
   !> Teter amplitude prefactor search grid
   real(dp), allocatable :: teter_scale_prefacs(:)
   !> Teter scale parameter search grid
   real(dp), allocatable :: teter_scale_params(:)
   !> Objective function for Teter parameter optimization
   character(len=1024) :: teter_objective_name
   !> d2Exc RMSE values from Teter grid search
   real(dp), allocatable :: grid_objective(:,:)
   !> Optimal Teter amplitude parameter from grid search
   real(dp) :: grid_opt_amp_param
   !> Optimal Teter scale parameter from grid search
   real(dp) :: grid_opt_scale_param
   real(dp) :: grid_opt_objective
   real(dp) :: nm_atol
   !> Optimal Teter amplitude parameter from Nelder-Mead search
   real(dp) :: nm_opt_amp_param
   !> Optimal Teter scale parameter from Nelder-Mead search
   real(dp) :: nm_opt_scale_param
   real(dp) :: nm_opt_objective
   !> Number of iterations from Nelder-Mead search
   integer :: nm_iter

   ! Wavefunction variables for postprocessing
   !> All-electron wavefunctions
   real(dp), allocatable :: uu_ae(:,:,:)
   !> Pseudo wavefunctions
   real(dp), allocatable :: uu_ps(:,:,:)
   !> All-electron wavefunction radial derivatives
   real(dp), allocatable :: up_ae(:,:,:)
   !> Pseudo wavefunction radial derivatives
   real(dp), allocatable :: up_ps(:,:,:)
   !> Matching points for all-electron wavefunctions
   integer :: mch_ae(mxprj,4)
   !> Matching points for pseudo wavefunctions
   integer :: mch_ps(mxprj,4)
   !> All-electron state energies
   real(dp) :: e_ae(mxprj,4)
   !> Pseudo state energies
   real(dp) :: e_ps(mxprj,4)
   !> Signs of all-electron wavefunctions at matching points
   real(dp) :: sign_ae(mxprj,4)
   !> Signs of pseudo wavefunctions at matching points
   real(dp) :: sign_ps(mxprj,4)
   !> .true. for scattering states, .false. for bound states
   logical :: is_scattering(mxprj,4)

   !> Phase shift variables for log derivative postprocessing
   !> Log radial grid mesh index at which log derivative is calculated
   integer :: irpsh(4)
   !> Radius at which log derivative is calculated
   real(dp) :: rpsh(4)
   !> Size of phase shift energies
   integer :: npsh
   !> Phase shift (log derivative) energies
   real(dp), allocatable :: epsh(:)
   !> All-electron phase shifts (log derivatives)
   real(dp), allocatable :: pshf(:,:)
   !> Pseudo phase shifts (log derivatives)
   real(dp), allocatable :: pshp(:,:)

   ! Test configuration variables
   !> Number of valence states per test
   integer :: nvt(5)
   integer :: nat(30,5)
   integer :: lat(30,5)
   real(dp) :: fat(30,3,5)
   real(dp) :: fat3(30,5)
   real(dp) :: eat(30,3,5)
   real(dp) :: eat3(30,5)
   real(dp) :: eatp(30,5)
   real(dp) :: eaetst(5), etsttot(5)

   character*2 :: atsym
   character*4 :: psfile

   logical :: cset

   integer, parameter :: INPUT_STDIN=1, INPUT_TEXT=2, INPUT_TOML=3
   integer :: input_mode
   integer :: unit
   character(len=1024) :: input_filename

#if (defined WITH_HDF5)
   character(len=1024) :: hdf5_filename
#endif

   input_mode = INPUT_STDIN
   input_filename = ''
#if (defined WITH_HDF5)
   hdf5_filename = ''
#endif

   write(6,'(a/a//)') &
      &      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
      &      'scalar-relativistic version 4.0.1 03/01/2019'

   write(6,'(a/a/a//)') &
      &      'While it is not required under the terms of the GNU GPL, it is',&
      &      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
      &      'in any publication utilizing these pseudopotentials.'

   parse_args: block
      integer :: i
      character(len=256) :: arg
      do i = 1, command_argument_count()
         call get_command_argument(i, arg)
         select case(arg)
          case('-h', '--help')
            write (stdout, '(a)') 'Usage: oncvpsp.x [options]'
            write (stdout, '(a)') 'Options:'
            write (stdout, '(a)') '  -h, --help       Show this help message and exit'
            write (stdout, '(a)') '  -v, --version    Show version information and exit'
            write (stdout, '(a)') '  -i, --input      Specify input file (text). For legacy input, use stdin.'
            write (stdout, '(a)') '  -t, --toml-input Specify input file (toml)'

            stop 0
          case('-v', '--version')
            write (stdout, '(a)') 'ONCVPSP (scalar-realtivistic) version 4.0.1'
            stop 0
          case('-i', '--input')
            if (i + 1 > command_argument_count()) then
               write (stderr, '(a)') 'Error: --input requires a filename argument'
               stop 1
            end if
            call get_command_argument(i + 1, input_filename)
            input_mode = INPUT_TEXT
#if (defined WITH_TOML)
          case('-t', '--toml-input')
            if (i + 1 > command_argument_count()) then
               write (stderr, '(a)') 'Error: --toml-input requires a filename argument'
               stop 1
            end if
            call get_command_argument(i + 1, input_filename)
            input_mode = INPUT_TOML
#else
          case('-t', '--toml-input')
            error stop 'Error: TOML input support not enabled in this build.'
#endif
#if (defined WITH_HDF5)
          case('-h5', '--hdf5-output')
            if (i + 1 > command_argument_count()) then
               write (stderr, '(a)') 'Error: --hdf5-output requires a filename argument'
               stop 1
            end if
            call get_command_argument(i + 1, hdf5_filename)
#endif
          case default
         end select
      end do
   end block parse_args

   ios = 0
   select case(input_mode)
    case(INPUT_STDIN, INPUT_TEXT)
      if (input_mode == INPUT_STDIN) then
         unit = stdin
      else
         open(newunit=unit, file=input_filename, status='old', action='read', iostat=ios)
      end if
      call read_input_text(unit, inline, &
                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                           lmax, rc, ep, ncon, nbas, qcut, &
                           lloc, lpopt, dvloc0, nproj, debl, &
                           icmod, fcfact, rcfact, &
                           teter_amp_min, teter_amp_max, teter_amp_step, & ! teter grid search
                           teter_scale_min, teter_scale_max, teter_scale_step, & ! teter grid search
                           epsh1, epsh2, depsh, rxpsh, &
                           rlmax, drl, &
                           ncnf, nvcnf, nacnf, lacnf, facnf)
      teter_relative = .true.
      teter_amp = fcfact
      teter_scale = rcfact
      teter_objective_name = 'd2exc_rmse'
      if (input_mode == INPUT_TEXT) close(unit)
#if (defined WITH_TOML)
    case(INPUT_TOML)
      open(newunit=unit, file=input_filename, status='old', action='read', iostat=ios)
      call read_input_toml(unit, &
                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, & ! atom
                           lmax, rc, ep, ncon, nbas, qcut, & ! pseudopotential
                           lloc, lpopt, dvloc0, nproj, debl, & ! local potential
                           icmod, fcfact, rcfact, & ! core charge
                           teter_relative, teter_amp, teter_scale, & ! teter core charge
                           teter_amp_min, teter_amp_max, teter_amp_step, & ! teter grid search
                           teter_scale_min, teter_scale_max, teter_scale_step, & ! teter grid search
                           teter_objective_name, & ! teter optimization objective
                           epsh1, epsh2, depsh, rxpsh, & ! logarithmic derivative / phase shift analysis
                           rlmax, drl, & ! linear output grid
                           ncnf, nvcnf, nacnf, lacnf, facnf) ! test configurations
      close(unit)
#endif
    case default
      write (stderr, '(a,i2)') 'Error: Unknown input mode =', input_mode
      stop 1
   end select

#if (defined WITH_HDF5)
   if (trim(hdf5_filename) /= '') then
      do_hdf5 = .true.
      call hdf_open_file(hdf5_file_id, hdf5_filename, 'REPLACE', 'WRITE')
   else
      do_hdf5 = .false.
   end if
#endif

   if(ios /= 0) then
      write(6,'(a)') 'oncvpsp: ERROR reading input file'
      stop
   end if

   nvcnf(1)=nv
   do ii=1,nc+nv
      nacnf(ii,1)=na(ii)
      lacnf(ii,1)=la(ii)
      facnf(ii,1)=fa(ii)
   end do

   do jj=2,ncnf+1
      do ii=1,nc
         nacnf(ii,jj)=na(ii)
         lacnf(ii,jj)=la(ii)
         facnf(ii,jj)=fa(ii)
      end do
   end do

   call check_data(atsym,zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl,fa,facnf, &
      &                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
      &                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

   nrl=int((rlmax/drl)-0.5d0)+1

   ! PWSCF wants an even number of mesh points
   ! if(trim(psfile)=='upf') then
   if(mod(nrl,2)/=0) nrl=nrl+1
   ! end if

   al=dlog(amesh)
   rr1=0.0005d0/zz
   rr1=dmin1(rr1,0.0005d0/10)
   mmax=dlog(45.0d0 /rr1)/al

   ! calculate zion for output
   zion=zz
   do ii=1,nc
      zion=zion-fa(ii)
   end do


   allocate(rr(mmax))
   allocate(rho(mmax),rhoc(mmax),rhot(mmax))
   allocate(uu(mmax),up(mmax),uupsa(mmax,30))
   allocate(evkb(mxprj,4), cvgplt(2,7,mxprj,4),qq(mxprj,mxprj))
   allocate(vp(mmax,5),vfull(mmax),vkb(mmax,mxprj,4),pswf(mmax,mxprj,4))
   allocate(vwell(mmax))
   allocate(vpuns(mmax,5))
   allocate(vo(mmax),vxc(mmax))
   allocate(rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax))
   allocate(npa(mxprj,6))
   allocate(epa(mxprj,6),fpa(mxprj,6))
   allocate(uua(mmax,mxprj),upa(mmax,mxprj))
   allocate(vr(mmax,mxprj,6))
   allocate(uu_ae(mmax,mxprj,4), up_ae(mmax,mxprj,4), uu_ps(mmax,mxprj,4), up_ps(mmax,mxprj,4))
   allocate(usratom(mmax,nc+nv), upsratom(mmax,nc+nv))

   vr(:,:,:)=0.0d0
   vp(:,:)=0.0d0
   vkb(:,:,:)=0.0d0
   epa(:,:)=0.0d0

   do ii=1,mmax
      rr(ii)=rr1*exp(al*(ii-1))
   end do
   call write_log_mesh_hdf5(hdf5_file_id, mmax, rr)
   !
   ! full potential atom solution
   !
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
      &           rr,vfull,vxc,zz,mmax,iexc,etot,ierr,srel,usratom,upsratom)
   call write_sratom_hdf5(hdf5_file_id, &
                          nc, nv, &
                          na, la, fa, &
                          mmax, rr, &
                          vfull, vxc, &
                          rhoc, rho, &
                          usratom, upsratom, &
                          rpk, ea, etot)

   ! Drop digits beyond 5 decimals for input rcs before making any use of them
   do l1=1,max(lmax+1,lloc+1)
      jj=int(rc(l1)*10.0d5)
      rc(l1)=jj/10.0d5
   end do

   rcmax=0.0d0
   do l1=1,lmax+1
      rcmax=dmax1(rcmax,rc(l1))
   end do
   do l1=1,lmax+1
      if(rc(l1)==0.0d0) then
         rc(l1)=rcmax
      end if
   end do

   nproj(lloc+1)=0
   rc0(:)=rc(:)

   ! output printing (echos input data, with all-electron eigenvalues added)
   call write_input_text(stdout, &
                         atsym, zz, nc, nv, iexc, psfile, &
                         na, la, fa, &
                         lmax, &
                         rc, ep, ncon, nbas, qcut, &
                         lloc, lpopt, dvloc0, &
                         nproj, debl, &
                         icmod, fcfact, rcfact, &
                         epsh1, epsh2, depsh, &
                         rlmax, drl, &
                         ncnf, nvcnf, nacnf, lacnf, facnf, &
                         ea)
   call write_reference_configuration_results_text(stdout, it, 100, etot)
#if (defined WITH_HDF5)
   if (do_hdf5) then
      ! These need to be written now because some of them are modified later
      ! (notably fcfact and rcfact)
      call write_input_hdf5(hdf5_file_id, &
                            atsym, zz, nc, nv, iexc, psfile, &
                            na, la, fa, &
                            lmax, &
                            rc, ep, ncon, nbas, qcut, &
                            lloc, lpopt, dvloc0, &
                            nproj, debl, &
                            icmod, fcfact, rcfact, &
                            epsh1, epsh2, depsh, &
                            rlmax, drl, &
                            ncnf, nvcnf, nacnf, lacnf, facnf)
      ! Ditto here; etot and ea are likely overwritten later
      call write_reference_configuration_results_hdf5(hdf5_file_id, nc + nv, it, 100, etot, ea)
   end if
#endif

   ! find log mesh point nearest input rc
   rcmax=0.0d0
   irc(:)=0
   do l1=1,max(lmax+1,lloc+1)
      rct=rc(l1)
      irc(l1)=0
      do ii=2,mmax
         if(rr(ii)>rct) then
            irc(l1)=ii
            rc(l1)=rr(ii)
            exit
         end if
      end do
      rcmax=dmax1(rcmax,rc(l1))
   end do

   cvgplt(:,:,:,:)=0.0d0
   uua(:, :)=0.0d0
   upa(:, :)=0.0d0
   !
   ! loop to construct pseudopotentials for all angular momenta
   !
   write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
      &      'and semi-local pseudopoentials for all angular momenta'
   ! temporarily set this to 1 so that the pseudo wave function needed for the
   ! local potential will be generated.  Reset after run_vkb.
   nproj(lloc+1)=1
   do l1=1,lmax+1
      ll=l1-1
      uu(:)=0.0d0; up(:)=0.d0; qq(:,:)=0.0d0; isboundpa(:)=.false.
      iprj=0

      ! get principal quantum number for the highest core state for this l
      npa(1,l1)=l1
      do kk=1,nc
         if(la(kk)==l1-1) npa(1,l1)=na(kk)+1
      end do !kk

      ! get all-electron bound states for projectors
      if(nv/=0) then
         do kk=nc+1,nc+nv
            if(la(kk)==l1-1) then
               iprj=iprj+1
               et=ea(kk)
               call lschfb(na(kk),la(kk),ierr,et, &
                  &                      rr,vfull,uu,up,zz,mmax,mch,srel)
               if(ierr /= 0) then
                  write(6,'(/a,3i4)') 'oncvpsp: lschfb convergence ERROR n,l,iter=', &
                     &           na(ii),la(ii),it
                  stop
               end if
               epa(iprj,l1)=ea(kk)
               npa(iprj,l1)=na(kk)
               uua(:,iprj)=uu(:)
               upa(:,iprj)=up(:)
               isboundpa(iprj)=.true.
            end if !la(kk)==l1-1
            if(iprj==nproj(l1)) exit
         end do !kk
      end if !nv/=0

      ! get all-electron well states for projectors
      ! if there were no valence states, use ep from input data for 1st well state
      ! otherwise shift up by input debl
      if(iprj==0) epa(1,l1)=ep(l1)
      if(iprj<nproj(l1))then
         do kk=1,nproj(l1)-iprj
            iprj=iprj+1
            if(iprj>1 .and. debl(l1)<=0.0d0) then
               write(6,'(a,f8.3,a/a)') 'oncvpsp: ERROR debl =',debl, 'for l=', &
                  &              ' ERROR not allowed with 2 or more scattering states', &
                  &              'program will stop'
               stop
            end if
            if(iprj>1) then
               epa(iprj,l1)=epa(iprj-1,l1)+debl(l1)
               npa(iprj,l1)=npa(iprj-1,l1)+1
            end if

            call wellstate(npa(iprj,l1),ll,irc(l1),epa(iprj,l1),rr, &
               &                     vfull,uu,up,zz,mmax,mch,srel)
            uua(:,iprj)=uu(:)
            upa(:,iprj)=up(:)
            isboundpa(iprj)=.false.
         end do !kk
      end if !iprj<nproj(l1)

      call write_optimize_inputs_hdf5(hdf5_file_id, &
                                      mxprj, ll, nproj(l1), &
                                      npa(:, l1), epa(:, l1), &
                                      mmax, &
                                      uua, upa, &
                                      vr(:, :, l1), &
                                      isboundpa)

      do iprj=1,nproj(l1)

         ! calculate relativistic correction to potential to force projectors to 0 at rc
         call vrel(ll,epa(iprj,l1),rr,vfull,vr(1,iprj,l1),uua(1,iprj),upa(1,iprj), &
                   zz,mmax,irc(l1),srel)

      end do

      ! get all-electron overlap matrix
      do jj=1,nproj(l1)
         do ii=1,jj
            call fpovlp(uua(1,ii),uua(1,jj),irc(l1),ll,zz,qq(ii,jj),rr,srel)
            qq(jj,ii)=qq(ii,jj)
         end do
      end do

      call run_optimize(epa(1,l1),ll,mmax,mxprj,rr,uua,qq, &
                        irc(l1),qcut(l1),qmsbf(l1),ncon(l1),nbas(l1),nproj(l1), &
                        pswf(1,1,l1),vp(1,l1),vkb(1,1,l1),vfull,cvgplt(1,1,1,l1))

   end do !l1

   ! construct Vanderbilt / Kleinman-Bylander projectors
   write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'
   call run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj,pswf,vfull,vp, &
                evkb,vkb,nlim,vr)

   ! restore this to its proper value
   nproj(lloc+1)=0

   deallocate(uua,upa)

   ! accumulate charge and eigenvalues
   ! pseudo wave functions are calculated with VKB projectors for
   ! maximum consistency of unscreening
   ! get all-electron and pseudopotential valence-state by valence-state
   ! charge densities

   ! null charge and eigenvalue accumulators
   uupsa(:,:)=0.0d0
   eeig=0.0d0
   zval=0.0d0
   rho(:)=0.0d0
   nodes(:)=0
   rhotae(:)=0.0d0
   irps=0
   do kk=1,nv

      ! Get all-electron valence state
      et=ea(nc+kk)
      ll=la(nc+kk)
      l1=ll+1
      call lschfb(na(nc+kk),ll,ierr,et, &
                  rr,vfull,uu,up,zz,mmax,mch,srel)
      if(ierr /= 0) then
         write(6,'(/a,3i4)') 'oncvpsp: lschfb convergence ERROR n,l,iter=', &
            &     na(ii),la(ii),it
         stop
      end if

      ! Accumulate all-electron charge density
      rhoae(:,kk)=(uu(:)/rr(:))**2
      rhotae(:)=rhotae(:) + fa(nc+kk)*rhoae(:,kk)

      ! Get pseudo valence state
      emax=0.75d0*et
      emin=1.25d0*et
      call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
                    rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
                    uu,up,mmax,mch)
      if(ierr/=0) then
         write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
         flush(6)
         stop
      end if

      ! save valence pseudo wave functions for upfout
      uupsa(:,kk)=uu(:)

      ! Accumulate pseudo charge density
      rhops(:,kk)=(uu(:)/rr(:))**2
      rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)

      eeig=eeig+fa(nc+kk)*et
      zval=zval+fa(nc+kk)
      nodes(l1)=nodes(l1)+1
      irps=max(irps,irc(l1))
   end do !kk

   ! Construct model core charge
   allocate(rhomod(mmax,5), source=0.0_dp)
   allocate(rhozero(mmax), source=0.0_dp)
   allocate(d2exc_dummy(nv, nv))
   allocate(d2excae(nv, nv))
   allocate(d2excps_no_rhom(nv, nv))
   allocate(d2excps_rhom(nv, nv))
   select case (icmod)
    case (1)
      call get_modcore1_match(irps, fcfact, mmax, rhoc, rho, modcore1_ircc, irmod)
      call modcore1(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
                    fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc, &
                    modcore1_iter)
    case (2)
      irmod = mmax
      call get_modcore2_match(mmax, rr, rhoc, rho, fcfact, modcore2_ircc, modcore2_a0, modcore2_b0)
      call get_modcore3_match(mmax, rr, rhoc, rho, modcore3_ircc, modcore3_rmatch, modcore3_rhocmatch)
      call modcore2(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
                    fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)
    case (3)
      irmod = mmax
      call get_modcore3_match(mmax, rr, rhoc, rho, modcore3_ircc, modcore3_rmatch, modcore3_rhocmatch)
      if (teter_relative) then
         teter_amp = fcfact * modcore3_rhocmatch
         teter_scale = rcfact * modcore3_rmatch
      else
         teter_amp = fcfact
         teter_scale = rcfact
      end if
      call modcore3(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
                    teter_amp,teter_scale,irps,mmax,rr,nc,nv,la,zion,iexc)
    case (4)
      irmod = mmax
      call get_modcore3_match(mmax, rr, rhoc, rho, modcore3_ircc, modcore3_rmatch, modcore3_rhocmatch)
      if (teter_relative) then
         call linspace(teter_amp_min, teter_amp_max, teter_amp_step, n_teter_amp, teter_amp_prefacs)
         call linspace(teter_scale_min, teter_scale_max, teter_scale_step, n_teter_scale, teter_scale_prefacs)
         allocate(teter_amp_params(n_teter_amp), teter_scale_params(n_teter_scale))
         teter_amp_params(:) = teter_amp_prefacs(:) * modcore3_rhocmatch
         teter_scale_params(:) = teter_scale_prefacs(:) * modcore3_rmatch
      else
         call linspace(teter_amp_min, teter_amp_max, teter_amp_step, n_teter_amp, teter_amp_params)
         call linspace(teter_scale_min, teter_scale_max, teter_scale_step, n_teter_scale, teter_scale_params)
         allocate(teter_amp_prefacs(n_teter_amp), teter_scale_prefacs(n_teter_scale))
         teter_amp_prefacs(:) = teter_amp_params(:) / modcore3_rhocmatch
         teter_scale_prefacs(:) = teter_scale_params(:) / modcore3_rmatch
      end if
      allocate(grid_objective(n_teter_amp, n_teter_scale))
      call modcore4(mmax, rr, nc, nv, la, zion, irps, iexc, &
                    rhops, rho, rhoc, rhoae, rhotae, &
                    teter_objective_name, &
                    n_teter_amp, n_teter_scale, teter_amp_params, teter_scale_params, &
                    grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, &
                    nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, nm_iter, &
                    rhomod)
      ! Expected values for psp8 output
      fcfact = 1.0_dp
      rcfact = 0.0_dp
    case default
   end select

   if (icmod > 0) then
      ! Compute d2Exc with core charge = true AE core charge
      call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2exc_dummy, d2exc_rmse_dummy, &
                   zion, iexc, nc, nv, la, irmod, mmax)
      ! Compute d2Exc with core charge = 0
      call der2exc(rho, rhozero, rhops, rr, d2excps_no_rhom, d2excae, d2exc_rmse_no_rhom, &
                   zion, iexc, nc, nv, la, irmod, mmax)
      ! Compute d2Exc with core charge = model core charge
      call der2exc(rho, rhomod(:, 1), rhops, rr, d2excps_rhom, d2excae, d2exc_rmse_rhom, &
                   zion, iexc, nc, nv, la, irmod, mmax)
   end if

   ! Write output
   select case (icmod)
    case (1)
      call write_modcore1_text(stdout, nv, d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, &
                               modcore1_iter, d2excps_rhom, d2exc_rmse_rhom)
    case (2)
      call write_modcore2_text(stdout, nv, d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, &
                               modcore3_rmatch, modcore3_rhocmatch, &
                               d2excps_rhom, d2exc_rmse_rhom, &
                               modcore2_a0, modcore2_b0)
    case (3)
      call write_modcore3_text(stdout, nv, d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, &
                               modcore3_rmatch, modcore3_rhocmatch, d2excps_rhom, d2exc_rmse_rhom)
    case (4)
      call write_modcore4_text(stdout, nv, d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, &
                               modcore3_rmatch, modcore3_rhocmatch, d2excps_rhom, d2exc_rmse_rhom, &
                               n_teter_amp, teter_amp_prefacs, n_teter_scale, teter_scale_prefacs, &
                               grid_objective, nm_opt_amp_param, nm_opt_scale_param, nm_iter)
    case default
   end select

! screening potential for pseudocharge

   call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! total energy output

   epstot= eeig + eexc - 0.5d0*eeel
   write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

   call run_diag(lmax,npa,epa,lloc,irc, &
      &                    vkb,evkb,nproj,rr,vfull,vp,zz,mmax,mxprj,srel)

   call run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
      &                    vkb,evkb,nproj,rr,vp,mmax,mxprj)

! unscreen semi-local potentials

   do l1=1,max(lmax+1,lloc+1)
      vpuns(:,l1)=vp(:,l1)-vo(:)
   end do

!fix unscreening error due to greater range of all-electron charge
   do ii=mmax,1,-1
      if(rho(ii)==0.0d0) then
         do l1=1,max(lmax+1,lloc+1)
            vpuns(ii,l1)=-zion/rr(ii)
         end do
      else
         exit
      end if
   end do

   ! loop over reference plus test atom configurations
   call run_test_configurations(ncnf,nacnf,lacnf,facnf,nc,nvcnf,rho,rhomod,rr,zz, &
                                rcmax,mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
                                lloc,vkb,evkb,srel, &
                                nvt,nat,lat,fat,eat,eatp,eaetst,etsttot)
   fat3(:,:)=fat(:,3,:)
   eat3(:,:)=eat(:,3,:)
   call write_test_configs_text(stdout,ncnf,nc,nvt,nat,lat,fat3,eat3,eatp,etot,eaetst,epstot,etsttot)

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
      &                                       n1, n2, n3, n4)
   write(stdout, '(/a)') 'DATA FOR PLOTTING'
   call write_rho_vpuns_text(stdout, mmax, lmax, drl, nrl, rr, irc, &
      &                          rho, vpuns)
   if (lloc == 4) then
      call write_vloc_text(stdout, mmax, rr, lmax, irc, drl, nrl, &
         &                        vpuns(:, lloc + 1))
   end if
   call write_rho_rhoc_rhom_text(stdout, mmax, rr, lmax, irc, drl, nrl, &
      &                              rho, rhoc, rhomod(:, 1))
   call get_wavefunctions(zz, srel, mmax, rr, vfull, lloc, vp, lmax, &
      &                       irc, nproj, drl, nrl, mxprj, epa, npa, vkb, evkb, &
      &                       uu_ae, up_ae, uu_ps, up_ps, mch_ae, mch_ps, e_ae, e_ps, &
      &                       sign_ae, sign_ps, is_scattering)
   call write_wavefunctions_vkb_text(stdout, mmax, rr, lmax, irc, drl, nrl, &
      &                                  lloc, mxprj, npa, nproj, sign_ae, sign_ps, uu_ae, uu_ps, is_scattering, &
      &                                  vkb)
   call write_convergence_profile_text(stdout, lmax, mxprj, &
      &                                    cvgplt)

   npsh = int(((epsh2 - epsh1) / depsh) - 0.5_dp) + 1
   allocate(epsh(npsh), pshf(npsh, 4), pshp(npsh, 4))
   call run_phsft(lmax,lloc,nproj,epa,epsh1,epsh2,depsh,rxpsh,npsh,vkb,evkb, &
      &               rr,vfull,vp,zz,mmax,mxprj,irc,srel, &
      &               irpsh,rpsh,epsh,pshf,pshp)
   call write_phsft_text(stdout,rpsh,npsh,epsh,pshf,pshp)

   call gnu_script(epa,evkb,lmax,lloc,mxprj,nproj)

#if (defined WITH_HDF5)
   if (do_hdf5) then
      call write_output_hdf5(hdf5_file_id, &
                             zz, nc, nv, mxprj, lmax, lloc, npa, epa, irc, nproj, &
                             mmax, rr, &  ! log radial mesh
                             drl, nrl, &  ! linear radial mesh
                             ncnf, nvt, nat, lat, fat3, eat3, eatp, etot, eaetst, epstot, etsttot, & ! test configurations
                             vfull, vp, vpuns, &  ! potentials
                             rhotae, rho, rhoc, &  ! charge densities
                             icmod, fcfact, rcfact, rhomod, &  ! model core charge density
                             modcore1_ircc, modcore1_iter, &
                             modcore2_ircc, modcore2_a0, modcore2_b0, &
                             modcore3_ircc, modcore3_rmatch, modcore3_rhocmatch, &
                             n_teter_amp, teter_amp_prefacs, teter_amp_params, &
                             n_teter_scale, teter_scale_prefacs, teter_scale_params, &
                             grid_objective, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, &
                             nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, nm_iter, &
                             d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, d2excps_rhom, d2exc_rmse_rhom, &
                             sign_ae, uu_ae, up_ae, mch_ae, e_ae, &  ! all-electron wavefunctions
                             sign_ps, uu_ps, up_ps, mch_ps, e_ps, &  ! pseudo wavefunctions
                             is_scattering, &  ! wavefunction scattering flags
                             vkb, evkb, &  ! KB projectors
                             rpsh, npsh, epsh1, epsh2, depsh, epsh, pshf, pshp, & ! phase shift / log derivative
                             cvgplt)  ! convergence profiles
      call hdf_close_file(hdf5_file_id)
   end if
#endif

   if(trim(psfile)=='psp8' .or. trim(psfile)=='both') then
      call linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
         &             rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
         &             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
         &             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
         &             epsh1,epsh2,depsh,rlmax,psfile)
   end if

   if(trim(psfile)=='upf' .or. trim(psfile)=='both') then
      call upfout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
         &             zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
         &             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
         &             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
         &             epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea)
   end if

   if(trim(psfile)=='psml' .or. trim(psfile)=='both') then
      write(stdout, '(a)') ' calling psmlout'
      call psmlout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
         &             irct, srel, &
         &             zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
         &             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
         &             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact, &
         &             epsh1,epsh2,depsh,rlmax,psfile)
   end if

   stop
end program oncvpsp

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
