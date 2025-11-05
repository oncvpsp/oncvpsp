module output_hdf5_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use constants_m, only: MAX_NUM_PROJ, MAX_NUM_ELL
   ! High-level HDF5 utilities
   use hdf5_utils_m, only: HID_T, &
      hdf_open_file, hdf_close_file, &
      hdf_create_group, hdf_open_group, hdf_close_group, &
      hdf_create_dataset, hdf_write_dataset, &
      hdf_write_attribute, &
      hdf_set_dimension_scale, hdf_attach_dimension_scale, hdf_label_dim, &
      hdf_exists
   ! HDF5 datasets
   use hdf5, only: h5dopen_f, h5dclose_f
   use hdf5, only: h5lexists_f
   implicit none
   private
   public :: write_input_hdf5, &
      write_log_mesh_hdf5, &
      write_sratom_hdf5, &
      write_optimize_inputs_hdf5, &
      write_optimize_results_hdf5, &
      write_bound_valence_states_hdf5, &
      write_reference_configuration_results_hdf5, &
      write_output_hdf5, &
      write_teter_optimization_hdf5

   !> HDF5 output file identifier
   integer(HID_T), public :: hdf5_file_id
   !> .true. to enable HDF5 output
   logical, public :: do_hdf5 = .false.
contains

   subroutine write_input_hdf5(file_id, &
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
      implicit none
      ! Input variables
      !> Output unit
      integer(HID_T), intent(in) :: file_id
      !> Atomic symbol
      character(len=2), intent(in) :: atsym
      !> Atomic number
      real(dp), intent(in) :: zz
      !> Number of core states
      integer, intent(in) :: nc
      !> Number of valence states
      integer, intent(in) :: nv
      !> Exchange-correlation functional ID
      integer, intent(in) :: iexc
      !> Pseudopotential output type
      character(len=4), intent(in) :: psfile
      !> Reference configuration rincipal quantum number array
      integer, intent(in) :: na(30)
      !> Reference configuration ngular momentum array
      integer, intent(in) :: la(30)
      !> Reference configuration occupation number array
      real(dp), intent(in) :: fa(30)
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> Core radii for each l
      real(dp), intent(in) :: rc(6)
      !> Reference energies for each l
      real(dp), intent(in) :: ep(6)
      !> Number of constraints for each l
      integer, intent(in) :: ncon(6)
      !> Number of basis functions for each l
      integer, intent(in) :: nbas(6)
      !> qcut for each l
      real(dp), intent(in) :: qcut(6)
      !> l for local potential
      integer, intent(in) :: lloc
      !>
      integer, intent(in) :: lpopt
      !> Local potential shift at origin
      real(dp), intent(in) :: dvloc0
      !> Number of projectors for each l
      integer, intent(in) :: nproj(6)
      !> Energy shift for each l
      real(dp), intent(in) :: debl(6)
      !> Model core charge type
      integer, intent(in) :: icmod
      !> Model core charge scaling factor
      real(dp), intent(in) :: fcfact
      !> Model core charge radial scaling factor
      real(dp), intent(in) :: rcfact
      !> Log derivative analysis energy minimum
      real(dp), intent(in) :: epsh1
      !> Log derivative analysis energy maximum
      real(dp), intent(in) :: epsh2
      !> Log derivative analysis energy step
      real(dp), intent(in) :: depsh
      !> Maximum radius for output grid
      real(dp), intent(in) :: rlmax
      !> Spacing of linear radial mesh
      real(dp), intent(in) :: drl
      !> Number of test configurations
      integer, intent(in) :: ncnf
      !> Number of valence states for each test configuration
      integer, intent(in) :: nvcnf(5)
      !> Principal quantum number array for test configurations
      integer, intent(in) :: nacnf(30,5)
      !> Angular momentum array for test configurations
      integer, intent(in) :: lacnf(30,5)
      !> Occupation number array for test configurations
      real(dp), intent(in) :: facnf(30,5)

      ! Local variables
      integer(HID_T) :: group_id
      integer(HID_T) :: subgroup_id
      integer(HID_T) :: subsubgroup_id
      integer :: psp_l(lmax + 1)
      integer :: ii, jj, l1
      character(len=8) :: name

      call hdf_create_group(file_id, 'input_parameters')
      call hdf_open_group(file_id, 'input_parameters', group_id)

      call hdf_create_group(group_id, 'oncvpsp')
      call hdf_write_attribute(group_id, 'oncvpsp', 'atsym', atsym)
      call hdf_write_attribute(group_id, 'oncvpsp', 'z', zz)
      call hdf_write_attribute(group_id, 'oncvpsp', 'nc', nc)
      call hdf_write_attribute(group_id, 'oncvpsp', 'nv', nv)
      call hdf_write_attribute(group_id, 'oncvpsp', 'iexc', iexc)

      call hdf_create_group(group_id, 'linear_mesh')
      call hdf_write_attribute(group_id, 'linear_mesh', 'rmax', rlmax)
      call hdf_write_attribute(group_id, 'linear_mesh', 'a', drl)

      call hdf_create_group(group_id, 'reference_configuration')
      call hdf_open_group(group_id, 'reference_configuration', subgroup_id)
      call hdf_write_dataset(subgroup_id, 'n', na(1:nc + nv))
      call hdf_write_dataset(subgroup_id, 'l', la(1:nc + nv))
      call hdf_write_dataset(subgroup_id, 'f', fa(1:nc + nv))
      call hdf_close_group(subgroup_id)

      call hdf_create_group(group_id, 'pseudopotentials')
      call hdf_open_group(group_id, 'pseudopotentials', subgroup_id)
      call hdf_write_attribute(group_id, 'pseudopotentials', 'lmax', lmax)
      do ii = 1, lmax + 1
         psp_l(ii) = ii - 1
      end do
      call hdf_write_dataset(subgroup_id, 'l', psp_l)
      call hdf_write_dataset(subgroup_id, 'rc', rc(1:lmax + 1))
      call hdf_write_dataset(subgroup_id, 'ep', ep(1:lmax + 1))
      call hdf_write_dataset(subgroup_id, 'ncon', ncon(1:lmax + 1))
      call hdf_write_dataset(subgroup_id, 'nbas', nbas(1:lmax + 1))
      call hdf_write_dataset(subgroup_id, 'qcut', qcut(1:lmax + 1))
      call hdf_close_group(subgroup_id)

      call hdf_create_group(group_id, 'local_potential')
      call hdf_write_attribute(group_id, 'local_potential', 'lloc', lloc)
      call hdf_write_attribute(group_id, 'local_potential', 'lpopt', lpopt)
      call hdf_write_attribute(group_id, 'local_potential', 'dvloc0', dvloc0)
      call hdf_write_attribute(group_id, 'local_potential', 'rcloc', rc(5))

      call hdf_create_group(group_id, 'vkb_projectors')
      call hdf_open_group(group_id, 'vkb_projectors', subgroup_id)
      call hdf_write_dataset(subgroup_id, 'l', psp_l)
      call hdf_write_dataset(subgroup_id, 'nproj', nproj(1:lmax + 1))
      call hdf_write_dataset(subgroup_id, 'debl', debl(1:lmax + 1))
      call hdf_close_group(subgroup_id)

      call hdf_create_group(group_id, 'model_core_charge')
      call hdf_write_attribute(group_id, 'model_core_charge', 'icmod', icmod)

      call hdf_create_group(group_id, 'log_derivative_analysis')
      call hdf_write_attribute(group_id, 'log_derivative_analysis', 'epsh1', epsh1)
      call hdf_write_attribute(group_id, 'log_derivative_analysis', 'epsh2', epsh2)
      call hdf_write_attribute(group_id, 'log_derivative_analysis', 'depsh', depsh)

      call hdf_create_group(group_id, 'pp_output')
      call hdf_write_attribute(group_id, 'pp_output', 'psfile', psfile)

      call hdf_create_group(group_id, 'test_configurations')
      call hdf_open_group(group_id, 'test_configurations', subgroup_id)
      do ii = 2, ncnf + 1
         write(name, '(i0)') ii - 1
         call hdf_create_group(subgroup_id, trim(name))
         call hdf_open_group(subgroup_id, trim(name), subsubgroup_id)
         call hdf_write_dataset(subsubgroup_id, 'n', nacnf(nc + 1:nc + nvcnf(ii), ii))
         call hdf_write_dataset(subsubgroup_id, 'l', lacnf(nc + 1:nc + nvcnf(ii), ii))
         call hdf_write_dataset(subsubgroup_id, 'f', facnf(nc + 1:nc + nvcnf(ii), ii))
         call hdf_close_group(subsubgroup_id)
      end do
      call hdf_close_group(subgroup_id)

      call hdf_close_group(group_id)
   end subroutine write_input_hdf5

   subroutine write_log_mesh_hdf5(file_id, mmax, rr)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: mmax
      real(dp), intent(in) :: rr(mmax)

      ! Local variables
      real(dp) :: a, b

      a = 0.01_dp * log(rr(101) / rr(1))
      b = rr(1)

      call hdf_write_dataset(file_id, 'log_mesh', rr)
      call hdf_set_dimension_scale(file_id, 'log_mesh', 'r (a.u.)')
      call hdf_write_attribute(file_id, 'log_mesh', 'a', a)
      call hdf_write_attribute(file_id, 'log_mesh', 'b', b)
      call hdf_write_attribute(file_id, 'log_mesh', 'description', 'Logarithmic radial mesh r(i) = b * exp(a * (i - 1))')
      call hdf_write_attribute(file_id, 'log_mesh', 'units', 'Bohr')

   end subroutine write_log_mesh_hdf5

   subroutine write_output_hdf5(file_id, &
      zz, nc, nv, mxprj, lmax, lloc, npa, epa, irc, nproj, &
      mmax, rr, &  ! log radial mesh
      drl, nrl, &  ! linear radial mesh
      ncnf, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot, & ! test configurations
      vfull, vp, vpuns, &  ! potentials
      rhotae, rho, rhoc, &  ! charge densities
      icmod, fcfact, rcfact, rhomod, &  ! model core charge density
      modcore1_ircc, modcore1_iter, &
      modcore2_ircc, modcore2_a0, modcore2_b0, &
      modcore3_ircc, modcore3_rmatch, modcore3_rhocmatch, &
      n_teter_amp, teter_amp_prefacs, teter_amp_params, &
      n_teter_scale, teter_scale_prefacs, teter_scale_params, &
      teter_objective_grid, grid_opt_amp_param, grid_opt_scale_param, grid_opt_objective, &
      nm_opt_amp_param, nm_opt_scale_param, nm_opt_objective, nm_iter, &
      d2excae, d2excps_no_rhom, d2exc_rmse_no_rhom, d2excps_rhom, d2exc_rmse_rhom, &
      sign_ae, uu_ae, up_ae, mch_ae, e_ae, &  ! all-electron wavefunctions
      sign_ps, uu_ps, up_ps, mch_ps, e_ps, &  ! pseudo wavefunctions
      is_scattering, &  ! wavefunction scattering flags
      vkb, evkb, &  ! KB projectors
      rpsh, npsh, epsh1, epsh2, depsh, epsh, pshf, pshp, & ! phase shift / log derivative
      cvgplt)  ! convergence profiles
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Atomic number
      real(dp), intent(in) :: zz
      !> Number of core states
      integer, intent(in) :: nc
      !> Number of valence states
      integer, intent(in) :: nv
      !> Maximum number of projectors
      integer, intent(in) :: mxprj
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> Local potential angular momentum
      integer, intent(in) :: lloc
      !> Principal quantum number for corresponding all-electron state
      integer, intent(in) :: npa(mxprj, 6)
      !> Bound-state or scattering state reference energies for VKB potentials
      real(dp), intent(in) :: epa(mxprj, 6)
      !> Indices of core radii
      integer, intent(in) :: irc(6)
      !> Number of VKB projectors for each l
      integer, intent(in) :: nproj(6)
      !> Number of points in logarithmic radial mesh
      integer, intent(in) :: mmax
      !> Logarithmic radial grid
      real(dp), intent(in) :: rr(mmax)
      !> Spacing of linear radial mesh
      real(dp), intent(in) :: drl
      !> Number of points in linear radial mesh
      integer, intent(in) :: nrl
      !> Number of test configurations
      integer, intent(in) :: ncnf
      !> Number of valence states for each test configuration
      integer, intent(in) :: nvt(5)
      !> Principal quantum numbers of states for each test configuration
      integer, intent(in) :: nat(30,5)
      !> Angular momenta of states for each test configuration
      integer, intent(in) :: lat(30,5)
      !> Occupation numbers of states for each test configuration
      real(dp), intent(in) :: fat(30,5)
      !> Eigenvalues of all-electron states for each test configuration
      real(dp), intent(in) :: eat(30,5)
      !> Eigenvalues of pseudo states for each test configuration
      real(dp), intent(in) :: eatp(30,5)
      !> Reference all-electron total energy
      real(dp), intent(in) :: etot
      !> All-electron total energies for each test configuration
      real(dp), intent(in) :: eaetst(5)
      !> Reference pseudo total energy
      real(dp), intent(in) :: epstot
      !> Pseudopotential total energies for each test configuration
      real(dp), intent(in) :: etsttot(5)
      !> All-electron potential
      real(dp), intent(in) :: vfull(mmax)
      !> Semi-local pseudopotentials
      real(dp), intent(in) :: vp(mmax, 5)
      !> Unscreened pseudopotentials
      real(dp), intent(in) :: vpuns(mmax, 5)
      !> All-electron total charge density
      real(dp), intent(in) :: rhotae(mmax)
      !> Valence pseudocharge
      real(dp), intent(in) :: rho(mmax)
      !> Core charge
      real(dp), intent(in) :: rhoc(mmax)
      !> Model core charge type
      integer, intent(in) :: icmod
      !> Model core charge amplitude factor
      real(dp), intent(in) :: fcfact
      !> Model core charge scale factor
      real(dp), intent(in) :: rcfact
      !> Model core charge
      real(dp), intent(in) :: rhomod(mmax, 5)
      !> icmod = 1 crossover index
      integer, intent(in) :: modcore1_ircc
      !> icmod = 1 number of iterations
      integer, intent(in) :: modcore1_iter
      !> icmod = 2 crossover index
      integer, intent(in) :: modcore2_ircc
      !> icmod = 2 a0 parameter
      real(dp), intent(in) :: modcore2_a0
      !> icmod = 2 b0 parameter
      real(dp), intent(in) :: modcore2_b0
      !> icmod = 3 crossover index
      integer, intent(in) :: modcore3_ircc
      !> icmod = 3 crossover radius
      real(dp), intent(in) :: modcore3_rmatch
      !> icmod = 3 core charge at matching radius
      real(dp), intent(in) :: modcore3_rhocmatch
      !> Number of Teter amplitudes for grid search
      integer, intent(in) :: n_teter_amp
      !> Teter amplitude prefactors for grid search
      real(dp), allocatable, intent(in) :: teter_amp_prefacs(:)
      !> Teter amplitude parameters for grid search
      real(dp), allocatable, intent(in) :: teter_amp_params(:)
      !> Number of Teter scales for grid search
      integer, intent(in) :: n_teter_scale
      !> Teter scale prefactors for grid search
      real(dp), allocatable, intent(in) :: teter_scale_prefacs(:)
      !> Teter scale parameters for grid search
      real(dp), allocatable, intent(in) :: teter_scale_params(:)
      !> d2Exc RMSE on grid of Teter parameters
      real(dp), allocatable, intent(in) :: teter_objective_grid(:, :)
      !> Optimal Teter amplitude parameter from grid search
      real(dp), intent(in) :: grid_opt_amp_param
      !> Optimal Teter scale parameter from grid search
      real(dp), intent(in) :: grid_opt_scale_param
      !> Minimum d2Exc RMSE from grid search
      real(dp), intent(in) :: grid_opt_objective
      !> Optimal Teter amplitude parameter from Nelder-Mead optimization
      real(dp), intent(in) :: nm_opt_amp_param
      !> Optimal Teter scale parameter from Nelder-Mead optimization
      real(dp), intent(in) :: nm_opt_scale_param
      real(dp), intent(in) :: nm_opt_objective
      !> Number of Nelder-Mead iterations
      integer, intent(in) :: nm_iter
      !> All-electron d2Exc
      real(dp), intent(in) :: d2excae(nv, nv)
      !> Pseudopotential d2Exc without model core charge
      real(dp), intent(in) :: d2excps_no_rhom(nv, nv)
      !> d2Exc RMSE AE vs pseudopotential without model core charge
      real(dp), intent(in) :: d2exc_rmse_no_rhom
      !> Pseudopotential d2Exc with model core charge
      real(dp), intent(in) :: d2excps_rhom(nv, nv)
      !> d2Exc RMSE AE vs pseudopotential with model core charge
      real(dp), intent(in) :: d2exc_rmse_rhom
      !> Sign of AE wavefunction at matching point
      real(dp), intent(in) :: sign_ae(mxprj, 4)
      !> AE wavefunction
      real(dp), intent(in) :: uu_ae(mmax, mxprj, 4)
      !> AE wavefunction derivative
      real(dp), intent(in) :: up_ae(mmax, mxprj, 4)
      !> AE wavefunction matching mesh index
      integer, intent(in) :: mch_ae(mxprj, 4)
      !> AE wavefunction energy
      real(dp), intent(in) :: e_ae(mxprj, 4)
      !> Sign of PS wavefunction at matching point
      real(dp), intent(in) :: sign_ps(mxprj, 4)
      !> PS wavefunction
      real(dp), intent(in) :: uu_ps(mmax, mxprj, 4)
      !> PS wavefunction derivative
      real(dp), intent(in) :: up_ps(mmax, mxprj, 4)
      !> PS wavefunction matching mesh index
      integer, intent(in) :: mch_ps(mxprj, 4)
      !> PS wavefunction energy
      real(dp), intent(in) :: e_ps(mxprj, 4)
      !> .true. for scattering states, .false. for bound states
      logical, intent(in) :: is_scattering(mxprj, 4)
      !> VKB projectors
      real(dp), intent(in) :: vkb(mmax, mxprj, 4)
      !> Coefficients of VKB projectors
      real(dp), intent(in) :: evkb(mxprj, 4)
      !> Radius at which phase shift is calculated
      real(dp), intent(in) :: rpsh(4)
      !> Number of phase shift energies
      integer, intent(in) :: npsh
      !> Phase shift energies
      real(dp), intent(in) :: epsh(npsh)
      !> Low energy limit for phase shift calculation
      real(dp), intent(in) :: epsh1
      !> High energy limit for phase shift calculation
      real(dp), intent(in) :: epsh2
      !> Energy increment for phase shift calculation
      real(dp), intent(in) :: depsh
      !> All-electron phase shifts
      real(dp), intent(in) :: pshf(npsh, 4)
      !> Pseudopotential phase shifts
      real(dp), intent(in) :: pshp(npsh, 4)
      !> Energy per electron error vs. cutoff
      !> cvgplt(1, error-value, projector, ang-mom) = Cutoff energies (Ha)
      !> cvgplt(2, error-value, projector, ang-mom) = [error-value] Energy error per electron (Ha)
      real(dp), intent(in) :: cvgplt(2, 7, mxprj, 4)

      ! Local variables
      !> Angular momentum
      integer :: ll
      !> Loop index
      integer :: ii
      !> Angular momentum index (l + 1)
      integer :: l1
      !> Temporary array for radial functions on the logarithmic mesh
      real(dp) :: f_tmp(mmax)

      ! HDF5 variables
      !> Group identifier
      integer(HID_T) :: group_id
      !> Group/dataset name
      character(len=1024) :: name

      ! Test results
      ! call write_test_results_hdf5(file_id, nc, ncnf, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot)
      ! Angular momenta
      call hdf_write_dataset(file_id, 'angular_momentum', [(ll, ll = 0, lmax)])
      call hdf_set_dimension_scale(file_id, 'angular_momentum', 'angular_momentum')
      call hdf_write_attribute(file_id, 'angular_momentum', 'description', 'Angular momentum quantum numbers')
      ! Derivative orders (for model core charge)
      call hdf_write_dataset(file_id, 'derivative_order', [(ii - 1, ii = 1, 5)])
      call hdf_set_dimension_scale(file_id, 'derivative_order', 'derivative_order')
      call hdf_write_attribute(file_id, 'derivative_order', 'description', &
         'Derivative orders for model core charge density')
      ! All-electron core charge density
      f_tmp(:) = rr(:)**2 * rhoc(:)
      call hdf_write_dataset(file_id, 'ae_core_charge_density', f_tmp)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'ae_core_charge_density')
      call hdf_write_attribute(file_id, 'ae_core_charge_density', 'description', &
         'All-electron core charge density')
      ! All-electron valence charge density
      f_tmp(:) = rr(:)**2 * rhotae(:)
      call hdf_write_dataset(file_id, 'ae_valence_charge_density', f_tmp)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'ae_valence_charge_density')
      call hdf_write_attribute(file_id, 'ae_valence_charge_density', 'description', &
         'All-electron valence charge density')
      ! Model core charge density
      call hdf_write_dataset(file_id, 'model_core_charge_density', rhomod(:, 1))
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'model_core_charge_density')
      call hdf_write_attribute(file_id, 'model_core_charge_density', 'description', &
         'Model core charge density')
      ! Pseudo valence charge density
      f_tmp(:) = rr(:)**2 * rho(:)
      call hdf_write_dataset(file_id, 'ps_valence_charge_density', f_tmp)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'ps_valence_charge_density')
      call hdf_write_attribute(file_id, 'ps_valence_charge_density', 'description', &
         'Pseudo valence charge density')
      select case(icmod)
       case(1)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod1_crossover_index', modcore1_ircc)
       case(2)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod2_crossover_index', modcore2_ircc)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod2_a0', modcore2_a0)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod2_b0', modcore2_b0)
       case(3)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_crossover_index', modcore3_ircc)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_rmatch', modcore3_rmatch)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_rhocmatch', modcore3_rhocmatch)
       case(4)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_crossover_index', modcore3_ircc)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_rmatch', modcore3_rmatch)
         call hdf_write_attribute(file_id, 'model_core_charge_density', 'icmod3_rhocmatch', modcore3_rhocmatch)
         call write_teter_optimization_hdf5(file_id, &
            modcore3_rhocmatch, n_teter_amp, teter_amp_prefacs, &
            modcore3_rmatch, n_teter_scale, teter_scale_prefacs, &
            teter_objective_grid, &
            nm_iter, 101, nm_opt_amp_param, nm_opt_scale_param)
       case default
      end select
      ! Local potential
      call hdf_write_dataset(file_id, 'local_potential', vpuns(:, lloc + 1))
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'local_potential')
      ! Screened local potential
      call hdf_write_dataset(file_id, 'screened_local_potential', vp(:, lloc + 1))
      call hdf_attach_dimension_scale(file_id, 'log_mesh', file_id, 'screened_local_potential')
      ! Unscreened pseudopotentials
      call hdf_create_group(file_id, 'unscreened_pseudopotentials')
      call hdf_open_group(file_id, 'unscreened_pseudopotentials', group_id)
      do l1 = 1, lmax + 1
         write(name, '(a,i0)') 'l_', l1 - 1
         call hdf_write_dataset(group_id, name, vpuns(:, l1))
         call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, name)
      end do
      call hdf_close_group(group_id)
      ! Semi-local pseudopotentials
      call hdf_create_group(file_id, 'semilocal_pseudopotentials')
      call hdf_open_group(file_id, 'semilocal_pseudopotentials', group_id)
      do l1 = 1, lmax + 1
         write(name, '(a,i0)') 'l_', l1 - 1
         call hdf_write_dataset(group_id, name, vp(:, l1))
         call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, name)
      end do
      call hdf_close_group(group_id)
      ! Wavefunctions
      call write_wavefunctions_hdf5(file_id, mmax, rr, lmax, &
         lloc, mxprj, npa, nproj, &
         sign_ae, uu_ae, up_ae, mch_ae, e_ae, &
         sign_ps, uu_ps, up_ps, mch_ps, e_ps, &
         is_scattering)
      ! VKB projectors
      call write_vkb_projectors_hdf5(file_id, mmax, mxprj, lmax, vkb, nproj)
      ! Convergence profiles
      call write_convergence_profiles_hdf5(file_id, cvgplt, mxprj, lmax, nproj)
      ! Log derivative phase shift analysis
      call write_phase_shift_hdf5(file_id, npsh, epsh1, epsh2, depsh, epsh, rpsh, pshf, pshp)
   end subroutine write_output_hdf5

   subroutine write_test_results_hdf5(file_id, nc, ncnf, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Number of core states
      integer, intent(in) :: nc
      !> Number of test configurations
      integer, intent(in) :: ncnf
      !> Number of valence states for each test configuration
      integer, intent(in) :: nvt(5)
      !> Principal quantum numbers of states for each test configuration
      integer, intent(in) :: nat(30,5)
      !> Angular momenta of states for each test configuration
      integer, intent(in) :: lat(30,5)
      !> Occupation numbers of states for each test configuration
      real(dp), intent(in) :: fat(30,5)
      !> Eigenvalues of all-electron states for each test configuration
      real(dp), intent(in) :: eat(30,5)
      !> Eigenvalues of pseudo states for each test configuration
      real(dp), intent(in) :: eatp(30,5)
      !> Reference all-electron total energy
      real(dp), intent(in) :: etot
      !> All-electron total energies for each test configuration
      real(dp), intent(in) :: eaetst(5)
      !> Reference pseudo total energy
      real(dp), intent(in) :: epstot
      !> Pseudopotential total energies for each test configuration
      real(dp), intent(in) :: etsttot(5)

      ! Local variables
      integer, allocatable :: ia(:)
      integer, allocatable :: na(:)
      integer, allocatable :: la(:)
      real(dp), allocatable :: fa(:)
      real(dp), allocatable :: ea(:)
      real(dp), allocatable :: eap(:)
      integer(HID_T) :: group_id
      integer(HID_T) :: test_group_id
      integer :: ii, jj
      character(len=1024) :: name

      call hdf_create_group(file_id, 'test_results')
      call hdf_open_group(file_id, 'test_results', group_id)
      do ii = 1, ncnf + 1
         if (ii == 1) then
            write(name, '(a)') 'reference_configuration'
         else
            write(name, '(a,i0)') 'test_configuration_', ii - 1
         end if
         allocate(ia(nvt(ii)), na(nvt(ii)), la(nvt(ii)), fa(nvt(ii)), ea(nvt(ii)), eap(nvt(ii)))
         do jj = 1, nvt(ii)
            ia(jj) = jj
         end do
         na(:) = nat(1:nvt(ii), ii)
         la(:) = lat(1:nvt(ii), ii)
         fa(:) = fat(1:nvt(ii), ii)
         ea(:) = eat(1:nvt(ii), ii)
         eap(:) = eatp(1:nvt(ii), ii)
         call hdf_create_group(group_id, trim(name))
         call hdf_open_group(group_id, trim(name), test_group_id)
         call hdf_write_dataset(test_group_id, 'state_index', ia)
         call hdf_write_dataset(test_group_id, 'principal_quantum_number', na)
         call hdf_write_dataset(test_group_id, 'angular_momentum', la)
         call hdf_write_dataset(test_group_id, 'occupation_number', fa)
         call hdf_write_dataset(test_group_id, 'all_electron_eigenvalue', ea)
         call hdf_write_dataset(test_group_id, 'pseudo_eigenvalue', eap)
         call hdf_write_attribute(test_group_id, '', 'number_of_core_states', nc)
         call hdf_write_attribute(test_group_id, '', 'number_of_valence_states', nvt(ii))
         call hdf_write_attribute(test_group_id, '', 'all_electron_total_energy_diff', eaetst - etot)
         call hdf_write_attribute(test_group_id, '', 'pseudo_total_energy_diff', etsttot - epstot)
         call hdf_write_attribute(test_group_id, '', 'psp_excitation_error', eaetst - etot - etsttot + epstot)
         call hdf_close_group(test_group_id)
         deallocate(na, la, fa, ea, eap)
      end do
      call hdf_close_group(group_id)
   end subroutine write_test_results_hdf5

   subroutine write_reference_configuration_results_hdf5(file_id, ncv, it, itmax, etot, ea)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Number of core + valence states
      integer, intent(in) :: ncv
      !> Iteration number
      integer, intent(in) :: it
      !> Maximum number of iterations
      integer, intent(in) :: itmax
      !> Total energy
      real(dp), intent(in) :: etot
      !> All-electron eigenvalues
      real(dp), intent(in) :: ea(ncv)

      ! Local variables
      integer(HID_T) :: group_id

      call hdf_create_group(file_id, 'reference_configuration_results')
      call hdf_open_group(file_id, 'reference_configuration_results', group_id)
      call hdf_write_attribute(file_id, 'reference_configuration_results', 'iteration', it)
      call hdf_write_attribute(file_id, 'reference_configuration_results', 'max_iterations', itmax)
      if (it < itmax) then
         call hdf_write_attribute(file_id, 'reference_configuration_results', 'converged', 'true')
      else
         call hdf_write_attribute(file_id, 'reference_configuration_results', 'converged', 'false')
      end if
      call hdf_write_attribute(file_id, 'reference_configuration_results', 'total_energy', etot)
      call hdf_write_dataset(group_id, 'eigenvalues', ea)
      call hdf_close_group(group_id)
   end subroutine write_reference_configuration_results_hdf5

   subroutine write_sratom_hdf5(file_id, nc, nv, na, la, fa, mmax, rr, &
      vfull, vxc, rhoc, rhov, uu, up, rpk, ea, etot)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: nc
      integer, intent(in) :: nv
      integer, intent(in) :: na(nc+nv)
      integer, intent(in) :: la(nc+nv)
      real(dp), intent(in) :: fa(nc+nv)
      integer, intent(in) :: mmax
      real(dp), intent(in) :: rr(mmax)
      real(dp), intent(in) :: vfull(mmax)
      real(dp), intent(in) :: vxc(mmax)
      real(dp), intent(in) :: rhoc(mmax)
      real(dp), intent(in) :: rhov(mmax)
      real(dp), intent(in) :: uu(mmax, nc+nv)
      real(dp), intent(in) :: up(mmax, nc+nv)
      real(dp), intent(in) :: rpk(nc+nv)
      real(dp), intent(in) :: ea(nc+nv)
      real(dp), intent(in) :: etot
      ! Local variables
      integer(HID_T) :: group_id
      integer :: i
      integer :: istate(nc+nv)
      real(dp) :: ftmp(mmax)

      call hdf_create_group(file_id, 'ae_atom')
      call hdf_open_group(file_id, 'ae_atom', group_id)
      call hdf_write_attribute(group_id, '', 'nc', nc)
      call hdf_write_attribute(group_id, '', 'nv', nv)
      call hdf_write_attribute(group_id, '', 'etot', etot)

      do i = 1, nc + nv
         istate(i) = i
      end do
      call hdf_write_dataset(group_id, 'state_index', istate)
      call hdf_set_dimension_scale(group_id, 'state_index', 'state')

      call hdf_write_dataset(group_id, 'n', na)
      call hdf_attach_dimension_scale(group_id, 'state_index', group_id, 'n')

      call hdf_write_dataset(group_id, 'l', la)
      call hdf_attach_dimension_scale(group_id, 'state_index', group_id, 'l')

      call hdf_write_dataset(group_id, 'f', fa)
      call hdf_attach_dimension_scale(group_id, 'state_index', group_id, 'f')

      call hdf_write_dataset(group_id, 'e', ea)
      call hdf_attach_dimension_scale(group_id, 'state_index', group_id, 'e')

      call hdf_write_dataset(group_id, 'rpk', rpk)
      call hdf_attach_dimension_scale(group_id, 'state_index', group_id, 'rpk')

      call hdf_write_dataset(group_id, 'vfull', vfull)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, 'vfull')

      call hdf_write_dataset(group_id, 'vxc', vxc)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, 'vxc')

      ftmp(:) = rr(:)**2 * rhoc(:)
      call hdf_write_dataset(group_id, 'rhoc', ftmp)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, 'rhoc')

      ftmp(:) = rr(:)**2 * rhov(:)
      call hdf_write_dataset(group_id, 'rhov', ftmp)
      call hdf_attach_dimension_scale(file_id, 'log_mesh', group_id, 'rhov')

      ! TODO: not sure how to properly attach dimension scales to 2D datasets
      call hdf_write_dataset(group_id, 'uu', uu)
      call hdf_write_dataset(group_id, 'up', up)

      call hdf_close_group(group_id)

   end subroutine write_sratom_hdf5

   subroutine write_optimize_inputs_hdf5(file_id, mxprj, ll, nproj_l, npa, epa, mmax, uua, upa, vr, isboundpa)
      implicit none
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: mxprj
      integer, intent(in) :: ll
      integer, intent(in) :: nproj_l
      integer, intent(in) :: npa(mxprj)
      real(dp), intent(in) :: epa(mxprj)
      integer, intent(in) :: mmax
      real(dp), intent(in) :: uua(mmax, mxprj)
      real(dp), intent(in) :: upa(mmax, mxprj)
      real(dp), intent(in) :: vr(mmax, mxprj)
      logical, intent(in) :: isboundpa(mxprj)

      integer(HID_T) :: group_id
      integer(HID_T) :: subgroup_id
      logical :: group_exists
      integer :: ierr
      character(len=1024) :: subgroup_name
      integer :: i
      integer :: iproj_l(mxprj)
      integer :: isboundpa_int(mxprj)

      if (.not. hdf_exists(file_id, 'optimize_inputs')) then
         call hdf_create_group(file_id, 'optimize_inputs')
      end if
      call hdf_open_group(file_id, 'optimize_inputs', group_id)

      write(subgroup_name, '(a,i0)') 'l_', ll
      call hdf_create_group(group_id, trim(subgroup_name))
      call hdf_open_group(group_id, trim(subgroup_name), subgroup_id)

      call hdf_write_attribute(subgroup_id, '', 'l', ll)
      call hdf_write_attribute(subgroup_id, '', 'nproj_l', nproj_l)

      do i = 1, mxprj
         iproj_l(i) = i
      end do
      call hdf_write_dataset(subgroup_id, 'projector_index', iproj_l)
      call hdf_set_dimension_scale(subgroup_id, 'projector_index', 'projector')

      call hdf_write_dataset(subgroup_id, 'npa', npa)
      call hdf_attach_dimension_scale(subgroup_id, 'projector_index', subgroup_id, 'npa')

      call hdf_write_dataset(subgroup_id, 'epa', epa)
      call hdf_attach_dimension_scale(subgroup_id, 'projector_index', subgroup_id, 'epa')

      ! TODO: not sure how to properly attach dimension scales to 2D datasets
      call hdf_write_dataset(subgroup_id, 'uu', uua)
      call hdf_write_dataset(subgroup_id, 'up', upa)
      call hdf_write_dataset(subgroup_id, 'vr', vr)

      isboundpa_int = 0
      do i = 1, mxprj
         if (isboundpa(i)) isboundpa_int(i) = 1
      end do
      call hdf_write_dataset(subgroup_id, 'is_bound', isboundpa_int)
      call hdf_attach_dimension_scale(subgroup_id, 'projector_index', subgroup_id, 'is_bound')

      call hdf_close_group(subgroup_id)
      call hdf_close_group(group_id)

   end subroutine write_optimize_inputs_hdf5

   subroutine write_optimize_results_hdf5(file_id, lmax, nproj, mmax, rr, ps_rpsi)
      implicit none
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: lmax
      integer, intent(in) :: nproj(6)
      integer, intent(in) :: mmax
      real(dp), intent(in) :: rr(mmax)
      real(dp), intent(in) :: ps_rpsi(mmax, MAX_NUM_PROJ, MAX_NUM_ELL)

      integer(HID_T) :: group_id
      integer(HID_T) :: subgroup_id
      integer :: ll
      integer :: l1
      integer :: iproj
      character(len=1024) :: name

      call hdf_create_group(file_id, 'optimize_results')
      call hdf_open_group(file_id, 'optimize_results', group_id)

      do ll = 0, lmax
         l1 = ll + 1
         write(name, '(a,i0)') 'l_', ll
         call hdf_create_group(group_id, trim(name))
         call hdf_open_group(group_id, trim(name), subgroup_id)

         call hdf_write_attribute(subgroup_id, '', 'l', ll)
         call hdf_write_attribute(subgroup_id, '', 'nproj_l', nproj(l1))
         do iproj = 1, nproj(l1)
            write(name, '(a,i0)') 'i_', iproj
            call hdf_write_dataset(subgroup_id, trim(name), ps_rpsi(:, iproj, l1))
         end do

         call hdf_close_group(subgroup_id)
      end do

      call hdf_close_group(group_id)

   end subroutine write_optimize_results_hdf5

   subroutine write_bound_valence_states_hdf5(file_id, nv, na, la, fa, mmax, &
                                              ae_uu, ae_up, &
                                              ps_uu, ps_up)
      implicit none
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: nv
      integer, intent(in) :: na(nv)
      integer, intent(in) :: la(nv)
      real(dp), intent(in) :: fa(nv)
      integer, intent(in) :: mmax
      real(dp), intent(in) :: ae_uu(mmax, nv)
      real(dp), intent(in) :: ae_up(mmax, nv)
      real(dp), intent(in) :: ps_uu(mmax, nv)
      real(dp), intent(in) :: ps_up(mmax, nv)

      integer(HID_T) :: group_id
      integer(HID_T) :: subgroup_id
      integer :: i
      integer :: istate(nv)

      call hdf_create_group(file_id, 'bound_valence_states')
      call hdf_open_group(file_id, 'bound_valence_states', group_id)
      call hdf_write_attribute(group_id, '', 'nv', nv)

      do i = 1, nv
         istate(i) = i
      end do
      call hdf_write_dataset(group_id, 'state_index', istate)
      call hdf_set_dimension_scale(group_id, 'state_index', 'state')

      call hdf_create_group(group_id, 'all_electron')
      call hdf_open_group(group_id, 'all_electron', subgroup_id)
      call hdf_write_dataset(subgroup_id, 'n', na)
      call hdf_attach_dimension_scale(group_id, 'state_index', subgroup_id, 'n')
      call hdf_write_dataset(subgroup_id, 'l', la)
      call hdf_attach_dimension_scale(group_id, 'state_index', subgroup_id, 'l')
      call hdf_write_dataset(subgroup_id, 'f', fa)
      call hdf_attach_dimension_scale(group_id, 'state_index', subgroup_id, 'f')
      call hdf_write_dataset(subgroup_id, 'uu', ae_uu)
      call hdf_write_dataset(subgroup_id, 'up', ae_up)
      call hdf_close_group(subgroup_id)

      call hdf_create_group(group_id, 'pseudo')
      call hdf_open_group(group_id, 'pseudo', subgroup_id)
      call hdf_write_dataset(subgroup_id, 'l', la)
      call hdf_attach_dimension_scale(group_id, 'state_index', subgroup_id, 'l')
      call hdf_write_dataset(subgroup_id, 'uu', ps_uu)
      call hdf_write_dataset(subgroup_id, 'up', ps_up)
      call hdf_close_group(subgroup_id)

   end subroutine write_bound_valence_states_hdf5

   subroutine write_wavefunctions_hdf5(file_id, mmax, rr, lmax, &
      lloc, mxprj, npa, nproj, &
      sgn_ae, uu_ae, up_ae, mch_ae, e_ae, &
      sgn_ps, uu_ps, up_ps, mch_ps, e_ps, &
      is_scattering)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Number of points in logarithmic radial mesh
      integer, intent(in) :: mmax
      !> Logarithmic radial grid
      real(dp), intent(in) :: rr(mmax)
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> Local potential angular momentum
      integer, intent(in) :: lloc
      !> Maximum number of projectors
      integer, intent(in) :: mxprj
      !> Principal quantum number for corresponding all-electron state
      integer, intent(in) :: npa(mxprj, 6)
      !> Number of VKB projectors for each l
      integer, intent(in) :: nproj(6)
      !> Sign of wavefunction at matching point
      real(dp), intent(in) :: sgn_ae(mxprj, 4)
      !> Wavefunction (r*psi(r))
      real(dp), intent(in) :: uu_ae(mmax, mxprj, 4)
      !> Wavefunction derivative (d(r*psi(r))/dr)
      real(dp), intent(in) :: up_ae(mmax, mxprj, 4)
      !> Wavefunction matching mesh index
      integer, intent(in) :: mch_ae(mxprj, 4)
      !> Eigenvalue
      real(dp), intent(in) :: e_ae(mxprj, 4)
      !> Sign of wavefunction at matching point
      real(dp), intent(in) :: sgn_ps(mxprj, 4)
      !> Wavefunction (r*psi(r))
      real(dp), intent(in) :: uu_ps(mmax, mxprj, 4)
      !> Wavefunction derivative (d(r*psi(r))/dr)
      real(dp), intent(in) :: up_ps(mmax, mxprj, 4)
      !> Wavefunction matching mesh index
      integer, intent(in) :: mch_ps(mxprj, 4)
      !> Eigenvalue
      real(dp), intent(in) :: e_ps(mxprj, 4)
      !> .true. for scattering states, .false. for bound states
      logical, intent(in) :: is_scattering(mxprj, 4)

      ! Local variables
      !> Group identifier
      integer(HID_T) :: group_id

      call hdf_create_group(file_id, 'wavefunctions')
      call hdf_open_group(file_id, 'wavefunctions', group_id)
      call write_wavefunctions_set_hdf5(file_id, group_id, mmax, rr, lmax, &
         lloc, mxprj, npa, nproj, sgn_ae, uu_ae, up_ae, &
         mch_ae, e_ae, is_scattering, .true.)
      call write_wavefunctions_set_hdf5(file_id, group_id, mmax, rr, lmax, &
         lloc, mxprj, npa, nproj, sgn_ps, uu_ps, up_ps, &
         mch_ps, e_ps, is_scattering, .false.)
      call hdf_close_group(group_id)
   end subroutine write_wavefunctions_hdf5

   subroutine write_wavefunctions_set_hdf5(file_id, group_id, mmax, rr, lmax, &
      lloc, mxprj, npa, nproj, sgn, uu, up, &
      mch, ee, is_scattering, is_ae)
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      integer(HID_T), intent(in) :: group_id
      !> Number of points in logarithmic radial mesh
      integer, intent(in) :: mmax
      !> Logarithmic radial grid
      real(dp), intent(in) :: rr(mmax)
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> Local potential angular momentum
      integer, intent(in) :: lloc
      !> Maximum number of projectors
      integer, intent(in) :: mxprj
      !> Principal quantum number for corresponding all-electron state
      integer, intent(in) :: npa(mxprj, 6)
      !> Number of VKB projectors for each l
      integer, intent(in) :: nproj(6)
      !> Sign of wavefunction at matching point
      real(dp), intent(in) :: sgn(mxprj, 4)
      !> Wavefunction (r*psi(r))
      real(dp), intent(in) :: uu(mmax, mxprj, 4)
      !> Wavefunction derivative (d(r*psi(r))/dr)
      real(dp), intent(in) :: up(mmax, mxprj, 4)
      !> Wavefunction matching mesh index
      integer, intent(in) :: mch(mxprj, 4)
      !> Eigenvalue
      real(dp), intent(in) :: ee(mxprj, 4)
      !> .true. for scattering states, .false. for bound states
      logical, intent(in) :: is_scattering(mxprj, 4)
      !> .true. for all-electron, .false. for pseudo wavefunctions
      logical, intent(in) :: is_ae

      ! Local variables
      !> Angular momentum
      integer :: ll
      !> Angular momentum index (l + 1)
      integer :: l1
      !> Projector index
      integer :: iproj
      !>
      real(dp) :: tmp(mmax)
      ! HDF5 variables
      integer(HID_T) :: ae_ps_group_id
      integer(HID_T) :: l_group_id
      integer(HID_T) :: iproj_group_id
      character(len=12) :: ae_ps_name
      character(len=2) :: ae_ps_short_name
      character(len=10) :: scattering_name
      character(len=1024) :: l_group_name
      character(len=1024) :: iproj_group_name

      if (is_ae) then
         ae_ps_name = 'all_electron'
         ae_ps_short_name = 'ae'
      else
         ae_ps_name = 'pseudo'
         ae_ps_short_name = 'ps'
      end if

      call hdf_create_group(group_id, ae_ps_name)  ! '[all_electron,pseudo]/'
      call hdf_open_group(group_id, ae_ps_name, ae_ps_group_id)  ! '[all_electron,pseudo]/'
      do l1 = 1, lmax + 1
         ll = l1 - 1
         l_group_name = ''
         write(l_group_name, '(a,i0)') 'l_', ll
         call hdf_create_group(ae_ps_group_id, l_group_name) ! '[all_electron,pseudo]/l_*/'
         call hdf_open_group(ae_ps_group_id, l_group_name, l_group_id) ! '[all_electron,pseudo]/l_*/'
         do iproj = 1, nproj(l1)
            if (is_scattering(iproj, l1)) then
               scattering_name = 'scattering'
            else
               scattering_name = 'bound'
            end if
            write(iproj_group_name, '(a,i0)') 'i_', iproj
            call hdf_create_group(l_group_id, iproj_group_name) ! '[all_electron,pseudo]/l_*/i_*/'
            call hdf_open_group(l_group_id, iproj_group_name, iproj_group_id) ! '[all_electron,pseudo]/l_*/i_*/'
            call hdf_write_attribute(iproj_group_id, '', 'principal_quantum_number', npa(iproj, l1))
            call hdf_write_attribute(iproj_group_id, '', 'angular_momentum', ll)
            call hdf_write_attribute(iproj_group_id, '', 'projector_index', iproj)
            call hdf_write_attribute(iproj_group_id, '', 'bound_scattering', scattering_name)
            call hdf_write_attribute(iproj_group_id, '', 'ae_ps', ae_ps_short_name)
            tmp(:) = sgn(iproj, l1) * uu(:, iproj, l1)
            call hdf_write_dataset(iproj_group_id, 'rpsi', tmp)
            call hdf_attach_dimension_scale(file_id, 'log_mesh', iproj_group_id, 'rpsi')
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'description', 'sign * r * psi(r)')
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'units', 'Bohr * Ha^(-1/2)')
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'matching_radius_index', mch(iproj, l1))
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'matching_radius', rr(mch(iproj, l1)))
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'sign', sgn(iproj, l1))
            call hdf_write_attribute(iproj_group_id, 'rpsi', 'eigenvalue', ee(iproj, l1))
            call hdf_write_dataset(iproj_group_id, 'drpsi_dr', up(:, iproj, l1))
            call hdf_attach_dimension_scale(file_id, 'log_mesh', iproj_group_id, 'drpsi_dr')
            call hdf_write_attribute(iproj_group_id, 'drpsi_dr', 'description', 'd(r * psi(r))/dr')
            call hdf_write_attribute(iproj_group_id, 'drpsi_dr', 'units', 'Ha^(-1/2)')
            call hdf_close_group(iproj_group_id) ! '[all_electron,pseudo]/l_*/i_*/'
         end do  ! iprj
         call hdf_close_group(l_group_id) !'[all_electron,pseudo]/l_*/'
      end do  ! l1
      call hdf_close_group(ae_ps_group_id) ! '[all_electron,pseudo]/'
   end subroutine write_wavefunctions_set_hdf5

   subroutine write_vkb_projectors_hdf5(file_id, mmax, mxprj, lmax, vkb, nproj)
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Number of points in logarithmic radial mesh
      integer, intent(in) :: mmax
      !> Maximum number of projectors
      integer, intent(in) :: mxprj
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> VKB projectors
      real(dp), intent(in) :: vkb(mmax, mxprj, 4)
      !> Number of VKB projectors for each l
      integer, intent(in) :: nproj(6)

      ! Local variables
      !> Angular momentum
      integer :: ll
      !> Angular momentum index (l + 1)
      integer :: l1
      !> Projector index
      integer :: iproj
      !> VKB projector group identifier
      integer(HID_T) :: group_id
      !> Angular momentum group identifier
      integer(HID_T) :: l_group_id
      ! Group/dataset names
      character(len=1024) :: l_name
      character(len=1024) :: i_name

      call hdf_create_group(file_id, 'vkb_projectors')
      call hdf_open_group(file_id, 'vkb_projectors', group_id)
      do l1 = 1, lmax + 1
         ll = l1 - 1
         write(l_name, '(a,i0)') 'l_', ll
         call hdf_create_group(group_id, l_name)
         call hdf_open_group(group_id, l_name, l_group_id)
         do iproj = 1, nproj(l1)
            write(i_name, '(a,i0)') 'i_', iproj
            call hdf_write_dataset(l_group_id, i_name, vkb(:, iproj, l1))
            call hdf_attach_dimension_scale(file_id, 'log_mesh', l_group_id, i_name)
         end do
         call hdf_close_group(l_group_id)
      end do
      call hdf_close_group(group_id)

   end subroutine write_vkb_projectors_hdf5

   subroutine write_convergence_profiles_hdf5(file_id, cvgplt, mxprj, lmax, nproj)
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Maximum angular momentum
      integer, intent(in) :: lmax
      !> Maximum number of projectors
      integer, intent(in) :: mxprj
      !> Number of projectors at each angular momentum
      integer, intent(in) :: nproj(6)
      !> Energy per electron error vs. cutoff
      !> cvgplt(1, error-value, projector, ang-mom) = Cutoff energies (Ha)
      !> cvgplt(2, error-value, projector, ang-mom) = [error-value] Energy error per electron (Ha)
      real(dp), intent(in) :: cvgplt(2, 7, mxprj, 4)

      ! Local variables
      !> Angular momentum index
      integer :: l1
      !> Angular momentum
      integer :: ll
      !> Projector index
      integer :: iproj
      !> Group identifiers
      integer(HID_T) :: group_id, l_group_id, i_group_id
      !> Group names
      character(len=4) :: l_name, i_name
      !> Temporary buffer for writing
      real(dp) :: buffer(7)

      call hdf_create_group(file_id, 'convergence_profiles')
      call hdf_open_group(file_id, 'convergence_profiles', group_id)

      do l1 = 1,lmax + 1
         ll = l1 - 1
         write (l_name, '(a,i0)') 'l_', ll
         call hdf_create_group(group_id, l_name)
         call hdf_open_group(group_id, l_name, l_group_id)
         do iproj = 1, nproj(l1)
            write (i_name, '(a,i0)') 'i_', iproj
            call hdf_create_group(l_group_id, i_name)
            call hdf_open_group(l_group_id, i_name, i_group_id)
            buffer(:) = cvgplt(1, 1:7, iproj, l1)
            call hdf_write_dataset(i_group_id, 'cutoff_energy', buffer)
            call hdf_write_attribute(i_group_id, 'cutoff_energy', 'units', 'Ha')
            buffer(:) = cvgplt(2, 1:7, iproj, l1)
            call hdf_write_dataset(i_group_id, 'energy_error_per_electron', buffer)
            call hdf_write_attribute(i_group_id, 'energy_error_per_electron', 'units', 'Ha')
            call hdf_close_group(i_group_id)
         end do
         call hdf_close_group(l_group_id)
      end do
      call hdf_close_group(group_id)
   end subroutine write_convergence_profiles_hdf5

   subroutine write_phase_shift_hdf5(file_id, npsh, epsh1, epsh2, depsh, epsh, rpsh, pshf, pshp)
      !Input variables
      !> HDF5 File ID
      integer(HID_T), intent(in) :: file_id
      !> Number of phase shift energies
      integer, intent(in) :: npsh
      !> Start of phase shift energy mesh
      real(dp), intent(in) :: epsh1
      !> End of phase shift energy mesh
      real(dp), intent(in) :: epsh2
      !> Step of phase shift energy mesh
      real(dp), intent(in) :: depsh
      !> Phase shift energies
      real(dp), intent(in) :: epsh(npsh)
      !> Radius at which phase shift is calculated
      real(dp), intent(in) :: rpsh(4)
      !> All-electron phase shifts
      real(dp), intent(in) :: pshf(npsh,4)
      !> Pseudopotential phase shifts
      real(dp), intent(in) :: pshp(npsh,4)

      ! Local variables
      integer :: ii
      integer :: ll
      integer :: l1
      character(len=6) :: name
      character(len=1024) :: error_msg
      integer(HID_T) :: psh_group_id, ae_group_id, ps_group_id

      call hdf_create_group(file_id, "log_derivative_phase_shift")
      call hdf_open_group(file_id, "log_derivative_phase_shift", psh_group_id)
      ! Write energy mesh with metadata and make it a data scale
      call hdf_write_dataset(file_id, "log_derivative_energy_mesh", epsh)
      call hdf_write_attribute(file_id, "log_derivative_energy_mesh", "epsh1", epsh1)
      call hdf_write_attribute(file_id, "log_derivative_energy_mesh", "epsh2", epsh2)
      call hdf_write_attribute(file_id, "log_derivative_energy_mesh", "depsh", depsh)
      call hdf_write_attribute(file_id, "log_derivative_energy_mesh", "description", "Phase shift energy mesh")
      call hdf_write_attribute(file_id, "log_derivative_energy_mesh", "units", "Ha")
      call hdf_set_dimension_scale(file_id, "log_derivative_energy_mesh", 'E (Ha)')
      ! Create and open PS and AE groups
      call hdf_create_group(psh_group_id, "all_electron")
      call hdf_open_group(psh_group_id, "all_electron", ae_group_id)
      call hdf_create_group(psh_group_id, "pseudo")
      call hdf_open_group(psh_group_id, "pseudo", ps_group_id)
      do l1 = 1, 4
         ll = l1 - 1
         write(name, '(a,i0)') 'l_', ll
         ! AE
         call hdf_write_dataset(ae_group_id, name, pshf(:, l1))
         call hdf_write_attribute(ae_group_id, name, 'r', rpsh(l1))
         call hdf_write_attribute(ae_group_id, name, 'angular_momentum', ll)
         call hdf_write_attribute(ae_group_id, name, 'ae_ps', 'ae')
         call hdf_attach_dimension_scale(file_id, "log_derivative_energy_mesh", ae_group_id, name)
         ! PS
         call hdf_write_dataset(ps_group_id, name, pshp(:, l1))
         call hdf_write_attribute(ps_group_id, name, 'r', rpsh(l1))
         call hdf_write_attribute(ps_group_id, name, 'angular_momentum', ll)
         call hdf_write_attribute(ps_group_id, name, 'ae_ps', 'ps')
         call hdf_attach_dimension_scale(file_id, "log_derivative_energy_mesh", ps_group_id, name)
      end do  ! l1

      call hdf_close_group(ps_group_id)
      call hdf_close_group(ae_group_id)
      call hdf_close_group(psh_group_id)
   end subroutine write_phase_shift_hdf5

   subroutine write_teter_optimization_hdf5(file_id, &
      rhocmatch, n_amp, amp_prefacs, &
      rmatch, n_scale, scale_prefacs, &
      teter_objective_grid, &
      iter, max_iter, amp_param, scale_param)
      implicit none
      ! Input variables
      integer(HID_T), intent(in) :: file_id
      !> Core radius used for matching
      real(dp), intent(in) :: rhocmatch
      !> Number of Teter amplitude prefactors tested
      integer, intent(in) :: n_amp
      !> Teter amplitude prefactors tested
      real(dp), allocatable, intent(in) :: amp_prefacs(:)
      !> Matching radius used for testing
      real(dp), intent(in) :: rmatch
      !> Number of Teter scale prefactors tested
      integer, intent(in) :: n_scale
      !> Teter scale prefactors tested
      real(dp), allocatable, intent(in) :: scale_prefacs(:)
      !> RMSE of d2Exc/drho2 at matching radius for each combination of
      !> Teter amplitude and scale prefactors
      real(dp), allocatable, intent(in) :: teter_objective_grid(:, :)
      !> Nelder-Mead iteration number
      integer, intent(in) :: iter
      !> Maximum number of Nelder-Mead iterations
      integer, intent(in) :: max_iter
      !> Optimal amplitude parameter
      real(dp), intent(in) :: amp_param
      !> Optimal scaling parameter
      real(dp), intent(in) :: scale_param

      ! Local variables
      integer(HID_T) :: group_id
      integer(HID_T) :: subgroup_id
      real(dp) :: amp_params(n_amp)
      real(dp) :: scale_params(n_scale)
      integer :: i

      if (.not. allocated(amp_prefacs)) error stop &
         'write_teter_optimization_hdf5: ERROR amp_prefacs not allocated'
      if (.not. allocated(scale_prefacs)) error stop &
         'write_teter_optimization_hdf5: ERROR scale_prefacs not allocated'
      if (.not. allocated(teter_objective_grid)) error stop &
         'write_teter_optimization_hdf5: ERROR teter_objective_grid not allocated'

      do i = 1, n_amp
         amp_params(i) = amp_prefacs(i) * rhocmatch
      end do
      do i = 1, n_scale
         scale_params(i) = scale_prefacs(i) * rmatch
      end do

      call hdf_create_group(file_id, 'teter_parameter_optimization')
      call hdf_write_attribute(file_id, 'teter_parameter_optimization', 'rhocmatch', rhocmatch)
      call hdf_write_attribute(file_id, 'teter_parameter_optimization', 'rmatch', rmatch)
      call hdf_open_group(file_id, 'teter_parameter_optimization', group_id)

      call hdf_create_group(group_id, 'grid_search')
      call hdf_open_group(group_id, 'grid_search', subgroup_id)
      call hdf_write_dataset(subgroup_id, 'amplitude_parameter', amp_params)
      call hdf_set_dimension_scale(subgroup_id, 'amplitude_parameter', 'amplitude_parameter')
      call hdf_write_dataset(subgroup_id, 'amplitude_prefactor', amp_prefacs)
      call hdf_write_dataset(subgroup_id, 'scale_parameter', scale_params)
      call hdf_set_dimension_scale(subgroup_id, 'scale_parameter', 'scale_parameter')
      call hdf_write_dataset(subgroup_id, 'scale_prefactor', scale_prefacs)
      call hdf_write_dataset(subgroup_id, 'teter_objective_grid', teter_objective_grid)
      call hdf_attach_dimension_scale(subgroup_id, 'amplitude_parameter', subgroup_id, 'teter_objective_grid', 0)
      call hdf_attach_dimension_scale(subgroup_id, 'scale_parameter', subgroup_id, 'teter_objective_grid', 1)
      call hdf_close_group(subgroup_id)

      call hdf_create_group(group_id, 'nelder_mead')
      call hdf_write_attribute(group_id, 'nelder_mead', 'iteration', iter)
      call hdf_write_attribute(group_id, 'nelder_mead', 'max_iterations', max_iter)
      if (iter > max_iter) then
         call hdf_write_attribute(group_id, 'nelder_mead', 'converged', 'false')
      else
         call hdf_write_attribute(group_id, 'nelder_mead', 'converged', 'true')
      end if
      call hdf_write_attribute(group_id, 'nelder_mead', 'amplitude_parameter', amp_param)
      call hdf_write_attribute(group_id, 'nelder_mead', 'amplitude_prefactor', amp_param / rhocmatch)
      call hdf_write_attribute(group_id, 'nelder_mead', 'scale_parameter', scale_param)
      call hdf_write_attribute(group_id, 'nelder_mead', 'scale_prefactor', scale_param / rmatch)

      call hdf_close_group(group_id)
   end subroutine write_teter_optimization_hdf5

end module output_hdf5_m
