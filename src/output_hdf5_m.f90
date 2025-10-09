module output_hdf5_m
   use, intrinsic :: iso_fortran_env, only: dp => real64
   ! High-level HDF5 utilities
   use hdf5_utils_m, only: HID_T, &
      hdf_open_file, hdf_close_file, &
      hdf_create_group, hdf_open_group, hdf_close_group, &
      hdf_create_dataset, hdf_write_dataset, &
      hdf_write_attribute, &
      hdf_set_data_scale, hdf_attach_data_scale, hdf_label_dim
   ! HDF5 datasets
   use hdf5, only: h5dopen_f, h5dclose_f
   implicit none
   private
   public :: write_output_hdf5
contains

subroutine write_output_hdf5(filename, &
                             zz, nc, mxprj, lmax, lloc, npa, epa, irc, nproj, &
                             mmax, rr, &  ! log radial mesh
                             drl, nrl, &  ! linear radial mesh
                             ncnf, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot, & ! test configurations
                             vfull, vp, vpuns, &  ! potentials
                             rho, rhoc, rhomod, &  ! charge densities
                             sign_ae, uu_ae, up_ae, mch_ae, e_ae, &  ! all-electron wavefunctions
                             sign_ps, uu_ps, up_ps, mch_ps, e_ps, &  ! pseudo wavefunctions
                             is_scattering, &  ! wavefunction scattering flags
                             vkb, evkb, &  ! KB projectors
                             rpsh, npsh, epsh1, epsh2, depsh, epsh, pshf, pshp, & ! phase shift / log derivative
                             cvgplt)  ! convergence profiles
   ! Input variables
   !> HDF5 filename
   character(*), intent(in) :: filename
   !> Atomic number
   real(dp), intent(in) :: zz
   !> Number of core states
   integer, intent(in) :: nc
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
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax, 5)
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

   ! HDF5 variables
   !> HDF5 file identifier
   integer(HID_T) :: file_id
   !> Group identifier
   integer(HID_T) :: group_id
   !> Group/dataset name
   character(len=1024) :: name

   ! Create HDF5 file
   call hdf_open_file(file_id, filename, 'REPLACE', 'WRITE')
   ! Test results
   ! write(*,*) 'Writing test results to HDF5 file: ', trim(filename)
   call write_test_results_hdf5(file_id, nc, ncnf, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot)
   ! Logarithmic radial grid
   ! write(*,*) 'Writing pseudopotential data to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'logarithmic_radial_mesh', rr)
   call hdf_set_data_scale(file_id, 'logarithmic_radial_mesh', 'r (a.u.)')
   call hdf_write_attribute(file_id, 'logarithmic_radial_mesh', 'description', 'Logarithmic radial mesh')
   call hdf_write_attribute(file_id, 'logarithmic_radial_mesh', 'units', 'Bohr')
   ! Angular momenta
   ! write(*,*) 'Writing angular momenta to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'angular_momentum', [(ll, ll = 0, lmax)])
   call hdf_set_data_scale(file_id, 'angular_momentum', 'angular_momentum')
   call hdf_write_attribute(file_id, 'angular_momentum', 'description', 'Angular momentum quantum numbers')
   ! Derivative orders (for model core charge)
   ! write(*,*) 'Writing derivative orders to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'derivative_order', [(ii - 1, ii = 1, 5)])
   call hdf_set_data_scale(file_id, 'derivative_order', 'derivative_order')
   call hdf_write_attribute(file_id, 'derivative_order', 'description', &
                            'Derivative orders for model core charge density')
   ! Pseudo valence charge density
   ! write(*,*) 'Writing charge densities to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'ps_valence_charge_density', rho)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'ps_valence_charge_density')
   call hdf_write_attribute(file_id, 'ps_valence_charge_density', 'description', &
                            'Pseudo valence charge density')
   ! Core charge density
   ! write(*,*) 'Writing charge densities to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'ae_core_charge_density', rhoc)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'ae_core_charge_density')
   call hdf_write_attribute(file_id, 'ae_core_charge_density', 'description', &
                            'All-electron core charge density')
   ! Model core charge density
   ! write(*,*) 'Writing charge densities to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'model_core_charge_density', rhomod)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'model_core_charge_density', 1)
   call hdf_attach_data_scale(file_id, 'derivative_order', file_id, 'model_core_charge_density', 2)
   call hdf_write_attribute(file_id, 'model_core_charge_density', 'description', &
                            'Model core charge density and its derivatives')
   ! Unscreened pseudopotentials
   ! write(*,*) 'Writing unscreened pseudopotentials to HDF5 file: ', trim(filename)
   call hdf_create_group(file_id, 'unscreened_pseudotentials')
   call hdf_open_group(file_id, 'unscreened_pseudotentials', group_id)
   do l1 = 1, lmax + 1
      write(name, '(a,i0)') 'l_', l1 - 1
      call hdf_write_dataset(group_id, name, vpuns(:, l1))
      call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', group_id, name)
   end do
   call hdf_close_group(group_id)
   ! Local potential
   ! write(*,*) 'Writing local potential to HDF5 file: ', trim(filename)
   call hdf_write_dataset(file_id, 'local_potential', vpuns(:, lloc + 1))
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'local_potential')
   ! Semi-local pseudopotentials
   ! write(*,*) 'Writing semi-local pseudopotentials to HDF5 file: ', trim(filename)
   call hdf_create_group(file_id, 'semilocal_pseudopotentials')
   call hdf_open_group(file_id, 'semilocal_pseudopotentials', group_id)
   do l1 = 1, lmax + 1
      write(name, '(a,i0)') 'l_', l1 - 1
      call hdf_write_dataset(group_id, name, vp(:, l1))
      call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', group_id, name)
   end do
   call hdf_close_group(group_id)
   ! Wavefunctions
   ! write(*,*) 'Writing wavefunctions to HDF5 file: ', trim(filename)
   call write_wavefunctions_hdf5(file_id, mmax, rr, lmax, &
                                 lloc, mxprj, npa, nproj, &
                                 sign_ae, uu_ae, up_ae, mch_ae, e_ae, &
                                 sign_ps, uu_ps, up_ps, mch_ps, e_ps, &
                                 is_scattering)
   ! VKB projectors
   ! write(*,*) 'Writing Vanderbilt-Kleinman-Bylander projectors to HDF5 file: ', trim(filename)
   call write_vkb_projectors_hdf5(file_id, mmax, mxprj, lmax, vkb, nproj)
   ! Convergence profiles
   ! write(*,*) 'Writing convergence profiles to HDF5 file: ', trim(filename)
   call write_convergence_profiles_hdf5(file_id, cvgplt, mxprj, lmax, nproj)
   ! Log derivative phase shift analysis
   ! write(*,*) 'Writing phase shifts to HDF5 file: ', trim(filename)
   call write_phase_shift_hdf5(file_id, npsh, epsh1, epsh2, depsh, epsh, rpsh, pshf, pshp)
   ! Close HDF5 file
   call hdf_close_file(file_id)

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
   integer, allocatable :: na(:)
   integer, allocatable :: la(:)
   real(dp), allocatable :: fa(:)
   real(dp), allocatable :: ea(:)
   real(dp), allocatable :: eap(:)
   integer(HID_T) :: group_id
   integer(HID_T) :: test_group_id
   integer :: ii
   character(len=1024) :: name

   call hdf_create_group(file_id, 'test_results')
   call hdf_open_group(file_id, 'test_results', group_id)

   do ii = 1, ncnf + 1
      if (ii == 1) then
         write(name, '(a)') 'reference_configuration'
      else
         write(name, '(a,i0)') 'test_configuration_', ii - 1
      end if
      allocate(na(nvt(ii)), la(nvt(ii)), fa(nvt(ii)), ea(nvt(ii)), eap(nvt(ii)))
      na = nat(1:nvt(ii), ii)
      la = lat(1:nvt(ii), ii)
      fa = fat(1:nvt(ii), ii)
      ea = eat(1:nvt(ii), ii)
      eap = eatp(1:nvt(ii), ii)
      call hdf_create_group(group_id, trim(name))
      call hdf_open_group(group_id, trim(name), test_group_id)
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
end subroutine write_test_results_hdf5

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
         call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', iproj_group_id, 'rpsi')
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'description', 'sign * r * psi(r)')
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'units', 'Bohr * Ha^(-1/2)')
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'matching_radius_index', mch(iproj, l1))
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'matching_radius', rr(mch(iproj, l1)))
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'sign', sgn(iproj, l1))
         call hdf_write_attribute(iproj_group_id, 'rpsi', 'eigenvalue', ee(iproj, l1))
         call hdf_write_dataset(iproj_group_id, 'drpsi_dr', up(:, iproj, l1))
         call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', iproj_group_id, 'drpsi_dr')
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
         call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', l_group_id, i_name)
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
   ! Write energy mesh wiith metadata and make it a data scale
   call hdf_write_dataset(psh_group_id, "energy_mesh", epsh)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "epsh1", epsh1)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "epsh2", epsh2)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "depsh", depsh)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "description", "Phase shift energy mesh")
   call hdf_write_attribute(psh_group_id, "energy_mesh", "units", "Ha")
   call hdf_set_data_scale(psh_group_id, "energy_mesh", 'E (Ha)')
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
      call hdf_attach_data_scale(psh_group_id, "energy_mesh", ae_group_id, name)
      ! PS
      call hdf_write_dataset(ps_group_id, name, pshp(:, l1))
      call hdf_write_attribute(ps_group_id, name, 'r', rpsh(l1))
      call hdf_write_attribute(ps_group_id, name, 'angular_momentum', ll)
      call hdf_write_attribute(ps_group_id, name, 'ae_ps', 'ps')
      call hdf_attach_data_scale(psh_group_id, "energy_mesh", ps_group_id, name)
   end do  ! l1

   call hdf_close_group(ps_group_id)
   call hdf_close_group(ae_group_id)
   call hdf_close_group(psh_group_id)
end subroutine write_phase_shift_hdf5

end module output_hdf5_m
