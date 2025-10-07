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

subroutine write_output_hdf5(filename, lmax, npa, epa, lloc, irc, &
                             vkb, evkb, nproj, rr, vfull, vp, vpuns, zz, mmax, mxprj, drl, nrl, &
                             rho, rhoc, rhomod, srel, cvgplt, epsh1, epsh2, depsh, rxpsh)
   ! Input variables
   !> HDF5 filename
   character(*), intent(in) :: filename
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Local potential angular momentum
   integer, intent(in) :: lloc
   !> Number of points in logarithmic radial mesh
   integer, intent(in) :: mmax
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Principal quantum number for corresponding all-electron state
   integer, intent(in) :: npa(mxprj, 6)
   !> Indices of core radii
   integer, intent(in) :: irc(6)
   !> Number of VKB projectors for each l
   integer, intent(in) :: nproj(6)
   !> Atomic number
   real(dp), intent(in) :: zz
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Logarithmic radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Semi-local pseudopotentials
   real(dp), intent(in) :: vp(mmax, 5)
   !> Unscreened pseudopotentials
   real(dp), intent(in) :: vpuns(mmax, 5)
   !> All-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax, 5)
   !> Bound-state or scattering state reference energies for VKB potentials
   real(dp), intent(in) :: epa(mxprj, 6)
   !> Coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4)
   !> Energy per electron error vs. cutoff
   !> cvgplt(1, error-value, projector, ang-mom) = Cutoff energies (Ha)
   !> cvgplt(2, error-value, projector, ang-mom) = [error-value] Energy error per electron (Ha)
   real(dp), intent(in) :: cvgplt(2, 7, mxprj, 4)
   !> Scalar-relativistic flag
   logical, intent(in) :: srel
   !> Minimum energy for log-der analysis
   real(dp), intent(in) :: epsh1
   !> Maximum energy for log-der analysis
   real(dp), intent(in) :: epsh2
   !> Energy spacing for log-der analysis
   real(dp), intent(in) :: depsh
   !> Radius for log-der analysis (if <= 0, use default value)
   real(dp), intent(in) :: rxpsh

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
   ! Logarithmic radial grid
   write(6, '(a)') 'write_output_hdf5: writing dataset logarithmic_radial_mesh'
   call hdf_write_dataset(file_id, 'logarithmic_radial_mesh', rr)
   call hdf_set_data_scale(file_id, 'logarithmic_radial_mesh', 'r (a.u.)')
   ! Angular momenta
   write(6, '(a)') 'write_output_hdf5: writing dataset angular_momentum'
   call hdf_write_dataset(file_id, 'angular_momentum', [(ll, ll = 0, lmax)])
   call hdf_set_data_scale(file_id, 'angular_momentum', 'angular_momentum')
   ! Derivative orders (for model core charge)
   write(6, '(a)') 'write_output_hdf5: writing dataset derivative_order'
   call hdf_write_dataset(file_id, 'derivative_order', [(ii - 1, ii = 1, 5)])
   call hdf_set_data_scale(file_id, 'derivative_order', 'derivative_order')
   ! Valence pseudo charge density
   write(6, '(a)') 'write_output_hdf5: writing dataset ps_valence_charge_density'
   call hdf_write_dataset(file_id, 'ps_valence_charge_density', rho)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'ps_valence_charge_density')
   ! Unscreened pseudopotentials
   write(6, '(a)') 'write_output_hdf5: writing dataset unscreened_pseudopotentials'
   call hdf_create_group(file_id, 'unscreened_pseudotentials')
   call hdf_open_group(file_id, 'unscreened_pseudotentials', group_id)
   do l1 = 1, lmax + 1
      write(name, '(a,i0)') 'l_', l1 - 1
      call hdf_write_dataset(group_id, name, vpuns(:, l1))
      call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', group_id, name)
   end do
   call hdf_close_group(group_id)
   ! Local potential
   write(6, '(a)') 'write_output_hdf5: writing dataset local_potential'
   call hdf_write_dataset(file_id, 'local_potential', vpuns(:, lloc + 1))
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'local_potential')
   ! Semi-local pseudopotentials
   write(6, '(a)') 'write_output_hdf5: writing dataset semilocal_pseudopotentials'
   call hdf_create_group(file_id, 'semilocal_pseudopotentials')
   call hdf_open_group(file_id, 'semilocal_pseudopotentials', group_id)
   do l1 = 1, lmax + 1
      write(name, '(a,i0)') 'l_', l1 - 1
      call hdf_write_dataset(group_id, name, vp(:, l1))
      call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', group_id, name)
   end do
   call hdf_close_group(group_id)
   ! Core charge density
   write(6, '(a)') 'write_output_hdf5: writing dataset ae_core_charge_density'
   call hdf_write_dataset(file_id, 'ae_core_charge_density', rhoc)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'ae_core_charge_density')
   ! Model core charge density
   write(6, '(a)') 'write_output_hdf5: writing dataset model_core_charge_density'
   call hdf_write_dataset(file_id, 'model_core_charge_density', rhomod)
   call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', file_id, 'model_core_charge_density', 1)
   call hdf_attach_data_scale(file_id, 'derivative_order', file_id, 'model_core_charge_density', 2)
   ! Wavefunctions
   call write_wavefunctions_hdf5(file_id, lmax, npa, epa, lloc, vkb, evkb, &
                                 nproj, rr, vfull, vp, zz, mmax, mxprj, drl, nrl, srel)
   ! VKB projectors
   call write_vkb_projectors_hdf5(file_id, mmax, mxprj, lmax, vkb, nproj)
   ! Convergence profiles:
   call write_convergence_profiles_hdf5(file_id, cvgplt, mxprj, lmax, nproj)
   ! Log derivative phase shift analysis
   call write_phase_shift_hdf5(file_id, lmax, lloc, nproj, epa, epsh1, epsh2, &
                               depsh, vkb, evkb, rr, vfull, vp, zz, mmax, mxprj, &
                               irc, rxpsh, srel)
   ! Close HDF5 file
   call hdf_close_file(file_id)

end subroutine write_output_hdf5

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

subroutine write_wavefunctions_hdf5(file_id, lmax, npa, epa, lloc, &
                                    vkb, evkb, nproj, rr, vfull, vp, zz, mmax, mxprj, drl, nrl, &
                                    srel)
   ! Input variables
   integer(HID_T), intent(in) :: file_id
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Local potential angular momentum
   integer, intent(in) :: lloc
   !> Number of points in logarithmic radial mesh
   integer, intent(in) :: mmax
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Principal quantum number for corresponding all-electron state
   integer, intent(in) :: npa(mxprj, 6)
   !> Number of VKB projectors for each l
   integer, intent(in) :: nproj(6)
   !> Atomic number
   real(dp), intent(in) :: zz
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Logarithmic radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Semi-local pseudopotentials
   real(dp), intent(in) :: vp(mmax, 5)
   !> All-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> Bound-state or scattering state reference energies for VKB potentials
   real(dp), intent(in) :: epa(mxprj, 6)
   !> Coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4)
   !> Scalar-relativistic flag
   logical, intent(in) :: srel

   ! Local variables
   !> Angular momentum
   integer :: ll
   !> Angular momentum index (l + 1)
   integer :: l1
   !> Loop index
   integer :: ii
   !> Error flag for lsch* routines
   integer :: ierr
   !> Matching index (dummy)
   integer :: mch
   !> ????
   integer :: n2
   !> Effective principal quantum number for scattering states
   integer :: nn
   !> Projector index
   integer :: iproj
   !> Number of projectors for given l
   integer :: nproj_l
   !> Logarithmic radial mesh spacing (log(amesh))
   real(dp) :: al
   !> Maximum energy for finding pseudo wavefunctions
   real(dp) :: emax
   !> Minimum energy for finding pseudo wavefunctions
   real(dp) :: emin
   !> Test energy for bound or scattering state
   real(dp) :: etest
   !> Sign of all-electron wavefunction at some radial point
   real(dp) :: sgnae
   !> Sign of pseudo wavefunction at some radial point
   real(dp) :: sgnps
   !> All-electron wavefunction
   real(dp), allocatable :: uae(:)
   !> Pseudo wavefunction
   real(dp), allocatable :: ups(:)
   !> Wavefunction derivative (dummy)
   real(dp), allocatable :: up(:)
   logical :: is_bound
   character(len=1024) :: error_message
   ! HDF5 variables
   character(len=1024) :: l_group_name
   character(len=1024) :: iproj_group_name
   integer(HID_T) :: dset_id
   integer(HID_T) :: wfn_group_id
   integer(HID_T) :: ae_wfn_group_id
   integer(HID_T) :: ps_wfn_group_id

   allocate (uae(mmax), ups(mmax), up(mmax))

   al = 0.01_dp * log(rr(101) / rr(1))
   n2 = int(log(dble(nrl) * drl / rr(1)) / al + 1.0_dp)

   call hdf_create_group(file_id, 'wavefunctions')
   call hdf_open_group(file_id, 'wavefunctions', wfn_group_id)
   call hdf_create_group(wfn_group_id, 'all_electron')
   call hdf_open_group(wfn_group_id, 'all_electron', ae_wfn_group_id)
   call hdf_create_group(wfn_group_id, 'pseudo')
   call hdf_open_group(wfn_group_id, 'pseudo', ps_wfn_group_id)
   do l1 = 1, lmax + 1
      ll = l1 - 1
      nproj_l = nproj(l1)
      if (ll == lloc) nproj_l = 0
      l_group_name = ''
      write(l_group_name, '(a,i0)') 'l_', ll
      call hdf_create_group(ae_wfn_group_id, l_group_name)
      call hdf_create_group(ps_wfn_group_id, l_group_name)
      do iproj = 1, nproj(l1)
         if (epa(iproj, l1) < 0.0_dp) then
            is_bound = .true.
            etest = epa(iproj, l1)
            call lschfb(npa(iproj, l1), ll, ierr, etest, rr, vfull, uae, up, zz, mmax, mch, srel)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschfb for l=', ll, ' n=', npa(iproj, l1)
               call hdf_close_group(ps_wfn_group_id)
               call hdf_close_group(ae_wfn_group_id)
               call hdf_close_group(wfn_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            emax = 0.9_dp * etest
            emin = 1.1_dp * etest
            call lschvkbb(ll + iproj, ll, nproj_l, ierr, etest, emin, emax, rr, vp(:, lloc + 1), vkb(:, :, l1), &
                          evkb(1, l1), ups, up, mmax, mch)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschvkbb for l=', ll, ' iprj=', iproj
               call hdf_close_group(ps_wfn_group_id)
               call hdf_close_group(ae_wfn_group_id)
               call hdf_close_group(wfn_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
         else
            is_bound = .false.
            etest = epa(iproj, l1)
            ! AE scattering state
            call lschfs(nn, ll, ierr, etest, rr, vfull, uae, up, zz, mmax, n2, srel)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschfs for l=', ll, ' iprj=', iproj
               call hdf_close_group(ps_wfn_group_id)
               call hdf_close_group(ae_wfn_group_id)
               call hdf_close_group(wfn_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            ! PS scattering state
            call lschvkbs(ll, nproj_l, etest, rr, vp(:, lloc + 1), vkb(:, :, l1), evkb(:, l1), ups, up, mmax, n2)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschvkbs for l=', ll, ' iprj=', iproj
               call hdf_close_group(ps_wfn_group_id)
               call hdf_close_group(ae_wfn_group_id)
               call hdf_close_group(wfn_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            sgnae = sign(1.0_dp, uae(n2))
            sgnps = sign(1.0_dp, ups(n2))
            do ii = 1, mmax
               uae(ii) = sgnae * uae(ii)
               ups(ii) = sgnps * ups(ii)
            end do
         end if
         write(iproj_group_name, '(a,i0)') 'n_', iproj
         ! AE wavefunction
         call hdf_open_group(ae_wfn_group_id, l_group_name, dset_id)
         call hdf_write_dataset(dset_id, iproj_group_name, uae)
         call hdf_write_attribute(dset_id, iproj_group_name, 'angular_momentum', ll)
         call hdf_write_attribute(dset_id, iproj_group_name, 'projector_number', iproj)
         call hdf_write_attribute(dset_id, iproj_group_name, 'ae_ps', 'ae')
         if (is_bound) then
            call hdf_write_attribute(dset_id, iproj_group_name, 'scattering_bound', 'bound')
            call hdf_write_attribute(dset_id, iproj_group_name, 'principal_quantum_number', npa(iproj, l1))
            call hdf_write_attribute(dset_id, iproj_group_name, 'energy', epa(iproj, l1))
         else
            call hdf_write_attribute(dset_id, iproj_group_name, 'scattering_bound', 'scattering')
            call hdf_write_attribute(dset_id, iproj_group_name, 'principal_quantum_number', nn)
            call hdf_write_attribute(dset_id, iproj_group_name, 'energy', etest)
         end if
         call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', dset_id, iproj_group_name)
         call hdf_close_group(dset_id)
         ! PS wavefunction
         call hdf_open_group(ps_wfn_group_id, l_group_name, dset_id)
         call hdf_write_dataset(dset_id, iproj_group_name, ups)
         call hdf_write_attribute(dset_id, iproj_group_name, 'angular_momentum', ll)
         call hdf_write_attribute(dset_id, iproj_group_name, 'projector_number', iproj)
         call hdf_write_attribute(dset_id, iproj_group_name, 'ae_ps', 'ps')
         if (is_bound) then
            call hdf_write_attribute(dset_id, iproj_group_name, 'scattering_bound', 'bound')
            call hdf_write_attribute(dset_id, iproj_group_name, 'principal_quantum_number', npa(iproj, l1))
            call hdf_write_attribute(dset_id, iproj_group_name, 'energy', epa(iproj, l1))
         else
            call hdf_write_attribute(dset_id, iproj_group_name, 'scattering_bound', 'scattering')
            call hdf_write_attribute(dset_id, iproj_group_name, 'principal_quantum_number', nn)
            call hdf_write_attribute(dset_id, iproj_group_name, 'energy', etest)
         end if
         call hdf_attach_data_scale(file_id, 'logarithmic_radial_mesh', dset_id, iproj_group_name)
         call hdf_close_group(dset_id)
      end do  ! iprj
   end do  ! l1
end subroutine write_wavefunctions_hdf5

subroutine write_phase_shift_hdf5(file_id, lmax, lloc, nproj, epa, epsh1, epsh2, depsh, vkb, evkb, &
                                  rr, vfull, vp, zz, mmax, mxprj, irc, rxpsh, srel)
   !Input variables
   !> HDF5 File ID
   integer(HID_T), intent(in) :: file_id
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Local potential angular momentum
   integer, intent(in) :: lloc
   !> Number of points in logarithmic radial mesh
   integer, intent(in) :: mmax
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Indices of core radii
   integer, intent(in) :: irc(6)
   !> Number of VKB projectors for each l
   integer, intent(in) :: nproj(6)
   !> Atomic number
   real(dp), intent(in) :: zz
   !> Logarithmic radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Semi-local pseudopotentials
   real(dp), intent(in) :: vp(mmax, 5)
   !> Bound-state or scattering state reference energies for VKB potentials
   real(dp), intent(in) :: epa(mxprj, 6)
   !> All-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> Coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4)
   !> Scalar-relativistic flag
   logical, intent(in) :: srel
   !> Minimum energy for log-der analysis
   real(dp) :: epsh1
   !> Maximum energy for log-der analysis
   real(dp) :: epsh2
   !> Energy spacing for log-der analysis
   real(dp) :: depsh
   !> Radius for log-der analysis
   real(dp) :: rxpsh

   ! Local variables
   integer :: ii
   integer :: irphs
   integer :: ll
   integer :: l1
   integer :: npsh
   integer :: xirphs
   real(dp), allocatable :: epsh(:), pshf(:), pshp(:)
   character(len=6) :: name
   character(len=1024) :: error_msg
   integer(HID_T) :: psh_group_id, ae_group_id, ps_group_id

   npsh = int(((epsh2 - epsh1) / depsh) - 0.5d0) + 1
   allocate (epsh(npsh), pshf(npsh), pshp(npsh))
   do ii = 1, npsh
      epsh(ii) = epsh2 - depsh * dble(ii - 1)
   end do

   call hdf_create_group(file_id, "log_derivative_phase_shift")
   call hdf_open_group(file_id, "log_derivative_phase_shift", psh_group_id)
   ! Write energy mesh wiith metadata and make it a data scale
   call hdf_write_dataset(psh_group_id, "energy_mesh", epsh)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "start", epsh1)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "stop", epsh2)
   call hdf_write_attribute(psh_group_id, "energy_mesh", "step", depsh)
   call hdf_set_data_scale(psh_group_id, "energy_mesh", 'E (Ha)')
   ! Create and open PS and AE groups
   call hdf_create_group(psh_group_id, "all_electron")
   call hdf_open_group(psh_group_id, "all_electron", ae_group_id)
   call hdf_create_group(psh_group_id, "pseudo")
   call hdf_open_group(psh_group_id, "pseudo", ps_group_id)
   do l1 = 1, 4
      ll = l1 - 1
      if (ll <= lmax) then
         irphs = irc(l1) + 2
      else
         irphs = irc(lloc + 1)
      end if

      ! Check that user-specified radius is at least as large as the default value
      ! If it's too small, error out; if it's larger, use it
      if (rxpsh > 0.d0) then
         xirphs = minloc(abs(rr - rxpsh), dim=1)
         if (xirphs < irphs) then
            write (error_msg, '(a,f6.2)') 'run_phsft: ERROR rxpsh for logder analysis too small (~< rc(l)) : ', rxpsh
            call hdf_close_group(ae_group_id)
            call hdf_close_group(ps_group_id)
            call hdf_close_file(file_id)
            error stop error_msg
         else
            irphs = xirphs
         end if
      end if

      ! Compute the phase shift atan(d(r*ψ(r))/dr / (rψ(r)))
      call fphsft(ll, epsh2, depsh, pshf, rr, vfull, zz, mmax, irphs, npsh, srel)
      if (ll == lloc) then
         call vkbphsft(ll, 0, epsh2, depsh, epa(:, l1), pshf, pshp, &
         &                   rr, vp(:, lloc + 1), vkb(:, :, l1), evkb(:, l1), &
         &                   mmax, irphs, npsh)
      else
         call vkbphsft(ll, nproj(l1), epsh2, depsh, epa(:, l1), pshf, pshp, &
         &                   rr, vp(:, lloc + 1), vkb(:, :, l1), evkb(:, l1), &
         &                   mmax, irphs, npsh)
      end if

      write(name, '(a,i0)') 'l_', ll
      ! AE
      call hdf_write_dataset(ae_group_id, name, pshf)
      call hdf_write_attribute(ae_group_id, name, 'r', rr(irphs))
      call hdf_write_attribute(ae_group_id, name, 'angular_momentum', ll)
      call hdf_write_attribute(ae_group_id, name, 'ae_ps', 'ae')
      call hdf_attach_data_scale(psh_group_id, "energy_mesh", ae_group_id, name)
      ! PS
      call hdf_write_dataset(ps_group_id, name, pshp)
      call hdf_write_attribute(ps_group_id, name, 'r', rr(irphs))
      call hdf_write_attribute(ps_group_id, name, 'angular_momentum', ll)
      call hdf_write_attribute(ps_group_id, name, 'ae_ps', 'ps')
      call hdf_attach_data_scale(psh_group_id, "energy_mesh", ps_group_id, name)
   end do  ! l1

   call hdf_close_group(ps_group_id)
   call hdf_close_group(ae_group_id)
   call hdf_close_group(psh_group_id)
   deallocate (epsh, pshf, pshp)

end subroutine write_phase_shift_hdf5

end module output_hdf5_m
