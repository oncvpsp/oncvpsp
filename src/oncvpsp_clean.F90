program oncvpsp
   ! Intrinsic modules
   use, intrinsic :: iso_fortran_env, only: stdin => input_unit, stdout => output_unit, stderr => error_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   ! External modules
   use hdf5_utils_m, only: HID_T, hdf_open_file, hdf_close_file
   use output_hdf5_m, only: do_hdf5, hdf5_file_id
   use output_hdf5_m, only: write_input_hdf5
   use output_hdf5_m, only: write_reference_configuration_results_hdf5
   use output_hdf5_m, only: write_output_hdf5
   ! Internal modules
   use input_toml_m, only: read_input_toml
   implicit none

#if (RELATIVISTIC == 0)
   logical, parameter :: scalar_relativistic = .false.
#else
   logical, parameter :: scalar_relativistic = .true.
#endif

   ! Atom definition variables
   character(len=2) :: atomic_symbol
   !> [input] Atomic charge
   real(dp) :: z
   !> [input] Number of core states
   integer :: n_core
   !> [input] Number of valence states
   integer :: n_valence
   !> Total number of states (core + valence)
   integer :: n_states

   ! Exchange-correlation variables
   !> [input] Exchange-correlation functional type
   integer :: exc_type

   ! Pseudopotential output variables
   !> [input] Pseudopotential file output format
   character(len=4) :: psp_format
   !> [input] Linear radial mesh maximum
   real(dp) :: lin_mesh_r_max
   !> [input] Linear radial mesh step
   real(dp) :: lin_mesh_dr
   !> Linear radial mesh size
   integer :: lin_mesh_size
   !> Linear radial mesh array
   real(dp), allocatable :: lin_mesh_r(:)

   ! All-electron reference configuration variables
   !> [input] Reference configuration principal quantum numbers
   integer, allocatable :: ae_n(:)
   !> [input] Reference configuration angular momenta
   integer, allocatable :: ae_l(:)
   !> [input] Reference configuration occupation numbers
   real(dp), allocatable :: ae_f(:)
   !> Reference configuration all-electron total potential
   real(dp), allocatable :: ae_v_tot(:)
   !> Reference configuration all-electron exchange-correlation potential
   real(dp), allocatable :: ae_v_xc(:)
   !> Reference configuration core charge density
   real(dp), allocatable :: ae_rho_core(:)
   !> Reference configuration valence charge density
   real(dp), allocatable :: ae_rho_val(:)
   !> Reference configuration wavefunctions (r*psi(r))
   real(dp), allocatable :: ae_rpsi(:, :)
   !> Reference configuration wavefunction derivatives (d(r*psi)/dr)
   real(dp), allocatable :: ae_drpsi(:, :)
   !> Reference configuration eigenvalues
   real(dp), allocatable :: ae_e(:)
   !> Reference configuration wavefunction outermost peak radii
   real(dp), allocatable :: ae_wfn_r_peak(:)
   !> Reference configuration total energy
   real(dp) :: ae_e_tot
   !> Number of iterations taken to solve the all-electron atom
   integer :: ae_niter
   !> Overlaps between all-electron wavefunctions
   real(dp), allocatable :: ae_ovlp(:, :)

   ! Logarithmic radial mesh variables [r(i) = b * exp(a * (i - 1))]
   !> Logarithmic radial mesh parameter a
   real(dp), parameter :: log_mesh_a = log(1.006_dp)
   !> Logarithmic radial mesh parameter b
   real(dp) :: log_mesh_b
   !> Logarithmic radial mesh maximum (approximate)
   real(dp), parameter :: log_mesh_r_max = 45.0_dp
   !> Logarithmic radial mesh size
   integer :: log_mesh_size
   !> Logarithmic radial mesh array
   real(dp), allocatable :: log_mesh_r(:)

   ! Pseudopotential variables
   !> [input] Maximum angular momentum for pseudopotential generation
   integer :: psp_lmax
   !> [input] Core radii for each angular momentum channel
   real(dp), allocatable :: psp_r_core(:)
   !> [input] Input (guess) energies for each angular momentum channel
   real(dp), allocatable :: psp_e_in(:)
   !> [input] Number of basis functions for each angular momentum channel
   integer, allocatable :: psp_n_basis(:)
   !> [input] Number of matching constraints for each angular momentum channel
   integer, allocatable :: psp_n_constr(:)
   !> [input] Cutoff wavevectors for each angular momentum channel
   real(dp), allocatable :: psp_q_cut(:)
   !> Principal quantum numbers of the highest core states for each angular momentum channel
   integer, allocatable :: psp_n(:, :)
   !> Semi-local pseudopotentials for each angular momentum channel
   real(dp), allocatable :: psp_v_sl(:, :)
   !> Local part of the KB pseudopotential
   real(dp), allocatable :: psp_v_loc(:)
   !> KB projectors
   real(dp), allocatable :: psp_p_kb(:, :)
   !> KB coupling coefficients (energies)
   real(dp), allocatable :: psp_e_kb(:)

   ! Kleinman-Bylander variables
   !> [input] Number of KB projectors for each angular momentum channel
   integer, allocatable :: kb_n_proj(:)
   !> [input] Energy shifts for each angular momentum channel
   real(dp), allocatable :: db_de(:)
   !> KB state principal quantum numbers
   integer, allocatable :: kb_state_n(:, :)
   !> KB state energies
   real(dp), allocatable :: kb_state_e(:, :)
   !> KB state inward-outward matching mesh indices
   integer, allocatable :: kb_state_match_idx(:, :)
   !> KB state wavefunctions (r*psi(r))
   real(dp), allocatable :: kb_state_rpsi(:, :, :)
   !> KB state wavefunction derivatives (d(r*psi)/dr)
   real(dp), allocatable :: kb_state_drpsi(:, :, :)
   !> Relativistic corrections to the potential for each angular momentum channel
   real(dp), allocatable :: kb_rel_corr(:, :)
   !> Kleinman-Bylander projectors
   real(dp), allocatable :: kb_p_kb(:, :)
   !> Kleinman-Bylander coupling coefficients (energies)
   real(dp), allocatable :: kb_e_kb(:)

   ! Local variables
   !> Principal quantum number
   integer :: enn
   !> Angular momentum
   integer :: ell

   ! Iteration indices
   !> Iteration index
   integer :: i
   !> Logarithmic radial mesh index
   integer :: i_r
   !> Angular momentum index
   integer :: i_l
   !> State index
   integer :: i_state
   !> Projector index
   integer :: i_proj
   !> Projector index
   integer :: j_proj
   !> Number of projectors at angular momentum l
   integer :: n_proj_l
   !> Error flag
   integer :: err_flag
   !> Error message
   character(:), allocatable :: err_msg

   ! Read input

   ! Set up logarithmic radial mesh
   log_mesh_b = min(5.0e-4 / z, 5.0e-4 / 10)
   log_mesh_size = int(log(log_mesh_r_max / log_mesh_b) / log_mesh_a)
   allocate(log_mesh_r(log_mesh_size))
   do i_r = 1, log_mesh_size
      log_mesh_r(i_r) = log_mesh_b * exp(log_mesh_a * real(i_r - 1, dp))
   end do

   ! Solve all-electron atom
   allocate(ae_v_tot(log_mesh_size))
   allocate(ae_v_xc(log_mesh_size))
   allocate(ae_rho_core(log_mesh_size))
   allocate(ae_rho_val(log_mesh_size))
   allocate(ae_rpsi(log_mesh_size, n_core + n_valence))
   allocate(ae_drpsi(log_mesh_size, n_core + n_valence))
   allocate(ae_e(n_core + n_valence))
   allocate(ae_wfn_r_peak(n_core + n_valence))
   call solve_ae_atom( &
   ! Inputs
                       z, scalar_relativistic, exc_type, &
                       n_core, n_valence, &
                       ae_n, ae_l, ae_f, &
                       log_mesh_size, log_mesh_r, &
   ! Outputs
                       ae_v_tot, ae_v_xc, &
                       ae_rho_core, ae_rho_val, &
                       ae_rpsi, ae_drpsi, ae_e, ae_wfn_r_peak, &
                       ae_e_tot, &
                       ae_niter, err_flag, err_msg
   )

   ! Construct optimized semi-local pseudopotentials
   do i_l = 1, lmax + 1
      ell = i_l - 1

      ! Set number of projectors to 1 for the local potential channel
      ! otherwise, use the user-specified number of projectors
      if (ell == lloc) then
         n_proj_l = 1
      else
         n_proj_l = kb_n_proj(i_l)
      end if

      ! Get principal quantum number of the highest-n core state for this l
      enn = i_l
      do i_state = 1, n_core
         if (ae_l(i_state) == ell) then
            enn = ae_n(i_state) + 1
         end if
      end do

      ! Get all-electron bound states for KB projectors
      i_proj = 0
      do i_state = n_core + 1, n_core + n_valence
         if (ae_l(i_state) == ell) then
            i_proj = i_proj + 1
            kb_state_e(i_proj, i_l) = ae_e(i_state)
            kb_state_n(i_proj, i_l) = ae_n(i_state)
            kb_state_rpsi(:, i_proj, i_l) = ae_rpsi(:, i_state)
            kb_state_drpsi(:, i_proj, i_l) = ae_drpsi(:, i_state)
         end if
         if (i_proj == kb_n_proj(i_l)) exit
      end do

      ! Find all-electron well states for KB projectors
      if (i_proj == 0) then
         ! If there were no valence states,
         ! use psp_e_in(l_idx) for the first well state
         kb_state_e(1, i_l) = psp_e_in(i_l)
      else
         ! otherwise, shift the valence state energy by kb_de(l_idx)
         ! and increment the principal quantum number
         do i = i_proj + 1, kb_n_proj(i_l)
            kb_state_e(i, i_l) = kb_state_e(i - 1, i_l) + db_de(i_l)
            kb_state_n(i, i_l) = kb_state_n(i - 1, i_l) + 1
         end do
      end if
      do i = 1, kb_n_proj(i_l) - i_proj
         call find_well_state( &
         ! Inputs
                              z, scalar_relativistic, &
                              log_mesh_size, log_mesh_r, &
                              ae_v_tot, &
                              kb_state_n(i_proj, i_l), ell, &
         ! Outputs
                              kb_state_e(i_proj, i_l), kb_state_match_idx(i_proj, i_l), &
                              kb_state_rpsi(:, i_proj, i_l), kb_state_drpsi(:, i_proj, i_l)
         )
      end do

      ! Calculate relativistic corrections to the potential
      ! to force the projectors to 0 at their core radii
      do i_proj = 1, kb_n_proj(i_l)
         call calculate_relativistic_correction( &
         ! Inputs
                                                z, scalar_relativistic, &
                                                log_mesh_size, log_mesh_r, &
                                                ae_v_tot, &
                                                ell, &
                                                psp_r_core_idx(i_l), &
                                                kb_state_e(i_proj, i_l), &
                                                kb_state_rpsi(:, i_proj, i_l), &
                                                kb_state_drpsi(:, i_proj, i_l), &
         ! Outputs
                                                kb_v_rel_corr(:, i_proj, i_l)
         )
      end do

      ! Calculate all-electron overlaps
      do j_proj = 1, kb_n_proj(i_l)
         do i_proj = 1, i_proj
            ae_ovlp(i_proj, j_proj, i_l) = compute_overlap(z, scalar_relativistic, &
                                                        log_mesh_size, log_mesh_r, &
                                                        ell, &
                                                        psp_r_core_idx(i_l), &
                                                        kb_state_rpsi(:, i_proj, i_l), &
                                                        kb_state_rpsi(:, j_proj, i_l))
            ae_ovlp(j_proj, i_proj, i_l) = ae_ovlp(i_proj, j_proj, i_l)
         end do
      end do

      ! Construct and optimized the Kleinman-Bylander form of the pseudopotential
      call construct_optimized_kleinman_bylander( &
      ! Inputs
      ! Outputs
      kb_q_max(i_l),
      ps_opt_wfn(:, :, i_l),
      psp_opt_v_sl(:, i_l),
      kb_p_kb(:, :, i_l),
      )

   end do

   ! Perform Kleinman-Bylander decomposition

   ! Compute pseudo-atomic valence states and charge density

   ! Construct model core charge
end program oncvpsp
