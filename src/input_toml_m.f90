module input_toml_m
   use, intrinsic :: iso_fortran_env, only: dp => real64, stdout => output_unit, stderr => error_unit
   use tomlf  ! TOML Fortran library: documentation at toml-f.readthedocs.io
   implicit none
   private
   public :: read_input_toml

   interface check_associated
      module procedure toml_table_check_associated
      module procedure toml_array_check_associated
   end interface check_associated
contains

subroutine read_input_toml(unit, &
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
   implicit none
   ! Input variables
   !> The Fortran unit number to read from
   integer, intent(in) :: unit

   ! Output variables
   ! [Atom and reference configuration]
   !> Atomic symbol
   character(len=2), intent(out) :: atsym
   !> Atomic number
   real(dp), intent(out) :: zz
   !> Number of core electrons
   integer, intent(out) :: nc
   !> Number of valence electrons
   integer, intent(out) :: nv
   !> Exchange-correlation functional code
   integer, intent(out) :: iexc
   !> Pseudopotential output file format flag
   character(len=4), intent(out) :: psfile
   !> Principal quantum number array
   integer, intent(out) :: na(30)
   !> Orbital angular momentum quantum number array
   integer, intent(out) :: la(30)
   !> Occupation number array
   real(dp), intent(out) :: fa(30)
   ! [Pseudopotential and optimization]
   !> Maximum angular momentum
   integer, intent(out) :: lmax
   !> Core radii for pseudopotentials
   real(dp), intent(out) :: rc(6)
   !> Pseudopotential energies
   real(dp), intent(out) :: ep(6)
   !> Number of matching constraints for each pseudopotential
   integer, intent(out) :: ncon(6)
   !> Number of basis functions for each pseudopotential
   integer, intent(out) :: nbas(6)
   !> Maximum wave number for each pseudopotential
   real(dp), intent(out) :: qcut(6)
   ! [Local potential]
   !> Angular momentum used for local potential
   integer, intent(out) :: lloc
   !>
   integer, intent(out) :: lpopt
   !> Local potential offset at origin
   real(dp), intent(out) :: dvloc0
   ! [Vanderbilt-Kleinman-Bylander projectors]
   !> Number of projectors for each angular momentum
   integer, intent(out) :: nproj(6)
   !> Energy shift for unbound states
   real(dp), intent(out) :: debl(6)
   ! [Model core charge]
   !> Model core charge flag
   integer, intent(out) :: icmod
   !> Scaling factor for core charge
   real(dp), intent(out) :: fcfact
   !> Core charge width factor
   real(dp), intent(out) :: rcfact
   !> Whether to interpret teter_amp and teter_scale as values relative
   !> to modcore3's `rhocmatch` and `rcmatch`
   logical, intent(out) :: teter_relative
   real(dp), intent(out) :: teter_amp
   real(dp), intent(out) :: teter_scale
   real(dp), intent(out) :: teter_amp_min
   real(dp), intent(out) :: teter_amp_max
   real(dp), intent(out) :: teter_amp_step
   real(dp), intent(out) :: teter_scale_min
   real(dp), intent(out) :: teter_scale_max
   real(dp), intent(out) :: teter_scale_step
   character(len=1024), intent(out) :: teter_objective_name
   ! [Log derivative analysis]
   !> Lower bound of energy range for log derivative analysis
   real(dp), intent(out) :: epsh1
   !> Upper bound of energy range for log derivative analysis
   real(dp), intent(out) :: epsh2
   !> Energy step size for log derivative analysis
   real(dp), intent(out) :: depsh
   !> Radius for log derivative analysis
   real(dp), intent(out) :: rxpsh
   ! [Output grid]
   !> Maximum radius for output grid
   real(dp), intent(out) :: rlmax
   !> Grid spacing for output grid
   real(dp), intent(out) :: drl
   ! [Test configurations]
   !> Number of test configurations
   integer, intent(out) :: ncnf
   !> Number of valence states in test configurations
   integer, intent(out) :: nvcnf(5)
   !> Principal quantum number array for test configurations
   integer, intent(out) :: nacnf(30, 5)
   !> Orbital angular momentum quantum number array for test configurations
   integer, intent(out) :: lacnf(30, 5)
   !> Occupation number array for test configurations
   real(dp), intent(out) :: facnf(30, 5)

   ! Local variables
   !> The TOML table object
   type(toml_table), allocatable :: table
   !> A pointer to a child table
   type(toml_table), pointer :: child
   !> A pointer to a toml array
   type(toml_array), pointer :: arr
   !> Status variable for optional sections
   integer :: stat
   !> Iteration index
   integer :: i
   !> Temporary string
   character(len=:), allocatable :: tmp_str

   nproj(:) = 0
   rc(:) = 0.0_dp
   ep(:) = 0.0_dp
   atsym(:) = '  '
   psfile(:) = '    '
   teter_objective_name(:) = ' '

   call toml_load(table, unit)
   if (.not. allocated(table)) then
      write (stderr, '(A)') 'Error: Failed to load TOML input file.'
      stop 1
   end if

   !! oncvpsp section (required)
   call get_value(table, "oncvpsp", child)
   call check_associated(child, "[oncvpsp] section")
   call get_value(child, "atsym", tmp_str)
   if (len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 2) then
      write (stderr, '(A)') 'Error: Invalid atomic symbol in [oncvpsp] section.'
      stop 1
   else
      atsym(1:len_trim(tmp_str)) = tmp_str(1:len_trim(tmp_str))
   end if
   deallocate(tmp_str)
   call get_value(child, "z", zz)
   call get_value(child, "nc", nc)
   call get_value(child, "nv", nv)
   call get_value(child, "iexc", iexc)

   !! Linear mesh section (required)
   call get_value(table, "linear_mesh", child)
   call check_associated(child, "[linear_mesh] section")
   call get_value(child, "a", drl)
   call get_value(child, "rmax", rlmax)

   !! Reference configuration section (required)
   call get_value(table, "reference_configuration", child)
   call check_associated(child, '[reference_configuration] section')
   ! n (required, length = nc + nv)
   call get_value(child, "n", arr)
   call check_associated(arr, 'n in [reference_configuration] section')
   call check_len(arr, nc + nv, 'n')
   do i = 1, len(arr)
      call get_value(arr, i, na(i))
   end do
   ! l (required, length = nc + nv)
   call get_value(child, "l", arr)
   call check_associated(arr, 'l in [reference_configuration] section')
   call check_len(arr, nc + nv, 'l')
   do i = 1, len(arr)
      call get_value(arr, i, la(i))
   end do
   ! f (required, length = nc + nv)
   call get_value(child, "f", arr)
   call check_associated(arr, 'f in [reference_configuration] section')
   call check_len(arr, nc + nv, 'f')
   do i = 1, len(arr)
      call get_value(arr, i, fa(i))
   end do

   !! Pseudopotential section (required)
   call get_value(table, "pseudopotentials", child)
   call check_associated(child, "[pseudopotentials] section")
   ! lmax (required)
   call get_value(child, "lmax", lmax)
   ! l (required, length = lmax + 1, must be in order from 0 to lmax)
   call get_value(child, "l", arr)
   call check_associated(arr, 'l in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'l in [pseudopotentials] section')
   pseudopotentials_l: block ! pseudopotentials.l
      integer :: l
      do i = 1, len(arr)
         call get_value(arr, i, l)
         if (i /= l + 1) then
            write (stderr, '(A,I0,A,I0)') 'Error: l values in [pseudopotentials] section must be in order ', &
               'starting from 0 to lmax without skipping any values. Expected ', i-1, ' but got ', l, '.'
            stop 1
         end if
      end do
   end block pseudopotentials_l
   ! rc (required, length = lmax + 1)
   call get_value(child, "rc", arr)
   call check_associated(arr, 'rc in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'rc in [pseudopotentials] section')
   do i = 1, lmax + 1
      call get_value(arr, i, rc(i))
   end do
   ! ep (required, length = lmax + 1)
   call get_value(child, "ep", arr)
   call check_associated(arr, 'ep in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'ep in [pseudopotentials] section')
   do i = 1, lmax + 1
      call get_value(arr, i, ep(i))
   end do
   ! ncon (required, length = lmax + 1)
   call get_value(child, "ncon", arr)
   call check_associated(arr, 'ncon in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'ncon in [pseudopotentials] section')
   do i = 1, lmax + 1
      call get_value(arr, i, ncon(i))
   end do
   ! nbas (required, length = lmax + 1)
   call get_value(child, "nbas", arr)
   call check_associated(arr, 'nbas in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'nbas in [pseudopotentials] section')
   do i = 1, lmax + 1
      call get_value(arr, i, nbas(i))
   end do
   ! qcut (required, length = lmax + 1)
   call get_value(child, "qcut", arr)
   call check_associated(arr, 'qcut in [pseudopotentials] section')
   call check_len(arr, lmax + 1, 'qcut in [pseudopotentials] section')
   do i = 1, lmax + 1
      call get_value(arr, i, qcut(i))
   end do

   !! Local potential section (required)
   call get_value(table, "local_potential", child)
   call check_associated(child, "[local_potential] section")
   call get_value(child, "lloc", lloc)
   call get_value(child, "lpopt", lpopt)
   call get_value(child, "rcloc", rc(5))
   call get_value(child, "dvloc0", dvloc0)

   !! VKB projectors section (required)
   call get_value(table, "vkb_projectors", child)
   call check_associated(child, "[vkb_projectors] section")
   ! l (required, length = lmax + 1, must be in order from 0 to lmax)
   call get_value(child, "l", arr)
   call check_associated(arr, 'l in [vkb_projectors] section')
   call check_len(arr, lmax + 1, 'l in [vkb_projectors] section')
   vkb_projectors_l: block ! vkb_projectors.l
      integer :: l
      do i = 1, len(arr)
         call get_value(arr, i, l)
         if (i /= l + 1) then
            write (stderr, '(A,I0,A,I0)') 'Error: l values in [vkb_projectors] section must be in order ', &
               'starting from 0 to lmax without skipping any values. Expected ', i-1, ' but got ', l, '.'
            stop 1
         end if
      end do
   end block vkb_projectors_l
   ! nproj (required, length = lmax + 1)
   call get_value(child, "nproj", arr)
   call check_associated(arr, 'nproj in [vkb_projectors] section')
   call check_len(arr, lmax + 1, 'nproj in [vkb_projectors] section')
   do i = 1, lmax + 1
      call get_value(arr, i, nproj(i))
   end do
   ! debl (required, length = lmax + 1)
   call get_value(child, "debl", arr)
   call check_associated(arr, 'debl in [vkb_projectors] section')
   call check_len(arr, lmax + 1, 'debl in [vkb_projectors] section')
   do i = 1, lmax + 1
      call get_value(arr, i, debl(i))
   end do

   !! Model core charge (nlcc) section (required)
   call get_value(table, "model_core_charge", child)
   call check_associated(child, "[model_core_charge] section")
   call get_value(child, "icmod", icmod)
   call get_value(child, "fcfact", fcfact)
   if (icmod >= 0 .and. icmod <= 2) then
      rcfact = 0.0_dp
   end if
   if (icmod >= 3 .and. icmod <= 4) then
      ! Ignore rcfact for icmod 3 and 4; use teter_amp and teter_scale instead
      ! call get_value(child, "rcfact", rcfact)
      call get_value(child, "teter_amp", fcfact)
      call get_value(child, "teter_scale", rcfact)
      call get_value(child, "teter_relative", teter_relative, default=.true.)
   end if
   if (icmod == 4) then
      call get_value(child, "teter_amp_min", teter_amp_min, default=1.5_dp)
      call get_value(child, "teter_amp_max", teter_amp_max, default=6.0_dp)
      call get_value(child, "teter_amp_step", teter_amp_step, default=0.5_dp)
      call get_value(child, "teter_scale_min", teter_scale_min, default=1.0_dp)
      call get_value(child, "teter_scale_max", teter_scale_max, default=1.9_dp)
      call get_value(child, "teter_scale_step", teter_scale_step, default=0.1_dp)
      call get_value(child, "teter_objective_name", tmp_str, default="d2exc_rmse")
      if (len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 1024) then
         write (stderr, '(A)') 'Error: Invalid teter_objective_name in [model] section.'
         stop 1
      else
         teter_objective_name(1:len_trim(tmp_str)) = tmp_str(1:len_trim(tmp_str))
      end if
      deallocate(tmp_str)
   end if
   if (icmod < 0 .or. icmod > 4) then
      write (stderr, '(A)') 'Error: Invalid icmod in [model_core_charge] section. Must be 0, 1, 2, 3, or 4.'
      stop 1
   end if

   !! Log-derivative analysis section (required)
   call get_value(table, "log_derivative_analysis", child)
   call check_associated(child, "[log_derivative_analysis] section")
   call get_value(child, "epsh1", epsh1)
   call get_value(child, "epsh2", epsh2)
   call get_value(child, "depsh", depsh)
   call get_value(child, "rxpsh", rxpsh, -1.0_dp)

   !! PP output section (required)
   call get_value(table, "pp_output", child)
   call check_associated(child, "[pp_output] section")
   call get_value(child, "psfile", tmp_str)
   if (len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 4) then
      write (stderr, '(A)') 'Error: Invalid psfile in [pp_out] section.'
      stop 1
   else
      psfile(1:len_trim(tmp_str)) = tmp_str(1:len_trim(tmp_str))
   end if
   deallocate(tmp_str)

   !! Test configurations section (optional)
   call get_value(table, "test_configurations", arr)
   if (associated(arr)) then
      test_configurations: block
         ! Test configuration index
         integer :: icnf
         ! Valence state index within test configuration
         integer :: iv
         ! Core state index
         integer :: ic
         ! Core + valence state index
         integer :: icv
         ! Configuration sub-table
         type(toml_table), pointer :: cnf_tbl
         ! Array within test configuration (e.g. n, l, f)
         type(toml_array), pointer :: cnf_arr
         ncnf = len(arr)
         nvcnf(1) = nv  ! First configuration is the reference configuration
         ! Determine maximum number of valence states in test configurations
         do icnf = 1, ncnf
            ! Get the configuration sub-table
            call get_value(arr, icnf, cnf_tbl)
            ! Get the n array to determine number of valence states in this configuration
            call get_value(cnf_tbl, "n", cnf_arr)
            call check_associated(cnf_arr, 'n in [[test_configurations]] section')
            nvcnf(icnf + 1) = len(cnf_arr)
         end do
         ! Read the test configurations
         do icnf = 1, ncnf
            ! Get the configuration sub-table
            call get_value(arr, icnf, cnf_tbl)
            ! n
            call get_value(cnf_tbl, "n", cnf_arr)
            call check_associated(cnf_arr, 'n in [[test_configurations]] section')
            do iv = 1, nvcnf(icnf + 1)
               call get_value(cnf_arr, iv, nacnf(nc + iv, icnf + 1))
            end do
            ! l
            call get_value(cnf_tbl, "l", cnf_arr)
            call check_associated(cnf_arr, 'l in [[test_configurations]] section')
            do iv = 1, nvcnf(icnf + 1)
               call get_value(cnf_arr, iv, lacnf(nc + iv, icnf + 1))
            end do
            ! f
            call get_value(cnf_tbl, "f", cnf_arr)
            call check_associated(cnf_arr, 'f in [[test_configurations]] section')
            do iv = 1, nvcnf(icnf + 1)
               call get_value(cnf_arr, iv, facnf(nc + iv, icnf + 1))
            end do
         end do
         ! Add the reference configuration as the first configuration
         do icv = 1, nc + nv
            nacnf(icv, 1) = na(icv)
            lacnf(icv, 1) = la(icv)
            facnf(icv, 1) = fa(icv)
         end do
         ! Add the reference core configuration to the test configurations
         do icnf = 1, ncnf
            do ic = 1, nc
               nacnf(ic, icnf + 1) = na(ic)
               lacnf(ic, icnf + 1) = la(ic)
               facnf(ic, icnf + 1) = fa(ic)
            end do
         end do
      end block test_configurations
   else
      ncnf = 0
      nvcnf(1) = nv  ! First configuration
      do i = 1, nc + nv
         nacnf(i, 1) = na(i)
         lacnf(i, 1) = la(i)
         facnf(i, 1) = fa(i)
      end do
   end if
end subroutine read_input_toml

subroutine toml_table_check_associated(ptr, name)
   implicit none
   type(toml_table), pointer :: ptr
   character(len=*), intent(in) :: name
   if (.not. associated(ptr)) then
      write (stderr, '(A,A)') 'Error: Missing or invalid table', trim(name)
      stop 1
   end if
end subroutine toml_table_check_associated

subroutine toml_array_check_associated(ptr, name)
   implicit none
   type(toml_array), pointer :: ptr
   character(len=*), intent(in) :: name
   if (.not. associated(ptr)) then
      write (stderr, '(A,A)') 'Error: Missing or invalid array', trim(name)
      stop 1
   end if
end subroutine toml_array_check_associated

subroutine check_len(arr, expected_len, name)
   implicit none
   type(toml_array), pointer :: arr
   integer, intent(in) :: expected_len
   character(len=*), intent(in) :: name
   if (len(arr) /= expected_len) then
      write (stderr, '(A,A,A,A,I0,A,I0)') 'Error: Length of', ' ', trim(name), '(len=', len(arr), &
         ') must be equal to ', expected_len
      stop 1
   end if
end subroutine check_len

end module input_toml_m
