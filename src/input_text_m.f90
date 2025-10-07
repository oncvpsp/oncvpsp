module input_text_m
   use, intrinsic :: iso_fortran_env, only: dp => real64, stdout => output_unit, stderr => error_unit
   implicit none
   private
   public :: read_input_text, read_input_text_r
contains

   subroutine read_input_text(unit, inline, &
                              atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                              lmax, rc, ep, ncon, nbas, qcut, &
                              lloc, lpopt, dvloc0, nproj, debl, &
                              icmod, fcfact, rcfact, &
                              fcfact_min, fcfact_max, fcfact_step, &
                              rcfact_min, rcfact_max, rcfact_step, &
                              epsh1, epsh2, depsh, rxpsh, &
                              rlmax, drl, &
                              ncnf, nvcnf, nacnf, lacnf, facnf)
      !> File unit number for input
      integer, intent(in) :: unit
      !> Current line number in input file
      integer, intent(out) :: inline
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
      !> Scaling factor range minimum for optimization
      real(dp), intent(out) :: fcfact_min
      !> Scaling factor range maximum for optimization
      real(dp), intent(out) :: fcfact_max
      !> Scaling factor step size for optimization
      real(dp), intent(out) :: fcfact_step
      !> Core charge width factor
      real(dp), intent(out) :: rcfact
      !> Core charge width factor minimum for optimization
      real(dp), intent(out) :: rcfact_min
      !> Core charge width factor maximum for optimization
      real(dp), intent(out) :: rcfact_max
      !> Core charge width factor step size for optimization
      real(dp), intent(out) :: rcfact_step
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

      !Local variables
      !> I/O status variable
      integer :: ios
      !> Loop indices
      integer :: ii, jj
      !> Angular momentum index (l+1)
      integer :: l1
      !> Temporary variable for reading angular momentum
      integer :: lt
      !> Temporary line string for reading model core charge, log derivative analysis parameters
      character(len=2048) :: line

      nproj(:) = 0
      fcfact = 0.d0
      rcfact = 0.d0
      rc(:) = 0.d0
      ep(:) = 0.d0
      rxpsh = -1.d0  ! negative value for auto setting of radius
      ! default values for icmod=4 optimization ranges
      ! overriden by input values if icmod=5
      fcfact_min = 1.5d0
      fcfact_max = 6.0d0
      fcfact_step = 0.5d0
      rcfact_min = 1.0d0
      rcfact_max = 1.9d0
      rcfact_step = 0.1d0

      ! read input data
      inline = 0

      ! atom and reference configuration
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
      call read_error(ios, inline)

      call cmtskp(unit, inline)
      do ii = 1, nc + nv
         read (unit, *, iostat=ios) na(ii), la(ii), fa(ii)
         call read_error(ios, inline)
      end do

      ! pseudopotential and optimization
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) lmax
      call read_error(ios, inline)

      call cmtskp(unit, inline)
      do l1 = 1, lmax + 1
         read (unit, *, iostat=ios) lt, rc(l1), ep(l1), ncon(l1), nbas(l1), qcut(l1)
         if (lt /= l1 - 1) ios = 999
         call read_error(ios, inline)
      end do

      ! local potential
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
      call read_error(ios, inline)

      ! Vanderbilt-Kleinman-Bylander projectors
      call cmtskp(unit, inline)
      do l1 = 1, lmax + 1
         read (unit, *, iostat=ios) lt, nproj(l1), debl(l1)
         if (lt /= l1 - 1) ios = 999
         call read_error(ios, inline)
      end do

      ! model core charge
      call cmtskp(unit, inline)
      read (unit, '(a)', iostat=ios) line
      ! icmod=0,1,2,4[default grid] : icmod, fcfact
      read (line, *, iostat=ios) icmod, fcfact
      ! icmod==3 : icmod, fcfact, rcfact
      if (ios == 0 .and. icmod == 3) then
         read (line, *, iostat=ios) icmod, fcfact, rcfact
      end if
      ! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
      !   fcfact_min, fcfact_max, fcfact_step,
      !   rcfact_min, rcfact_max, rcfact_step
      if (ios == 0 .and. icmod >= 4) then
         read (line, *, iostat=ios) icmod, fcfact, rcfact, &
            fcfact_min, fcfact_max, fcfact_step, &
            rcfact_min, rcfact_max, rcfact_step
         ! icmod==4[implicit grid]: *_min, *_max, *_step not provided, fall back to defaults
         if (ios /= 0) then
            read (line, *, iostat=ios) icmod, fcfact
            rcfact = 0.d0
            fcfact_min = 1.5d0
            fcfact_max = 6.0d0
            fcfact_step = 0.5d0
            rcfact_min = 1.0d0
            rcfact_max = 1.9d0
            rcfact_step = 0.1d0
         end if
      end if
      call read_error(ios, inline)

      ! log derivative analysis
      call cmtskp(unit, inline)
      read (unit, '(a)', iostat=ios) line
      read (line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
      if (ios /= 0) then
         read (line, *, iostat=ios) epsh1, epsh2, depsh
         rxpsh = -1.d0  ! negative value for auto setting of radius
      end if
      call read_error(ios, inline)

      ! output grid
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) rlmax, drl
      call read_error(ios, inline)

      ! test configurations
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) ncnf
      call read_error(ios, inline)

      do jj = 2, ncnf + 1
         call cmtskp(unit, inline)
         read (unit, *, iostat=ios) nvcnf(jj)
         call read_error(ios, inline)
         do ii = nc + 1, nc + nvcnf(jj)
            call cmtskp(unit, inline)
            read (unit, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
            call read_error(ios, inline)
         end do
      end do

      ! end of reading input data
   end subroutine read_input_text

   subroutine read_input_text_r(unit, inline, &
                                atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                                lmax, rc, ep, ncon, nbas, qcut, &
                                lloc, lpopt, dvloc0, nproj, debl, &
                                icmod, fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step, &
                                rcfact_min, rcfact_max, rcfact_step, &
                                epsh1, epsh2, depsh, rxpsh, &
                                rlmax, drl, &
                                ncnf, nvcnf, nacnf, lacnf, facnf)
      !> File unit number for input
      integer, intent(in) :: unit
      !> Current line number in input file
      integer, intent(out) :: inline
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
      !> Maximum angular momentum
      integer, intent(out) :: lmax
      !> Pseudopotential cutoff radii
      real(dp), intent(out) :: rc(6)
      !> Pseudopotential energies
      real(dp), intent(out) :: ep(6, 2)
      !> Number of matching constraints for each pseudopotential
      integer, intent(out) :: ncon(6)
      !> Number of basis functions for each pseudopotential
      integer, intent(out) :: nbas(6)
      !> Maximum wave number for each pseudopotential
      real(dp), intent(out) :: qcut(6)
      !> Angular momentum used for local potential
      integer, intent(out) :: lloc
      !>
      integer, intent(out) :: lpopt
      !> Local potential offset at origin
      real(dp), intent(out) :: dvloc0
      !> Number of projectors for each angular momentum
      integer, intent(out) :: nproj(6)
      !> Energy shift for unbound states
      real(dp), intent(out) :: debl(6)
      !> Model core charge flag
      integer, intent(out) :: icmod
      !> Scaling factor for core charge
      real(dp), intent(out) :: fcfact
      !> Scaling factor range minimum for optimization
      real(dp), intent(out) :: fcfact_min
      !> Scaling factor range maximum for optimization
      real(dp), intent(out) :: fcfact_max
      !> Scaling factor step size for optimization
      real(dp), intent(out) :: fcfact_step
      !> Core charge width factor
      real(dp), intent(out) :: rcfact
      !> Core charge width factor minimum for optimization
      real(dp), intent(out) :: rcfact_min
      !> Core charge width factor maximum for optimization
      real(dp), intent(out) :: rcfact_max
      !> Core charge width factor step size for optimization
      real(dp), intent(out) :: rcfact_step
      !> Lower bound of energy range for log derivative analysis
      real(dp), intent(out) :: epsh1
      !> Upper bound of energy range for log derivative analysis
      real(dp), intent(out) :: epsh2
      !> Energy step size for log derivative analysis
      real(dp), intent(out) :: depsh
      !> Radius for log derivative analysis
      real(dp), intent(out) :: rxpsh
      !> Maximum radius for output grid
      real(dp), intent(out) :: rlmax
      !> Grid spacing for output grid
      real(dp), intent(out) :: drl
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

      !Local variables
      !> I/O status variable
      integer :: ios
      !> Loop indices
      integer :: ii, jj
      !> Angular momentum index (l+1)
      integer :: l1
      !> Temporary variable for reading angular momentum
      integer :: lt
      !> Temporary line string for reading log derivative analysis
      character(len=1024) :: line

      nproj(:) = 0
      fcfact = 0.d0
      rcfact = 0.d0
      rc(:) = 0.d0
      ep(:, :) = 0.d0
      rxpsh = -1.d0  ! negative value for auto setting of radius
      ! default values for icmod=4 optimization ranges
      fcfact_min = 1.5d0
      fcfact_max = 6.0d0
      fcfact_step = 0.5d0
      rcfact_min = 1.0d0
      rcfact_max = 1.9d0
      rcfact_step = 0.1d0

      ! read input data
      inline = 0

      ! atom and reference configuration
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
      call read_error(ios, inline)

      call cmtskp(unit, inline)
      do ii = 1, nc + nv
         read (unit, *, iostat=ios) na(ii), la(ii), fa(ii)
         call read_error(ios, inline)
      end do

      ! pseudopotential and optimization
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) lmax
      call read_error(ios, inline)

      call cmtskp(unit, inline)
      do l1 = 1, lmax + 1
         read (unit, *, iostat=ios) lt, rc(l1), ep(l1, 1), ncon(l1), nbas(l1), qcut(l1)
         if (lt /= l1 - 1) ios = 999
         call read_error(ios, inline)
         ep(l1, 2) = ep(l1, 1)
      end do

      ! local potential
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
      call read_error(ios, inline)

      ! Vanderbilt-Kleinman-Bylander projectors
      call cmtskp(unit, inline)
      do l1 = 1, lmax + 1
         read (unit, *, iostat=ios) lt, nproj(l1), debl(l1)
         if (lt /= l1 - 1) ios = 999
         call read_error(ios, inline)
      end do

      ! model core charge
      call cmtskp(unit, inline)
      read (unit, '(a)', iostat=ios) line
      ! icmod=0,1,2,4[default grid] : icmod, fcfact
      read (line, *, iostat=ios) icmod, fcfact
      ! icmod==3 : icmod, fcfact, rcfact
      if (ios == 0 .and. icmod == 3) then
         read (line, *, iostat=ios) icmod, fcfact, rcfact
      end if
      ! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
      !   fcfact_min, fcfact_max, fcfact_step,
      !   rcfact_min, rcfact_max, rcfact_step
      if (ios == 0 .and. icmod >= 4) then
         read (line, *, iostat=ios) icmod, fcfact, rcfact, &
            fcfact_min, fcfact_max, fcfact_step, &
            rcfact_min, rcfact_max, rcfact_step
         ! icmod==4[implicit grid]: *_min, *_max, *_step not provided, fall back to defaults
         if (ios /= 0) then
            read (line, *, iostat=ios) icmod, fcfact
            rcfact = 0.d0
            fcfact_min = 1.5d0
            fcfact_max = 6.0d0
            fcfact_step = 0.5d0
            rcfact_min = 1.0d0
            rcfact_max = 1.9d0
            rcfact_step = 0.1d0
         end if
      end if
      call read_error(ios, inline)

      ! log derivative analysis
      call cmtskp(unit, inline)
      read (unit, '(a)', iostat=ios) line
      read (line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
      if (ios /= 0) then
         read (line, *, iostat=ios) epsh1, epsh2, depsh
         rxpsh = -1.d0  ! negative value for auto setting of radius
      end if
      call read_error(ios, inline)

      ! output grid
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) rlmax, drl
      call read_error(ios, inline)

      ! test configurations
      call cmtskp(unit, inline)
      read (unit, *, iostat=ios) ncnf
      call read_error(ios, inline)

      do jj = 2, ncnf + 1
         call cmtskp(unit, inline)
         read (unit, *, iostat=ios) nvcnf(jj)
         call read_error(ios, inline)
         do ii = nc + 1, nc + nvcnf(jj)
            call cmtskp(unit, inline)
            read (unit, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
            call read_error(ios, inline)
         end do
      end do

      ! end of reading input data
   end subroutine read_input_text_r

   subroutine cmtskp(unit, inline)
      ! skips lines of standard input (file 5) whose first character is #

      !In/Out variable
      integer :: inline
      integer, intent(in) :: unit

      !Local variable
      character(len=1) :: tst

      tst = '#'
      do while (tst == '#')
         read (unit, *) tst
         inline = inline + 1
      end do
      backspace (unit)
      inline = inline - 1

      return
   end subroutine cmtskp

   subroutine read_error(ios, inline)
      ! report data read error and stop

      !Input variables
      integer :: ios, inline

      inline = inline + 1
      if (ios /= 0) then
         write (6, '(a,i4)') 'Read ERROR, input data file line', inline
         write (6, '(a)') 'Program will stop'
         stop
      end if

      return
   end subroutine read_error

end module input_text_m
