module output_text_m
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use postprocess_m, only: get_pseudo_linear_mesh_parameters
   implicit none
   private
   public :: write_config_text, &
      write_rho_vpuns_text, &
      write_vloc_text, &
      write_rho_rhoc_rhom_text, &
      write_wavefunctions_text, &
      write_vkb_projectors_text, &
      write_wavefunctions_vkb_text, &
      write_convergence_profile_text, &
      write_phsft_text
contains

subroutine write_input_with_ae_eig_text(unit, &
                                        atsym, zz, nc, nv, iexc, psfile, &
                                        na, la, fa, ea, &
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
   integer, intent(in) :: unit
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
   !> Reference configuration all-electron eigenvalue array
   real(dp), intent(in) :: ea(30)
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
   integer :: ii, jj, l1

   write(unit, '(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write(unit, '(a)') '# atsym  z   nc   nv     iexc    psfile'
   write(unit, '(a,a,f6.2,2i5,i8,2a)') '  ', trim(atsym), zz, nc, nv, iexc, '      ', psfile
   write(unit, '(a/a)') '#','#   n    l    f        energy (Ha)'
   do ii = 1, nc + nv
      write(unit, '(2i5,f8.2,1pe18.7)') na(ii), la(ii), fa(ii), ea(ii)
   end do

   write(unit, '(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
   write(unit, '(i5)')  lmax
   write(unit, '(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
   do l1=1,lmax+1
      write(unit, '(i5,2f10.5,2i5,f10.5)') l1-1, rc(l1), ep(l1), ncon(l1), nbas(l1), qcut(l1)
   end do

   write(unit, '(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', '   dvloc0'
   write(unit, '(2i5,f10.5,a,f10.5)') lloc, lpopt, rc(5), '   ', dvloc0

   write(unit, '(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', '# l, nproj, debl'
   do l1 = 1, lmax + 1
      write(unit, '(2i5,f10.5)') l1-1, nproj(l1), debl(l1)
   end do

   write(unit, '(a/a/a)') '#','# MODEL CORE CHARGE', '# icmod, fcfact, rcfact'
   write(unit, '(i5,2f10.5)') icmod, fcfact, rcfact

   write(unit, '(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', '# epsh1, epsh2, depsh'
   write(unit, '(3f8.2)') epsh1, epsh2, depsh

   write(unit, '(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
   write(unit, '(2f8.4)') rlmax, drl

   write(unit, '(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
   write(unit, '(i5)') ncnf
   write(unit, '(a/a)') '# nvcnf','#   n    l    f'
   do jj = 2, ncnf + 1
      write(unit, '(i5)') nvcnf(jj)
      do ii = nc + 1, nc + nvcnf(jj)
         write(unit, '(2i5,f8.2)') nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
      end do
      write(unit, '(a)') '#'
   end do

end subroutine write_input_with_ae_eig_text

subroutine write_config_text(unit, nc, nvt, nat, lat, fat, eat, eatp, etot, eaetst, epstot, etsttot)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states for this test configuration
   integer, intent(in) :: nvt
   !> Array of state principal quantum numbers
   integer, intent(in) :: nat(30)
   !> Array of state angular momenta
   integer, intent(in) :: lat(30)
   !> Array of occupation numbers and energies for test configuration
   real(dp), intent(in) :: fat(30)
   !> All-electron energies for test configuration
   real(dp), intent(in) :: eat(30)
   !> Pseudo energies for test configuration
   real(dp), intent(in) :: eatp(30)
   !> Reference all-electron total energy
   real(dp), intent(in) :: etot
   !> Test all-electron total energy
   real(dp), intent(in) :: eaetst
   !> Reference pseudo total energy
   real(dp), intent(in) :: epstot
   !> Test pseudo total energy
   real(dp), intent(in) :: etsttot

   ! Local variables
   integer :: ii

   write(unit,'(/a)') '   n   l     f        eae           eps        diff'
   do ii = 1, nc
      write(unit,'(2i4,f8.4,f14.8)') nat(ii), lat(ii), fat(ii), eat(ii)
   end do
   do ii = 1, nvt
      if (abs(fat(ii+nc)) < tiny(0.0_dp)) cycle
      write(unit,'(2i4,f8.4,2f14.8,1p,d12.2)') nat(ii + nc), lat(ii + nc), fat(ii + nc), eat(ii + nc), &
         eatp(ii), eatp(ii) - eat(ii + nc)
   end do

   write(unit,'(/a)') '    Total energies and differences'
   write(unit,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=', etot, '  AE_tst=', eaetst, '  dif=', eaetst - etot
   write(unit,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=', epstot, '  PS_tst=', etsttot, '  dif=', etsttot - epstot
   write(unit,'(a,1p,d10.2)') '      PSP excitation error=', eaetst - etot - etsttot + epstot
end subroutine write_config_text

subroutine write_rho_vpuns_text(unit, mmax, lmax, drl, nrl, rr, irc, rho, vpuns)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Unscreened pseudopotentials
   real(dp), intent(in) :: vpuns(mmax,5)

   ! Local variables
   integer :: ii, l1, n1, n2, n3, n4
   real(dp) :: dr, r0

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)
   dr = min(drl, rr(n4) / 200.0_dp)
   r0 = rr(n1) - dr
   write(unit,'(/2a/)') ' radii, charge,',' pseudopotentials (ll=0, 1, lmax)'
   do ii = n1, n4
      if (rr(ii) > r0) then
         write(unit, '(a,6(f12.7,1x))') '!p', rr(ii), rho(ii), (vpuns(ii, l1), l1 = 1, lmax + 1)
         r0 = rr(ii) + dr
      end if
   end do
end subroutine write_rho_vpuns_text

subroutine write_vloc_text(unit, mmax, rr, lmax, irc, drl, nrl, vloc)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Local pseudopotential
   real(dp), intent(in) :: vloc(mmax)

   ! Local variables
   integer :: ii, n1, n2, n3, n4
   real(dp) :: dr, r0

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)
   dr = min(drl, rr(n4) / 200.0_dp)
   r0 = rr(n1) - dr
   do ii=n1,n4
      if (rr(ii) > r0) then
         write(unit, '(a, 2(f12.7,1x))') ' !L', rr(ii), vloc(ii)
         r0 = rr(ii) + dr
      end if
   end do
end subroutine write_vloc_text

subroutine write_rho_rhoc_rhom_text(unit, mmax, rr, lmax, irc, drl, nrl, rho, rhoc, rhomod)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax)

   ! Local variables
   integer :: ii, n1, n2, n3, n4
   real(dp) :: dr, r0, rmx

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)

   rmx = 0.0_dp
   do ii = n1, n2
      rmx = max(rmx, 5.0_dp * rho(ii), rhomod(ii))
   end do

   write(unit, '(//2a/)') ' radii, charge,',' core charge, model core charge'
   dr = min(drl, rr(n2) / 200.0_dp)
   r0 = rr(n1) - dr
   do ii = n1, n2
      if (rr(ii) > r0) then
         if (rhoc(ii) < rmx)then
            write(unit, '(a,8(f12.7,1x))') '!r', rr(ii), rho(ii), rhoc(ii), rhomod(ii)
         else if (rr(ii) > 0.01_dp) then
            write(unit, '(a,8(f12.7,1x))') '!r', rr(ii), rho(ii), rmx, rhomod(ii)
         end if
         r0 = rr(ii) + dr
      end if
   end do

end subroutine write_rho_rhoc_rhom_text

subroutine write_wavefunctions_text(unit, mmax, rr, lmax, irc, drl, nrl, &
                                    n, l, iproj, sign_ae, sign_ps, uu_ae, uu_ps, is_scattering)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Principal quantum number of state (for bound states)
   integer, intent(in) :: n
   !> Angular momentum of state
   integer, intent(in) :: l
   !> Index of projector at this l
   integer, intent(in) :: iproj
   !> Sign of AE wavefunction at matching point
   real(dp), intent(in) :: sign_ae
   !> Sign of PS wavefunction at matching point
   real(dp), intent(in) :: sign_ps
   !> AE wavefunction
   real(dp), intent(in) :: uu_ae(mmax)
   !> PS wavefunction
   real(dp), intent(in) :: uu_ps(mmax)
   !> .true. for scattering states, .false. for bound states
   logical, intent(in) :: is_scattering

   ! Local variables
   integer :: ii, n1, n2, n3, n4
   real(dp) :: dr, r0

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)

   if (is_scattering) then
      write(unit,'(//a,i2,a,i2,a,a/)') 'scattering, iprj=',iproj,',  l=', l, ', all-electron wave function', ', pseudo w-f'
   else
      write(unit,'(//a,i2,a,i2,a,a/)') 'n=', n, ',  l=', l, ', all-electron wave function', ', pseudo w-f'
   end if

   dr = min(drl, rr(n2) / 200.0_dp)
   r0 = rr(n1) - dr
   do ii = n1, n2
      if (rr(ii) > r0) then
         write(unit,'(a,i5,i1,3(f12.6,1x))') '&', iproj, l, rr(ii), sign_ae * uu_ae(ii), sign_ps * uu_ps(ii)
         r0 = rr(ii) + dr
      end if
   end do
end subroutine write_wavefunctions_text

subroutine write_vkb_projectors_text(unit, mmax, rr, lmax, irc, drl, nrl, mxprj, &
                                     ll, nproj, vkb)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Angular momentum of state
   integer, intent(in) :: ll
   !> Number of projectors for this l
   integer, intent(in) :: nproj
   !> VKB projector
   real(dp), intent(in) :: vkb(mmax, mxprj)

   ! Local variables
   integer :: ii, iproj, n1, n2, n3, n4
   real(dp) :: dr, r0

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)

   write(unit,'(/)')
   dr = min(drl, rr(n3) / 200.0_dp)
   r0 = rr(n1) - dr
   do ii = n1, n3
      if (rr(ii) > r0) then
         write(unit,'(a,i6,6(f12.6,1x))') '!J', ll, rr(ii), (vkb(ii, iproj), iproj = 1, nproj)
         r0 = rr(ii) + dr
      end if
   end do


end subroutine write_vkb_projectors_text

subroutine write_wavefunctions_vkb_text(unit, mmax, rr, lmax, irc, drl, nrl, &
                                        lloc, mxprj, npa, nproj, sign_ae, sign_ps, uu_ae, uu_ps, is_scattering, &
                                        vkb)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Size of radial grid
   integer, intent(in) :: mmax
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Indices of core radii on the logarithmic radial mesh
   integer, intent(in) :: irc(6)
   !> Spacing of linear radial mesh
   real(dp), intent(in) :: drl
   !> Number of points in linear radial mesh
   integer, intent(in) :: nrl
   !> Angular momentum of local potential
   integer, intent(in) :: lloc
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Principal quantum number of state (for bound states)
   integer, intent(in) :: npa(mxprj, 6)
   !> Number of projectors for each l
   integer, intent(in) :: nproj(6)
   !> Sign of AE wavefunction at matching point
   real(dp), intent(in) :: sign_ae(mxprj, 4)
   !> Sign of PS wavefunction at matching point
   real(dp), intent(in) :: sign_ps(mxprj, 4)
   !> AE wavefunction
   real(dp), intent(in) :: uu_ae(mmax, mxprj, 4)
   !> PS wavefunction
   real(dp), intent(in) :: uu_ps(mmax, mxprj, 4)
   !> .true. for scattering states, .false. for bound states
   logical, intent(in) :: is_scattering(mxprj, 4)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 6)

   ! Local variables
   integer :: ii, iproj, l1, ll, n1, n2, n3, n4

   call get_pseudo_linear_mesh_parameters(mmax, rr, lmax, irc, drl, nrl, &
                                          n1, n2, n3, n4)
   do l1 = 1, lmax + 1
      ll = l1 - 1
      do iproj = 1, nproj(l1)
         call write_wavefunctions_text(stdout, mmax, rr, lmax, irc, drl, nrl, &
                                       npa(iproj, l1), ll, iproj, &
                                       sign_ae(iproj, l1), sign_ps(iproj, l1), &
                                       uu_ae(:, iproj, l1), uu_ps(:, iproj, l1), &
                                       is_scattering(iproj, l1))
      end do
      call write_vkb_projectors_text(stdout, mmax, rr, lmax, irc, drl, nrl, &
                                     mxprj, ll, nproj(l1), vkb(:, :, l1))
   end do

end subroutine write_wavefunctions_vkb_text

subroutine write_convergence_profile_text(unit, lmax, mxprj, cvgplt)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Convergence profile array
   real(dp), intent(in) :: cvgplt(2,7,mxprj,4)

   ! Local variables
   integer :: ll, l1, ii, jj

   write(stdout, '(/a)') 'convergence profiles, (ll=0,lmax)'
   write(stdout, *) 'lmax', lmax
   do l1 = 1, lmax + 1
      ll = l1 - 1
      do jj = 1, 7
         if (abs(cvgplt(1, jj, 1, l1)) > 1e-22_dp) then
            write(stdout, '(a,i6,3(f12.6,1x))') '!C', ll, cvgplt(1, jj, 1, l1), cvgplt(2, jj, 1, l1)
         end if
      end do
   end do
end subroutine write_convergence_profile_text

subroutine write_phsft_text(unit, rpsh, npsh, epsh, pshf, pshp)
   implicit none
   ! Input variables
   !> Output unit
   integer, intent(in) :: unit
   !> Radius at which phase shift is calculated
   real(dp), intent(in) :: rpsh(4)
   !> Number of phase shift energies
   integer, intent(in) :: npsh
   !> Phase shift energies
   real(dp), intent(in) :: epsh(npsh)
   !> All-electron phase shifts
   real(dp), intent(in) :: pshf(npsh,4)
   !> Pseudopotential phase shifts
   real(dp), intent(in) :: pshp(npsh,4)

   ! Local variables
   integer :: l1, ll, ii

   do l1 = 1, 4
      ll = l1 - 1
      write(unit,'(/a,i2)') 'log derivativve data for plotting, l=',ll
      write(unit,'(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=',rpsh(l1)
      write(unit,'(a/)') 'l, energy, all-electron, pseudopotential'
      do ii = 1, npsh
         write(unit,'(a, i6, 3f12.6)') '! ', ll, epsh(ii), pshf(ii,l1), pshp(ii,l1)
      end do
   end do
end subroutine write_phsft_text

end module output_text_m
