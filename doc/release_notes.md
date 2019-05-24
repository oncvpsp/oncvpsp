## Version 4.0.1 - beta release, please report issues

1) The wellstate approach for relativistic calculations has been modified
   in light of the wellstate changes described below.  The confining well
   potential and node counts found for J=L+1/2 are now carried over and 
   used for J=L-1/2 to avoid possible inconsistencies in the corresponding
   projectors.

2) Several minor bugs (uninitialized variables, single-precision constants,
   potential array overwrites, overflowing buffers, etc.) have been
   corrected.  Thanks to  J. E. Pask and A. C. Medina

3) modcore2.f90 and modcore3.f90 have been modified to remove numerical noise
   in the 3rd and 4th derivatives that could influence the output on the linear
   grid at very small radii.

## Version 4.0.0 - beta release, please report issues

1) The code has been reorganized to permit 3 or more projectors per angular
   momentum L.  Bound states in successive shells and scattering states are
   processed together in loops over [1:nproj(l1)] throughout the code. This
   replaces the separate treatment of the first and second projectors used
   in earlier versions.  The maximum allowed nproj is set by mxprj, 
   presently 5, in oncvpsp.f90. (Going beyond 5 is probably useless, and 
   would require modifying some format statements by hand.)

  The following routines have been reorganized in this manner:

	gnu_script.f90		run_diag_sr_so_r.f90
	gnu_script_r.f90	run_ghosts.f90
	linout.f90		run_optimize.f90
	linout_r.f90		run_phsft.f90
	oncvpsp.f90		run_phsft_r.f90
	oncvpsp_nr.f90		run_plot.f90
	oncvpsp_r.f90		run_plot_r.f90
	psatom.f90		run_vkb.f90
	psatom_r.f90		run_vkb_r.f90
	run_config.f90		sbf_basis_con.f90
	run_config_r.f90	sr_so_r.f90
	run_diag.f90		upfout.f90
	run_diag_r.f90		upfout_r.f90


2) The overall optimization approach is retained, in the sense that each
   successive projector pseudo wave function (ppwf) is calculated independently
   enforcing off-diagonal overlap constraints with lower-index (presumably
   lower energy) ppwfs.  The underlying spherical Bessel function basis is 
   calculated based on the first ppwf's value and slope constraints at rc 
   as in the two-projector code, but is expanded at each stage to retain
   the same number of "null" functions for consistent variational flexibily
   in convergence optimization.  Convergence is typically found to be
   rather uniform, improving on that of previous two-projector results.

3) The wellstate algorithm for simulating scattering states with normalizable
   bound states for convergence optimization has changed.  A "soft step"
   potential is added to the all-electron potential as before, but the
   asymptotic potential is a fixed energy (presently 0.5 Ha) above the
   specified scattering-state energy.  Only the width of the step is
   adjusted to match that energy.  In the previous version, the number of
   nodes of this simulated bound state was fixed at one more than the
   next-lower ppwf.  The new algorithm explores this possibility  first,
   but then allows successively more nodes.  The new algorithm is much
   more robust, and will converge for considerably higher energies and
   larger intervals.  Applied to input data that were developed for the old
   algorithm, however, the results are very similar.  At present, the
   input data ep is used for the first well state if there are no bound
   ppwfs.  Otherwise the first is debl above the last bound ppwf as in
   older versions, and higher well states are debl above each other.  If
   a highly-localized bound state is found when searching for the first
   specified scattering state, it is substituted for that state since ghosts
   and poor performance would otherwise invariably result.

4) More projectors are an invitation to ghost states.  The ghost-detection
   routine from version 3.3.1 has been modified to ensure convergence when
   deeper core shells are treated as valence.  Ghost-detection has been
   incorporated in the fully-relativistic code.  Note that a legitimate
   low-energy scattering resonance may be reported as a positive-energy
   ghost, especially for L=3.  

5) The graphics output has been modified.  The usual semi-local plus
  local potential and charge density plots are displayed first.  Then,
  for each L, the all-electron and pseudo wave functions are compared
  for successive projector energies, either bound or scattering, and
  labeled by energy in the title.  Finally, the orthogonal projectors
  are displayed, with their coefficients in the title.  These plots
  span only radii 0 to 1.1*rc, and give a better idea of the real-space
  scales required if projectors are to be represented accurately on a real-
  space mesh.  The radial log mesh, which has been made denser to deal
  with deeper cores, is now "sampled" for these plots to reduce the total
  output size.  The final plot is the total kinetic energy error versus
  plane-wave cutoff for the first projector as in previous releases.
  The new graphics is illustrated in 73_Ta_plots.pdf in this directory.


6) The excited-configuration tests have been significantly modified.
   The solution if the radial Schroedinger equation for the pseudo wave
   functions is more challenging with more projectors, and must be started
   with tight energy bounds around the expected solution.  This is achieved
   by starting all excited-atom runs from the reference configuration,
   and "adiabatically" varying the occupation numbers concomitantly with
   self-consistency convergence.  This is in fact done in two stages,
   initially from the reference to the most highly ionized intermediate
   configuration, and then from that to the final one.  A more detailed
   description is in psatom.f90.

7) The fact that the scalar-relativistic all-electron radial Schroedinger
   equation contains a first-order derivative operator violates the
   conditions under which generalized norm conservation is exact.  This
   is seen in the "B matrix Hermiticity error" in the *.out file, and
   treated by an ad-hoc symmetrization.  Inner shells of heavy atoms
   reveal a related problem - projectors do not go to zero at rc as they
   do in a non-relativistic calculation.  A plausible ad-hoc fix for
   this is to smoothly blend an effective potential into the projector
   calculation over the last 5% of [0:rc].  This potential is calculated
   in vrel.f90.  These issues eventually limit the accuracy of the psp
   at sufficiently small rc.

8) Thanks to a contribution from A. C. Medina, the libxc capabilities are now
   extended to libxc versions 4.0.0 (just a coincidence) and above.  To use
   earlier versions, make.inc must be modified as indicated. The test/data
   files 38_Sr_lxc.dat and 79_Au_lxc.dat require libxc.

9) The tests/refs *.out files are computed with the code in this release,
   and show small differences from previous releases.  While the psps
   have not all been tested in benchmark calculations, enough testing
   has been done to expect very comparable performance.  Several datasets
   are new or have been modified to demonstrate new capabilities. Input data
   from the two presently published oncvpsp periodic tables runs without error.

     http://www.pseudo-dojo.org/ 
     http://www.quantum-simulation.org/potentials/sg15_oncv/ 

10) The unscreened local potetial at the origin of the linear output mesh
    is now set to the value on the first point of the internal log mesh.
    This suggestion of F. Gygi partially mitigates the erratic behavior of 
    GGA potentials arising from very small values of the psedo-chargre 
    density near the origin.  (Even a small model core charge completely
    eliminates this problem.)

11) For convenience in testing  many sets of input data, tests/multirun.sh
    and tests/multirun_r.sh will launch groups of NPROC=4 jobs at once with
    a command like <path>/multirun.sh *.dat for all the data files in the
    current directory.  NPROC cah be edited to suit your system.

12) I thank John Pask for stimulating my interest in this extension of
    oncvpsp, and for a great deal of testing and useful feedback in the
    course of the its development.

## Version 3.3.1

1) The recently-introduced ghost detector has been upgraded to scan positive
   energies as well.  Adjustments have also been made to overcome
   a few seldom-encountered errors reported to me.  See doc/ghosts.txt for
   details, including advice on fixing ghosts found using published
   oncvpsp input data.

2) The log-derivative routines have been slightly modified since ghosts
   found by the detector would occasionally fail to show up as steps in the
   pseudopotential log-derivative plots.  This primarily involved reducing
   the radii at which these derivatives are compared, so the plots
   will differ slightly from those of previous releases.  The energy ranges
   have been increased for some of the tests/data examples.

3) An indexing bug in the diagnostic tests which could change the test radii 
   from those intended has been fixed.  Some numbers in the Diagnosics
   will change as a result.

4) vkboutwf was modified to remedy a node-counting problem at very small
   core redii (eg., treating 1s as valence for a first-row atom)

5)  Perdew-Wang LDA using libxc with iexc==-001012 is now accepted for
    upf output for pwscf.

6)  Errors that cause stops will now write "ERROR" with some identification
    on either the standard output or standard error (both go to the *.out
    file in my scripts).  This will facilitate examination of results
    in which many sets of data are run by automated scripts with a simple
    "grep" command.  WARNINGs remain in caps as well.

7)  The optimize to old_optimize comparison has been deactivated in
    run_optimize.f90, given extensive experience with the corrected
    routine as discussed under 3.3.0 below.  I will continue to carry
    old_optimize.f90 for a while and it can be simply substituted
    by moving the comment "!"s following lines 171 and 346 in
    run_optimize.  This will reproduce the results of 3.2.3 while
    including the new ghost diagnostics and other improvements above.

8)  While the tests/refs/*.out files differ from 3.3.0 reflecting these
    changes, the embedded pseudopotential output files have no significant
    differences.

  Thanks to J. E. Pask and M. van Setten for their contributions.

## Version 3.3.0

1) An unfortunate math error made its way from my notes into the oncvpsp
   paper and into optimize.f90.  Fortunately, the only effect of this has
   been to keep this routine from reaching the exact minimum of the 
   residual kinetic energy above the cutoff qc, which has proven a very
   good proxy for plane-wave code convergence.  The corrected routine also
   significantly simplifies the algorithm to find the minimum.  Details
   are explained in Erratum_PhysRevB.88.085117.pdf in this directory which
   has been submitted.  The new version of run_optimize.f90 calls the
   renamed old_optimize.f90 and reports the new minimum and the typically
   extremely small difference.  Benchmark tests against all-electron results 
   using input designed with the old code show insignificant differences.
   The 3.2.3 results can be reproduced by commenting out lines 177-180 and
   293-296 of run_optimize.  The preceding old_optimize calls will be
   dropped in some later release.

2) Deep-lying ghosts far below the energy range of the usual log-derivative
   plots have surprised a few users with bizarre results in applications.
   A new ghost detector has been designed which should reliably detect
   these.  It is described in ghosts.txt in this directory.

## Version 3.2.3

1) An error common to modcore2 and modcore3 has been corrected.  The intended
   "crossover blend" from the Teter function to the all-electron core
   charge tail was incorrectly coded.  The result was that the blend
   occurred far out in the tails of both functions, so that for all
   practical purposes, the model core was "pure Teter."  Correcting this
   to the intended blending region (from the model-core / AE-core crossing
   radius to the first zero of the Teter function) produced significant 
   undesirable structure in the derivatives of the blended model charge,
   likely to cause convergence errors, and especially acoustic sum rule 
   violations in phonon calculations.  The new versions of these routines
   simply use the Teter function for all radii, which essentially reproduces
   the results of 3.2.2 and 3.2.1.

2) The core radii in the input data echo in *.upf files have been corrected.
   Previous versions erroneously reported rcs as shifted to the next-larger
   log mesh point for internal use, which could create inconsistencies
   due to rounding when used to reproduce these potentials.  The original
   input rcs are now truncated to 5 decimals before use, and reported
   to 5 decimals in the *.out file and the *.upf echo, ensuring consistency
   in re-use.

3) A similar echo of the input data is now appended to the end of the
   psp8 files, and is ignored by abinit.  (The actual log-grid-point
   rcs are given as "r_core=" on the first psp8 line, and the input rcs 
   in the trailing echo.)

4) The psfile input datum may now have the values "psp8," "upf," and "both."
   The extract.sh script will now extract both psp8 and upf potential files
   when both are present in the *.out file.  Two tests/data examples have
   been changed (Ge and Hg) to use this option.

5) In keeping with the conventions of ABINIT version 7.10.5 and later,
   the extension_switch in line 6 of psp8 files is now 1 (non-relativistic)
   or 3 (relativistic) to indicate the presence of the valence pseudo-charge
   density (added in 3.1.1).  This must be changed to  0 (nr) or 2 (r) for
   earlier ABINIT versions.

6) The valence charge listing in psp8 files is now supplemented by the
   all-electron valence charge and all-electron core charge in two
   additional columns.  This is intended for auxiliary programs which
   do Bader and related analyses of charge densities.  It is signaled
   by a 2nd integer following the "extension_switch" integer on line
   6 of the psp8 file which is ignored by abinit.

7) Numerous small corrections to enhance inter-platform transferability
   have been incorporated (uninitialized variables, single-precision 
   constants, overreaching array bounds, etc.)
   (Thanks to J. E. Pask)

8) A simple modification of FFLAGS in make.inc allows the use of newer
   libxc releases.  This has been tested with libxc-2.2.1, used in
   abinit-8.0.8.  The list of XC functions and their identifying 
   3-digit numbers in libxc_use.txt has been updated to 2.2.1.
   Be sure to uncomment the appropriate "LIBS +=" line for your
   make.inc.  upfout.f90 has not been updated to include more
   functionals for pwscf (volunteers?).

9) Be aware that two sets of well-tested oncvpsp potentials have been
   published for most of the periodic table and are continuing to be
   expanded:

   http://www.abinit.org/downloads/pseudodojo/pseudodojo

   http://www.quantum-simulation.org/potentials/sg15_oncv/

   The input data is included on these web sites, and in general will
   work fine to drive this release to get an alternative output format,
   a different Exc functional, or relativistic potentials with spin-orbit.

   While I expect the sg15 potentials will be updated to reflect the error
   corrected in item 2 above, be sure you are using the actual input data
   files rather than the echo in their upf files (at least until their rcs
   appear in 5-decimal format).
   
(Thanks for suggestions and feedback to M. van Setten, M. Giantomassi,
 M. Verstraete and F. Gygi.)


## Version 3.2.2

1) The new core features crashed for the relativistic treatment of a 
   partiularly challenging data set one of you came up with.  The solution
   was to use multiplicity-weighted averages of the fully-relativistic 
   state-by-state all-electron charges and valence pseudocharges in the 
   Teter metric analysis.  There is essentially no difference in the output.  

   Efficient coding of this scheme required a significant reorgainzation
   of oncvpsp*.f90 and modcore*.f90, but the underlying math is unchanged
   except as noted above.  Scalar-relativistic output from the built-in
   tests is identical, and fully-relativistic has very insignificant
   differences in the Teter analysis and hence icmod==4.

2) Minor bug fixes.

## Version 3.2.1


1) The shell script tests/data/TEST.sh has been introduced.  This will run
   all the tests with plotting suppressed, and summarize the results
   in TEST.report. The new build procedure (see INSTALL) automatically runs
   this.  (Thanks to M. Verstraete)

2) Several new options and tests have been added for model core charges.
   These are discussed in detail in core_correction.txt in this directory.
   Some tests/data and tests/ref files have been modified to test these.
   (Thanks to M. van Setten and M. Giantomassi for motivating and testing
   these additions)

3) The upf output now includes pseudo wave funtions for the occupied states.
   This will allow LDA+U calculations in some codes.  Please note that 
   ONCVPSP potentials have not been tested or benchmarked for this use.
   Also be aware that LDA+U is a semi-empirical mean field theory, and
   does NOT represent many-body phyysics.  (Thanks to R. Sundararaman
   for providing and testing these code modifications)

4) The ability to run totally non-relativistic calculations was inadverently
   disabled in 3.0.0 and has been restored.  This can be used to verify
   that the small errors discussed following Eq.(24) in my paper arising
   from the scalar-relativistic (and relativistic) issues with general
   norm conservation are corrected.  

   This also applies to the accuracy of convergence of the unscreened 
   potentials to  Coulombic beyond r_c discussed in item 6 below.

5) The automated plots now include log-derivative comparisons for angular
   momenta from lmax+1 to 3.  These are calculated using the local
   potential, and may be used to decide if lmax should be increased.  They
   are placed before the convergence plot.  This motivated changes in
   two tests/data files.

6) An inconsistency in tne unscreening of the pseudopotentials due to the
   slightly longer tail of the all-electron charge density has been fixed
   in an appropriate ad-hoc manner.  This was totally inconsequential and
   only showed up when rlmax was set to a large value. This should never be
   necessary since the (corrected) pseudopotentials are very well converged
   to Coulombic beyond the largest r_c.  This change is in oncvpsp*.f90,
   following the comment "fix unscreening error ..."

   The underlying reason for the non-overlapping tail region arises from
   the fact that the bound-state radial eqution solvers begin their inward
   integration at a radius which is 10X that of the classical turning
   point.  When this point is inside r_c, this will be different for
   the all-electron and pseudopotentials.  The charge densities in
   this region are negligible (~1e-12) but their one-third power in the
   exchange-correlation potential is large enough to notice. (Thanks to 
   C. Fortmann)

7) A typo in a format statement in oncvpsp_r.f90 which bothered some
   compilers has been fixed. (Thanks to G. Marco)

## Version 3.1.1

1) The routine dp3int which performed cubic polynomial interpolation from
   the internal log radial grid to our preferred linear output grid has
   been upgraded to enable nth-order interpolation, with the internal
   parameter npoly presently set to 7.  This allows the polymial model 
   core charge for the non-linear core correction to be exact on the
   linear grid, and eliminates noise when the application code take
   derivatives internally.  dpnint is also used in eresid.  The changes
   in the output are negligible.  The output routines linout, linout_r,
   upfout, and upfout_r have been "cleaned up" to take further advantage
   of dpnint.  (Thanks to Matthieu Verstrate and Alberto Garcia for 
   calling my attention to the external-derivative issue.)

2) The gnu_script* routines have been modified to take advantage of
   features in gnuplot-4.6 (and presumably later) and to improve the
   visibility of the graphics.  With older versions, changing the
   text 'set termoption dash' to '#set termoption dash', line 137
   in gnu_script.f90 and 181 in gnu_script_r.f90 should suffice to
   run the graphics without dashed lines.  Several intermediate points 
   have been added to the "Energy Error" convergence plots. The file
   docs/40_Zr_plots.pdf illustrates the new graphics.

3) The valence pseudo-charge density has been appended to the end of
   the psp8 output for potential future use in abinit or other codes.

## Version 3.0.0

1) The major advance of this version is the contribution by Matthieu
   Verstraete of an interface to the libxc library of exchange-correlation
   functionals.  Instructions for the option of installing this are in the
   INSTALL file.  For consistency with libxc, the built-in exc*.f90
   routines have been redone, basically for more significant figures
   in the constants.  excpwca was updated in this manner and renamed
   excpzca to reflect the fact that it actually computed the Perdew-
   Zunger-Ceperly-Alder lda referenced in its comments.  Changes
   in the ref/*.out files are in the 5-6th decimal place.  Details
   are in the new libxc_use.txt file in this directory.

2) All system-dependent editing in now confined to the make.inc file

3) The scripts directory in the main directory has been introduced.  This
   has copies of the shell scripts in tests, and fldiff.pl, borrowed from
   Abinit for the new compare.sh script which allows easy testing of
   your installation.  The shell script set_path must be run in the
   main directory to set the correct paths in these scripts.

4) Minor bug fixes, none of which detectably influence the output have been
   included, again with thanks to Matthieu Verstraete.

5) Some of the coding standards described in the file in this directory
   have been relaxed to expidite the libxc interface (multiple routines
   in one file, modules, *.F file to be pre-processed), but I would
   prefer that this be an exception rather than a model for future
   contributions.

6) As you must already know, the directory into which the tarball now
   unpacks is tagged with the version number.

## Version 2.1.2

1) A bug has been fixed so that the "First projector wave function 
   outermost peak radius" diagnostic result is calculated correctly for
   positive-energy well-bound states. Nothing else in the output changes. 

## Version 2.1.1

1) The routine sbf8 has been replaced with an improved version that gives
   an overall speedup of program execution by ~X5!

   The original sbf8 was written decades ago (in connection with my development
   of FLAPW) using a recursive algorithm for spherical Bessel functions which
   retains accuracy at large l values encountered in this application.  For
   large values of the argument (not usual there) it becomes very slow.  The
   routine eresid has dominated the run time of this whole code, calling
   sbf8 inside a large double loop with small l and mostly large arguments.
   The small-l cases (l = 1-3) are now treated with the sin-cos expressions
   for spherical Bessel functions, which are much faster for large arguments.
   The new and old versions agree to essentially machine accuracy.

2) I have been informed that abinit-7.6.3 will incorporate the changes to
   use psp8 files with spin-orbit, and this and subsequent releases will no
   longer need a patch.

3) I'm now encouraging users to submit input data for new elements.

4) The non-relativistic script run_nr.sh now creates <prefix>_nr.out.

5) A minor change to the run*.sh and replot.sh scripts eliminated a
   possible problem with gnuplot.

## Version 2.0.2

1) Minor bug fix in pspot: In line 125, finite small eps in if statement
   caused screened semi-local psps to be set to zero at very small radii.

2) Minor bug fix in run_diag*: For some choices of rc and positive 
   second-projector energies, scattering states can be nodeless inside rc. 
   If lschvkbbe returns an error searching for a state with a node
   that matches the all-electron log derivative, it is called again to
   search for a nodeless one before run_diag* reports an error.

Neither of these bugs had any effect on the pseudopotentials in the
examples, but there were a few spurious diagnostic error reports.

3) A few small parameter adjustments were made to the examples to take
   advantage of the new convergende graphics and make the convergence
   more uniform in a few cases

## Version 2.0.1

1) A fully relativistic version giving psps including spin-orbit effects
   has been implemented (executable src/oncvpspr.x, script tests/run_r.sh).
   The same input files drive the non- scalar- and fully-relativistic codes.
   upf-format psp files have been tested with quantumespresso-5.0.4 PWSCF.
   psp8 files for ABINIT-7.4.3 require the supplied patch be applied.
   Hopefully these changes can be incorporated in later abinit releases.
   Details are in doc/relativistic.

2) Ghost detection using the arctan(log derivative) plots has been made more
   robust.  In particular, the occasional occurrence of false-positive
   with the relativistic psp8 potentials, and has no effect on the results
   for solids.  (The OPTION local variable in run_vkb can be reset to restore
   the behavior of earlier releases.)

5) Minor changes in the main routines (oncvpsp*) correct their behavior when
   a negative-energy is specified for the second projector and this is not
   just an occupied valence state over a semi-core being treated as valence.  
   A WARNING will be issued since this is probably inadvertent, and it is
   probably better to adjust debl for a small positive energy.


## Version 1.1.1

No real changes to the code

1) Phys. Rev. B reference now printed in output header and UPF files

2) Minor cosmetic bug-fix in linout.f90

## Version 1.1.0

1) Pseudopotential output in the UPF format has been implemented so that
   ONCVPSP psps can be used in PWSCF calculations.  The input variable
   psfile must now be appended to the first line of the input data with
   values psp8 for ABINIT and upf for PWSCF (www.quantumespresso.org).
   the script tests/extract.sh detects which format is in the *.out file
   and extracts it as before.  Tests comparing PWSCF and ABINIT results
   give excellent agreement.

2) Some related changes are made in the documentation.

## Version 1.0.2

No real changes to code

1) arXiv reference now printed in output header

2) Revised preprint in doc

3) Some labeling corrected in tests/data/<prefix>.dat files

## Version 1.0.1

Very minor changes

1) Main program now prints the radius of outermost peak of all-electron wave
   function for the first projector as a guide for choosing rc for each l.

2) Inconsistencies in the comments in some of the tests/data/*.dat files
   have been corrected.  None of the data has been changed.

3) tests/data/03_Li.dat has been added, but not tested in any calculations 
   for solids.

3) Item 2 under "Potential Issues" in doc/users_guide.txt has been expanded.
