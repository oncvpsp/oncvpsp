﻿## OVERVIEW OF ONCVPSP

References to "the paper" are to oncvpsp6.pdf, which is reproduced in this
directory (now published as Phys. Rev. B 88, 085117 (2013)).

The main program is oncvpsp.f90, and has the following major sections:

1) Data readin: Oncvpsp reads an input file, say *32_Ge.dat*, from its standard
   input.  While this file may contain comment lines with the initial character
   '#', its format is fixed, and input data must appear in the proper order,
   line by line.  Examples in tests/data as well as *32_Ge_annotated.dat* in
   this directory should clarify this.  The program will stop, identifying
   the first bad line encountered in the input file if it does not find what
   it expects.  After a successful readin, the data is scanned by the
   routine *check_data* for violations of bounds, etc., and offending data will
   be identified prior to halting.  There are, quite deliberately, no defaults.
   Everything must be thought about and specified (but there is not that much
   input data).

2) The routine sratom will then run to perform a self-consistent scalar-
   relativistic calculation for the reference configuration of the atom,
   as described in the input.  I strongly believe that neutral ground-state
   reference configurations provide the best pseudopotentials for solids,
   and oncvpsp incorporates a number of features to facilitate their use.
   After some bookkeeping, the main program echoes the input data in a format
   similar to a typical input file, with the addition of the reference-
   configuration eigenvalues to the configuration description.  All output
   printing is to the standard output.

3) Next the main program loops through angular momenta from 0 to lmax to
   calculate one or two optimized pseudo wave functions.  For each angular
   momentum, there are four major sub-sections:

  A) The all-electron wave function is calculated for the lowest occupied
     state declared as valence or for a barrier-confined "well state" for a
     specified positive energy corresponding to a scattering state
     (subroutines lschfb and wellstate).  The barrier potential is given
     in Eq.(17) of the paper.  To find a bound state at the desired
     energy ep, wellstate starts with a large r_b = 8*r_c and searches
     for a suitable v_inf within the range [ep, ep + 1 Ha].

  B) If two projectors using the norm-conserving Vanderbilt method are
     specified, the all-electron wave function with one additional node is
     calculated for either a second bound valence state (when the first
     is a shallow core) or for a second "well state" at a specified energy
     above the lowest state at the current angular momentum.

  C) One or two pseudo wave functions are calculated using the RRKJ
     convergence optimization approach.  Both will satisfy continuity
     and norm-conserving conditions with their corresponding all-electron
     wave functions.  The second will have one node and will also satisfy
     the off-diagonal Vanderbilt generalized norm-conserving condition,
     which is dependent of the first pseudo wave function as well as the
     overlap of the two all-electron wave functions.  Details of the
     routine run_optimize are expanded below. Diagnostic information is
     printed.

  D) The semi-local pseudopotentials are calculated.  More diagnostics
     and the convergence profile (Energy error per electron vs. cutoff)
     are printed.

4) The routine *run_vkb* then is then called.  It constructs the local
   potential from a polynomial extrapolation of the all-electron potential
   to zero (the recommended option, especially to take advantage of the
   Vanderbilt two-projector method), or uses one of the semi-local
   potentials, which may be satisfactory when using a single Kleinman=
   Bylander projector.  It then loops over angular momenta, calculating
   the "raw" projectors (chi in Eq.(18) of the paper) and the B matrix
   (Eq.(21)).  For two projectors, any residual asymmetry  of B is
   removed (see the discussion following Eq.(24)), it is diagonalized,
   and the final diagonal and normalized separable non-local potential
   of Eq.(24) is produced.  A small amount of diagnostic information is
   printed.  As of release 2.0, orthonormal projectors are subsequently
   constructed, and overlap information is omitted.  The "OPTION"
   variable in *run_vkb* can be reset to restore the old non-orthogonal
   output, although only for non- and scalar-relativistic cases.

5) The valence pseudocharge density is calculated based on the separable
   potential, and then that of any "extra" valence state (ie., when cores
   are included), using subroutine lschvkbb.  If specified, a monotonic
   polynomial model core charge is calculated which joins smoothly to
   the all-electron core charge at the radius at which this charge density
   equals a specified fraction of the valence pseudocharge density
   (subroutine modcore).  The Hartree and exchange-correlation potentials
   of the pseudocharge are computed (subroutine vout)and the pseudopotentials
   are "unscreened."

6) The routine *run_diag* is next called to test the semi-local and separable
   pseudopotentials.  For bound states, the normalized bound-state
   eigenvalues and pseudo wave functions are calculated.  The eigenvalues
   and amplitudes and first derivatives at `r_c` are compared to the all-
   electron results (as ratios).  For positive-energy states, the all-
   electron log derivative is computed at *r_c*, and the "eigenvalue" of
   the pseudo wave function which matches this is calculated.  Value and
   slope ratios (equal by construction in this case) are printed.
   For Vanderbilt two-projector calculations, these calculations are
   repeated at the energies used to construct the second projector.
   "Extra" valence states are treated as bound states here.

7) The routine *run_config* performs all-electron and VKB pseudopotential
   self-consistent atom calculations on several configurations.
   Configuration 0 is always the reference configuration, and minimal
   errors here are primarily a consistency check.  Up to 4 additional
   valence configurations may be specified in the input data.
   Eigenvalues are compared to all-electron results, and energy differences
   relative to the reference configuration indicate the excitation energy
   error of the pseudopotential.  The VKB calculations are performed with
   the specified number of projectors and non-linear core correction.
   This routine duplicates a capability of the OPIUM code, but its utility
   is not clear (see paragraph 3 of Sec. V of the paper).

8) The routine *run_plot* generates output that will be automatically
   plotted by gnuplot when running the main shell script.  These include the
   unscreened semi-local potentials, with the local potential added in
   the more general case of a polynomial interpolation.  The pseudo-valence
   charge, core charge (up to a cutoff), and model core charge if calculated
   are plotted next.  Then, all-electron and pseudo wave functions are
   compared, looping from l=0 to l = lmax.  In the shallow core with occupied
   valence case, both wave functions are computed, the single-node pseudo-
   valence wave function being calculated with the specified number of VKB
   projectors.  For two projectors, the first- and second-projector pseudo
   wave functions are printed for comparison.  The cutoff energy convergence
   predictions based on the residual kinetic energy are output in form
   suitable for plotting as of version 2.0.  The *run_plot* output is also
   printed on the standard output below the heading "DATA FOR PLOTTING."

9) The routine *run_phsft* compares the log derivatives of the all-electron
   and pseudo wave functions (bound or scattering) at the maximum *r_c*.
   (When shallow cores are included, this value shifts to max(rc)+1.)
   Since simple plots of log derivatives contain multiple divergences, and
   are difficult to interpret visually, I plot atan(r * ((d psi(r)/dr)/psi(r)))
   at this radius, which is akin but not equal to a scattering phase shift.
   These are plotted using the specified number of vkb projectors.  If the
   local potential is that of a single l, no projectors are used for this l.
   These results are printed to the standard output, and plotted after
   the wave functions plots for the corresponding l.  A smoothed downward step
   ~pi in these plots at negative energies indicates a bound state.  The energy
   range of these plots is specified in the input data, and starting the
   range at an energy somewhat below the lowest relevant core state energy
   serves as a ghost test.  A step in the pseudopotential plot which is not
   matched by one in the all-electron plot indicates a ghost.  The
   pseudopotential results are shifted by an appropriate multiple of
   pi to align the plots at the first reference energy.

10) The routine gnuscript produces an input script for gnuplot, which is
    extracted and appropriately handled for plotting by the main shell
    script tests/run.sh.

11) The routine linout creates a linear radial mesh, whose size and spacing
    are specified in the input data.  All relevant potentials and projectors
    are interpolated onto this mesh, which is printed to the output file in
    a format that can be read by ABINIT as its psp format 8.  The linear
    mesh has advantages of economy for pseudopotentials at small radii, and
    of greater accuracy for the Fourier transforms that will by performed by
    ABINIT at larger radii.  In addition, format 8 accepts numerical
    projectors and their coefficients, which allows the multi-projector
    flexibility that ONCVPSP achieves to be used.  The script extract.sh
    transfers this part of the single ONCVPSP output file to an ABINIT
    pseudopotential file in a specified directory.  The routine upfout
    uses the same linear mesh but produces a PWSCF-conpatible pseudopotential
    file in the UPF format.  Which one is called depends on the "psfile"
    input variable.  While all the psp files presently linked to PWSCF use
    a log mesh, the linear mesh works fine and produces results in excellent
    agreement with ABINIT results.


The auxiliary main program *oncvpsp_nr.f90* is identical to oncvpsp except it
sets a switch to do non-relativistic all-electron calculations.  The only
reason for including this feature is to support the discussion following
Eq.(24) of the paper. The script *tests/run_nr.sh* calls it.

The subroutine *run_optimize* executes the math in Sec. II of the paper for
each l, and has the following steps:

1) The routine *wf_rc_der* calculates the value and 4 derivatives of the
   input all-electron wave function at *r_c*.

2) The routine qroots selects the wave vectors *q_i* for the xi^B basis,
   Eq.(3) according to the principles discussed following Eq.(10).

3) The routine *sbf_basis* creates the xi^O basis using Eq.(4).

4) The routine *Const_basis* executes Eqs.(6-10), producing the xi^N basis
   and the *phi_0* component of the pseudo wave function.  Diagnostics
   confirming that the constraints are satisfied are printed.

5) The routine eresid calculates the matrix elements in Eq.(11) of the
   E^r "operator" defined in Eqs.(1-2).  The q integration is performed
   inward from a large value ("infinity"), saving snapshots of the matrix
   elements at regular intervals as well as the values at the specific
   specified *q_c* cutoff.  These snapshots allow the convergence profile
   of the optimized pseudo wave function to be calculated very economically.
   This is the single most time-consuming step in the code, and setting
   the integration increments dq and dr to larger values will save time
   at some expense in accuracy.  These are set in the first executable
   statements of *run_optimize*.  With the new sbf8 routine of release 2.1.1
   you shouldn't need to  bother.

6) The routine optimize diagonalizes the E^r_ij matrix, calculates the
   xi_R basis functions, and the coefficients in Eq.(14).  It then
   executes the minimization procedure discussed in the paragraph
   surrounding Eqs.(15-16), giving the optimum x_i in Eq.(13).
   It finds the coefficients of the optimized pseudo wave function
   in the xi^N basis and uses the "snapshot" values of the terms in
   Eq.(11) to compute the convergence profile.  Finally, it finds
   the coefficients in the xi^O and xi^B basis sets.

7) The routine pspot calculates the semi-local potential and the
   (epsilon_i - T)|ph_i> component of the |chi_i> projector, Eq.(18).
   It also calculates the total kinetic energy from the Laplacian to
   be compared with the Fourier-transform result (the q_c = 0 end
   of the "snapshot" convergence profile).

At this point, run_optimize exits if only one projector is wanted.
If a second projector is wanted, it repeats step 2 above, computing
derivatives of the second all-electron wave function.  The existing
xi^B basis appears to work well for the second projector, so step 3
(choosing a set of q_i) is skipped.  The constraint vector and matrix,
Eqs.(6-7) is then incremented following the discussion surrounding
Eq.(25), and steps 4-7 are completed to provide an optimized second
projector.

*run_optimize* prints some additional diagnostic comparisons following
several of the subroutine calls, and prints the plane-wave energy
cutoffs corresponding to kinetic-energy-per-electron errors of 10^-2 to 10^-5 Ha.
