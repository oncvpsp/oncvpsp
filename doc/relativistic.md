﻿## OVERVIEW OF RELATIVISTIC CHANGES

The relativistic code utilizes a new main program, oncvpsp_r, rather than
complicate oncvpsp with lots of branching.  Many arrays have acquired an
additional final index, and most loops over angular momenta include inner
loops over this "ikap=1,2" index referring to the additional quantum number of
the radial Dirac equations kappa (kap) =l, -(l+1) for j=l -/+ 1/2.

The routines ldiracfb and ldiracfs are introduced for bound and scattering
solutions of the radial Dirac equations. The following routines are
straightforward modifications of the old versions, mostly just incorporating
the ikap loops and Dirac equations for full-potential calculations:

* fphsft_r
* fpovlp_r
* gnu_script_r (plot lines are labeled S+, P+, P-, ... indicting j=l +/- 1/2)
* linout_r
* oncvpsp_r
* psatom_r
* renorm_r
* run_config_r
* run_diag_r
* run_phsft_r
* run_plot_r (the Ion potential plots are scalar averages)
* run_vkb_r
* upfout_r
* wellstate_r

The following routines are introduced to serve new functions:

* renorm_r:

Recall the discussion in the paper about the scalar-relativistic Schroedinger
equation being incompatible, strictly speaking, with the derivation of
Vanderbilt's generalized norm-conservation theorem.  The same is obviously
true of the radial Dirac equations.  Two approaches were investigated to
compute the all-electron norms and overlaps used to constrain the pseudo
wave functions (Eq. 23) -- using the small and large components of the
Dirac wave function, and using the renormalized large component only.  The
second approach gave errors in the B matrix non-Hermiticity of the same order
found for the scalar-relativistic case, ~10^-4, with ~10^-5 - 10^-4 errors
in computed quantities after symmetrization.  The first approach gave larger
errors by approximately an order of magnitude, so the second was adapted.
renorm_r performs the large-component renormalization.  All comparisons, wave
functions, etc. are made with straight, unrenormalized Dirac-solution
large components.

    sr_so_r & run_diag_sr_so_r:

The j = l +/- 1/2 non-local pseusopotentials V(r,r') can be used directly
for input to PWSCF in the UPF format.  ABINIT is structured to work with
separate scalar and spin-orbit projectors.  While the decomposition of
the two j potentials into the weighted-average and simple difference for
scalar (sr) and so is simple for semi-local potentials, this is cumbersome and
not advisable for two-projector fully-nonlocal potentials.  It results
in 4 projectors each for sr and so, with large subtractions forced upon
the electronic structure code.  The routine sr-so_r re-diagonalizes each
of these operators, creating different sr and so projectors.  Some of
the coefficients of these orthonormal projectors are very small and can
be neglected.  The threshold eprmin in oncvpsp_r is presently set to
2E-5 Ha, but can be easily changed if desired.  Another advantage of the
decomposition is that the strength of the so term can be artificially
scaled.  While no additional input variable is provided for this, to
maintain compatibility  with the scalar-relativistic data structure, the
variable soscale, currently 1.0, can be reset or read in after the
standard data.  The routine sr_sr_r performs the decomposition, and
run_diag_sr_so re-runs the diagnostics using a multi-projector version
of vkboutwf.f90 and the truncated versions of the sr-so projectors.

The psp8 format (introduced originally into ABINIT by me ten years or so
ago) was intended for increased flexibility, along with some other
features I considered desirable.  It incorporated a header line consisting
of the integer zero with the label "extension switch."  For relativistic
psp8 files now produced by ONCVPSP, the integer is two. To read these,
the supplied patch must be applied in the main abinit-7.4.3 directory
and the code recompiled.  With the patch, the code will go on to read
a line giving the numbers of spin-orbit projectors per angular momentum
channel.  It will then read the scalar projectors and local potential,
and finally the spin-orbit projectors.  ABINIT will run with 8 projectors
(4 for scalar and 4 for spin-orbit) per angular momentum.  It's been
tested with band structure and phonon calculations.  The patched ABINIT
code is fully compatible with older psp8 files.  Abinit versions 7.6.3
and beyond have incorporated these changes and don't need patching.
