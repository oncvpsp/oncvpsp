## ONCVPSP USERS' GUIDE

### Introduction

There is no detailed users' manual, and this is primarily an introductory 
tutorial for new users.  With the included Phys. Rev. paper, the
program overview, the simple structure of the input data, the examples
provided, and the automated graphics, it should be relatively easy to get
to know and use this code.  As of release 4.0.1, most of the potental
problems disussed in the last few paragraphs are highly improbable.

You should start by reading the paper, and some of the references if
necessary.  The input and output files in tests/data and tests/refs include
all the "OV" pseudopotentials used in the Sec.IV tests in the paper. Next, 
I recommend reading ../INSTALL, and [program_overview](program_overview.md).

The gnuplot graphics are very important, and a good place to go next.
The file 73_Ta_plots.pdf in this directory is a "slide show" which walks you
through a sample of what you will see whenever you run the program.  After
viewing the pdf, go to directory ../tests.refs and type

    ../replot.sh 73_Ta

If this doesn't produce the same sequence of plots, stop here and try to fix 
your gnuplot installation and accessibility.  Once this works,  it will pay to 
run replot for all the examples in this directory.  This should illustrate what 
your own output plots should look like.  If your own data give results that
look odd by these standards, they probably are and should not be used until
adjustments of your input data can produce "nice" plots as exemplified by
these cases.

The input files for oncvpsp, examples of which you'll find in
../tests/data, are NOT THE KEYWORD FILES YOU MAY BE USED TO.  The # comment
lines define a strict template in which ALL input parameters must be specified,
in the correct order.  There are no defaults.  (On the other hand,  there are
not so many input parameters, and you need to understand them all to get good results anyway.)

The file `32_Ge_annotated.dat` adds more comments (###) to the usual
`*.dat` file describing each input variable.  Reading this should be your next
assignment.  The code doesn't mind extra comments, and you can in fact
use this as a combined template and abbreviated users' manual for your
own calculations instead of the shorter versions in ../tests/data.

Next, look at `tests/refs/32_Ge.out`.  The output is pretty-much
self-explanatory.  First, the input data is echoed, adding the all-electron
eigenvalues to the configuration information.  This is formatted so that it
can be copied and re-used as an input file (the eigenvalues will be ignored).
The code then loops over angular momenta and reports some information and
consistency test results for each pseudo wave function. For unbound states,
the 'half-point radius' is that at which the added barrier has reached half
its asymptotic value.  The qq values are the generalized norm-conserving
conditions, Eq.(23).  The constraint-consistency test should be very well
satisfied, or there is an error in the code. The following output is just
information ("Fraction of norm ...").  For the kinetic energy and potential
tests, the ratio should be close to one, and the difference should be small.
Last comes the cutoff convergence profile for the current wave function.

After the pseudo-wave-function loop is completed, the projectors are
constructed.  The Hermiticity error is discussed in the paragraph following
Eq.(24).  It will be greatly reduced in a non-relativistic run. The data on
the projectors is not of much interest, but if the coefficients are huge or
the overlap is very close to one, it probably pays to change some parameters.
As of release 2.0.1, orthonormal projectors are constructed by default,
so the overlap is no longer reported.

The "Diagnostic tests" results are explained in more detail in
program_overview.txt.  Basically, e in and e test should closely match, and
the norm and slope values should be close to one.  For two projectors and a
non-relativistic calculation, the agreement should be within numerical
convergence error for all these quantities.  Finally,  the total energies
and eigenvalues of self-consistent pseudo atoms are compared with all-electron
results for the reference and test valence configurations.  You should
always look at these diagnostics, even if the automated graphics all look
good.

Everything below the "DATA FOR PLOTTING"line is for the scripts to read,
not you (plotting data, template for a custom gnuplot input script, and
the psp file).

After you've gotten through the above and have compiled the code
following INSTALL, copy one of the `<prefix>.dat` files in tests/data to a 
directory of your own, run it using

    ../tests/run.sh <prefix>

The sample `<prefix>.dat` files in tests/data should provide a starting 
point to produce good psps for new elements. Reading the `<prefix>.out` files,
you should be sure that the various consistency tests look good compared to
the examples. Then, pay attention to the plots.

### Potential issues

1) The code will stop if it doesn't find an input data entry it is looking for.
   It will also run some checks to be sure the data is within bounds. The first
   stop will just give you a line number, but the error might be above that.
   The second kind should give more explanation.

2) The code will stop for a number of other error conditions you shouldn't
   encounter if it is built properly and your parameters are reasonable.  
   All stops should give some information.  One stop you probably will 
   encounter is 

        "pspot:  first pseudo wave function has node."

   While in principle the Vanderbilt method can deal with this, it isn't a 
   good idea and it's better to have a semi-local potential to look at to make
   sure it looks reasonable.  The fix in this case is usually to increase qcut
   a bit, although changes in nbas, rc or ep for a scattering state may also
   work.

   Another you might encounter is

        "optimize: norm constraint violated (line 260)"

   which will probably only occur if you have an unphysically small rc,
   indicated by a very small value for "Fraction of norm inside rc."  The
   code now gives "First projector wave function outermost peak radius," 
   and rc should be > that value (usually >1.5 * r_peak, see examples).

   Another is 

       "wellstate: well potential iteration failed"

   which should be fixable by changing ep (for the first projector) or
   a combination of ep and debl for the second.  rc can also influence
   this.

3) A related problem is an unreasonably large semi-local 1=0 potential at r=0.
   Once again, this bad behavior won't carry over into the projector, and
   might not do any real damage (unless lloc=0,  which is always a bad idea
   anyway).  However, it is best to eliminate it by parameter changes like
   those for the nodes.

4) The next-to-last paragraph of Sec. II of the paper remarks on the
   "strong short-wavelength oscillations" discussed by Troullier and
   Martins as an issue with the Rappe et al. optimization scheme.  While
   I have not observed these for the usual ranges of parameters in the
   potentials I've generated with the present optimization approach,
   it can be made to happen.  It is certainly a good idea to eliminate them
   by adjusting the usual suspects (qcut, nbas, rc, ep).

5) The ncon continuity-constraints may be less important than I initially
   thought, and don't have that much influence on the results (see
   Fig. 3).  I prefer to use ncon = 4 or 5.  When looking at pseudo
   wave functions, a higher the degree of continuity does allows a larger 
   rc for comparable agreement with the all-electron wave function.

6) Choosing qcut is a compromise.  While a smaller qcut will typically give
   better convergence up to the corresponding energy, the "foot" of the
   convergence curve  (see Fig. 4) will have an earlier onset and be flatter.
   The cutoff convergence profiles in the output (and in the graphics
   as of release 2.0) can guide you.  Increasing qcut by a small amount,
   say ~0.5, can usually give more uniform convergence at only a small
   sacrifice of low-cutoff convergence.

7) I have had trouble with pseudo atom calculations failing to converge in
   the configuration-testing section.  The code won't stop, and will skip
   the "bad" configuration.  These problems seem to arise from the pseudo
   wave functions acquiring nodes at small radii, which can happen for the
   fully-non-local radial Schroedinger equation.  My present solution for
   this is to only start counting nodes for r>0.5 a_B (in the routine vkboutwf).

8) The shell scripts run.sh run_nr.sh etc. assume everything has been left
   in in place the directroy tree in which is was unpacked and built.  
   If you move the executables and scripts elsewhere, you can probbly
   use set_path to fix things.
