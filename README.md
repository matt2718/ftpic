# Particle-in-Fourier plasma simulation

This is a proof of concept for efficient Fourier-space particle simultion. For
more information about the algorithm, see <https://arxiv.org/abs/1808.03742>.
Note that, although this repository and several of the code samples are titled
"FT-PIC", the name of the algorithm has been changed to "particle-in-Fourier"
(PIF).

This code is compiled using the USFFT library, which is not open source. In
order to build it, you'll need a copy of the compiled USFFT library, or you'll
need to rewrite the code to use another library that provides support for fast
nonuniform FFTs (which we might do in the future).

You'll also need the QDSP library to view phase plots, which can be obtained
from <https://github.com/matt2718/qdsp>. A dummy QDSP header is provided in the
`include` directory for building on headless systems.

Two source files are provided: `ftpic.c` contains our fast gridless algorithm,
while `oldpic.c` contains a reference 1D PIC code with first-order deposition
and interpolation for the purposes of comparison with PIF. They use the same
options, which you can view with the `-h` option.

## Analysis scripts

- `eplot.py` displays a plot of kinetic, field, and total energy over time based
  on the output file of a single simulation.

- `modeplot.py` displays a logarithmic plot of the growth or damping of one or
  more modes from one or more output files.
  
- `ecomp.py` plots the total energies over time for several simulations.
  
- `npng.sh` runs a scaling study varying the number of particles and the number
  of grid cells/modes.

- `threadscale.sh` runs a scaling study varying the number of threads.

- `dtscale.sh` and `dtscale2d.sh` run convergence studies by varying the
  timestep size, for 1D and 2D programs respectively. The same scripts can be
  used to analyze the output of convergence studies with the `-a` flag.

- `getmax.py` prints the maximum difference in total between two timesteps for
  some output file.
  
- `convplot.py` displays a convergence plot based on the output of `dtscale.sh`
