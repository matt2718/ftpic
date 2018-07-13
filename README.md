# FT-PIC

This is a proof of concept for efficient Fourier-space particle simultion.
FT-PIC is probably a bad name, since anything with a spectral field solve
involves a Fourier transform.

This code is compiled using the USFFT library, which is not open source. In
order to build it, you'll need a copy of the compiled USFFT library, or you'll
need to rewrite the code to use another library that provides support for fast
nonuniform FFTs (which we might do in the future).

You'll also need the QDSP library to view phase plots, which can be obtained
from <https://github.com/matt2718/qdsp>. A dummy QDSP header is provided in the
`include` directory for building on headless systems.

Two source files are provided: `ftpic.c` contains our fast gridless algorithm,
while `oldpic.c` contains a reference 1-D PIC code with first-order deposition
and interpolation for the purposes of comparison with FT-PIC. They use the same
options, which you can view with the `-h` option.

You can view a plot of the kinetic and field energy of either program by piping
it into `eplot.py`. You can also use `modeplot.py` to look at various modes over
time. More documentation is coming soon.
