// this file contains common functions used by both the reference PIC
// implementation and Particle in Fourier
#include <math.h>
#include <omp.h>

#include "common.h"

// helper function for initRhoSin
static double bisect(double x1, double x2, double y, double XMAX) {
	double xmid = (x1 + x2) / 2;
	double ymid = xmid + 0.25 * XMAX/(2*M_PI) * (1 - cos(2 * M_PI * xmid/XMAX));
	if (ymid - y > 1e-9) return bisect(x1, xmid, y, XMAX);
	else if (ymid - y < -1e-9) return bisect(xmid, x2, y, XMAX);
	else return xmid;
}

// 2-stream instability, standard test case
void init2Stream(double *x, double *v, int *color, int PART_NUM,
                 double XMAX, double BEAM_SPEED) {
	
	double stddev = sqrt(500 / (5.1e5));
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;

		if (i % 2) {
			v[i] = BEAM_SPEED;
			color[i] = 0xff0000;
		} else {
			v[i] = -BEAM_SPEED;
			color[i] = 0x0000ff;
		}

		// box-mueller
		double r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		double r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		v[i] += stddev * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
	}
}

// standing wave, displaced charges. this should oscillate as a single mode, but
// this doesn't happen in standard PIC. see section 4 of:
// Huang, et al, 2016. Finite grid instability and spectral fidelity of the
// electrostatic Particle-In-Cell algorithm. Computer Physics Communications
// 207, 123–135.
void initDisplace(double *x, double *v, int *color, int PART_NUM,
                  double XMAX, double OMEGA_P) {
	int mode = 9;
	double ampl = 0.1;
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;
		x[i] += ampl * XMAX  / (2 * M_PI * mode)
			* cos(2 * M_PI * mode * x[i] / XMAX);

		v[i] = 0.01;
		v[i] = ampl * OMEGA_P * XMAX / (2 * M_PI * mode)
			* sin(2 * M_PI * mode * x[i] / XMAX);

		color[i] = 0x0000ff;
	}
}

// sinusoidal charge dist with no initial velocity
// not super useful as a test case, but can make sure the code isn't broken
void initRhoSin(double *x, double *v, int *color, int PART_NUM, double XMAX) {
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = bisect(0, XMAX, i * XMAX / PART_NUM, XMAX);
		v[i] = 0;
		color[i] = 0x0000ff;
	}
}
