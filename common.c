// this file contains common functions used by both the reference PIC
// implementation and Particle in Fourier
#include <math.h>
#include <omp.h>

#include "common.h"

const double XMAX = 16.0; // system length
const int NGRID = 256; // grid size
double DX;

// particle number and properties
const int PART_NUM = 5000;
const double PART_MASS = 0.005;
const double PART_CHARGE = -0.01;
const double EPS_0 = 1.0;

const double BEAM_SPEED = 8.0;

double OMEGA_P;

// wave periods per system length for standing wave and landau damping
const int WAVE_MODE = 3;

// 2-stream instability, standard test case
void init2Stream(double *x, double *v, int *color) {
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
void initStanding(double *x, double *v, int *color) {
	OMEGA_P = sqrt((PART_NUM/XMAX) * (PART_CHARGE*PART_CHARGE)
	               / (PART_MASS * EPS_0));
	
	double ampl = 0.3;
	
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;
		x[i] += ampl * XMAX  / (2 * M_PI * WAVE_MODE)
			* cos(2 * M_PI * WAVE_MODE * x[i] / XMAX);

		v[i] = ampl * OMEGA_P * XMAX / (2 * M_PI * WAVE_MODE)
			* sin(2 * M_PI * WAVE_MODE * x[i] / XMAX);

		color[i] = 0x0000ff;
	}
}

// helper function for initLandau
// use newton's method to find inverse of position CDF
static double newton(double u) {
	double k = WAVE_MODE * 2*M_PI/XMAX;
	double ampl = 0.4;
	
	// find root of x/L + a/(k L) sin(k x) - u = 0
	// derivative is 1/L + a/L cos(k x)
	double x = u * XMAX;
	double f = x/XMAX + ampl / (k * XMAX) * sin(k * x) - u;

	while (fabs(f) > 1e-9) {
		double fprime = 1/XMAX + ampl / XMAX * cos(k * x);
		x -= f / fprime;
		f = x/XMAX + ampl / (k * XMAX) * sin(k * x) - u;
	}

	return x;
}

// wave in maxwellian, used to test Landau damping
// distribution function: f(x,v) = exp(-1/2 (v/v_th)^2) / (sqrt(2 Pi) * v_th)
//                                 * (1 + a cos(k x)) / L
void initLandau(double *x, double *v, int *color) {
	double vth = 1; // sqrt(kt/m)

	for (int i = 0; i < PART_NUM; i++) {
		// charge distribution produced via inverse transform sampling
		x[i] = newton(i / (double)PART_NUM);

		// init maxwellian distro
		// box-mueller
		double r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		double r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		v[i] = vth * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
		
		color[i] = 0x0000ff;
	}
}
