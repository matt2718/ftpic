// this file contains common functions used by both the reference PIC
// implementation and Particle in Fourier
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common2d.h"

// time info
double DT = 0.001;
double TMAX = 10;

// system size
const double XMAX = 32.0;
const double YMAX = 32.0;

// grid size and spacing
int NGRIDX = 16;
int NGRIDY = 16;
double DX, DY; // spacing

// particle number and properties
int PART_NUM = 10000;
const double PART_MASS = 0.5;
const double PART_CHARGE = -2.0;
const double EPS_0 = 1.0;

double OMEGA_P;

const double BEAM_SPEED = 16.0;
const double V_TH = 4.0; // sqrt(kt/m)

// wave periods per system length for landau damping
const int LANDAU_MODE = 2;
const double PERTURB_AMPL = 0.25;

const int WAVE_MODE = 9;

// logging
const int MODELOG_MAX = 32;
FILE *modeLog = NULL;
int printTime = 0;

int commonInit(int argc, char **argv,
               double **x, double **y, double **vx, double **vy, int **color,
               QDSPplot **xyPlot) {

	FILE *paramLog = NULL;
	
	// whether to plot
	int quiet = 0;

	// which test case?
	int initType = 0;
	
	// parse arguments
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-h")) {
			// help message
			printf("Usage: { oldpic | ftpic } -c CASE [OPTIONS]\n"
			       "  -c CASE   Specific test case to run. Choices are:\n"
			       "              2stream : 2-stream instability\n"
			       "              landau : Landau damping\n"
			       "              wave : standing wave in cold plasma\n"
			       "  -np NUM   Number of particles to simulate.\n"
			       "  -ng SIZE  Number of grid cells (for PIC) or modes (for FT-PIC).\n"
			       "  -t DT,MAX Time step size and length of simulation.\n"
			       "  -q        Quiet, do not render plots.\n"
			       "  -pt       Print the wall time per step.\n"
			       "  -p FILE   Dump parameters to file FILE.\n"
			       "  -m FILE   Dump E-field modes at each step to file FILE.\n"
			       );
			return 1;
		} else if (!strcmp(argv[i], "-c")) {
			// test case type
			if (++i == argc) return 1;

			if (!strcmp(argv[i], "2stream")) initType = 1;
			else if (!strcmp(argv[i], "landau")) initType = 2;
			else if (!strcmp(argv[i], "wave")) initType = 3;

		} else if (!strcmp(argv[i], "-np")) {
			// number of particles
			if (++i == argc) return 1;
			sscanf(argv[i], "%d", &PART_NUM);

		} else if (!strcmp(argv[i], "-ng")) {
			// number of grid cells/modes
			if (++i == argc) return 1;
			sscanf(argv[i], "%d", &NGRIDX);
			sscanf(argv[i], "%d", &NGRIDY);

		} else if (!strcmp(argv[i], "-t")) {
			// time step and limit
			if (++i == argc) return 1;
			sscanf(argv[i], "%lf,%lf", &DT, &TMAX);

		} else if (!strcmp(argv[i], "-q")) {
			// quiet, do not render plots
			quiet = 1;
			
		} else if (!strcmp(argv[i], "-pt")) {
			// print the time per step
			printTime = 1;
			
		} else if (!strcmp(argv[i], "-p")) {
			// log file for parameters
			if (++i == argc) return 1;
			paramLog = fopen(argv[i], "w");

		} else if (!strcmp(argv[i], "-m")) {
			// log file for modes
			if (++i == argc) return 1;
			modeLog = fopen(argv[i], "w");
		}
	}

	if (initType == 0 || NGRIDX <= 0 || DT <= 0 || PART_NUM <= 0) {
		if (paramLog != NULL) fclose(paramLog);
		if (modeLog != NULL) fclose(modeLog);
		return 1;
	}

	DX = XMAX / NGRIDX;
	DY = YMAX / NGRIDY;

	*x = malloc(PART_NUM * sizeof(double));
	*y = malloc(PART_NUM * sizeof(double));
	*vx = malloc(PART_NUM * sizeof(double));
	*vy = malloc(PART_NUM * sizeof(double));
	*color = malloc(PART_NUM * sizeof(int));

	// initialize particles
	if (initType == 1) {
		init2Stream(*x, *y, *vx, *vy, *color);
	} else if (initType == 2) {
		initLandau(*x, *y, *vx, *vy, *color);
	} else if (initType == 3) {
		initLandau(*x, *y, *vx, *vy, *color);
		//initWave(*x, *v, *color);
	}
	
	// dump parameters
	if (paramLog) {
		// basic
		fprintf(paramLog, " particles: %i\n", PART_NUM);
		fprintf(paramLog, "  timestep: %e\n", DT);
		fprintf(paramLog, "       L_x: %e\n", XMAX);
		fprintf(paramLog, "       L_y: %e\n", YMAX);
		fprintf(paramLog, "    v_beam: %e\n", BEAM_SPEED);
		fprintf(paramLog, "      mass: %e\n", PART_MASS);
		fprintf(paramLog, "    charge: %e\n", PART_CHARGE);
		fprintf(paramLog, "     eps_0: %e\n", EPS_0);
		fprintf(paramLog, "\n");

		// Debye length
		double kt;
		if (initType == 1) kt = PART_MASS * BEAM_SPEED * BEAM_SPEED;
		else kt = PART_MASS * V_TH * V_TH;
		
		double dens = PART_NUM / XMAX;
		double ne2 = (dens * PART_CHARGE * PART_CHARGE);
		fprintf(paramLog, "    lambda: %e\n", sqrt(EPS_0 * kt / ne2));

		// plasma frequency
		OMEGA_P = sqrt(ne2 / (PART_MASS * EPS_0));
		fprintf(paramLog, " frequency: %e\n", OMEGA_P);

		fclose(paramLog);
	}

	// header for modes
	if (modeLog) {
		fprintf(modeLog, "time");
		for (int i = 1; i <= MODELOG_MAX; i++)
			fprintf(modeLog, ",m%d", i);
		fprintf(modeLog, "\n");
	}

	// set up plots
	if (!quiet) {
		*xyPlot = qdspInit("Particles");
		qdspSetBounds(*xyPlot, 0, XMAX, 0, YMAX);
		qdspSetGridX(*xyPlot, 0, 2, 0x888888);
		qdspSetGridY(*xyPlot, 0, 2, 0x888888);
		qdspSetPointColor(*xyPlot, 0x000000);
		qdspSetBGColor(*xyPlot, 0xffffff);
	}
	
	return 0;
}

// 2-stream instability, standard test case
void init2Stream(double *x, double *y, double *vx, double *vy, int *color) {
	double stddevx = sqrt(500 / (5.1e5));
	double stddevy = sqrt(0 / (5.1e5));
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;
		y[i] = YMAX * (rand() + 1) / ((double)RAND_MAX + 1);

		//x[i] = XMAX/2 + 0.1 * x[i] - 0.05 * XMAX;
		//y[i] = YMAX/2 + 0.1 * y[i] - 0.05 * XMAX;
		/*
		double xoff = XMAX * ((rand() + 1) / ((double)RAND_MAX + 1) - 0.5) / 4.0;
		double yoff = YMAX * ((rand() + 1) / ((double)RAND_MAX + 1) - 0.5);
		x[i] = XMAX/2 + xoff;
		y[i] = YMAX/2 + yoff;
		*/
		if (i % 2) {
			vx[i] = BEAM_SPEED;
			color[i] = 0xff0000;
		} else {
			vx[i] = -BEAM_SPEED;
			color[i] = 0x0000ff;
		}
		vy[i] = 0;
		/*
		vy[i] = 0;
		vx[i] = 0;
		*/
		//
		// box-mueller

		double r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		double r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		vx[i] += stddevx * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
		r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		vx[i] += stddevy * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);

	}
}

// helper function for initLandau
// use newton's method to find inverse of position CDF
static double newton(double u) {
	double k = LANDAU_MODE * 2*M_PI/XMAX;
	
	// find root of x/L + a/(k L) sin(k x) - u = 0
	// derivative is 1/L + a/L cos(k x)
	double x = u * XMAX;
	double f = x/XMAX + PERTURB_AMPL / (k * XMAX) * sin(k * x) - u;

	while (fabs(f) > 1e-9) {
		double fprime = 1/XMAX + PERTURB_AMPL / XMAX * cos(k * x);
		x -= f / fprime;
		f = x/XMAX + PERTURB_AMPL / (k * XMAX) * sin(k * x) - u;
	}

	return x;
}

// wave in maxwellian, used to test Landau damping
// distribution function: f(x,v) = exp(-1/2 (v/v_th)^2) / (sqrt(2 Pi) * v_th)
//                                 * (1 + a cos(k x)) / L
void initLandau(double *x, double *y, double *vx, double *vy, int *color) {
	for (int i = 0; i < PART_NUM; i++) {
		// charge distribution produced via inverse transform sampling
		x[i] = newton(i / (double)PART_NUM);

		// uniform in y
		y[i] = YMAX * (rand() + 1) / ((double)RAND_MAX + 1);

		// init maxwellian distro
		// box-mueller
		double r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		double r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		vx[i] = V_TH * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);

		r1 = (rand() + 1) / ((double)RAND_MAX + 1);
		r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		vy[i] = V_TH * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
		
		color[i] = 0x0000ff;
	}
}

// traveling wave, displaced charges. this should oscillate as a single mode, but
// this doesn't happen in standard PIC. see section 4 of:
// Huang, et al, 2016. Finite grid instability and spectral fidelity of the
// electrostatic Particle-In-Cell algorithm. Computer Physics Communications
// 207, 123–135.
void initWave(double *x, double *v, int *color) {
	OMEGA_P = sqrt((PART_NUM/XMAX) * (PART_CHARGE*PART_CHARGE)
	               / (PART_MASS * EPS_0));
	
	double ampl = 0.1;
	
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;
		x[i] += ampl * XMAX  / (2 * M_PI * WAVE_MODE)
			* cos(2 * M_PI * WAVE_MODE * x[i] / XMAX);

		x[i] = fmod(x[i], XMAX);
		
		v[i] = ampl * OMEGA_P * XMAX / (2 * M_PI * WAVE_MODE)
			* sin(2 * M_PI * WAVE_MODE * x[i] / XMAX);

		color[i] = 0x0000ff;
	}
}
