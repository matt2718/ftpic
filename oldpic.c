#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <omp.h>

#include <qdsp.h>

#include "common.h"

const int MODELOG_MAX = 32;

// time info
double DT = 0.001;
double TMAX = 20;

// plans and buffers for fft
fftw_plan rhoFFT;
fftw_plan phiIFFT;
double *rhoxBuf, *phixBuf;
fftw_complex *rhokBuf, *phikBuf;

FILE *modeLog = NULL;
FILE *paramLog = NULL;

void deposit(double *x, double *rho);
void fields(double *rho, double *e, double *phi, double *potential);
void xPush(double *x, double *v);
void vHalfPush(double *x, double *v, double *e, int forward);

double kineticEnergy(double *v);
double momentum(double *v);
void vDist();

int main(int argc, char **argv) {
	DX = XMAX / NGRID;

	// whether to plot
	int phasePlotOn = 1;
	int phiPlotOn = 1;
	int rhoPlotOn = 1;
	
	// parse arguments
	for (int i = 1; i < argc; i++) {
		// log file for parameters
		if (!strcmp(argv[i], "-p")) {
			if (++i == argc) return 1;
			paramLog = fopen(argv[i], "w");
		}

		// log file for modes
		if (!strcmp(argv[i], "-m")) {
			if (++i == argc) return 1;
			modeLog = fopen(argv[i], "w");
		}

		// time step and limit
		if (!strcmp(argv[i], "-t")) {
			if (++i == argc) return 1;
			sscanf(argv[i], "%lf,%lf", &DT, &TMAX);
			if (DT <= 0) return 1;
		}

		// quiet, do not render plots
		if (!strcmp(argv[i], "-q")) {
			phasePlotOn = 0;
			phiPlotOn = 0;
			rhoPlotOn = 0;
		}
	}

	// dump parameters
	if (paramLog) {
		// basic
		fprintf(paramLog, " particles: %i\n", PART_NUM);
		fprintf(paramLog, "  timestep: %e\n", DT);
		fprintf(paramLog, "    length: %e\n", XMAX);
		fprintf(paramLog, "    v_beam: %e\n", BEAM_SPEED);
		fprintf(paramLog, "      mass: %e\n", PART_MASS);
		fprintf(paramLog, "    charge: %e\n", PART_CHARGE);
		fprintf(paramLog, "     eps_0: %e\n", EPS_0);
		fprintf(paramLog, "\n");

		// Debye length
		double kt = PART_MASS * BEAM_SPEED * BEAM_SPEED;
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
	
	// allocate memory
	double *x = malloc(PART_NUM * sizeof(double));
	double *v = malloc(PART_NUM * sizeof(double));
	int *color = malloc(PART_NUM * sizeof(int));

	double *rho = malloc(NGRID * sizeof(double));
	double *eField = malloc(NGRID * sizeof(double));
	double *phi = malloc(NGRID * sizeof(double));

	// transform buffers
	rhoxBuf = fftw_malloc(NGRID * sizeof(double));
	phixBuf = fftw_malloc(NGRID * sizeof(double));
	rhokBuf = fftw_malloc(NGRID/2 * sizeof(fftw_complex));
	phikBuf = fftw_malloc(NGRID/2 * sizeof(fftw_complex));
	// plan transforms
	rhoFFT = fftw_plan_dft_r2c_1d(NGRID, rhoxBuf, rhokBuf, FFTW_MEASURE);
	phiIFFT = fftw_plan_dft_c2r_1d(NGRID, phikBuf, phixBuf, FFTW_MEASURE);

	// initialize particles
	//initLandau(x, v, color);
	init2Stream(x, v, color);

	QDSPplot *phasePlot = NULL;
	QDSPplot *phiPlot = NULL;
	QDSPplot *rhoPlot = NULL;

	if (phasePlotOn) {
		phasePlot = qdspInit("Phase plot");
		qdspSetBounds(phasePlot, 0, XMAX, -30, 30);
		qdspSetGridX(phasePlot, 0, 2, 0x888888);
		qdspSetGridY(phasePlot, 0, 10, 0x888888);
		qdspSetPointColor(phasePlot, 0x000000);
		qdspSetBGColor(phasePlot, 0xffffff);
	}
		
	if (phiPlotOn) {
		phiPlot = qdspInit("Phi(x)");
		qdspSetBounds(phiPlot, 0, XMAX, -100, 100);
		qdspSetGridX(phiPlot, 0, 2, 0x888888);
		qdspSetGridY(phiPlot, 0, 20, 0x888888);
		qdspSetConnected(phiPlot, 1);
		qdspSetPointColor(phiPlot, 0x000000);
		qdspSetBGColor(phiPlot, 0xffffff);
	}

	if (rhoPlotOn) {
		rhoPlot = qdspInit("Rho(x)");
		qdspSetBounds(rhoPlot, 0, XMAX, -100, 100);
		qdspSetGridX(rhoPlot, 0, 2, 0x888888);
		qdspSetGridY(rhoPlot, 0, 10, 0x888888);
		qdspSetConnected(rhoPlot, 1);
		qdspSetPointColor(rhoPlot, 0x000000);
		qdspSetBGColor(rhoPlot, 0xffffff);
	}
		
	double *xar = malloc(NGRID * sizeof(double));
	for (int j = 0; j < NGRID; j++) xar[j] = j * DX;
	
	double potential;
	
	deposit(x, rho);
	fields(rho, eField, phi, NULL);

	vHalfPush(x, v, eField, 0); // push backwards

	int open = 1;
	
	printf("time,potential,kinetic,total,momentum\n");

	for (int n = 0; open && n * DT < TMAX; n++) {
		if (modeLog) fprintf(modeLog, "%f", n * DT);
		
		deposit(x, rho);
		fields(rho, eField, phi, &potential);

		vHalfPush(x, v, eField, 1);

		if (phasePlotOn)
			open = qdspUpdateIfReady(phasePlot, x, v, color, PART_NUM);

		// logging
		if (n % 10 == 0) {
			double kinetic = kineticEnergy(v);
			printf("%f,%f,%f,%f,%f\n",
			       n * DT,
			       potential,
			       kinetic,
			       potential + kinetic,
			       momentum(v));
		}

		if (phiPlotOn)
			phiPlotOn = qdspUpdateIfReady(phiPlot, xar, phi, NULL, NGRID);

		if (rhoPlotOn)
			rhoPlotOn = qdspUpdateIfReady(rhoPlot, xar, rho, NULL, NGRID);
		
		vHalfPush(x, v, eField, 1);
		xPush(x, v);
	}

	// cleanup
	if(modeLog) fclose(modeLog);
	
	free(x);
	free(v);
	free(color);

	free(rho);
	free(eField);
	free(phi);

	fftw_free(rhoxBuf);
	fftw_free(rhokBuf);
	fftw_free(phixBuf);
	fftw_free(phikBuf);

	fftw_destroy_plan(rhoFFT);
	fftw_destroy_plan(phiIFFT);

	if (phasePlot) qdspDelete(phasePlot);
	if (phiPlot) qdspDelete(phiPlot);
	if (rhoPlot) qdspDelete(rhoPlot);

	return 0;
}

void deposit(double *x, double *rho) {
	// neutralizing bg
	for (int j = 0; j < NGRID; j++)
		rho[j] = -PART_NUM * PART_CHARGE / XMAX;

#pragma omp parallel
	{
		double *myRho = calloc(NGRID, sizeof(double));

#pragma omp for
		for (int i = 0; i < PART_NUM; i++) {
			int j = x[i] / DX;
			double xg = j * DX;
			myRho[j] += PART_CHARGE * (xg + DX - x[i]) / (DX * DX);
			myRho[(j+1) % NGRID] += PART_CHARGE * (x[i] - xg) / (DX * DX);
		}
				
#pragma omp critical
		{
			for (int j = 0; j < NGRID; j++)
				rho[j] += myRho[j];
		}
		free(myRho);
	}
}

// determines E and phi from rho
void fields(double *rho, double *e, double *phi, double *potential) {
	// rho(x) -> rho(k)
	for (int j = 0; j < NGRID; j++) {
		// nomalization
		rhoxBuf[j] = rho[j] / NGRID;
	}
	fftw_execute(rhoFFT);

	// rho(k) -> phi(k)
	phikBuf[0][0] = 0;
	phikBuf[0][1] = 0;
	for (int j = 1; j < NGRID/2; j++) {
		double k = 2 * M_PI * j / XMAX;
		double ksqi = k * pow(sin(k*DX/2) / (k*DX/2), 2);
		phikBuf[j][0] = rhokBuf[j][0] / (k * k * EPS_0);
		phikBuf[j][1] = rhokBuf[j][1] / (k * k * EPS_0);
	}
	
	// potential energy calculation
	if (potential != NULL) {
		double pot = 0;
		for (int j = 0; j < NGRID/2; j++) {
			pot += phikBuf[j][0] * rhokBuf[j][0] + phikBuf[j][1] * rhokBuf[j][1];
		}
		*potential = pot * XMAX;

		if (modeLog) {
			for (int j = 1; j <= MODELOG_MAX; j++) {
				double etmp = phikBuf[j][0] * rhokBuf[j][0]
					+ phikBuf[j][1] * rhokBuf[j][1];
				fprintf(modeLog, ",%e", etmp);
			}
			fprintf(modeLog, "\n");
		}
	}
	
	// phi(k) -> phi(x)
	fftw_execute(phiIFFT);
	memcpy(phi, phixBuf, NGRID * sizeof(double));

	// determine E(x) via finite difference
	e[0] = (phi[NGRID-1] - phi[1]) / (2 * DX);
	e[NGRID-1] = (phi[NGRID-2] - phi[0]) / (2 * DX);
	for (int j = 1; j < NGRID - 1; j++)
		e[j] = (phi[j-1] - phi[j+1]) / (2 * DX);
}

// moves particles given velocities
void xPush(double *x, double *v) {
#pragma omp parallel for
	for (int i = 0; i < PART_NUM; i++) {
		x[i] += DT * v[i];

		// periodicity
		// (not strictly correct, but if a particle is moving several grid
		// lengths in 1 timestep, something has gone horribly wrong)
		if (x[i] < 0) x[i] += XMAX;
		if (x[i] >= XMAX) x[i] -= XMAX;
	}
}

// pushes particles, electric field calculated via linear interpoation
void vHalfPush(double *x, double *v, double *e, int forward) {
#pragma omp parallel for
	for (int i = 0; i < PART_NUM; i++) {
		// calculate index and placement between grid points
		int jint = (int)(x[i] * NGRID / XMAX);
		double jfrac = (x[i] * NGRID / XMAX) - jint;

		// interpolate e(x_i)
		double e1 = e[jint];
		double e2 = e[(jint + 1) % NGRID];
		double ePart = e1 + (e2 - e1) * jfrac;

		// push
		if (forward)
			v[i] += DT/2 * (PART_CHARGE / PART_MASS) * ePart;
		else
			v[i] -= DT/2 * (PART_CHARGE / PART_MASS) * ePart;
	}
}

double kineticEnergy(double *v) {
	double kinetic = 0;
	
#pragma omp parallel for reduction(+:kinetic)
	for (int i = 0; i < PART_NUM; i++) {
		kinetic += v[i] * v[i] * PART_MASS / 2;
	}
	
	return kinetic;
}

double momentum(double *v) {
	double p = 0;

#pragma omp parallel for reduction(+:p)
	for (int i = 0; i < PART_NUM; i++) {
		p += PART_MASS * v[i];
	}

	return p;
}
