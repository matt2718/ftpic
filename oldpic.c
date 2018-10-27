#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <fftw3.h>

#include <qdsp.h>

#include "common.h"

// deposition buffer
double *tmpRho;

// plans and buffers for fft
fftw_plan rhoFFT;
fftw_plan phiIFFT;
double *rhoxBuf, *phixBuf;
fftw_complex *rhokBuf, *phikBuf;

//FILE *modeLog = NULL;
FILE *paramLog = NULL;

void deposit(double *x, double *rho);
void fields(double *rho, double *e, double *phi, double *potential);
void xPush(double *x, double *v);
void interpField(double *x, double *e, double *ePart);
void vHalfPush(double *v, double *ePart, int forward);

double kineticEnergy(double *v);
double momentum(double *v);
void vDist();

int main(int argc, char **argv) {
	double *x;
	double *v;
	int *color;

	QDSPplot *phasePlot = NULL;
	QDSPplot *phiPlot = NULL;
	QDSPplot *rhoPlot = NULL;

	// parse command line arguments, initialize simulation, set up logging
	// and plotting
	int ret = commonInit(argc, argv,
	                     &x, &v, &color,
	                     &phasePlot, &phiPlot, &rhoPlot);
	if (ret) return ret;
	
	double *rho = malloc(NGRID * sizeof(double));
	double *eField = malloc(NGRID * sizeof(double));
	double *phi = malloc(NGRID * sizeof(double));

	// e field at particle positions
	double *ePart = malloc(PART_NUM * sizeof(double));
		
	// we want to deposit in parallel
	tmpRho = malloc(2 * omp_get_max_threads() * NGRID * sizeof(double));

	// transform buffers
	rhoxBuf = fftw_malloc(NGRID * sizeof(double));
	phixBuf = fftw_malloc(NGRID * sizeof(double));
	rhokBuf = fftw_malloc(NGRID * sizeof(fftw_complex));
	phikBuf = fftw_malloc(NGRID * sizeof(fftw_complex));
	// plan transforms
	rhoFFT = fftw_plan_dft_r2c_1d(NGRID, rhoxBuf, rhokBuf, FFTW_MEASURE);
	phiIFFT = fftw_plan_dft_c2r_1d(NGRID, phikBuf, phixBuf, FFTW_MEASURE);

	// axis for plots
	double *xar = malloc(NGRID * sizeof(double));
	for (int j = 0; j < NGRID; j++) xar[j] = j * DX;

	double potential;

	deposit(x, rho);
	fields(rho, eField, phi, NULL);

	interpField(x, eField, ePart);
	vHalfPush(v, eField, 0); // push backwards

	int open = 1;

	printf("time,potential,kinetic,total,momentum\n");

	// start time logging
	struct timespec time1, time2;
	clock_gettime(CLOCK_MONOTONIC, &time1);

	int n;
	for (n = 0; open && n * DT < TMAX; n++) {
		if (modeLog) fprintf(modeLog, "%f", n * DT);

		deposit(x, rho);
		fields(rho, eField, phi, &potential);

		interpField(x, eField, ePart);
		vHalfPush(v, eField, 1);

		if (phasePlot)
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

		if (phiPlot) {
			int on = qdspUpdateIfReady(phiPlot, xar, phi, NULL, NGRID);
			if (!on) phiPlot = NULL;
		}

		if (rhoPlot) {
			int on = qdspUpdateIfReady(rhoPlot, xar, rho, NULL, NGRID);
			if (!on) rhoPlot = NULL;
		}

		//interpField(x, eField, ePart);
		vHalfPush(v, eField, 1);
		xPush(x, v);
	}
	clock_gettime(CLOCK_MONOTONIC, &time2);

	if (printTime) {
		double elapsed = (time2.tv_sec - time1.tv_sec) * 1000.0;
		elapsed += (time2.tv_nsec - time1.tv_nsec) / 1000000.0;
		fprintf(stderr, "%f ms per step\n", elapsed / n);
	}

	// cleanup
	if(modeLog) fclose(modeLog);

	free(xar);
	
	free(x);
	free(v);
	free(color);

	free(rho);
	free(eField);
	free(phi);
	
	free(ePart);
	free(tmpRho);

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
		double *myRho = tmpRho + omp_get_thread_num() * NGRID;
		for (int j = 0; j < NGRID; j++)
			myRho[j] = 0;

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

// calculated E field at each particle
void interpField(double *x, double *e, double *ePart) {
#pragma omp parallel for
	for (int i = 0; i < PART_NUM; i++) {
		// calculate index and placement between grid points
		int jint = (int)(x[i] * NGRID / XMAX);
		double jfrac = (x[i] * NGRID / XMAX) - jint;

		// interpolate e(x_i)
		double e1 = e[jint];
		double e2 = e[(jint + 1) % NGRID];
	      ePart[i] = e1 + (e2 - e1) * jfrac;
	}
}
	
// pushes particles
void vHalfPush(double *v, double *ePart, int forward) {
	double factor = DT/2 * (PART_CHARGE / PART_MASS);
	if (!forward) factor = -factor;

#pragma omp parallel for
	for (int m = 0; m < PART_NUM; m++) {
		v[m] += factor * ePart[m];
	}
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
