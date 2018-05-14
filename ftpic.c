#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>

#include <qdsp.h>

#include "common.h"

// fft plans and buffers
fftw_plan phiIFFT;
fftw_complex *phikBuf;
double *phixBuf;

// USFFT buffers
fftw_complex *zcBuf;
double *xpBuf;
double *fpBuf;

fftw_complex *ekBuf;
fftw_complex *epBuf;

extern void uf1t_(int*, double*, int*, double*, double*, int*, int*);
extern void uf1a_(int*, double*, int*, double*, double*, int*, int*);

double shape(double x);

void deposit(double *x, fftw_complex *rhok, fftw_complex *sk);
void fields(fftw_complex *rhok, fftw_complex *sk, double *phi, double *potential);
void xPush(double *x, double *v);
void vHalfPush(double *x, double *v, int forward);

double kineticEnergy(double *v);
double momentum(double *v);

int main(int argc, char **argv) {
	DX = XMAX / NGRID;

	// allocate memory
	double *x = malloc(PART_NUM * sizeof(double));
	double *v = malloc(PART_NUM * sizeof(double));
	int *color = malloc(PART_NUM * sizeof(int));
	
	fftw_complex *rhok = fftw_malloc(NGRID * sizeof(fftw_complex));
	double *rhox = fftw_malloc(NGRID * sizeof(double));
	double *phix = fftw_malloc(NGRID * sizeof(double));

	double *sx = fftw_malloc(NGRID * sizeof(double));
	fftw_complex *sk = fftw_malloc(NGRID * sizeof(fftw_complex));

	// transform buffers
	phikBuf = fftw_malloc(NGRID/2 * sizeof(fftw_complex));
	phixBuf = fftw_malloc(NGRID * sizeof(double));

	// USFFT buffers
	zcBuf = malloc(NGRID * sizeof(fftw_complex));
	xpBuf = malloc(PART_NUM * sizeof(double));
	fpBuf = malloc(2*PART_NUM * sizeof(double));

	// reverse USFFT
	ekBuf = malloc(2 * NGRID * sizeof(fftw_complex));
	epBuf = malloc(PART_NUM * sizeof(fftw_complex));
	
	// plan transforms
	fftw_plan rhoIFFT = fftw_plan_dft_c2r_1d(NGRID, rhok, rhox, FFTW_MEASURE);
	phiIFFT = fftw_plan_dft_c2r_1d(NGRID, phikBuf, phixBuf, FFTW_MEASURE);

	// determine s(k)
	fftw_plan sFFT = fftw_plan_dft_r2c_1d(NGRID, sx, sk, FFTW_ESTIMATE);

	double sxsum = 0;
	for (int j = 0; j < NGRID; j++) {
		double xcur = j * XMAX / NGRID;
		sx[j] = shape(xcur) + shape(XMAX - xcur);
		sxsum += sx[j];
	}

	// USFFT coefficients are 1 + 0i
	// change this so we can convolute in USFFT
	for (int m = 0; m < PART_NUM; m++) {
		fpBuf[2*m] = 1;
		fpBuf[2*m + 1] = 0;
	}

	sxsum *= DX;
	for (int j = 0; j < NGRID; j++) sx[j] /= sxsum;

	fftw_execute(sFFT);
	fftw_destroy_plan(sFFT);

	QDSPplot *phasePlot = NULL;
	QDSPplot *phiPlot = NULL;
	QDSPplot *rhoPlot = NULL;

	// parse command line arguments, initialize simulation, and set up logging
	int ret = commonInit(argc, argv, x, v, color, 
	                     &phasePlot, &phiPlot, &rhoPlot);
	
	if (ret) return ret;
		
	double *xar = malloc(NGRID * sizeof(double));
	for (int j = 0; j < NGRID; j++) xar[j] = j * DX;

	double potential;
	
	deposit(x, rhok, sk);
	fields(rhok, sk, phix, &potential);
	
	vHalfPush(x, v, 0);

	int open = 1;

	printf("time,potential,kinetic,total,momentum\n");

	// check momentum conservation (not currently used)
	double minp = 1/0.0;
	double maxp = 0.0;
	
	for (int n = 0; open && n * DT < TMAX; n++) {
		if (modeLog) fprintf(modeLog, "%f", n * DT);

		deposit(x, rhok, sk);
		fields(rhok, sk, phix, &potential);
		vHalfPush(x, v, 1);

		if (phasePlot)
			open = qdspUpdateIfReady(phasePlot, x, v, color, PART_NUM);
		
		// logging
		if (n % 10 == 0) {
			double kinetic = kineticEnergy(v);
			double curp = momentum(v);
			printf("%f,%f,%f,%f,%f\n",
			       n * DT,
			       potential,
			       kinetic,
			       potential + kinetic,
			       curp);
			if (curp < minp) minp = curp;
			if (curp > maxp) maxp = curp;
		}

		if (phiPlot) {
			int on = qdspUpdateIfReady(phiPlot, xar, phix, NULL, NGRID);
			if (!on) phiPlot = NULL;
		}
		
		if (rhoPlot) {
			fftw_execute(rhoIFFT);
			int on = qdspUpdateIfReady(rhoPlot, xar, rhox, NULL, NGRID);
			if (!on) rhoPlot = NULL;
		}

		vHalfPush(x, v, 1);
		xPush(x, v);
	}

	// cleanup
	if(modeLog) fclose(modeLog);

	free(x);
	free(v);
	free(color);

	free(zcBuf);
	free(xpBuf);
	free(fpBuf);
	
	fftw_free(rhok);
	fftw_free(rhox);
	fftw_free(phix);

	fftw_free(sx);
	fftw_free(sk);

	fftw_destroy_plan(phiIFFT);
	fftw_destroy_plan(rhoIFFT);

	if (phasePlot) qdspDelete(phasePlot);
	if (phiPlot) qdspDelete(phiPlot);
	if (rhoPlot) qdspDelete(rhoPlot);

	return 0;
}

// particle shape function, centered at 0, gaussian in this case
double shape(double x) {
	//const double sigma = 1;
	//return exp(-x*x / (2 * sigma * sigma)) / sqrt(2 * M_PI * sigma * sigma);
	return 1.0 * (x == 0);
	//return fmax(1 - fabs(x/DX), 0);
}

// determines rho(k) from list of particle positions
void deposit(double *x, fftw_complex *rhok, fftw_complex *sk) {

	int nc = NGRID;
	int np = PART_NUM;
	int isign = -1;
	int order = 5;

	#pragma omp parallel for
	for (int m = 0; m < PART_NUM; m++) {
		xpBuf[m] = x[m] / XMAX;
		//if (xpBuf[m] >= 0.5) xpBuf[m] -= 1.0;
	}

	uf1t_(&nc, (double*)zcBuf, &np, xpBuf, fpBuf, &isign, &order);

	#pragma omp parallel for
	for (int j = 0; j < NGRID/2; j++) {
		double real = PART_CHARGE * zcBuf[NGRID/2 + j][0] / NGRID;
		double imag = PART_CHARGE * zcBuf[NGRID/2 + j][1] / NGRID;
		//printf("%d\t%f,%f\t%f,%f\n", j, real, imag);
		rhok[j][0] = real * sk[j][0] - imag * sk[j][1];
		rhok[j][1] = real * sk[j][1] + imag * sk[j][0];
	}

	// background
	rhok[0][0] = 0;
}

// determine phi and e from rho(k) 
void fields(fftw_complex *rhok, fftw_complex *sk, double *phi, double *potential) {
	
	// rho(k) -> phi(k)
	phikBuf[0][0] = 0;
	phikBuf[0][1] = 0;
	for (int j = 1; j < NGRID/2; j++) {
		double k = 2 * M_PI * j / XMAX;
		
		double phikRe = rhok[j][0] / (k * k * EPS_0);
		double phikIm = rhok[j][1] / (k * k * EPS_0);

		phikBuf[j][0] = (sk[j][0] * phikRe - sk[j][1] * phikIm) * DX;
		phikBuf[j][1] = (sk[j][0] * phikIm + sk[j][1] * phikRe) * DX;

		ekBuf[NGRID/2 + j][0] = k * phikBuf[j][1];
		ekBuf[NGRID/2 + j][1] = -k * phikBuf[j][0];
	}

	// find PE
	if (potential != NULL) {
		double pot = 0;
		for (int j = 1; j < NGRID / 2; j++) {
			double k = 2 * M_PI * j / XMAX;
			// we don't use phi[j] directly because it's already been convoluted
			double tmp = rhok[j][0] * rhok[j][0] + rhok[j][1] * rhok[j][1];
			pot += tmp / (k * k * EPS_0);
		}
		*potential = pot * XMAX;

		if (modeLog) {
			for (int j = 1; j <= MODELOG_MAX; j++) {
				double etmp = phikBuf[j][0] * rhok[j][0]
					+ phikBuf[j][1] * rhok[j][1];
				fprintf(modeLog, ",%e", etmp);
			}
			fprintf(modeLog, "\n");
		}
	}

	// phi(k) -> phi(x)
	fftw_execute(phiIFFT);
	memcpy(phi, phixBuf, NGRID * sizeof(double));
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

// interpolates E field on particles and accelerates them 1/2 timestep
void vHalfPush(double *x, double *v, int forward) {
	// calculate forces from E(k)
	int nc = NGRID;
	int np = PART_NUM;
	int isign = 1;
	int order = 5;

	// first element of ekBuf must be 0
	ekBuf[0][0] = 0;
	ekBuf[0][1] = 0;
	for (int j = 0; j < NGRID/2; j++) {
		// array must be Hermitian
		ekBuf[NGRID/2 - j][0] = ekBuf[j + NGRID/2][0];
		ekBuf[NGRID/2 - j][1] = -ekBuf[j + NGRID/2][1];
	}
	
#pragma omp parallel for
	for (int m = 0; m < PART_NUM; m++) {
		xpBuf[m] = x[m] / XMAX;
	}
	
	uf1a_(&nc, (double*)ekBuf, &np, xpBuf, (double*)epBuf, &isign, &order);
	for (int m = 0; m < PART_NUM; m++) {
		// interpolated e(x_m)
		double ePart = epBuf[m][0];

		// push
		if (forward)
			v[m] += DT/2 * (PART_CHARGE / PART_MASS) * ePart;
		else
			v[m] -= DT/2 * (PART_CHARGE / PART_MASS) * ePart;
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
