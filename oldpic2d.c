#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <fftw3.h>

#include <qdsp.h>

#include "common2d.h"

double *rhoxBuf, *phixBuf, *exBuf, *eyBuf;
fftw_complex *rhokBuf, *phikBuf, *exkBuf, *eykBuf;

fftw_plan rhoFFT, exIFFT, eyIFFT, phiIFFT;

void deposit(double *x, double *y, double *rho);
void fields(double *rho, double *phi, double *potential);
void interpField(double *x, double *y, double *eBuf, double *ePart);
void vHalfPush(double *v, double *ePart, int forward);
void xPush(double *x, double *v);

double kineticEnergy(double *v);
double momentum(double *v);

int main(int argc, char **argv) {
	double *x;
	double *y;
	double *vx;
	double *vy;

	int *color;

	QDSPplot *xyPlot = NULL;

	// parse command line arguments, initialize simulation, set up logging
	// and plotting
	int ret = commonInit(argc, argv,
	                     &x, &y, &vx, &vy, &color,
	                     &xyPlot);
	if (ret) return ret;
	
	double *rho = malloc(NGRIDX * NGRIDY * sizeof(double));
	double *phi = malloc(NGRIDX * NGRIDY * sizeof(double));

	// transform buffers
	rhoxBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(double));
	rhokBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	phixBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(double));
	phikBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	exBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(double));
	exkBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));
	eyBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(double));
	eykBuf = fftw_malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	// field at each particle
	double *exPart = malloc(PART_NUM * sizeof(double));
	double *eyPart = malloc(PART_NUM * sizeof(double));
	
	// plan transforms
	rhoFFT = fftw_plan_dft_r2c_2d(NGRIDX, NGRIDY, rhoxBuf, rhokBuf, FFTW_MEASURE);
	exIFFT = fftw_plan_dft_c2r_2d(NGRIDX, NGRIDY, exkBuf, exBuf, FFTW_MEASURE);
	eyIFFT = fftw_plan_dft_c2r_2d(NGRIDX, NGRIDY, eykBuf, eyBuf, FFTW_MEASURE);
	phiIFFT = fftw_plan_dft_c2r_2d(NGRIDX, NGRIDY, phikBuf, phixBuf, FFTW_MEASURE);

	double potential;

	deposit(x, y, rho);
	fields(rho, phi, NULL);
	//vHalfPush(x, y, vx, vy, 0);
	interpField(x, y, exBuf, exPart);
	interpField(x, y, eyBuf, eyPart);
	// push backwards
	vHalfPush(vx, exPart, 0);
	vHalfPush(vy, eyPart, 0);
	
	int open = 1;

	printf("time,potential,kinetic,total,momentum\n");

	// start time logging
	struct timespec time1, time2;
	clock_gettime(CLOCK_MONOTONIC, &time1);

	int n;
	for (n = 0; open && n * DT < TMAX; n++) {
		if (modeLog) fprintf(modeLog, "%f", n * DT);

		deposit(x, y, rho);
		fields(rho, phi, &potential);

		//vHalfPush(x, y, vx, vy, 1);
		interpField(x, y, exBuf, exPart);
		interpField(x, y, eyBuf, eyPart);
		vHalfPush(vx, exPart, 1);
		vHalfPush(vy, eyPart, 1);

		if (xyPlot)
			open = qdspUpdateIfReady(xyPlot, x, y, color, PART_NUM);

		// logging
		if (n % 10 == 0) {
			double kinetic = kineticEnergy(vx) + kineticEnergy(vy);
			printf("%f,%f,%f,%f,%f\n",
			       n * DT,
			       potential,
			       kinetic,
			       potential + kinetic,
			       momentum(vx));
		}

		vHalfPush(vx, exPart, 1);
		vHalfPush(vy, eyPart, 1);
		xPush(x, vx);
		xPush(y, vy);
	}
	clock_gettime(CLOCK_MONOTONIC, &time2);

	if (printTime) {
		double elapsed = (time2.tv_sec - time1.tv_sec) * 1000.0;
		elapsed += (time2.tv_nsec - time1.tv_nsec) / 1000000.0;
		fprintf(stderr, "%f ms per step\n", elapsed / n);
	}

	////////////////////////////////
	// CLEANUP

	if(modeLog) fclose(modeLog);

	// allocated by commonInit
	free(x);
	free(y);
	free(vx);
	free(vy);
	free(color);

	free(rho);
	free(phi);

	free(exPart);
	free(eyPart);

	fftw_free(rhoxBuf);
	fftw_free(rhokBuf);
	fftw_free(phixBuf);
	fftw_free(phikBuf);

	fftw_free(exBuf);
	fftw_free(exkBuf);
	fftw_free(eyBuf);
	fftw_free(eykBuf);
	
	fftw_destroy_plan(rhoFFT);
	fftw_destroy_plan(exIFFT);
	fftw_destroy_plan(eyIFFT);
	fftw_destroy_plan(phiIFFT);

	if (xyPlot) qdspDelete(xyPlot);

	return 0;
}

void deposit(double *x, double *y, double *rho) {
	for (int j = 0; j < NGRIDX * NGRIDY; j++) {
		rho[j] = -PART_NUM * PART_CHARGE / (XMAX * YMAX);
	}

	for (int i = 0; i < PART_NUM; i++) {
		int jx1 = (int)(x[i] / DX);
		int jy1 = (int)(y[i] / DY);
		int jx2 = (jx1 + 1) % NGRIDX;
		int jy2 = (jy1 + 1) % NGRIDY;
		
		double xfrac = x[i] / DX - jx1;
		double yfrac = y[i] / DY - jy1;

		rho[jx1 + NGRIDX * jy1] += PART_CHARGE * (1 - xfrac) * (1 - yfrac);
		rho[jx2 + NGRIDX * jy1] += PART_CHARGE * xfrac       * (1 - yfrac);
		rho[jx1 + NGRIDX * jy2] += PART_CHARGE * (1 - xfrac) * yfrac;
		rho[jx2 + NGRIDX * jy2] += PART_CHARGE * xfrac       * yfrac;
	}
}

// determines E and phi from rho
void fields(double *rho, double *phi, double *potential) {
	// rho(x) -> rho(k)
	for (int j = 0; j < NGRIDX * NGRIDY; j++) {
		// nomalization
		rhoxBuf[j] = rho[j] / (XMAX * YMAX);
	}
	fftw_execute(rhoFFT);

	// rho(k) -> phi(k)
	phikBuf[0][0] = 0;
	phikBuf[0][1] = 0;
	exkBuf[0][0] = 0;
	exkBuf[0][1] = 0;
	eykBuf[0][0] = 0;
	eykBuf[0][1] = 0;
	for (int j = 1; j < NGRIDX * NGRIDY; j++) {
		double kx = 2 * M_PI * (j % (NGRIDX/2+1)) / XMAX;
		double ky = 2 * M_PI * (int)(j / (NGRIDX/2+1)) / YMAX;

		if (ky > M_PI * NGRIDY / YMAX)
			ky -= 2 * M_PI * NGRIDY / YMAX;

		double ksqr = kx * kx + ky * ky;
		phikBuf[j][0] = rhokBuf[j][0] / (ksqr * EPS_0);
		phikBuf[j][1] = rhokBuf[j][1] / (ksqr * EPS_0);
		
		double k = sqrt(ksqr);

		exkBuf[j][0] =  kx * phikBuf[j][1];
		exkBuf[j][1] = -kx * phikBuf[j][0];

		eykBuf[j][0] =  ky * phikBuf[j][1];
		eykBuf[j][1] = -ky * phikBuf[j][0];
	}
	
	// potential energy calculation
	if (potential != NULL) {
		double pot = 0;
		for (int j = 0; j < NGRIDX * NGRIDY / 2; j++) {
			pot += phikBuf[j][0] * rhokBuf[j][0] + phikBuf[j][1] * rhokBuf[j][1];
		}
		*potential = pot * XMAX * YMAX;
	}

	fftw_execute(exIFFT);
	fftw_execute(eyIFFT);

	/*
	// phi(k) -> phi(x)
	fftw_execute(phiIFFT);
	memcpy(phi, phixBuf, NGRIDX * NGRIDY * sizeof(double));

	// determine E(x) via finite difference
	e[0] = (phi[NGRID-1] - phi[1]) / (2 * DX);
	e[NGRID-1] = (phi[NGRID-2] - phi[0]) / (2 * DX);
	for (int j = 1; j < NGRID - 1; j++)
		e[j] = (phi[j-1] - phi[j+1]) / (2 * DX);
	*/
}

void interpField(double *x, double *y, double *eBuf, double *ePart) {
	for (int i = 0; i < PART_NUM; i++) {
		int jx1 = (int)(x[i] / DX);
		int jy1 = (int)(y[i] / DY);
		int jx2 = (jx1 + 1) % NGRIDX;
		int jy2 = (jy1 + 1) % NGRIDY;
		
		double xfrac = x[i] / DX - jx1;
		double yfrac = y[i] / DY - jy1;

		ePart[i] = 0; // field at particle
		
		ePart[i] += (1 - xfrac) * (1 - yfrac) * eBuf[jx1 + NGRIDX * jy1];
		ePart[i] += xfrac       * (1 - yfrac) * eBuf[jx2 + NGRIDX * jy1];
		ePart[i] += (1 - xfrac) * yfrac       * eBuf[jx1 + NGRIDX * jy2];
		ePart[i] += xfrac       * yfrac       * eBuf[jx2 + NGRIDX * jy2];
	}
}

void vHalfPush(double *v, double *ePart, int forward) {
	double factor = DT/2 * (PART_CHARGE / PART_MASS);
	if (!forward) factor = -factor;
	for (int m = 0; m < PART_NUM; m++)
		v[m] += factor * ePart[m];
}
/*
// pushes particles, electric field calculated via linear interpoation
void vHalfPush(double *x, double *y, double *vx, double *vy, int forward) {
	//#pragma omp parallel for
	for (int i = 0; i < PART_NUM; i++) {
		int jx1 = (int)(x[i] / DX);
		int jy1 = (int)(y[i] / DY);
		int jx2 = (jx1 + 1) % NGRIDX;
		int jy2 = (jy1 + 1) % NGRIDY;
		
		double xfrac = x[i] / DX - jx1;
		double yfrac = y[i] / DY - jy1;

		double exPart = 0, eyPart = 0; // x and y field at particle
		
		exPart += (1 - xfrac) * (1 - yfrac) * exBuf[jx1 + NGRIDX * jy1];
		exPart += xfrac       * (1 - yfrac) * exBuf[jx2 + NGRIDX * jy1];
		exPart += (1 - xfrac) * yfrac       * exBuf[jx1 + NGRIDX * jy2];
		exPart += xfrac       * yfrac       * exBuf[jx2 + NGRIDX * jy2];

		eyPart += (1 - xfrac) * (1 - yfrac) * eyBuf[jx1 + NGRIDX * jy1];
		eyPart += xfrac       * (1 - yfrac) * eyBuf[jx2 + NGRIDX * jy1];
		eyPart += (1 - xfrac) * yfrac       * eyBuf[jx1 + NGRIDX * jy2];
		eyPart += xfrac       * yfrac       * eyBuf[jx2 + NGRIDX * jy2];

		// push
		if (forward) {
			vx[i] += DT/2 * (PART_CHARGE / PART_MASS) * exPart;
			vy[i] += DT/2 * (PART_CHARGE / PART_MASS) * eyPart;
		} else {
			vx[i] -= DT/2 * (PART_CHARGE / PART_MASS) * exPart;
			vy[i] -= DT/2 * (PART_CHARGE / PART_MASS) * eyPart;
		}
	}
}
*/
// moves particles given velocities
void xPush(double *x, double *v) {
	//#pragma omp parallel for
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

	//#pragma omp parallel for reduction(+:kinetic)
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
