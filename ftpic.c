#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>

#include <qdsp.h>

const double XMAX = 16.0; // system length
const int NGRID = 256; // grid size
double DX;

// particle number and properties
const int PART_NUM = 1000;
const double PART_MASS = 1.0;
const double PART_CHARGE = -0.01;
const double EPS_0 = 1.0;

// time info
const double DT = 0.0005;
const double NMAX = 1000;

// fft plans and buffers
fftw_plan phiIFFT;
fftw_complex *phikBuf;
double *phixBuf;

// USFFT buffers
fftw_complex *zcBuf;
double *xpBuf;
double *fpBuf;

extern void uf1t_(int*, double*, int*, double*, double*, int*, int*);

double shape(double x);
void init(double *x, double *v, int *color);

void deposit(double *x, fftw_complex *rhok, fftw_complex *sk);
void fields(fftw_complex *rhok, fftw_complex *sk, double *e, double *phi,
            double *potential);
void xPush(double *x, double *v);
void vHalfPush(double *x, double *v, double *e, int forward);

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
	double *ex = fftw_malloc(NGRID * sizeof(double));

	double *sx = fftw_malloc(NGRID * sizeof(double));
	fftw_complex *sk = fftw_malloc(NGRID * sizeof(fftw_complex));

	// transform buffers
	phikBuf = fftw_malloc(NGRID * sizeof(fftw_complex));
	phixBuf = fftw_malloc(NGRID * sizeof(double));

	// USFFT buffers
	zcBuf = malloc(2 * NGRID * sizeof(fftw_complex));
	xpBuf = malloc(PART_NUM * sizeof(double));
	fpBuf = malloc(2*PART_NUM * sizeof(double));
	
	// plan transforms
	fftw_plan rhoIFFT = fftw_plan_dft_c2r_1d(NGRID, rhok, rhox, FFTW_MEASURE);
	phiIFFT = fftw_plan_dft_c2r_1d(NGRID, phikBuf, phixBuf, FFTW_MEASURE);

	// determine s(k)
	double sxsum = 0;
	for (int j = 0; j < NGRID; j++) {
		double xcur = j * XMAX / NGRID;
		sx[j] += shape(xcur) + shape(XMAX - xcur);
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

	fftw_plan sFFT = fftw_plan_dft_r2c_1d(NGRID, sx, sk, FFTW_ESTIMATE);
	fftw_execute(sFFT);
	fftw_destroy_plan(sFFT);
	
	// initialize particles
	init(x, v, color);

	QDSPplot *phasePlot = qdspInit("Phase plot");
	qdspSetBounds(phasePlot, 0, XMAX, -30, 30);
	qdspSetGridX(phasePlot, 0, 2, 0x888888);
	qdspSetGridY(phasePlot, 0, 10, 0x888888);
	qdspSetPointColor(phasePlot, 0x000000);
	qdspSetBGColor(phasePlot, 0xffffff);

	QDSPplot *phiPlot = qdspInit("Phi(x)");
	qdspSetBounds(phiPlot, 0, XMAX, -10000, 10000);
	qdspSetGridX(phiPlot, 0, 2, 0x888888);
	qdspSetGridY(phiPlot, 0, 2000, 0x888888);
	qdspSetConnected(phiPlot, 1);
	qdspSetPointColor(phiPlot, 0x000000);
	qdspSetBGColor(phiPlot, 0xffffff);

	QDSPplot *rhoPlot = qdspInit("Rho(x)");
	qdspSetBounds(rhoPlot, 0, XMAX, -1000, 1000);
	qdspSetGridX(rhoPlot, 0, 2, 0x888888);
	qdspSetGridY(rhoPlot, 0, 200, 0x888888);
	qdspSetConnected(rhoPlot, 1);
	qdspSetPointColor(rhoPlot, 0x000000);
	qdspSetBGColor(rhoPlot, 0xffffff);
	
	double *xar = malloc(NGRID * sizeof(double));
	for (int j = 0; j < NGRID; j++) xar[j] = j * DX;

	double potential;
	
	deposit(x, rhok, sk);
	fields(rhok, sk, ex, phix, &potential);
	vHalfPush(x, v, ex, 0);

	int open = 1;
	int phiOn = 1;
	int rhoOn = 1;

	printf("time,potential,kinetic,total,momentum\n");
	
	for (int n = 0; open; n++) {
		deposit(x, rhok, sk);
		fields(rhok, sk, ex, phix, &potential);

		vHalfPush(x, v, ex, 1);
		open = qdspUpdateIfReady(phasePlot, x, v, color, PART_NUM);
		// logging
		if (n % 50 == 0) {
			double kinetic = kineticEnergy(v);
			printf("%f,%f,%f,%f,%f\n",
			       n * DT,
			       potential,
			       kinetic,
			       potential + kinetic,
			       momentum(v));
		}

		if (phiOn) phiOn = qdspUpdateIfReady(phiPlot, xar, phix, NULL, NGRID);
		if (rhoOn) {
			fftw_execute(rhoIFFT);
			rhoOn = qdspUpdateIfReady(rhoPlot, xar, rhox, NULL, NGRID);
		}

		vHalfPush(x, v, ex, 1);
		xPush(x, v);
	}

	// cleanup
	free(x);
	free(v);
	free(color);

	free(zcBuf);
	free(xpBuf);
	free(fpBuf);
	
	fftw_free(rhok);
	fftw_free(rhox);
	fftw_free(phix);
	fftw_free(ex);

	fftw_free(sx);
	fftw_free(sk);

	fftw_destroy_plan(phiIFFT);
	fftw_destroy_plan(rhoIFFT);

	qdspDelete(phasePlot);
	qdspDelete(phiPlot);
	qdspDelete(rhoPlot);

	return 0;
}

// particle shape function, centered at 0, gaussian in this case
double shape(double x) {
	const double sigma = 0.05;
	return exp(-x*x / (2 * sigma * sigma)) / sqrt(2 * M_PI * sigma * sigma);
	//return 1.0 * (x == 0);
}

void init(double *x, double *v, int *color) {
	double stdev = sqrt(5000 / (5.1e5));
	for (int i = 0; i < PART_NUM; i++) {
		x[i] = i * XMAX / PART_NUM;

		if (i % 2) {
			v[i] = 8.0;
			color[i] = 0xff0000;
		} else {
			v[i] = -8.0;
			color[i] = 0x0000ff;
		}

		// box-mueller
		double r1 = (rand() + 1) / ((double)RAND_MAX + 1); // log(0) breaks stuff
		double r2 = (rand() + 1) / ((double)RAND_MAX + 1);
		v[i] += stdev * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
	}
}

// determines rho(k) from list of particle positions
void deposit(double *x, fftw_complex *rhok, fftw_complex *sk) {

	int nc = 2 * NGRID;
	int np = PART_NUM;
	int isign = -1;
	int order = 5;

	#pragma omp parallel for
	for (int m = 0; m < PART_NUM; m++) {
		xpBuf[m] = x[m] / XMAX;
	}

	uf1t_(&nc, (double*)zcBuf, &np, xpBuf, fpBuf, &isign, &order);
	
	#pragma omp parallel for
	for (int j = 0; j < NGRID; j++) {
		double real = PART_CHARGE * zcBuf[NGRID + j][0];
		double imag = PART_CHARGE * zcBuf[NGRID + j][1];
		rhok[j][0] = real * sk[j][0] - imag * sk[j][1];
		rhok[j][1] = real * sk[j][1] + imag * sk[j][0];
	}

	// background
	rhok[0][0] = 0;
}


void fields(fftw_complex *rhok, fftw_complex *sk, double *e, double *phi,
            double *potential) {
	
	// rho(k) -> phi(k)
	phikBuf[0][0] = 0;
	phikBuf[0][1] = 0;
	for (int j = 1; j < NGRID; j++) {
		double k = 2 * M_PI * j / XMAX;
		
		double phikRe = rhok[j][0] / (k * k * EPS_0);
		double phikIm = rhok[j][1] / (k * k * EPS_0);

		phikBuf[j][0] = sk[j][0] * phikRe - sk[j][1] * phikIm;
		phikBuf[j][1] = sk[j][0] * phikIm + sk[j][1] * phikRe;
	}


	// find PE
	if (potential != NULL) {
		double pot = 0;
		for (int j = 1; j < NGRID; j++) {
			pot += phikBuf[j][0] * rhok[j][0] + phikBuf[j][1] * rhok[j][1];
		}
		*potential = pot / XMAX;
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

// pushes particles, electric field calculated via linear interpoation
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
