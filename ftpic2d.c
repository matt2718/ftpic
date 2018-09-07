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

#include "/home/msm/Documents/fftw-test/hsv2rgb.h"

double absClamp(fftw_complex z, double scale) {
	return fmax(0.0, fmin(1.0, sqrt(z[0]*z[0] + z[1]*z[1]) / scale));
}

double argDeg(fftw_complex z) {
	return 180.0 / M_PI * fmod(atan2(z[1], z[0]) + 2*M_PI, 2*M_PI);
}

int complex2hex(fftw_complex z) {
	hsv colHSV = {argDeg(z), absClamp(z, 1), 1.0};
	rgb colRGB = hsv2rgb(colHSV);
	return (int)(255*colRGB.r) * 0x010000
           + (int)(255*colRGB.g) * 0x000100
	     + (int)(255*colRGB.b);
}

extern void uf2t_(int*, int*, double*, int*, double*, double*, double*, int*, int*);
extern void uf2a_(int*, int*, double*, int*, double*, double*, double*, int*, int*);

void deposit(double *x, double *y);
void fields(double *potential);
void interpField(double *x, double *y, fftw_complex *ek, fftw_complex *ePart);
void vHalfPush(double *v, fftw_complex *ePart, int forward);
void xPush(double *x, double *y, double *vx, double *vy);

double kineticEnergy(double *v);
double momentum(double *v);

fftw_complex *zcBuf;
double *xpBuf;
double *ypBuf;
fftw_complex *fpBuf;

void arraydump_double(double *array, const char* fname) {
	FILE *outfile = fopen(fname, "w");

	for (int r = 0; r < NGRIDY; r++) {
		fprintf(outfile, "%f", array[NGRIDX * r]);
		for (int c = 1; c < NGRIDX; c++) {
			fprintf(outfile, ",%f", array[NGRIDX * r + c]);
		}
		fprintf(outfile, "\n");
	}

	fclose(outfile);
}

void arraydump_complex(fftw_complex *array, const char* fname, int im) {
	FILE *outfile = fopen(fname, "w");

	for (int r = 0; r < NGRIDY; r++) {
		fprintf(outfile, "%f", array[NGRIDX * r][im]);
		for (int c = 1; c < NGRIDX; c++) {
			fprintf(outfile, ",%f", array[NGRIDX * r + c][im]);
		}
		fprintf(outfile, "\n");
	}

	fclose(outfile);
}

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

	
	//fftw_complex *rhok = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));
	//double *phi = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	// transformed quantity buffers
	rhoxBuf = malloc(NGRIDX * NGRIDY * sizeof(double));
	rhokBuf = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	phixBuf = malloc(NGRIDX * NGRIDY * sizeof(double));
	phikBuf = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	exBuf  = malloc(NGRIDX * NGRIDY * sizeof(double));
	eyBuf  = malloc(NGRIDX * NGRIDY * sizeof(double));
	exkBuf = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));
	eykBuf = malloc(NGRIDX * NGRIDY * sizeof(fftw_complex));

	// transform buffers
	zcBuf = malloc(PART_NUM * sizeof(fftw_complex));
	xpBuf = malloc(PART_NUM * sizeof(double));
	ypBuf = malloc(PART_NUM * sizeof(double));
	fpBuf = malloc(PART_NUM * sizeof(fftw_complex));
 
	fftw_complex *exPart = malloc(PART_NUM * sizeof(fftw_complex));
	fftw_complex *eyPart = malloc(PART_NUM * sizeof(fftw_complex));
	
	double potential;

	deposit(x, y);
	fields(NULL);
	interpField(x, y, exkBuf, exPart);
	interpField(x, y, eykBuf, eyPart);
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

		deposit(x, y);
		fields(&potential);
		interpField(x, y, exkBuf, exPart);
		interpField(x, y, eykBuf, eyPart);
		vHalfPush(vx, exPart, 0);
		vHalfPush(vy, eyPart, 0);

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

		vHalfPush(vx, exPart, 0);
		vHalfPush(vy, eyPart, 0);
		xPush(x, y, vx, vy);
	}
	clock_gettime(CLOCK_MONOTONIC, &time2);

	if (printTime) {
		double elapsed = (time2.tv_sec - time1.tv_sec) * 1000.0;
		elapsed += (time2.tv_nsec - time1.tv_nsec) / 1000000.0;
		fprintf(stderr, "%f ms per step\n", elapsed / n);
	}

	// cleanup
	if(modeLog) fclose(modeLog);

	free(x);
	free(y);
	free(vx);
	free(vy);
	free(color);

	//free(rho);
	//free(phi);

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

void deposit(double *x, double *y) {
	int ncx = NGRIDX, ncy = NGRIDY;
	int np = PART_NUM;
	int isign = -1;
	int order = 5;

	for (int m = 0; m < PART_NUM; m++) {
		xpBuf[m] = x[m] / XMAX;
		ypBuf[m] = y[m] / XMAX;
		fpBuf[m][0] = 1;
		fpBuf[m][1] = 0;
	}

	uf2t_(&ncx, &ncy, (double*)zcBuf,
	      &np, xpBuf, ypBuf, (double*)fpBuf,
	      &isign, &order);

	for (int jx = 0; jx < NGRIDX; jx++) {
		for (int jy = 0; jy < NGRIDX; jy++) {
			int j = jx + NGRIDX*jy;
			int j2 = NGRIDY*jx + jy;

			double real = PART_CHARGE * zcBuf[j2][0] / (NGRIDX*NGRIDY);
			double imag = PART_CHARGE * zcBuf[j2][1] / (NGRIDX*NGRIDY);
			//printf("%d\t%f,%f\n", j, real, imag);

			//rhok[j][0] = real * sk[j][0] - imag * sk[j][1];
			//rhok[j][1] = real * sk[j][1] + imag * sk[j][0];
			rhokBuf[j][0] = real;
			rhokBuf[j][1] = imag;
		}
	}

	// neutralizing background
	rhokBuf[NGRIDX/2 + (NGRIDY/2)*NGRIDX][0] = 0;
	rhokBuf[NGRIDX/2 + (NGRIDY/2)*NGRIDX][1] = 0;
}

// determines E and phi from rho
void fields(double *potential) {
	// rho(k) -> phi(k)
	for (int j = 0; j < NGRIDX * NGRIDY; j++) {
		double kx = 2 * M_PI * ((j % NGRIDX) - NGRIDX/2) / XMAX;
		double ky = 2 * M_PI * ((int)(j / NGRIDX) - NGRIDY/2) / YMAX;

		double ksqr = kx * kx + ky * ky;
		if (ksqr != 0) {
			phikBuf[j][0] = rhokBuf[j][0] / (ksqr * EPS_0);
			phikBuf[j][1] = rhokBuf[j][1] / (ksqr * EPS_0);
		} else {
			phikBuf[j][0] = 0;
			phikBuf[j][1] = 0;
		}
		
		double k = sqrt(ksqr);

		exkBuf[j][0] =  ky * phikBuf[j][1];
		exkBuf[j][1] = -ky * phikBuf[j][0];

		eykBuf[j][0] =  kx * phikBuf[j][1];
		eykBuf[j][1] = -kx * phikBuf[j][0];
	}
	/*
	arraydump_complex(rhokBuf, "rhok_re.csv", 0);
	arraydump_complex(rhokBuf, "rhok_im.csv", 1);
	arraydump_complex(phikBuf, "phik_re.csv", 0);
	arraydump_complex(phikBuf, "phik_im.csv", 1);
	exit(1);
	*/
	// potential energy calculation

	if (potential != NULL) {
		double pot = 0;
		for (int j = 0; j < NGRIDX * NGRIDY / 2; j++) {
			pot += phikBuf[j][0] * rhokBuf[j][0] + phikBuf[j][1] * rhokBuf[j][1];
		}
		*potential = pot * XMAX * YMAX;
	}
	/*
	fftw_execute(exIFFT);
	fftw_execute(eyIFFT);
	*/
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

// finds fields on particles from positions 
void interpField(double *x, double *y, fftw_complex *ek, fftw_complex *ePart) {
	int ncx = NGRIDX, ncy = NGRIDY;
	int np = PART_NUM;
	int isign = 1;
	int order = 5;

	for (int m = 0; m < PART_NUM; m++) {
		xpBuf[m] = x[m] / XMAX;
		ypBuf[m] = y[m] / YMAX;
	}

	// copy in E_k array
	for (int jx = 0; jx < NGRIDX; jx++) {
		for (int jy = 0; jy < NGRIDX; jy++) {
			int j = jx + NGRIDX*jy;
			int j2 = NGRIDY*jx + jy;

			zcBuf[j2][0] = ek[j][0];
			zcBuf[j2][1] = ek[j][1];
		}
	}

	// interpolate E_x on particles
	uf2a_(&ncx, &ncy, (double*)zcBuf,
	      &np, xpBuf, ypBuf, (double*)ePart,
	      &isign, &order);
}

// pushes particles forward by half a timestep
// ePart is complex because we're using USFFT output for the field
void vHalfPush(double *v, fftw_complex *ePart, int forward) {
	double factor = -DT/2 * (PART_CHARGE / PART_MASS);
	if (!forward) factor = -factor;
	for (int m = 0; m < PART_NUM; m++)
		v[m] += factor * ePart[m][0];
}

// moves particles given velocities
void xPush(double *x, double *y, double *vx, double *vy) {
	//#pragma omp parallel for
	for (int i = 0; i < PART_NUM; i++) {
		x[i] += DT * vx[i];
		y[i] += DT * vy[i];

		// periodicity
		// (not strictly correct, but if a particle is moving several grid
		// lengths in 1 timestep, something has gone horribly wrong)
		if (x[i] < 0) x[i] += XMAX;
		if (x[i] >= XMAX) x[i] -= XMAX;

		if (y[i] < 0) y[i] += YMAX;
		if (y[i] >= XMAX) y[i] -= YMAX;
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
