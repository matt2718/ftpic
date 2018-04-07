#ifndef _FTPIC_COMMON_H
#define _FTPIC_COMMON_H

extern const double XMAX; // system length
extern const int NGRID; // grid size
extern double DX;

// particle number and properties
extern const int PART_NUM;
extern const double PART_MASS;
extern const double PART_CHARGE;
extern const double EPS_0;

extern const double BEAM_SPEED;

extern double OMEGA_P;

void init2Stream(double *x, double *v, int *color);

void initStanding(double *x, double *v, int *color);

void initLandau(double *x, double *v, int *color);

void initRhoSin(double *x, double *v, int *color);

#endif
