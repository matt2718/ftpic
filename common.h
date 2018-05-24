#ifndef _FTPIC_COMMON_H
#define _FTPIC_COMMON_H

#include <qdsp.h>

extern double DT;
extern double TMAX;

extern const double XMAX; // system length
extern int NGRID; // grid size
extern double DX;

// particle number and properties
extern int PART_NUM;
extern const double PART_MASS;
extern const double PART_CHARGE;
extern const double EPS_0;

extern const double BEAM_SPEED;
extern const double V_TH;

extern double OMEGA_P;

extern const int MODELOG_MAX;
extern FILE *modeLog;
extern int printTime;

int commonInit(int argc, char **argv,
               double **x, double **v, int **color,
               QDSPplot **phasePlot, QDSPplot **phiPlot, QDSPplot **rhoPlot);

void init2Stream(double *x, double *v, int *color);

void initLandau(double *x, double *v, int *color);

void initWave(double *x, double *v, int *color);

#endif
